import os
import sys
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

os.environ['DEVITO_ARCH'] = 'gcc'

if not hasattr(os, 'getuid'):
    os.getuid = lambda: 1000

# Devito patch for Windows memory allocation
import devito.data.allocators as allocators
class Win32NumpyAllocator(allocators.MemoryAllocator):
    def alloc(self, shape, dtype, padding=0):
        data = np.zeros(shape, dtype=dtype)
        def memfree(): pass
        return data, {"memfree": memfree}
    def free(self, *args, **kwargs): pass

if sys.platform == 'win32':
    allocators.PosixAllocator.alloc = Win32NumpyAllocator.alloc
    allocators.PosixAllocator.free = Win32NumpyAllocator.free

# Patch codepy string compiler to fix Windows paths for GCC
import devito.arch.compiler as comp
orig_compile = comp.compile_from_string
def safe_win_compile(*args, **kwargs):
    args = [str(a).replace('\\', '/') if isinstance(a, (str, os.PathLike)) else a for a in args]
    return orig_compile(*args, **kwargs)
comp.compile_from_string = safe_win_compile

from devito import configuration, Grid, TimeFunction, Function, Eq, solve, Operator, SparseTimeFunction
configuration['log-level'] = 'WARNING'
configuration['language'] = 'C'

def main():
    # 1. Setup Grid and True Model
    shape = (100, 100)
    spacing = (0.01, 0.01) # 1m x 1m overall
    origin = (0., 0.)
    grid = Grid(shape=shape, extent=(shape[0]*spacing[0], shape[1]*spacing[1]), origin=origin)

    v_bg = 2.0 # km/s or similar units (scaled for devito)
    v_true = np.full(shape, v_bg)
    
    # Scatterer 1: fast anomaly
    v_true[20:40, 40:60] = 2.5
    # Scatterer 2: slow anomaly
    v_true[60:80, 50:70] = 1.5

    # 2. Source and Receiver Geometry
    nshots = 8
    nreceivers = 50
    
    tx_coords = np.zeros((nshots, 2))
    tx_coords[:, 0] = np.linspace(0.1, 0.9, nshots)
    tx_coords[:, 1] = 0.05 # Near the top

    rx_coords = np.zeros((nreceivers, 2))
    rx_coords[:, 0] = np.linspace(0.05, 0.95, nreceivers)
    rx_coords[:, 1] = 0.05 # Near the top

    # 3. Source Time Function (Ricker)
    f0 = 0.1 # 100 Hz center freq (in typical Devito kHz/ms scale = 0.1)
    dt = 0.001
    nt = 800
    time_range = np.arange(nt) * dt
    src_time = (1.0 - 2.0 * (np.pi * f0 * (time_range - 1.0/f0))**2) * np.exp(-(np.pi * f0 * (time_range - 1.0/f0))**2)

    # 4. Define Operators
    def create_operators(grid, nt):
        u = TimeFunction(name='u', grid=grid, time_order=2, space_order=4)
        m = Function(name='m', grid=grid, space_order=4)
        
        # Forward PDE
        pde = m * u.dt2 - u.laplace
        stencil = Eq(u.forward, solve(pde, u.forward))
        
        src = SparseTimeFunction(name='src', grid=grid, nt=nt, npoint=1)
        rec = SparseTimeFunction(name='rec', grid=grid, nt=nt, npoint=nreceivers)
        
        src_term = src.inject(field=u.forward, expr=src * dt**2 / m)
        rec_term = rec.interpolate(expr=u.forward)
        op_fwd = Operator([stencil] + src_term + rec_term, subs=grid.spacing_map, name="Forward")
        
        # Adjoint PDE
        v = TimeFunction(name='v', grid=grid, time_order=2, space_order=4)
        stencil_adj = Eq(v.backward, solve(m * v.dt2 - v.laplace, v.backward))
        rec_adj = SparseTimeFunction(name='rec_adj', grid=grid, nt=nt, npoint=nreceivers)
        rec_adj_term = rec_adj.inject(field=v.backward, expr=rec_adj * dt**2 / m)
        
        # Gradient
        grad = Function(name='grad', grid=grid)
        grad_update = Eq(grad, grad - u * v.dt2)
        op_adj = Operator([stencil_adj] + rec_adj_term + [grad_update], subs=grid.spacing_map, name="Adjoint")
        
        return m, u, v, src, rec, rec_adj, grad, op_fwd, op_adj

    print("Compiling FWI operators... this may take a moment.")
    m, u, v, src, rec, rec_adj, grad, op_fwd, op_adj = create_operators(grid, nt)

    # 5. Generate Synthetic "True" Data
    print("Generating synthetic observed data from true model...")
    m.data[:] = 1.0 / (v_true**2)
    d_obs = np.zeros((nshots, nreceivers, nt))

    for shot in range(nshots):
        u.data[:] = 0.0
        src.coordinates.data[0, :] = tx_coords[shot, :]
        src.data[:, 0] = src_time
        rec.coordinates.data[:] = rx_coords
        
        op_fwd.apply(time_M=nt-1, dt=dt)
        d_obs[shot, :, :] = rec.data.T
        
    print("True data generation completed.")

    # 6. FWI Wrapper Function
    def fwi_wrapper(m_flat):
        m.data[:] = m_flat.reshape(shape)
        grad.data[:] = 0.0
        F = 0.0
        
        for shot in range(nshots):
            # Reset fields
            u.data[:] = 0.0
            v.data[:] = 0.0
            
            # Forward pass
            src.coordinates.data[0, :] = tx_coords[shot, :]
            src.data[:, 0] = src_time
            rec.coordinates.data[:] = rx_coords
            op_fwd.apply(time_M=nt-1, dt=dt)
            
            # Calculate residual
            res = rec.data - d_obs[shot, :, :].T
            F += 0.5 * np.linalg.norm(res)**2
            
            # Adjoint pass
            rec_adj.coordinates.data[:] = rx_coords
            rec_adj.data[:] = res
            op_adj.apply(time_M=nt-1, dt=dt)
            
        print(f"Objective Function Value: {F:.4e}")
        return F, grad.data.flatten()

    # 7. Run FWI Optimization
    print("Starting FWI optimization loop...")
    m0 = np.full(shape, 1.0 / (v_bg**2)) # Start with homogeneous background
    bounds = [(1.0/(3.5**2), 1.0/(1.0**2)) for _ in range(np.prod(shape))]
    
    res = optimize.minimize(fwi_wrapper,
                            m0.flatten(),
                            method='L-BFGS-B',
                            jac=True,
                            bounds=bounds,
                            options={'maxiter': 20, 'disp': True})

    # Reshape optimized model back to physical velocity
    v_inv = 1.0 / np.sqrt(res.x.reshape(shape))

    # 8. Plot Setup
    print("Plotting results...")
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    
    vmin, vmax = 1.0, 3.0
    
    im1 = ax[0].imshow(v_true.T, cmap='viridis', extent=[0, 1.0, 1.0, 0], vmin=vmin, vmax=vmax)
    ax[0].set_title('True Model')
    ax[0].set_xlabel('X (km)')
    ax[0].set_ylabel('Depth (km)')
    plt.colorbar(im1, ax=ax[0])

    m0_v = 1.0 / np.sqrt(m0)
    im2 = ax[1].imshow(m0_v.T, cmap='viridis', extent=[0, 1.0, 1.0, 0], vmin=vmin, vmax=vmax)
    ax[1].set_title('Initial Background Model')
    ax[1].set_xlabel('X (km)')
    plt.colorbar(im2, ax=ax[1])

    im3 = ax[2].imshow(v_inv.T, cmap='viridis', extent=[0, 1.0, 1.0, 0], vmin=vmin, vmax=vmax)
    ax[2].set_title('FWI Result (Inverted Model)')
    ax[2].set_xlabel('X (km)')
    plt.colorbar(im3, ax=ax[2])

    plt.tight_layout()
    plt.savefig('Full_Waveform_Inversion/fwi_dummy_result.png')
    print("Done! Check 'Full_Waveform_Inversion/fwi_dummy_result.png' for the output diagram.")

if __name__ == '__main__':
    main()
