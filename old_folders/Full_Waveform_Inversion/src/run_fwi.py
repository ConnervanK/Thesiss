import os
import sys
import numpy as np

if not hasattr(os, 'getuid'):
    os.getuid = lambda: 1000

os.environ['DEVITO_ARCH'] = 'gcc'

# Patch Devito to use simple numpy arrays for memory allocation on Windows
import devito.data.allocators as allocators
class Win32NumpyAllocator(allocators.MemoryAllocator):
    def alloc(self, shape, dtype, padding=0):
        # We must return an ndarray and a dictionary with 'memfree'
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

from devito import configuration
configuration['log-level'] = 'INFO'
configuration['language'] = 'C'

from devito import Grid, TimeFunction, Function, Eq, solve, Operator
from devito import SparseTimeFunction
from scipy import optimize

from old_folders.Full_Waveform_Inversion.src.fwi_prepare_data import prepare_fwi_data

# Grid settings matching GprMax domain
SHAPE = (500, 250)      # 1.0/0.002 = 500, 0.5/0.002 = 250
EXTENT = (1.0, 0.5)     # Physical size in meters
SPACING = (0.002, 0.002)
ORIGIN = (0., 0.)
V_BG = 0.1224 # c / sqrt(6) (ice)

def create_fwi_operators(grid, m_initial, nt, rx_coords, tx_coords):
    """
    Creates and returns Devito operators for forward and adjoint passes 
    """
    # 1. Forward Wavefield
    u = TimeFunction(name='u', grid=grid, time_order=2, space_order=4)
    m = Function(name='m', grid=grid, space_order=4)
    m.data[:] = m_initial
    
    pde = m * u.dt2 - u.laplace
    stencil = Eq(u.forward, solve(pde, u.forward))
    
    # Define source and receiver geometries dynamically via names
    src = SparseTimeFunction(name='src', grid=grid, nt=nt, npoint=1)
    rec = SparseTimeFunction(name='rec', grid=grid, nt=nt, npoint=len(rx_coords))
    
    src_term = src.inject(field=u.forward, expr=src * grid.time_dim.spacing**2 / m)
    rec_term = rec.interpolate(expr=u.forward)
    
    op_fwd = Operator([stencil] + src_term + rec_term, subs=grid.spacing_map, name="Forward")
    
    # 2. Adjoint Wavefield and Gradient
    v = TimeFunction(name='v', grid=grid, time_order=2, space_order=4)
    
    # Adjoint PDE (backwards in time)
    pde_adj = m * v.dt2 - v.laplace
    stencil_adj = Eq(v.backward, solve(pde_adj, v.backward))
    
    rec_adj = SparseTimeFunction(name='rec_adj', grid=grid, nt=nt, npoint=len(rx_coords))
    rec_adj_term = rec_adj.inject(field=v.backward, expr=rec_adj * grid.time_dim.spacing**2 / m)
    
    # Gradient computation
    grad = Function(name='grad', grid=grid)
    grad_update = Eq(grad, grad - u * v.dt2)
    
    op_adj = Operator([stencil_adj] + rec_adj_term + [grad_update], subs=grid.spacing_map, name="Adjoint")
    
    return m, op_fwd, op_adj, u, v, src, rec, rec_adj, grad


def loss_and_gradient(m_update, grid, dat_obs, dt_val, rx_coords, tx_coords, f0,
                      m, op_fwd, op_adj, u, v, src, rec, rec_adj, grad):
    """
    Computes objective F and gradient grad_F for current model flattened m_update.
    """
    m_full = m_update.reshape(grid.shape)
    m.data[:] = m_full
    
    if hasattr(grad, 'data') and grad.data is not None:
        grad.data[:] = 0.0 
        
    F = 0.0
    nt = dat_obs.shape[2]
    time_range = np.arange(nt) * dt_val
    src_time = (1.0 - 2.0 * (np.pi * f0 * (time_range - 1.0/f0))**2) * \
               np.exp(-(np.pi * f0 * (time_range - 1.0/f0))**2)
               
    for shot in range(dat_obs.shape[0]):
        # Clear fields
        if hasattr(u, 'data'): u.data[:] = 0.0
        if hasattr(v, 'data'): v.data[:] = 0.0
        
        # Setup source
        src.coordinates.data[0, :] = tx_coords[shot, :]
        src.data[:, 0] = src_time
        rec.coordinates.data[:] = rx_coords
        
        # Run forward
        op_fwd(time_M=nt-1, dt=dt_val)
        
        # Data residual
        res = rec.data - dat_obs[shot, :, :].T
        F += .5 * np.linalg.norm(res)**2
        
        # Adjoint backwards
        rec_adj.coordinates.data[:] = rx_coords
        rec_adj.data[:] = res
        
        op_adj(time_M=nt-1, dt=dt_val)
        
    return F, grad.data.flatten()


def main():
    d_obs, dt, t_new = prepare_fwi_data("RTM/subwavelength_scatters*.out", 20, 20)
    print(f"Loaded {d_obs.shape[0]} shots, shape: {d_obs.shape}")
    
    tx_coords = np.zeros((20, 2))
    tx_coords[:, 0] = np.linspace(0.1, 0.9, 20)
    tx_coords[:, 1] = 0.45
    
    rx_coords = np.zeros((20, 2))
    rx_coords[:, 0] = np.linspace(0.1, 0.9, 20)
    rx_coords[:, 1] = 0.45
    
    def fwi_wrapper(m_flat):
        import devito
        return loss_and_gradient(m_flat, grid, d_obs, dt, rx_coords, tx_coords, 1.5,
                                 m, op_fwd, op_adj, u, v, src, rec, rec_adj, grad)
        
    m0 = np.full(SHAPE, 1.0 / (V_BG**2))
    grid = Grid(shape=SHAPE, extent=EXTENT, origin=ORIGIN)
    
    print("Compiling Devito operators... This might take a minute.")
    nt = d_obs.shape[2]
    m, op_fwd, op_adj, u, v, src, rec, rec_adj, grad = create_fwi_operators(grid, m0, nt, rx_coords, tx_coords)
    
    v_max = 0.3 # air
    v_min = 0.033 # water
    bounds = [(1.0/(v_max**2), 1.0/(v_min**2)) for _ in range(np.prod(SHAPE))]
    
    print("Starting L-BFGS-B Optimization...")
    res = optimize.minimize(fwi_wrapper,
                            m0.flatten(),
                            method='L-BFGS-B',
                            jac=True,
                            bounds=bounds,
                            options={'maxiter': 5, 'disp': True})
                            
    np.save('Full_Waveform_Inversion/fwi_result.npy', res.x.reshape(SHAPE))
    print("Optimization Completed!")

if __name__ == '__main__':
    main()
