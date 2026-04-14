import numpy as np
from devito import Grid, TimeFunction, Function, Eq, solve, Operator, SparseTimeFunction, configuration
import sys
import os
import matplotlib.pyplot as plt
from scipy import optimize

configuration['log-level'] = 'INFO'
configuration['language'] = 'C'

os.environ['DEVITO_ARCH'] = 'gcc'

if not hasattr(os, 'getuid'):
    os.getuid = lambda: 1000

import devito.data.allocators as allocators

class Win32MsvcrtAllocator(allocators.MemoryAllocator):
    @classmethod
    def initialize(cls):
        import ctypes
        try:
            cls.lib = ctypes.cdll.msvcrt
        except OSError:
            cls.lib = None

    def _alloc_C_libcall(self, size, ctype):
        import ctypes
        c_bytesize = size * ctypes.sizeof(ctype)
        alignment = 64
        
        lib = ctypes.cdll.msvcrt
        lib._aligned_malloc.restype = ctypes.c_void_p
        lib._aligned_malloc.argtypes = [ctypes.c_size_t, ctypes.c_size_t]
        
        ptr = lib._aligned_malloc(c_bytesize, alignment)
        if ptr:
            c_pointer = ctypes.cast(ptr, ctypes.c_void_p)
            return c_pointer, (c_pointer, )
        else:
            return None, None

    def free(self, c_pointer):
        import ctypes
        lib = ctypes.cdll.msvcrt
        lib._aligned_free.argtypes = [ctypes.c_void_p]
        lib._aligned_free(c_pointer)

if sys.platform == 'win32':
    allocators.PosixAllocator._alloc_C_libcall = Win32MsvcrtAllocator._alloc_C_libcall
    allocators.PosixAllocator.free = Win32MsvcrtAllocator.free
    allocators.PosixAllocator.initialize = Win32MsvcrtAllocator.initialize
    allocators.PosixAllocator._attempted_init = False

# Patch codepy string compiler to fix Windows paths for GCC
import devito.arch.compiler as comp
orig_compile = comp.compile_from_string
def safe_win_compile(*args, **kwargs):
    args = [str(a).replace('\\', '/') if isinstance(a, (str, os.PathLike)) else a for a in args]
    return orig_compile(*args, **kwargs)
comp.compile_from_string = safe_win_compile


def run_em_fwi_dummy():
    # 1. Setup Grid and True Model
    shape = (200, 200)
    spacing = (0.01, 0.01) # 2m x 2m overall to avoid boundary reflections entirely over 8ns
    origin = (0., 0.)
    grid = Grid(shape=shape, extent=(shape[0]*spacing[0], shape[1]*spacing[1]), origin=origin)

    # Note: Devito solves Acoustic, which behaves identically to 2D TE EM field when:
    # m = epsilon * mu. Let's use velocity v = c / sqrt(er) -> v_bg = 0.3 / sqrt(1) = 0.3 m/ns for Air
    # Let's say background is 0.15 m/ns
    v_bg = 0.15 
    v_true = np.full(shape, v_bg)
    
    # Scatterer 1: fast anomaly
    v_true[70:90, 90:110] = 0.25 # Air cavity maybe
    # Scatterer 2: slow anomaly
    v_true[110:130, 100:120] = 0.05 # Water/metal

    # 2. Source and Receiver Geometry
    nshots = 10
    nreceivers_side = 50
    nreceivers = nreceivers_side * 2
    
    tx_coords = np.zeros((nshots, 2))
    tx_coords[:, 0] = np.linspace(0.5, 1.5, nshots)
    tx_coords[:, 1] = 0.5 # Near the top

    rx_coords = np.zeros((nreceivers, 2))
    # Top array
    rx_coords[:nreceivers_side, 0] = np.linspace(0.4, 1.6, nreceivers_side)
    rx_coords[:nreceivers_side, 1] = 0.5
    # Bottom array (provides transmission data needed to perfectly resolve the full shape)
    rx_coords[nreceivers_side:, 0] = np.linspace(0.4, 1.6, nreceivers_side)
    rx_coords[nreceivers_side:, 1] = 1.5

    # 3. Source Time Function
    f0 = 1.0 # 1.0 GHz
    dt = 0.01
    nt = 800
    time_range = np.arange(nt) * dt
    src_time = (1.0 - 2.0 * (np.pi * f0 * (time_range - 1.0/f0))**2) * np.exp(-(np.pi * f0 * (time_range - 1.0/f0))**2)

    # 4. Operators
    def create_operators(grid, nt):
        # critical fix: save=nt is required to compute full wavefield gradients correctly!
        u = TimeFunction(name='u', grid=grid, time_order=2, space_order=4, save=nt)
        m = Function(name='m', grid=grid, space_order=4)
        damp = Function(name='damp', grid=grid)
        
        # Initialize simple sponge layer damping
        pad = 20
        damp_profile = np.linspace(0, 50.0, pad)**2
        for i in range(pad):
            damp.data[i, :] = np.maximum(damp.data[i, :], damp_profile[-1-i])
            damp.data[-1-i, :] = np.maximum(damp.data[-1-i, :], damp_profile[-1-i])
            damp.data[:, i] = np.maximum(damp.data[:, i], damp_profile[-1-i])
            damp.data[:, -1-i] = np.maximum(damp.data[:, -1-i], damp_profile[-1-i])
        
        # Forward PDE
        pde = m * u.dt2 + damp * u.dt - u.laplace
        stencil = Eq(u.forward, solve(pde, u.forward))
        
        src = SparseTimeFunction(name='src', grid=grid, nt=nt, npoint=1)
        rec = SparseTimeFunction(name='rec', grid=grid, nt=nt, npoint=nreceivers)
        
        src_term = src.inject(field=u.forward, expr=src * dt**2 / m)
        rec_term = rec.interpolate(expr=u.forward)
        op_fwd = Operator([stencil] + src_term + rec_term, subs=grid.spacing_map, name="Forward")
        
        # Adjoint PDE
        v = TimeFunction(name='v', grid=grid, time_order=2, space_order=4)
        # Note: adjoint PDE propagates backwards in time and sign of damp is reversed
        stencil_adj = Eq(v.backward, solve(m * v.dt2 - damp * v.dt - v.laplace, v.backward))
        rec_adj = SparseTimeFunction(name='rec_adj', grid=grid, nt=nt, npoint=nreceivers)
        rec_adj_term = rec_adj.inject(field=v.backward, expr=rec_adj * dt**2 / m)
        
        # Gradient
        grad = Function(name='grad', grid=grid)
        grad_update = Eq(grad, grad - u * v.dt2)
        op_adj = Operator([stencil_adj] + rec_adj_term + [grad_update], subs=grid.spacing_map, name="Adjoint")
        
        return m, u, v, damp, src, rec, rec_adj, grad, op_fwd, op_adj

    print("Compiling FWI operators...")
    m, u, v, damp, src, rec, rec_adj, grad, op_fwd, op_adj = create_operators(grid, nt)

    # 5. Generate True Data using Devito instead of gprMax
    print("Generating synthetic data exclusively through Devito...")
    
    # We must be careful because of Windows patch. Try native array copy:
    m.data[:] = 1.0 / (v_true**2)
    d_obs = np.zeros((nshots, nreceivers, nt))

    for shot in range(nshots):
        u.data[:] = 0.0
        src.coordinates.data[0, :] = tx_coords[shot, :]
        src.data[:, 0] = src_time
        rec.coordinates.data[:] = rx_coords
        
        op_fwd.apply(time_M=nt-2, dt=dt)
        d_obs[shot, :, :] = rec.data.copy().T
        
    print("Devito Forward Modeling Completed.")

    # 6. FWI Run Loop
    def fwi_wrapper(m_flat):
        # We have to be careful with the monkey patch memory updating logic. 
        # Devito operators on Window Numpy array patch might not register C pointer updates perfectly 
        # across scopes if replaced entirely. We insert explicitly:
        m.data[:] = m_flat.reshape(shape)
        grad.data[:] = 0.0
        F = 0.0
        
        # Keep a persistent numpy gradient array because devito func doesn't persist data when patched sometimes
        np_grad = np.zeros(shape)
        
        for shot in range(nshots):
            u.data[:] = 0.0
            v.data[:] = 0.0
            
            src.coordinates.data[0, :] = tx_coords[shot, :]
            src.data[:, 0] = src_time
            rec.coordinates.data[:] = rx_coords
            op_fwd.apply(time_M=nt-2, dt=dt)
            
            res = rec.data - d_obs[shot, :, :].T
            if shot == 0:
                print(f"Norm rec.data: {np.linalg.norm(rec.data)}")
                print(f"Norm d_obs: {np.linalg.norm(d_obs[shot, :, :].T)}")
                print(f"Norm res: {np.linalg.norm(res)}")
            F += 0.5 * np.linalg.norm(res)**2
            
            rec_adj.coordinates.data[:] = rx_coords
            rec_adj.data[:] = res
            op_adj.apply(time_M=nt-2, dt=dt)
            
            if shot == 0:
                print(f"Max abs u: {np.max(np.abs(u.data))}")
                print(f"Max abs v: {np.max(np.abs(v.data))}")
                print(f"Max abs grad: {np.max(np.abs(grad.data))}")

            np_grad += grad.data
            grad.data[:] = 0.0
            
        print(f"Objective Function Value: {F:.4e}")
        
        # Scale objective and gradient to avoid numerical underflow in L-BFGS-B
        scale = 1e15
        
        # Additional trick: FWI gradients at the source point explode (singularities).
        # We should mask out the top 20 rows of the gradient to prevent source artifacts!
        np_grad[0:20, :] = 0.0
        
        return F * scale, np_grad.flatten() * scale

    print("Starting FWI optimization loop...")
    m_initial = np.full(shape, 1.0 / (v_bg**2))
    m0 = m_initial.copy()
    
    from scipy.ndimage import gaussian_filter

    # Manual Gradient Descent for 50 iterations
    for i in range(50):
        print(f"--- Iteration {i+1} ---")
        F, grad_val = fwi_wrapper(m0.flatten())
        
        # Smooth the gradient to remove source/receiver footprints and artifacts
        grad_reshaped = grad_val.reshape(shape)
        grad_smooth = gaussian_filter(grad_reshaped, sigma=1.0)
        grad_val = grad_smooth.flatten()
        
        m0_flat = m0.flatten()
        max_grad = np.max(np.abs(grad_val)) + 1e-12
        
        # Step size (learning rate)
        alpha = 2.0 / max_grad 
        
        m0_new = m0_flat - alpha * grad_val
        m0_new = np.clip(m0_new, 1.0/(0.3**2), 1.0/(0.01**2))
        m0 = m0_new.reshape(shape)

    v_inv = 1.0 / np.sqrt(m0)

    print("Plotting results...")
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    
    vmin, vmax = 0.05, 0.25
    
    im1 = ax[0].imshow(v_true.T, cmap='viridis', extent=[0, 2.0, 2.0, 0], vmin=vmin, vmax=vmax)
    ax[0].set_title('True Model')
    plt.colorbar(im1, ax=ax[0])

    m0_v = 1.0 / np.sqrt(m_initial)
    im2 = ax[1].imshow(m0_v.T, cmap='viridis', extent=[0, 2.0, 2.0, 0], vmin=vmin, vmax=vmax)
    ax[1].set_title('Initial Background Model')
    plt.colorbar(im2, ax=ax[1])

    im3 = ax[2].imshow(v_inv.T, cmap='viridis', extent=[0, 2.0, 2.0, 0], vmin=vmin, vmax=vmax)
    ax[2].set_title('FWI Result')
    plt.colorbar(im3, ax=ax[2])

    plt.tight_layout()
    plt.savefig('Full_Waveform_Inversion/fwi_devito_only.png')
    print("Done! Check 'Full_Waveform_Inversion/fwi_devito_only.png' ")


if __name__ == '__main__':
    run_em_fwi_dummy()
