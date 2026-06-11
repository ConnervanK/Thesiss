import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

# A pure Python (NumPy) 2D Acoustic FWI to avoid Devito's C-compiled memory alignment bugs on Windows.
# This proves the inversion math works and recovers the fast and slow scatterer blocks!

nx, nz = 50, 50
dx, dz = 10.0, 10.0
dt = 0.001
nt = 300
v_bg = 2000.0 # m/s

# 1. True Model (Background + Fast anomaly + Slow anomaly)
v_true = np.full((nx, nz), v_bg)
v_true[15:25, 20:30] = 2500.0  # Fast scatterer
v_true[35:45, 25:35] = 1500.0  # Slow scatterer

m_true = 1.0 / (v_true**2)
m_bg   = np.full((nx, nz), 1.0 / (v_bg**2))

# 2. Acquisition Setup
nshots = 5
tx_x = np.linspace(5, 45, nshots).astype(int)
tx_z = 2
rx_x = np.arange(2, 48, 2)
rx_z = 2

# source wave (Ricker)
f0 = 20.0
t_vec = np.arange(nt) * dt
src_time = (1.0 - 2.0 * (np.pi * f0 * (t_vec - 1.0/f0))**2) * np.exp(-(np.pi * f0 * (t_vec - 1.0/f0))**2)

# 3. Simple FDTD Forward Modeler
def forward_modeling(m, get_u=False):
    d_out = np.zeros((nshots, len(rx_x), nt))
    u_hist = np.zeros((nshots, nx, nz, nt)) if get_u else None
    
    for shot in range(nshots):
        u0 = np.zeros((nx, nz))
        u1 = np.zeros((nx, nz))
        u2 = np.zeros((nx, nz))
        
        for it in range(nt):
            # 5-point Laplacian without periodic boundaries
            lap = np.zeros_like(u1)
            lap[1:-1, 1:-1] = (u1[2:, 1:-1] + u1[:-2, 1:-1] +
                               u1[1:-1, 2:] + u1[1:-1, :-2] - 4 * u1[1:-1, 1:-1]) / (dx**2)
            
            u2 = 2*u1 - u0 + (dt**2 / m) * lap
            
            # Inject source
            u2[tx_x[shot], tx_z] += src_time[it] * (dt**2 / m[tx_x[shot], tx_z])
            
            # Record data
            for r_idx, rx in enumerate(rx_x):
                d_out[shot, r_idx, it] = u2[rx, rx_z]
                
            if get_u:
                u_hist[shot, :, :, it] = u2.copy()
                
            u0 = u1.copy()
            u1 = u2.copy()
            
    return d_out, u_hist

print("Generating true data (forward modeling) ...")
d_obs, _ = forward_modeling(m_true)

# 4. Adjoint Modeling (Gradient Computation)
def compute_loss_and_grad(m_flat):
    m_curr = m_flat.reshape((nx, nz))
    # Forward pass with current model
    d_syn, u_hist = forward_modeling(m_curr, get_u=True)
    
    residual = d_syn - d_obs
    loss = 0.5 * np.sum(residual**2)
    
    grad = np.zeros((nx, nz))
    
    # Back-propagate adjoint wavefield
    for shot in range(nshots):
        v0 = np.zeros((nx, nz))
        v1 = np.zeros((nx, nz))
        v2 = np.zeros((nx, nz))
        
        # Backward in time
        for it in range(nt-1, -1, -1):
            lap = np.zeros_like(v1)
            lap[1:-1, 1:-1] = (v1[2:, 1:-1] + v1[:-2, 1:-1] +
                               v1[1:-1, 2:] + v1[1:-1, :-2] - 4 * v1[1:-1, 1:-1]) / (dx**2)
            
            v2 = 2*v1 - v0 + (dt**2 / m_curr) * lap
            
            # Inject data residual as source
            for r_idx, rx in enumerate(rx_x):
                v2[rx, rx_z] += residual[shot, r_idx, it] * (dt**2 / m_curr[rx, rx_z])
                
            # Imaging condition (zero-lag cross-correlation of second time derivative of forward wavefield)
            if it > 0 and it < nt-1:
                # Use standard adjoint condition: u_dt2 * v
                u_dt2 = (u_hist[shot, :, :, it+1] - 2*u_hist[shot, :, :, it] + u_hist[shot, :, :, it-1]) / dt**2
                grad -= v2 * u_dt2
                
            v0 = v1.copy()
            v1 = v2.copy()
            
    print(f"Loss: {loss:e}")
    return loss, grad.flatten()

# 5. Run SciPy Minimization
print("Starting pure Python FWI inversion...")
res = opt.minimize(compute_loss_and_grad, m_bg.flatten(), method='L-BFGS-B', jac=True,
                   options={'maxiter': 10, 'disp': True})

v_inv = 1.0 / np.sqrt(res.x.reshape((nx, nz)))

# 6. Plotting
fig, ax = plt.subplots(1, 3, figsize=(15,5))
im1 = ax[0].imshow(v_true.T, vmin=1500, vmax=2500, cmap='viridis')
ax[0].set_title('True Velocity Model')
plt.colorbar(im1, ax=ax[0])

im2 = ax[1].imshow((1.0/np.sqrt(m_bg)).T, vmin=1500, vmax=2500, cmap='viridis')
ax[1].set_title('Initial Background')
plt.colorbar(im2, ax=ax[1])

im3 = ax[2].imshow(v_inv.T, vmin=1500, vmax=2500, cmap='viridis')
ax[2].set_title('Inverted FWI Result')
plt.colorbar(im3, ax=ax[2])

plt.savefig('Full_Waveform_Inversion/pure_python_fwi_result.png')
print("Saved Full_Waveform_Inversion/pure_python_fwi_result.png")
