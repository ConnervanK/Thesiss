import h5py
import numpy as np
import matplotlib.pyplot as plt
import os

def trmusic_bscan_imaging():
    num_runs = 40
    num_rx = 11
    traces = None
    
    # Tx trajectory
    tx_start = 0.2
    tx_end = 0.8
    tx_x_arr = np.linspace(tx_start, tx_end, num_runs)
    
    # Relative Rx offsets
    offsets = np.linspace(-0.1, 0.1, num_rx)
    
    for tx_i in range(1, num_runs + 1):
        out_file = f'RTM/subwavelength_bscan_array{tx_i}.out'
        
        if not os.path.exists(out_file):
            print(f"File {out_file} not found. Run simulation first.")
            return

        with h5py.File(out_file, 'r') as f:
            if traces is None:
                dt = f.attrs['dt']
                num_time_steps = len(f['rxs']['rx1']['Ez'])
                traces = np.zeros((num_runs, num_rx, num_time_steps))
            
            rx_data = f['rxs']
            for rx_j in range(1, num_rx + 1):
                traces[tx_i-1, rx_j-1, :] = rx_data[f'rx{rx_j}']['Ez']
                
    # Direct wave removal via time-muting (assuming no baseline)
    t_arr = np.arange(traces.shape[2]) * dt
    mute_mask = t_arr < 2.5e-9
    for tx_i_idx in range(num_runs):
        for rx_j_idx in range(num_rx):
            traces[tx_i_idx, rx_j_idx, mute_mask] *= (0.5 * (1 - np.cos(np.pi * t_arr[mute_mask] / 2.5e-9)))

    print("Executing Time Reversal MUSIC for Multi-Receiver B-Scan...")
    
    c = 299792458.0
    eps_r = 6.0
    v = c / np.sqrt(eps_r)

    grid_x = np.linspace(0.3, 0.7, 100)
    grid_y = np.linspace(0.1, 0.4, 100)
    image = np.zeros((len(grid_y), len(grid_x)), dtype=float)

    D_f = np.fft.rfft(traces, axis=-1)
    freqs = np.fft.rfftfreq(traces.shape[2], dt)
    
    f_min, f_max = 0.5e9, 2.5e9
    valid_f_indices = np.where((freqs >= f_min) & (freqs <= f_max))[0]
    t_delay = 0.82e-9 
    num_scatterers = 2

    # Here, we treat the moving sub-array by forming a sliding sub-aperture MDM.
    # A single (1 x 11) data vector yields a rank-1 subspace, so 2 scatterers will smear.
    # By stacking neighboring Tx positions (e.g., 5 Tx positions), we get a 5x11 local MDM.
    # The local MDM has enough rank to separate the targets and preserves the near-paraxial
    # block factorization of the Green's function.
    sub_ap_size = 5
    half_sub = sub_ap_size // 2

    for f_idx in valid_f_indices:
        f = freqs[f_idx]
        k = 2 * np.pi * f / v 
        I_f = np.zeros((len(grid_y), len(grid_x)), dtype=float)
        
        # Evaluate Sub-Aperture TR-MUSIC 
        for tx_i_idx in range(half_sub, num_runs - half_sub):
            # Local Multistatic Data Matrix for Sub-aperture (sub_ap_size x num_rx)
            D_mat = D_f[tx_i_idx - half_sub : tx_i_idx + half_sub + 1, :, f_idx] * np.exp(1j * 2 * np.pi * f * t_delay)
            
            # SVD on local MDM
            # D_mat size: 5 x 11
            U, S, Vh = np.linalg.svd(D_mat, full_matrices=True)
            
            # Signal subspace rank up to num_scatterers=2
            # Tx noise subspace Un: 5 x (5-2) = 5 x 3
            Un = U[:, num_scatterers:]
            # Rx noise subspace Vn: 11 x (11-2) = 11 x 9
            Vn = Vh[num_scatterers:, :].T
            
            tx_x_local = tx_x_arr[tx_i_idx - half_sub : tx_i_idx + half_sub + 1]
            tx_x_center = tx_x_arr[tx_i_idx]
            rx_x_local = tx_x_center + offsets
            
            for y_idx, z in enumerate(grid_y):
                X_focal = grid_x
                
                # Tx steering vector (5 x Nx)
                d_tx = np.sqrt((X_focal[:, None] - tx_x_local[None, :])**2 + (z - 0.45)**2)
                g_tx = np.exp(-1j * k * d_tx) / np.sqrt(d_tx)
                g_tx = g_tx / (np.linalg.norm(g_tx, axis=1, keepdims=True) + 1e-15)
                
                # Rx steering vector (11 x Nx)
                d_rx = np.sqrt((X_focal[:, None] - rx_x_local[None, :])**2 + (z - 0.45)**2)
                g_rx = np.exp(-1j * k * d_rx) / np.sqrt(d_rx)
                g_rx = g_rx / (np.linalg.norm(g_rx, axis=1, keepdims=True) + 1e-15)
                
                # Projections onto Noise subspaces
                proj_tx = np.sum(np.abs(np.conj(g_tx) @ Un)**2, axis=1)
                proj_rx = np.sum(np.abs(np.conj(g_rx) @ Vn)**2, axis=1)
                
                # Double-sided projection pseudo-spectrum
                I_f[y_idx, :] += 1.0 / (proj_tx + proj_rx + 1e-15)
                
        image += I_f / (len(valid_f_indices) * (num_runs - 2*half_sub))

    # Convert to dB
    image = 10 * np.log10(image / np.max(image))

    plt.figure(figsize=(10, 6))
    plt.imshow(image, extent=[grid_x[0], grid_x[-1], grid_y[0], grid_y[-1]], 
               origin='lower', aspect='auto', cmap='plasma')
    
    circle1 = plt.Circle((0.47, 0.25), 0.01, color='lime', fill=False, lw=2)
    circle2 = plt.Circle((0.53, 0.25), 0.01, color='lime', fill=False, lw=2)
    plt.gca().add_patch(circle1)
    plt.gca().add_patch(circle2)
    
    plt.title('High-Resolution TR-MUSIC (Moving SIMO B-Scan)')
    plt.xlabel('X (m)')
    plt.ylabel('Y / Depth (m)')
    plt.tight_layout()
    
    out_img = 'RTM/trmusic_result_bscan_array.png'
    plt.savefig(out_img, dpi=300)
    print(f"Saved generated image to {out_img}")

if __name__ == '__main__':
    trmusic_bscan_imaging()