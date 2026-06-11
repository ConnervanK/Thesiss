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

    # Direct wave removal via Background Subtraction
    for rx_j_idx in range(num_rx):
        mean_trace = np.mean(traces[:, rx_j_idx, :], axis=0)
        for tx_i_idx in range(num_runs):
            traces[tx_i_idx, rx_j_idx, :] -= mean_trace


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
    
    # Use Broadband Covariance Matrix method.
    # Stack all Tx-Rx traces into a single measurement vector per frequency.
    # Over multiple frequencies, we build a spatial covariance matrix, which 
    # provides enough snapshots to resolve coherent scatterers and properly
    # handles arbitrary moving Tx/Rx geometries (unlike strict MDM SVD which requires fixed arrays).
    
    num_elements = num_runs * num_rx
    R_cov = np.zeros((num_elements, num_elements), dtype=complex)
    
    # Pre-build data vectors
    d_f_all = np.zeros((num_elements, len(valid_f_indices)), dtype=complex)
    
    for idx, f_idx in enumerate(valid_f_indices):
        f = freqs[f_idx]
        d_vec = D_f[:, :, f_idx].flatten() * np.exp(1j * 2 * np.pi * f * t_delay)
        d_f_all[:, idx] = d_vec
        R_cov += np.outer(d_vec, np.conj(d_vec))
        
    R_cov /= len(valid_f_indices)
    
    U, S, Vh = np.linalg.svd(R_cov, full_matrices=True)
    # broadband covariance MUSIC often needs more "signal" dimensions for two 
    # spatially distributed coherent pulses. Try 2 or 3 targets.
    # We revert the numbers back to 2 targets.
    num_scatterers = 2
    Un = U[:, num_scatterers:]
    
    # Pseudo-Spectrum computation over spatial grid
    # Steering vector will be size 440x1 at each point
    tx_x_flat = np.repeat(tx_x_arr, num_rx)
    rx_x_flat = np.zeros(num_elements)
    for i in range(num_runs):
        rx_x_flat[i*num_rx : (i+1)*num_rx] = tx_x_arr[i] + offsets

    I_map = np.zeros((len(grid_y), len(grid_x)), dtype=float)

    # Broadband incoherent sum
    for idx, f_idx in enumerate(valid_f_indices):
        f = freqs[f_idx]
        k = 2 * np.pi * f / v
        
        for y_idx, z in enumerate(grid_y):
            X_focal = grid_x
            
            d_tx = np.sqrt((X_focal[:, None] - tx_x_flat[None, :])**2 + (z - 0.45)**2)
            d_rx = np.sqrt((X_focal[:, None] - rx_x_flat[None, :])**2 + (z - 0.45)**2)
            
            # g size: (Nx, 440). We use Green's function for two-way propagation
            g = np.exp(-1j * k * (d_tx + d_rx)) / (np.sqrt(d_tx * d_rx) + 1e-15)
            
            # Normalize steering vector
            g_norm = np.linalg.norm(g, axis=1, keepdims=True)
            g = g / (g_norm + 1e-15)
            
            # Project onto noise subspace
            # g is (Nx, 440), Un is (440, 440-Ns)
            proj = np.sum(np.abs(np.conj(g) @ Un)**2, axis=1)
            
            I_map[y_idx, :] += 1.0 / (proj + 1e-15)
            
    image = I_map / len(valid_f_indices)

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