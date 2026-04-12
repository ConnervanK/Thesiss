import h5py
import numpy as np
import matplotlib.pyplot as plt
import os

def trmusic_imaging():
    num_runs = 20
    traces = None
    
    rx_x_arr = np.linspace(0.1, 0.9, 20)
    tx_x_arr = np.linspace(0.1, 0.9, 20)
    rx_y = 0.45
    tx_y = 0.45

    for tx_i in range(1, num_runs + 1):
        out_file = f'RTM/subwavelength_scatters{tx_i}.out'
        baseline_file = f'RTM/baseline{tx_i}.out'
        
        if not os.path.exists(out_file):
            print(f"File {out_file} not found. Run simulation first.")
            return

        with h5py.File(out_file, 'r') as f:
            if traces is None:
                dt = f.attrs['dt']
                num_time_steps = len(f['rxs']['rx1']['Ez'])
                traces = np.zeros((num_runs, 20, num_time_steps))
            
            rx_data = f['rxs']
            for rx_j in range(1, 21):
                traces[tx_i-1, rx_j-1, :] = rx_data[f'rx{rx_j}']['Ez']
                
        if os.path.exists(baseline_file):
            with h5py.File(baseline_file, 'r') as f_base:
                rx_data_base = f_base['rxs']
                for rx_j in range(1, 21):
                    traces[tx_i-1, rx_j-1, :] -= rx_data_base[f'rx{rx_j}']['Ez']
        else:
            t_arr = np.arange(num_time_steps) * dt
            mute_mask = t_arr < 2.5e-9
            for tx_i_idx in range(num_runs):
                for rx_j_idx in range(20):
                    traces[tx_i_idx, rx_j_idx, mute_mask] *= (0.5 * (1 - np.cos(np.pi * t_arr[mute_mask] / 2.5e-9)))

    print("Executing Time Reversal MUSIC (TR-MUSIC) imaging...")
    
    c = 299792458.0
    eps_r = 6.0
    v = c / np.sqrt(eps_r)

    grid_x = np.linspace(0.3, 0.7, 100)
    grid_y = np.linspace(0.1, 0.4, 100)
    image = np.zeros((len(grid_y), len(grid_x)), dtype=float)

    D_f = np.fft.rfft(traces, axis=-1)
    freqs = np.fft.rfftfreq(num_time_steps, dt)
    
    f_min, f_max = 0.5e9, 2.5e9
    valid_f_indices = np.where((freqs >= f_min) & (freqs <= f_max))[0]

    t_delay = 0.82e-9 
    
    # We expect 2 scatterers, so signal subspace dimension is 2
    num_scatterers = 6

    for f_idx in valid_f_indices:
        f = freqs[f_idx]
        
        # Multistatic Data Matrix K
        K = D_f[:, :, f_idx].T * np.exp(1j * 2 * np.pi * f * t_delay)
        
        # SVD of K
        U, S, Vh = np.linalg.svd(K, full_matrices=True)
        
        # Noise subspace U_n consists of singular vectors from index num_scatterers to end
        Un = U[:, num_scatterers:]
        
        k = 2 * np.pi * f / v # wave number
        
        # Precompute the pseudo-spectrum at this frequency
        I_f = np.zeros((len(grid_y), len(grid_x)), dtype=float)
        
        for y_idx, z in enumerate(grid_y):
            X_focal = grid_x
            
            d_to_array = np.sqrt((X_focal[:, None] - rx_x_arr[None, :])**2 + (z - rx_y)**2)
            
            # Steering vector g
            # For 2D, green function is proportional to exp(-jkd) / sqrt(d) (Hankel function approx)
            # In far-field or basic TR-MUSIC: g = exp(-jkd)
            g = np.exp(-1j * k * d_to_array) / np.sqrt(d_to_array)
            
            # Normalize steering vectors
            g_norm = np.linalg.norm(g, axis=1, keepdims=True)
            g = g / (g_norm + 1e-15)
            
            # Projection onto noise subspace
            # proj = sum_m |g^H * u_m|^2
            # g is (Nx, Nrx)
            # Un is (Nrx, Nrx - num_scatterers)
            proj = np.sum(np.abs(np.conj(g) @ Un)**2, axis=1)
            
            I_f[y_idx, :] = 1.0 / (proj + 1e-15)
            
        # Add to broadband image
        image += I_f / len(valid_f_indices)

    # Convert to log scale for better visualization as pseudo-spectrum peaks can be very sharp
    #image = 10 * np.log10(image / np.max(image))

    plt.figure(figsize=(10, 6))
    plt.imshow(image, extent=[grid_x[0], grid_x[-1], grid_y[0], grid_y[-1]], 
               origin='lower', aspect='auto', cmap='plasma')
    
    circle1 = plt.Circle((0.47, 0.25), 0.01, color='lime', fill=False, lw=2, label='True target positions')
    circle2 = plt.Circle((0.53, 0.25), 0.01, color='lime', fill=False, lw=2)
    plt.gca().add_patch(circle1)
    plt.gca().add_patch(circle2)
    
    plt.title('High-Resolution TR-MUSIC Pseudo-spectrum')
    plt.xlabel('X (m)')
    plt.ylabel('Y / Depth (m)')
    plt.legend(loc='lower left')
    plt.colorbar(label='MUSIC Profile (dB)')
    
    out_img = 'RTM/trmusic_result_fmc.png'
    plt.savefig(out_img, dpi=300)
    print(f"Saved generated image map to {out_img}")

if __name__ == '__main__':
    trmusic_imaging()
