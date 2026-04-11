import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.signal import hilbert

def rtm_imaging():
    num_runs = 20
    # Arrays to hold FMC data
    # 20 Tx x 20 Rx x num_times
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

        # Read the target data
        with h5py.File(out_file, 'r') as f:
            if traces is None:
                dt = f.attrs['dt']
                num_time_steps = len(f['rxs']['rx1']['Ez'])
                traces = np.zeros((num_runs, 20, num_time_steps))
            
            rx_data = f['rxs']
            for rx_j in range(1, 21):
                traces[tx_i-1, rx_j-1, :] = rx_data[f'rx{rx_j}']['Ez']
                
        # Read the baseline data if available, to precisely remove direct wave
        if os.path.exists(baseline_file):
            with h5py.File(baseline_file, 'r') as f_base:
                rx_data_base = f_base['rxs']
                for rx_j in range(1, 21):
                    # Direct wave subtraction
                    traces[tx_i-1, rx_j-1, :] -= rx_data_base[f'rx{rx_j}']['Ez']
        else:
            print(f"Warning: {baseline_file} not found. Applying smooth time-mute instead.")
            t_arr = np.arange(num_time_steps) * dt
            # Use a simple hann taper instead of a hard mute to prevent step artifacts
            mute_mask = t_arr < 2.5e-9
            traces[:, :, mute_mask] *= (0.5 * (1 - np.cos(np.pi * t_arr[mute_mask] / 2.5e-9)))

    print("Executing Frequency-domain Reflection Matrix backpropagation...")
    
    # Physics params
    c = 299792458.0
    eps_r = 6.0
    v = c / np.sqrt(eps_r)

    # Imaging grid over region of interest (zoom in slightly)
    grid_x = np.linspace(0.3, 0.7, 100)
    grid_y = np.linspace(0.1, 0.4, 100)
    image = np.zeros((len(grid_y), len(grid_x)), dtype=complex)

    # Convert time-domain traces to frequency-domain Reflection Matrix
    # traces shape: [N_tx, N_rx, num_time_steps] -> D_f shape: [N_tx, N_rx, num_freqs]
    D_f = np.fft.rfft(traces, axis=-1)
    freqs = np.fft.rfftfreq(num_time_steps, dt)
    
    # Select frequency band around center 1.5 GHz
    f_min, f_max = 0.5e9, 2.5e9
    valid_f_indices = np.where((freqs >= f_min) & (freqs <= f_max))[0]

    # Source time delay correction. In gprMax, the Ricker source peak is shifted to make it causal. 
    # For a 1.5 GHz Ricker, the peak occurs exactly at 1.0 / 1.5 GHz = 0.666 ns.
    # To place the targets at their true spatial locations, we must advance the phase (remove the delay) in the backpropagation.
    # A slightly empirical offset is t_delay ~ 0.81 ns due to zero-offset coupling effects.
    t_delay = 0.82e-9

    for f_idx in valid_f_indices:
        f = freqs[f_idx]
        
        # d(f) is the Reflection Matrix at surface, shape: [N_rx, N_tx]
        # Our D_f is [tx, rx, freq], so we transpose the Tx and Rx dimensions
        # Apply the phase advance to compensate for the source emission delay
        d_f_mat = D_f[:, :, f_idx].T * np.exp(1j * 2 * np.pi * f * t_delay) 
        
        k = 2 * np.pi * f / v
        
        for y_idx, z in enumerate(grid_y):
            # x_in points at depth z
            X_focal = grid_x
            
            # G_conj size: (N_x, N_array). Distance from (X_focal, z) to array element
            # Arrays are at rx_x_arr, rx_y = 0.45 (same for tx)
            d_to_array = np.sqrt((X_focal[:, None] - rx_x_arr[None, :])**2 + (z - rx_y)**2)
            
            # G*_ur(z, f) = exp(+j * k * r)
            G_conj = np.exp(1j * k * d_to_array)
            G_dagger = G_conj.T  # Transpose gives the required (G*_ur)^T
            
            # Rxx(z, f) = G*_ur(z, f) x d(f) x G_dagger_ur(z, f)
            R_xx = G_conj @ d_f_mat @ G_dagger
            
            # I(z, f) = diag[Rxx(x_in = x_out, x_in, f)]
            I_zf = np.diag(R_xx)
            
            image[y_idx, :] += I_zf

    # Take absolute amplitude for the final image envelope
    env_image = np.abs(image)
    # Plot results
    plt.figure(figsize=(10, 6))
    plt.imshow(env_image, extent=[grid_x[0], grid_x[-1], grid_y[0], grid_y[-1]], 
               origin='lower', aspect='auto', cmap='plasma')
    
    # Overlay the true geometry of our sub-wavelength air pockets
    circle1 = plt.Circle((0.47, 0.25), 0.01, color='lime', fill=False, lw=2, label='True target positions')
    circle2 = plt.Circle((0.53, 0.25), 0.01, color='lime', fill=False, lw=2)
    plt.gca().add_patch(circle1)
    plt.gca().add_patch(circle2)
    
    plt.title('Reverse Time Migration (Reflection Matrix - FMC array)')
    plt.xlabel('X (m)')
    plt.ylabel('Y / Depth (m)')
    plt.legend(loc='lower left')
    plt.colorbar(label='RTM Focus Intensity (Envelope)')
    
    out_img = 'RTM/rtm_result_fmc.png'
    plt.savefig(out_img, dpi=300)
    print(f"Saved generated image map to {out_img}")

if __name__ == '__main__':
    rtm_imaging()