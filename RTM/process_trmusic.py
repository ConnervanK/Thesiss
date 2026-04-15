import h5py
import numpy as np
import matplotlib.pyplot as plt
import os


def trmusic_imaging():
    num_runs = 20
    num_rx = 20
    num_tx = 20
    traces = None

    rx_x_arr = np.linspace(0.1, 0.9, num_rx)
    tx_x_arr = np.linspace(0.1, 0.9, num_tx)
    rx_y = 0.45
    tx_y = 0.45

    f_min, f_max = 0.0e9, 10e9
    mute_time_s = 2.5e-9
    use_delay_comp = False
    t_delay = 0.82e-9
    num_scatterers = 5
    eps = 1e-15

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
                traces = np.zeros((num_runs, num_rx, num_time_steps))
                t_arr = np.arange(num_time_steps) * dt
                mute_mask = t_arr < mute_time_s
                taper = 0.5 * (1 - np.cos(np.pi * t_arr[mute_mask] / mute_time_s))

            rx_data = f['rxs']
            for rx_j in range(1, num_rx + 1):
                traces[tx_i - 1, rx_j - 1, :] = rx_data[f'rx{rx_j}']['Ez']

        if os.path.exists(baseline_file):
            with h5py.File(baseline_file, 'r') as f_base:
                rx_data_base = f_base['rxs']
                for rx_j in range(1, num_rx + 1):
                    traces[tx_i - 1, rx_j - 1, :] -= rx_data_base[f'rx{rx_j}']['Ez']
        else:
            traces[tx_i - 1, :, mute_mask] *= taper

    print("Executing Time Reversal MUSIC (TR-MUSIC) imaging...")

    c = 299792458.0
    eps_r = 6.0
    v = c / np.sqrt(eps_r)

    grid_x = np.linspace(0.3, 0.7, 100)
    grid_y = np.linspace(0.1, 0.4, 100)
    image = np.zeros((len(grid_y), len(grid_x)), dtype=float)

    D_f = np.fft.rfft(traces, axis=-1)
    freqs = np.fft.rfftfreq(num_time_steps, dt)

    valid_f_indices = np.where((freqs >= f_min) & (freqs <= f_max))[0]
    if len(valid_f_indices) == 0:
        print("No valid frequency bins found in the requested band.")
        return

    signal_dim = min(num_scatterers, min(num_rx, num_tx) - 1)

    for f_idx in valid_f_indices:
        f = freqs[f_idx]

        # Multistatic matrix K(freq) with shape (Nrx, Ntx)
        K = D_f[:, :, f_idx].T
        if use_delay_comp:
            K = K * np.exp(1j * 2 * np.pi * f * t_delay)

        # TR-MUSIC uses both left and right singular vectors for receive/transmit spaces.
        U, _, Vh = np.linalg.svd(K, full_matrices=False)
        Un = U[:, signal_dim:]
        Vn = Vh.conj().T[:, signal_dim:]

        k = 2 * np.pi * f / v
        I_f = np.zeros((len(grid_y), len(grid_x)), dtype=float)

        for y_idx, z in enumerate(grid_y):
            d_rx = np.sqrt((grid_x[:, None] - rx_x_arr[None, :]) ** 2 + (z - rx_y) ** 2)
            d_tx = np.sqrt((grid_x[:, None] - tx_x_arr[None, :]) ** 2 + (z - tx_y) ** 2)

            g_rx = np.exp(-1j * k * d_rx) / np.sqrt(np.maximum(d_rx, eps))
            g_tx = np.exp(-1j * k * d_tx) / np.sqrt(np.maximum(d_tx, eps))

            g_rx = g_rx / (np.linalg.norm(g_rx, axis=1, keepdims=True) + eps)
            g_tx = g_tx / (np.linalg.norm(g_tx, axis=1, keepdims=True) + eps)

            proj_rx = np.sum(np.abs(np.conj(g_rx) @ Un) ** 2, axis=1)
            proj_tx = np.sum(np.abs(np.conj(g_tx) @ Vn) ** 2, axis=1)

            I_f[y_idx, :] = 1.0 / (proj_rx * proj_tx + eps)

        image += I_f / len(valid_f_indices)

    image = 10 * np.log10(image / (np.max(image) + eps) + eps)

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
