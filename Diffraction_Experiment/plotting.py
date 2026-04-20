import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import pyvista as pv
import h5py

from processing import apply_svd_filter, calculate_cross_correlation, envelope, kirchhoff_migration, compute_cwt_image

def get_snapshots_data(snapshot_folder, snapshot_prefix, snapshot_indices):
    loaded_snapshots = []
    for snap_num in snapshot_indices:
        snapshot_path = os.path.join(snapshot_folder, f'{snapshot_prefix}{snap_num}.vti')
        if not os.path.exists(snapshot_path):
            print(f"Warning: {snapshot_path} not found. Skipping plot.")
            return []

        mesh = pv.read(snapshot_path)
        dims = mesh.dimensions

        field_name = 'E-field'
        data = mesh[field_name]
        if data.shape[1] == 3:
            Ez_data = data[:, 2]
        else:
            Ez_data = data.flatten()

        data_2d = Ez_data.reshape(dims[1] - 1, dims[0] - 1)
        loaded_snapshots.append((snap_num, data_2d))
    return loaded_snapshots

def do_plot(loaded_snapshots, domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time, save_filename, title, n_blocks=None, block_width=None):
    if not loaded_snapshots:
        return

    length = domain_width
    depth  = domain_height
    upper_layer_thickness = air_thickness
    fracture_y0 = domain_height - fracture_bottom
    fracture_y1 = domain_height - fracture_top

    early_data = [d for s, d in loaded_snapshots if s <= 16]
    mid_data = [d for s, d in loaded_snapshots if 17 <= s <= 24]
    late_data = [d for s, d in loaded_snapshots if s > 24]

    early_abs_max = max(np.abs(d).max() for d in early_data) if early_data else 1.0
    mid_pct_ref = max(np.percentile(np.abs(d), 99.2) for d in mid_data) if mid_data else early_abs_max
    late_pct_ref = max(np.percentile(np.abs(d), 99.7) for d in late_data) if late_data else early_abs_max

    mid_abs_max = min(mid_pct_ref, early_abs_max * 0.22)
    mid_abs_max = max(mid_abs_max, early_abs_max * 0.035)

    late_abs_max = min(late_pct_ref, early_abs_max * 0.32)
    late_abs_max = max(late_abs_max, early_abs_max * 0.06)

    n_plots = len(loaded_snapshots)
    ncols = 2
    nrows = int(np.ceil(n_plots / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 3.3 * nrows), sharex=True, sharey=True)
    axes = np.array(axes).reshape(-1)

    for i, (snap_num, data_2d) in enumerate(loaded_snapshots):
        ax = axes[i]

        if snap_num <= 8:
            vlim = early_abs_max
            scale_tag = 'regular scale'
        elif snap_num <= 24:
            vlim = mid_abs_max
            scale_tag = 'reflection-enhanced'
        else:
            vlim = late_abs_max
            scale_tag = 'reflection-enhanced'

        im = ax.imshow(
            data_2d,
            aspect='auto',
            cmap='seismic',
            origin='lower',
            vmin=-vlim,
            vmax=vlim,
            extent=[0, length, 0, depth],
        )

        time_ns = snap_num * snapshot_time * 1e9
        ax.set_title(f'Snapshot {snap_num} (t = {time_ns:.2f} ns) | {scale_tag}', fontsize=10, weight='bold')

        ax.add_patch(Rectangle((0, depth - upper_layer_thickness), length, upper_layer_thickness, facecolor='#cfe8ff', edgecolor='blue', alpha=0.2))
        ax.add_patch(Rectangle((0, fracture_y0), length, fracture_y1 - fracture_y0, facecolor='#ffe6e6', edgecolor='red', linewidth=1.0, alpha=0.2, label='Fracture' if i == 0 else ''))

        if n_blocks is not None and block_width is not None:
            # Alternating Fracture Plot ('_alt' suffix -> diff or alt images)
            is_alternating = ('_alt' in save_filename) or ('_diff' in save_filename)
        
            for b in range(n_blocks):
                x_pos = b * block_width
                block_w = min(block_width, domain_width - x_pos)
                
                # Determine block material/color
                if is_alternating:
                    # Plus (air) is lightblue, Minus (water/high perm) is blue
                    block_color = '#ADD8E6' if b % 2 == 0 else '#00008B'
                    alpha_val = 0.15 if b % 2 == 0 else 0.25
                else:
                    # Homogenous is somewhere in between
                    block_color = '#4169E1'
                    alpha_val = 0.2
                    
                ax.add_patch(Rectangle((x_pos, fracture_y0), block_w, fracture_y1 - fracture_y0, 
                                       facecolor=block_color, edgecolor='none', alpha=alpha_val))
                                       
            # Illustrate tx/rx positions
            y_air_bottom = domain_height - air_thickness
            tx_x = domain_width / 2
            
            # Draw receivers
            for b in range(n_blocks):
                rx_x = (b + 0.5) * block_width
                if rx_x > 0 and rx_x < domain_width:
                    ax.plot(rx_x, y_air_bottom, 'g^', markersize=3, alpha=0.6, label='Receiver' if (b == 0 and i == 0) else '')
                    
            # Draw transmitter
            ax.plot(tx_x, y_air_bottom, 'r*', markersize=6, alpha=0.9, label='Transmitter' if i == 0 else '')
                                       
            # Just draw the separating lines
            for b in range(1, n_blocks):
                x_pos = b * block_width
                ax.vlines(x=x_pos, ymin=fracture_y0, ymax=fracture_y1, color='red', linestyle='--', linewidth=0.5, alpha=0.5)

        ax.set_xlabel('Length (m)', fontsize=10)
        ax.set_ylabel('Depth (m)', fontsize=10)

        if i == 0:
            ax.legend(loc='upper right', fontsize=9)

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cbar.set_label('Ez (V/m)', fontsize=8)
        cbar.ax.tick_params(labelsize=8)

    for j in range(n_plots, len(axes)):
        axes[j].axis('off')

    fig.suptitle(title, fontsize=14, weight='bold', y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig(save_filename, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Plot saved to: {save_filename}")

def plot_time_traces(out_alt, out_homo, n_blocks, block_width, rx_offset, save_filename, title, rx_per_block=1):
    fig, ax = plt.subplots(figsize=(10, 8))
    try:
        with h5py.File(out_alt, 'r') as fa, h5py.File(out_homo, 'r') as fh:
            iterations = fa.attrs['Iterations']
            dt = fa.attrs['dt']
            time = np.arange(iterations) * dt * 1e9 # ns
            
            diff_traces = []
            homo_traces = []
            rx_x_positions = []
            
            rx_idx = 1
            rx_dx = block_width / rx_per_block
            for i in range(n_blocks * rx_per_block):
                rx_x_init = (i + 0.5) * rx_dx
                dx_dy_dz = 0.002 # safe minimum
                domain_width = n_blocks * block_width
                if rx_x_init > dx_dy_dz and rx_x_init < domain_width - dx_dy_dz:
                    rx_name = f'rx{rx_idx}'
                    try:
                        ez_alt = np.array(fa['rxs'][rx_name]['Ez'])
                        ez_homo = np.array(fh['rxs'][rx_name]['Ez'])
                        homo_traces.append(ez_homo)
                        diff_traces.append(ez_alt - ez_homo)
                        rx_x_positions.append(rx_x_init)
                    except KeyError:
                        pass
                    rx_idx += 1
            
            if not diff_traces:
                return

            diff_traces = np.array(diff_traces)
            homo_traces = np.array(homo_traces)
            
            # --- Standard Difference Plot ---
            plot_traces_kernel(ax, diff_traces, time, n_blocks, block_width, title, rx_x_positions, rx_per_block=rx_per_block)
            plt.tight_layout()
            plt.savefig(save_filename, bbox_inches='tight', dpi=300)
            plt.close(fig)
            print(f"Plot saved to: {save_filename}")

            # --- SVD Filtering ---
            # Remove 1 (first) singular component to suppress strong reflectors/events
            n_comp = 0
            # svd_filtered_traces = apply_svd_filter(diff_traces, n_components_to_mute=n_comp)
            
            # fig_svd, ax_svd = plt.subplots(figsize=(10, 8))
            # plot_traces_kernel(ax_svd, svd_filtered_traces, time, n_blocks, block_width, f"{title} (SVD Muted: {n_comp} comp)")
            # svd_save_filename = save_filename.replace('.png', '_svd.png')
            # plt.tight_layout()
            # plt.savefig(svd_save_filename, bbox_inches='tight', dpi=300)
            # plt.close(fig_svd)
            # print(f"Plot saved to: {svd_save_filename}")

            # # --- Cross Correlation ---
            # cc_traces = calculate_cross_correlation(homo_traces, diff_traces)
            # fig_cc, ax_cc = plt.subplots(figsize=(10, 8))
            # # Just reusing plot_traces_kernel. Note that CC output is dimensionless or energy squared
            # # the time axis might represent lag instead depending on the goal, but 'same' mode centers it
            # plot_traces_kernel(ax_cc, cc_traces, time, n_blocks, block_width, f"Cross-Correlation (Avg Homo vs Diff)")
            # cc_save_filename = save_filename.replace('.png', '_cc.png')
            # plt.tight_layout()
            # plt.savefig(cc_save_filename, bbox_inches='tight', dpi=300)
            # plt.close(fig_cc)
            # print(f"Plot saved to: {cc_save_filename}")

            # # --- Depth Migration & Envelope ---
            # # Derive velocities and tx position for migration
            # v_ice = 3e8 / np.sqrt(6)
            # f_central = 1.5e9
            # wavelength_fracture = (3e8 / np.sqrt(6)) / f_central
            # rx_x_init_array = np.array([(i + 0.5) * block_width for i in range(n_blocks)])
            # source_receiver_steps = wavelength_fracture / 10
            # x_first_measurement = 1/2 * wavelength_fracture
            # tx_start_x = x_first_measurement + 54 * source_receiver_steps
            # max_depth = 0.5 # We know fracture is around 0.3-0.4m depth
            
            # # Step 1: Migrate the SVD filtered traces using Ice velocity
            # migrated_img, depths = kirchhoff_migration(
            #     traces=svd_filtered_traces, time_array=time, rx_x_array=rx_x_init_array, 
            #     tx_x=tx_start_x, velocity=v_ice, max_depth=max_depth, dz=0.002
            # )
            
            # # Step 2: Apply Hilbert envelope to get the Energy Blobs from wavelets
            # migrated_envelope = envelope(migrated_img)
            
            # # Plot the migrated image
            # fig_mig, ax_mig = plt.subplots(figsize=(10, 6))
            # im_m = ax_mig.imshow(migrated_envelope, aspect='auto', cmap='hot', 
            #                      extent=[0, n_blocks * block_width, np.max(depths), np.min(depths)])
            # ax_mig.set_title("Migrated Block Model Structure (SVD + Envelope)")
            # ax_mig.set_xlabel('Length (m)', fontsize=12)
            # ax_mig.set_ylabel('Depth (m)', fontsize=12)
            # cbar_m = fig_mig.colorbar(im_m, ax=ax_mig, pad=0.02)
            # cbar_m.set_label('Energy', fontsize=10)
            # print(f"Maximum migrated envelope value: {migrated_envelope.max()}")
            # mig_save_filename = save_filename.replace('.png', '_migrated.png')
            # plt.tight_layout()
            # plt.savefig(mig_save_filename, bbox_inches='tight', dpi=300)
            # plt.close(fig_mig)
            # print(f"Plot saved to: {mig_save_filename}")

            # --- Wavelet Transform (Time-Frequency Image) ---
            freqs, cwt_image = compute_cwt_image(diff_traces, dt, wavelet='cmor1.5-1.0')
            if cwt_image is not None:
                fig_cwt, ax_cwt = plt.subplots(figsize=(10, 6))
                
                # We plot Time (x-axis) vs Frequency (y-axis)
                # Note: cwt_image has shape (len(scales), len(time))
                # freqs are ordered from highest to lowest or lowest to highest depending on the creation
                # The extent must align properly. Time -> x from 0 to time[-1].
                # CWT is generally evaluated bottom up or top down, so we specify origin='lower' assuming
                # frequencies increase row-wise since we passed them generated from linspace (lowest to highest)
                
                # Freq to GHz for cleaner display
                freqs_ghz = freqs / 1e9
                
                im_cwt = ax_cwt.imshow(
                    cwt_image, aspect='auto', cmap='jet', origin='lower',
                    extent=[time[0], time[-1], freqs_ghz[0], freqs_ghz[-1]]
                )
                
                ax_cwt.set_title("Wavelet Transform of Global Average Difference Trace")
                ax_cwt.set_xlabel("Time (ns)", fontsize=12)
                ax_cwt.set_ylabel("Frequency (GHz)", fontsize=12)
                
                import matplotlib.ticker as ticker
                ax_cwt.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
                ax_cwt.tick_params(axis='x', which='minor', length=4, color='k')
                
                cbar_cwt = fig_cwt.colorbar(im_cwt, ax=ax_cwt, pad=0.02)
                cbar_cwt.set_label('CWT Magnitude', fontsize=10)
                
                cwt_save_filename = save_filename.replace('.png', '_cwt.png')
                plt.tight_layout()
                plt.savefig(cwt_save_filename, bbox_inches='tight', dpi=300)
                plt.close(fig_cwt)
                print(f"Plot saved to: {cwt_save_filename}")

    except Exception as e:
        print(f"Error plotting traces: {e}")

def plot_traces_kernel(ax, traces, time, n_blocks, block_width, title, rx_x_positions=None, rx_per_block=1):
    """Helper method to plot wiggle traces on a given axis."""
    max_val = np.max(np.abs(traces))
    if max_val == 0: max_val = 1
    rx_dx = block_width / rx_per_block
    scale = (rx_dx * 0.8) / max_val  # Amplify amplitudes
    
    for i, trace in enumerate(traces):
        if rx_x_positions is not None and i < len(rx_x_positions):
            x_base = rx_x_positions[i]
        else:
            x_base = (i + 0.5) * rx_dx
            
        scaled_trace = x_base + trace * scale
        
        ax.plot(scaled_trace, time, 'k-', linewidth=0.8)
        ax.fill_betweenx(time, x_base, scaled_trace, where=(scaled_trace > x_base), facecolor='k', alpha=0.5)
        ax.axvline(x_base, color='k', linestyle=':', linewidth=0.5, alpha=0.3)
        
    ax.set_ylim(time[-1], time[0])
    ax.set_xlim(0, n_blocks * block_width)
    ax.set_xlabel('Length (m)', fontsize=12)
    ax.set_ylabel('Time (ns)', fontsize=12)
    ax.set_title(title, fontsize=14, weight='bold')

    import matplotlib.ticker as ticker
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
    ax.tick_params(axis='y', which='minor', length=4, color='k')

def plot_snapshots(domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time, dx_dy_dz, n_blocks, block_width, rx_offset, rx_per_block=1):
    print("Plotting snapshots...")
    snapshot_prefix = 'snapshot_mid_x_'
    snapshot_indices = list(range(1, 37))

    # 1. Alternating snapshots
    folder_alt = r'horizontal_scattering_0p5lambda_snaps'
    alt_data = get_snapshots_data(folder_alt, snapshot_prefix, snapshot_indices)
    do_plot(
        alt_data, domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time,
        'gpr_snapshots_result_alt.png', 'GPR Forward Modeling Snapshots (Alternating Block Fracture)',
        n_blocks=n_blocks, block_width=block_width
    )

    # 2. Homogeneous snapshots
    folder_homo = r'horizontal_scattering_homogeneous_snaps'
    homo_data = get_snapshots_data(folder_homo, snapshot_prefix, snapshot_indices)
    do_plot(
        homo_data, domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time,
        'gpr_snapshots_result_homo.png', 'GPR Forward Modeling Snapshots (Homogeneous Fracture)',
        n_blocks=n_blocks, block_width=block_width
    )

    # 3. Difference snapshots (Alternating - Homogeneous)
    if alt_data and homo_data:
        diff_data = []
        for (snap_num, data_alt), (_, data_homo) in zip(alt_data, homo_data):
            diff_data.append((snap_num, data_alt - data_homo))
        do_plot(
            diff_data, domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time,
            'gpr_snapshots_result_diff.png', 'GPR Forward Modeling Snapshots (Diff: Alternating - Homogeneous)',
            n_blocks=n_blocks, block_width=block_width
        )
        
    # 4. Difference Time Traces (Alternating - Homogeneous) from .out files
    plot_time_traces(
        'horizontal_scattering_0p5lambda.out',
        'horizontal_scattering_homogeneous.out',
        n_blocks, block_width, rx_offset,
        'gpr_snapshots_result_diff_traces.png', 'GPR Difference Time Traces (Alternating - Homogeneous)',
        rx_per_block=rx_per_block
    )
