import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import glob
import numpy as np
from scipy.signal import csd as scipy_csd

def plot_wiggle_traces(model, traces, title, save_filename, domain_width=None):
    """Plots standard wiggle traces from the model data."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    max_val = np.max(np.abs(traces))
    if max_val == 0: max_val = 1
    
    rx_dx = model.rx_spacing if hasattr(model, 'rx_spacing') and model.rx_spacing > 0 else (model.block_width / model.rx_per_block)
    scale = (rx_dx * 0.8) / max_val
    rx_x_arr = model.get_rx_x_array()
    
    for i, trace in enumerate(traces):
        if hasattr(model, 'mode') and model.mode == 'bscan' and hasattr(model, 'actual_traces') and model.actual_traces == 1:
            x_base = rx_x_arr[0] if len(rx_x_arr) > 0 else 0
        else:
            x_base = rx_x_arr[i] if i < len(rx_x_arr) else (i + 0.5) * rx_dx

        scaled_trace = x_base + trace * scale
        
        ax.plot(scaled_trace, model.time, 'k-', linewidth=0.8)
        ax.fill_betweenx(model.time, x_base, scaled_trace, where=(scaled_trace > x_base), facecolor='k', alpha=0.5)
        ax.axvline(x_base, color='k', linestyle=':', linewidth=0.5, alpha=0.3)
        
    ax.set_ylim(model.time[-1], model.time[0])
    
    if domain_width:
        ax.set_xlim(0, domain_width)
    else:
        min_x = min(rx_x_arr[0], model.tx_start_x) if hasattr(model, 'tx_start_x') else rx_x_arr[0]
        max_x = max(rx_x_arr[-1], model.tx_start_x) if hasattr(model, 'tx_start_x') else rx_x_arr[-1]
        padding = (max_x - min_x) * 0.1 if max_x > min_x else rx_dx * 2
        ax.set_xlim(min_x - padding, max_x + padding)

    ax.set_xlabel('Length (m)', fontsize=12)
    ax.set_ylabel('Time (ns)', fontsize=12)
    ax.set_title(title, fontsize=14, weight='bold')

    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
    ax.tick_params(axis='y', which='minor', length=4, color='k')
    
    plt.tight_layout()
    plt.savefig(save_filename, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Plot saved to: {save_filename}")

def plot_fk_image(k, freqs, fk_mag, title, save_filename, block_width=None):
    """Plots Frequency-Wavenumber (F-K) spectra."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Scale to GHz for easier reading
    freqs_ghz = freqs / 1e9
    
    # Use log scale to visualize both bright and dim spectral lobes
    fk_log = np.log1p(fk_mag / np.max(fk_mag) * 1e4)

    im = ax.imshow(
        fk_log.T, aspect='auto', cmap='jet', origin='lower',
        extent=[k[0], k[-1], freqs_ghz[0], freqs_ghz[-1]]
    )

    # Add vertical lines for +/- 2*pi / d where d is the block_width
    if block_width is not None and block_width > 0:
        k_lateral = 2 * np.pi / block_width
        ax.axvline(x=k_lateral, color='w', linestyle='dashed', linewidth=1.5, alpha=0.8, label=r'$k = \pm 2\pi/d$')
        ax.axvline(x=-k_lateral, color='w', linestyle='dashed', linewidth=1.5, alpha=0.8)
        ax.legend(loc='upper right', fontsize=12)
    
    ax.set_title(title)
    ax.set_xlabel("Wavenumber $k_x$ [rad/m]", fontsize=12)
    ax.set_ylabel("Frequency [GHz]", fontsize=12)
    
    # Ensure standard quadrant display (positive freq is usually enough, but full bounds are safe)
    ax.set_ylim(0, 5) # focus up to 5 GHz
    
    # The max measured wavenumber is pi/d (Nyquist limit). Since we are drawing 2*pi/d, 
    # we need to artificially widen the x-axis to show the lines relative to the spectrum.
    if block_width is not None and block_width > 0:
        k_lateral = 2 * np.pi / block_width
        ax.set_xlim(-k_lateral * 1.1, k_lateral * 1.1)

    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label('Log Magnitude', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(save_filename, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Plot saved to: {save_filename}")

def plot_cwt_image(time, freqs, cwt_mag, title, save_filename, vmin=None, vmax=None, cmap='jet'):
    """Plots Continuous Wavelet Transform (Time-Frequency)."""
    fig, ax = plt.subplots(figsize=(10, 6))
    freqs_ghz = freqs / 1e9
    
    im = ax.imshow(
        cwt_mag, aspect='auto', cmap=cmap, origin='lower',
        extent=[time[0], time[-1], freqs_ghz[0], freqs_ghz[-1]],
        vmin=vmin, vmax=vmax
    )
    
    ax.set_title(title)
    ax.set_xlabel("Time (ns)", fontsize=12)
    ax.set_ylabel("Frequency (GHz)", fontsize=12)

    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
    ax.tick_params(axis='x', which='minor', length=4, color='k')

    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label('CWT Value', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(save_filename, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Plot saved to: {save_filename}")

def plot_cwt_cross_sections(time, freqs, cwt_data, f1, f2, title, save_filename):
    """Plots a 1D cross section of the CWT phase/magnitude at two distinct frequencies."""
    f1_idx = np.argmin(np.abs(freqs - f1))
    f2_idx = np.argmin(np.abs(freqs - f2))
    
    fig, ax = plt.subplots(figsize=(10, 4))
    
    ax.plot(time, cwt_data[f1_idx, :], 'b-', label=f'{freqs[f1_idx]/1e9:.2f} GHz')
    ax.plot(time, cwt_data[f2_idx, :], 'r-', label=f'{freqs[f2_idx]/1e9:.2f} GHz')
    
    ax.set_title(title)
    ax.set_xlabel("Time (ns)", fontsize=12)
    ax.set_ylabel("Phase (rad)", fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend(loc='best')
    
    plt.tight_layout()
    plt.savefig(save_filename, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Plot saved to: {save_filename}")


def plot_aligned_traces(time, trace1, trace2, title, save_filename):
    """
    Plots two single traces (one shifted or windowed) overlapping to visualize alignment.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(time * 1e9, trace1, label='Trace 1 (Windowed)', linewidth=1.5, alpha=0.9)
    ax.plot(time * 1e9, trace2, label='Trace 2 (Shifted & Windowed)', linewidth=1.5, alpha=0.9, linestyle='--')
    
    ax.set_title(title, fontsize=14)
    ax.set_xlabel("Time (ns)", fontsize=12)
    ax.set_ylabel("Amplitude", fontsize=12)
    ax.legend(loc='best')
    ax.grid(True, linestyle='--', alpha=0.6)
    
    plt.tight_layout()
    plt.savefig(save_filename, dpi=300)
    plt.close(fig)
    print(f"Plot saved to: {save_filename}")

def plot_xwt_phase_arrows(time, freqs, xwt_power, xwt_phase, title, save_filename, vmin=None, vmax=None, cmap='jet'):
    """
    Plots the Cross-Wavelet Transform power with phase arrows.
    Right-pointing (→): In-phase (0)
    Left-pointing (←): Out-of-phase (180 / π)
    Down-pointing (↓): T1 leads T2 by 90° (π/2)
    Up-pointing (↑): T2 leads T1 by 90° (-π/2)
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    freqs_ghz = freqs / 1e9
    
    im = ax.imshow(
        xwt_power, aspect='auto', cmap=cmap, origin='lower',
        extent=[time[0], time[-1], freqs_ghz[0], freqs_ghz[-1]],
        vmin=vmin, vmax=vmax
    )
    
    # Overlay phase arrows. Subsample to avoid clutter.
    dt = max(1, len(time) // 30)
    df = max(1, len(freqs) // 15)

    X, Y = np.meshgrid(time, freqs_ghz)

    # Subsampled phase matrix
    Ph = xwt_phase[::df, ::dt]
    # Unit direction vectors from phase
    U = np.cos(Ph)
    V = -np.sin(Ph)

    # Compensate for axis aspect: frequency axis is in GHz, time in ns
    time_span = time[-1] - time[0]
    freq_span = freqs_ghz[-1] - freqs_ghz[0]
    span_ratio = freq_span / time_span if time_span > 0 else 1.0

    # Scale V to match aspect so arrows visually represent phase direction
    V_scaled = V * span_ratio

    # Normalize vectors to unit length, then set a fixed arrow length in data units
    norm = np.sqrt(U**2 + V_scaled**2)
    norm[norm == 0] = 1.0
    U_unit = U / norm
    V_unit = V_scaled / norm

    # Arrow length chosen as a fraction of the smaller plot span to keep arrows readable
    arrow_len = max(0.02, 0.06 * min(time_span, freq_span))

    U_plot = U_unit * arrow_len
    V_plot = V_unit * arrow_len

    # Use scale_units='xy' and scale=1 to interpret U_plot/V_plot in data coordinates
    ax.quiver(
        X[::df, ::dt], Y[::df, ::dt], U_plot, V_plot,
        angles='xy', pivot='mid', color='black', scale=1, scale_units='xy',
        width=0.003, headwidth=3, headlength=4, alpha=0.8
    )
    
    ax.set_title(title)
    ax.set_xlabel("Time (ns)", fontsize=12)
    ax.set_ylabel("Frequency (GHz)", fontsize=12)

    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.tick_params(axis='x', which='minor', length=4, color='k')

    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label('XWT Power', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(save_filename, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Plot saved to: {save_filename}")

def plot_migrated_image(migrated_img, z_array, rx_x_array, tx_x, title, save_filename, fracture_depth=None):
    fig, ax = plt.subplots(figsize=(10, 6))

    vlim = np.percentile(np.abs(migrated_img), 99)
    if vlim == 0:
        vlim = 1.0

    im = ax.imshow(
        migrated_img, aspect='auto', cmap='seismic',
        extent=[rx_x_array[0], rx_x_array[-1], z_array[-1], z_array[0]],
        vmin=-vlim, vmax=vlim,
        origin='upper'
    )

    if fracture_depth is not None:
        ax.axhline(fracture_depth, color='lime', linestyle='--', linewidth=1.2,
                   label=f'Fracture depth ({fracture_depth:.2f} m)')

    ax.plot(tx_x, z_array[0], 'r*', markersize=10, zorder=5, label='Tx')
    ax.legend(loc='upper right', fontsize=9)

    ax.set_title(title, fontsize=14, weight='bold')
    ax.set_xlabel('x (m)', fontsize=12)
    ax.set_ylabel('Depth (m)', fontsize=12)

    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label('Amplitude', fontsize=10)

    plt.tight_layout()
    plt.savefig(save_filename, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Plot saved to: {save_filename}")


def plot_bscan_section(traces, time, x_array, title, save_filename,
                       fracture_depth=None, velocity=None, tx_rx_offset=0.0):
    fig, ax = plt.subplots(figsize=(12, 8))

    time_ns = time * 1e9
    vlim = np.percentile(np.abs(traces), 99)
    if vlim == 0:
        vlim = 1.0

    im = ax.imshow(
        traces.T,
        aspect='auto',
        cmap='seismic',
        origin='upper',
        vmin=-vlim,
        vmax=vlim,
        extent=[x_array[0], x_array[-1], time_ns[-1], time_ns[0]],
    )

    if fracture_depth is not None and velocity is not None:
        half = tx_rx_offset / 2.0
        t_frac_ns = np.sqrt(half**2 + fracture_depth**2) * 2 / velocity * 1e9
        ax.axhline(t_frac_ns, color='lime', linestyle='--', linewidth=1.2,
                   label=f'Expected fracture arrival ({t_frac_ns:.2f} ns)')
        ax.legend(loc='upper right', fontsize=9)

    ax.set_title(title, fontsize=14, weight='bold')
    ax.set_xlabel('Midpoint position (m)', fontsize=12)
    ax.set_ylabel('Two-way travel time (ns)', fontsize=12)

    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.25))
    ax.tick_params(axis='y', which='minor', length=4, color='k')

    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label('Amplitude (V/m)', fontsize=10)

    plt.tight_layout()
    plt.savefig(save_filename, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Plot saved to: {save_filename}")


import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import pyvista as pv
import h5py

from utils.processing import apply_svd_filter, calculate_cross_correlation, envelope, kirchhoff_migration, compute_cwt_image

def _resolve_snapshot_folder(snapshot_folder, snapshot_prefix, snapshot_indices):
    if os.path.isdir(snapshot_folder):
        candidate = os.path.join(snapshot_folder, f'{snapshot_prefix}{snapshot_indices[0]}.vti')
        if os.path.exists(candidate):
            return snapshot_folder

    parent_dir = os.path.dirname(snapshot_folder)
    base_name = os.path.basename(snapshot_folder)
    candidates = sorted(
        folder for folder in glob.glob(os.path.join(parent_dir, f"{base_name}*"))
        if os.path.isdir(folder)
    )

    for folder in candidates:
        candidate = os.path.join(folder, f'{snapshot_prefix}{snapshot_indices[0]}.vti')
        if os.path.exists(candidate):
            return folder

    # Fallback: climb to the project root and look in the sibling configs/ directory.
    current_dir = os.path.abspath(snapshot_folder)
    project_root = None
    for _ in range(6):
        current_dir = os.path.dirname(current_dir)
        if os.path.isdir(os.path.join(current_dir, 'configs')):
            project_root = current_dir
            break

    if project_root is not None:
        fallback_root = os.path.join(project_root, 'configs')
        fallback_candidates = sorted(
            folder for folder in glob.glob(os.path.join(fallback_root, f"{base_name}*"))
            if os.path.isdir(folder)
        )
        for folder in fallback_candidates:
            candidate = os.path.join(folder, f'{snapshot_prefix}{snapshot_indices[0]}.vti')
            if os.path.exists(candidate):
                return folder

    return snapshot_folder


def get_snapshots_data(snapshot_folder, snapshot_prefix, snapshot_indices):
    snapshot_folder = _resolve_snapshot_folder(snapshot_folder, snapshot_prefix, snapshot_indices)
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

def do_plot(loaded_snapshots, domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time, save_filename, title, n_blocks=None, block_width=None, rx_per_block=1, is_diff=False, tx_x=2.0, rx_start_x=1.0, rx_spacing=None, rx_count=None):
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

    if is_diff:
        # Relax constraints for diff data so it doesn't saturate
        mid_abs_max = mid_pct_ref
        late_abs_max = late_pct_ref
    else:
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
            
            # Draw receivers at actual positions (if provided)
            if rx_start_x is not None and rx_spacing is not None and rx_count is not None:
                for rx_idx in range(rx_count):
                    rx_x = rx_start_x + rx_idx * rx_spacing
                    if rx_x > 0 and rx_x < domain_width:
                        ax.plot(rx_x, y_air_bottom, 'g^', markersize=3, alpha=0.6, label='Receiver' if (rx_idx == 0 and i == 0) else '')
            else:
                # Fallback: use block-based distribution if receiver params not provided
                rx_dx = block_width / rx_per_block
                for rx_idx in range(n_blocks * rx_per_block):
                    rx_x = (rx_idx + 0.5) * rx_dx
                    if rx_x > 0 and rx_x < domain_width:
                        ax.plot(rx_x, y_air_bottom, 'g^', markersize=3, alpha=0.6, label='Receiver' if (rx_idx == 0 and i == 0) else '')
                    
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

def plot_snapshots(domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time, dx_dy_dz, n_blocks, block_width, rx_offset, rx_per_block=1, output_dir=None):
    if output_dir is None:
        output_dir = os.path.join('data', 'outputs')
    
    print("Plotting snapshots...")
    snapshot_prefix = 'snapshot_mid_x_'
    snapshot_indices = list(range(1, 37))

    # 1. Alternating snapshots
    folder_alt = os.path.join(output_dir, 'horizontal_scattering_0p5lambda_snaps')
    alt_data = get_snapshots_data(folder_alt, snapshot_prefix, snapshot_indices)
    do_plot(
        alt_data, domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time,
        os.path.join(output_dir, 'gpr_snapshots_result_alt.png'), 'GPR Forward Modeling Snapshots (Alternating Block Fracture)',
        n_blocks=n_blocks, block_width=block_width, rx_per_block=rx_per_block
    )

    # 2. Homogeneous snapshots
    folder_homo = os.path.join(output_dir, 'horizontal_scattering_homogeneous_snaps')
    homo_data = get_snapshots_data(folder_homo, snapshot_prefix, snapshot_indices)
    do_plot(
        homo_data, domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time,
        os.path.join(output_dir, 'gpr_snapshots_result_homo.png'), 'GPR Forward Modeling Snapshots (Homogeneous Background)',
        n_blocks=n_blocks, block_width=block_width, rx_per_block=rx_per_block
    )

    # 3. Difference snapshots (Alternating - Homogeneous)
    if alt_data and homo_data:
        diff_data = []
        for (snap_num, data_alt), (_, data_homo) in zip(alt_data, homo_data):
            diff_data.append((snap_num, data_alt - data_homo))
        do_plot(
            diff_data, domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time,
            os.path.join(output_dir, 'gpr_snapshots_result_diff.png'), 'GPR Forward Modeling Snapshots (Diff: Alternating - Homogeneous)',
            n_blocks=n_blocks, block_width=block_width, rx_per_block=rx_per_block, is_diff=True
        )
        
    # 4. Difference Time Traces (Alternating - Homogeneous) from .out files
    plot_time_traces(
        os.path.join(output_dir, 'horizontal_scattering_0p5lambda.out'),
        os.path.join(output_dir, 'horizontal_scattering_homogeneous.out'),
        n_blocks, block_width, rx_offset,
        os.path.join(output_dir, 'gpr_snapshots_result_diff_traces.png'), 'GPR Difference Time Traces (Alternating - Homogeneous)',
        rx_per_block=rx_per_block
    )

def plot_csd(trace1, trace2, fs_hz, title, save_filename, nperseg=256, source_freq_hz=None, x_limit_ghz=10.0):
    """Plots cross-spectral magnitude and phase difference for two traces."""
    freqs_hz, cross_psd = scipy_csd(
        trace1,
        trace2,
        fs=fs_hz,
        nperseg=nperseg,
        noverlap=nperseg // 2,
        scaling='density',
        return_onesided=True,
    )

    freqs_ghz = freqs_hz / 1e9
    magnitude_db = 10.0 * np.log10(np.maximum(np.abs(cross_psd), 1e-30))
    phase_rad = np.unwrap(np.angle(cross_psd))

    fig, (ax_mag, ax_phase) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    ax_mag.plot(freqs_ghz, magnitude_db, color='navy', linewidth=1.2)
    ax_mag.set_ylabel('Magnitude (dB/Hz)', fontsize=12)
    ax_mag.set_title(title, fontsize=14)
    ax_mag.grid(True, linestyle='--', alpha=0.6)
    if source_freq_hz is not None:
        ax_mag.axvline(source_freq_hz / 1e9, color='crimson', linestyle='--', linewidth=1.2, alpha=0.8)

    ax_phase.plot(freqs_ghz, phase_rad, color='darkgreen', linewidth=1.2)
    ax_phase.set_xlabel('Frequency (GHz)', fontsize=12)
    ax_phase.set_ylabel('Phase Difference (rad)', fontsize=12)
    ax_phase.grid(True, linestyle='--', alpha=0.6)
    if source_freq_hz is not None:
        ax_phase.axvline(source_freq_hz / 1e9, color='crimson', linestyle='--', linewidth=1.2, alpha=0.8)

    if x_limit_ghz is not None:
        upper_limit = min(x_limit_ghz, float(freqs_ghz[-1]))
        ax_phase.set_xlim(0, upper_limit)

        visible_mask = freqs_ghz <= upper_limit
        if np.any(visible_mask):
            visible_mag = magnitude_db[visible_mask]
            visible_phase = phase_rad[visible_mask]

            mag_max = float(np.max(visible_mag))
            mag_min = float(np.min(visible_mag))
            ax_mag.set_ylim(mag_min - 5.0, mag_max + 5.0)

            phase_min = float(np.min(visible_phase))
            phase_max = float(np.max(visible_phase))
            phase_pad = max(1.0, 0.05 * (phase_max - phase_min))
            ax_phase.set_ylim(phase_min - phase_pad, phase_max + phase_pad)

    plt.tight_layout()
    plt.savefig(save_filename, dpi=300)
    plt.close(fig)
