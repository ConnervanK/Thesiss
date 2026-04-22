import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

def plot_wiggle_traces(model, traces, title, save_filename):
    """Plots standard wiggle traces from the model data."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    max_val = np.max(np.abs(traces))
    if max_val == 0: max_val = 1
    
    rx_dx = model.block_width / model.rx_per_block
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
    ax.set_xlim(0, model.n_blocks * model.block_width)
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
    
    Ph = xwt_phase[::df, ::dt]
    # Ph = angle(W1) - angle(W2). If T1 leads T2 by 90deg, phase is pi/2.
    # sin(pi/2) = 1, cos = 0. To make it point down, we use V = -sin.
    U = np.cos(Ph)
    V = -np.sin(Ph)
    
    # Adjust scale to keep arrows proportional
    # Find a good scale heuristic based on plot bounds
    time_span = time[-1] - time[0]
    freq_span = freqs_ghz[-1] - freqs_ghz[0]
    span_ratio = freq_span / time_span if time_span > 0 else 1.0

    ax.quiver(X[::df, ::dt], Y[::df, ::dt], U, V * span_ratio, 
              angles='xy', pivot='mid', 
              color='black', scale=70, width=0.003, headwidth=3, headlength=4, alpha=0.7)
    
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

def plot_migrated_image(model, img, depths, title, save_filename):
    fig, ax = plt.subplots(figsize=(10, 6))
    
    im = ax.imshow(
        img, aspect='auto', cmap='hot', 
        extent=[0, model.n_blocks * model.block_width, np.max(depths), np.min(depths)]
    )
    
    ax.set_title(title)
    ax.set_xlabel('Length (m)', fontsize=12)
    ax.set_ylabel('Depth (m)', fontsize=12)
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label('Energy / Amplitude', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(save_filename, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Plot saved to: {save_filename}")