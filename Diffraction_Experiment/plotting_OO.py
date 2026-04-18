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
    
    for i, trace in enumerate(traces):
        x_base = (i + 0.5) * rx_dx
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

def plot_cwt_image(time, freqs, cwt_mag, title, save_filename):
    """Plots Continuous Wavelet Transform (Time-Frequency)."""
    fig, ax = plt.subplots(figsize=(10, 6))
    freqs_ghz = freqs / 1e9
    
    im = ax.imshow(
        cwt_mag, aspect='auto', cmap='jet', origin='lower',
        extent=[time[0], time[-1], freqs_ghz[0], freqs_ghz[-1]]
    )
    
    ax.set_title(title)
    ax.set_xlabel("Time (ns)", fontsize=12)
    ax.set_ylabel("Frequency (GHz)", fontsize=12)

    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
    ax.tick_params(axis='x', which='minor', length=4, color='k')

    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label('CWT Magnitude', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(save_filename, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Plot saved to: {save_filename}")

def plot_migrated_image(model, img, depths, title, save_filename):
    """Plots 2D depth migration."""
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