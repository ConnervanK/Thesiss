# ─── Section 8: Cross-Phase Spectrogram (TWT vs. Frequency) ──────────────────
# For the central trace (x = x_sc0, the scatterer's baseline position — the
# fixed reference column shared by every monitor scenario), use CLSSA (the
# same narrowband decomposition as Section 7) to get the baseline and monitor
# phase spectra φ_base(f,τ) and φ_mon(f,τ), then form the cross-phase
#     ΔΦ(τ,f) = angle[ exp(i·(φ_base − φ_mon)) ]   (radians, wrapped to ±π)
# which is the phase of the complex cross-spectrum XS = S_base · S_mon* without
# ever forming S explicitly. CLSSA samples frequency directly (no FFT
# bin-spacing limit like the earlier STFT version had), so the frequency axis
# is smooth regardless of analysis window length.
#
# No amplitude masking — the whole point is to see the noise-vs-coherence
# contrast across the full image, not just the signal-bearing band.
#
# For separations below 1/4 lambda the phase change at the scatterer is very
# subtle, so those scenarios additionally:
#   - zoom the TWT axis in on a window around the scatterer depth
#   - clip the colorbar to a much smaller range, to make the subtle phase
#     contrast visible

import sys as _sys
_sys.path.insert(0, r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks')
from TimeLapse_Notebooks.helper_functions.phase_decomposition import clssa_phase_decomposition

dt_twt = float(2 * dz_mig / v_ice)
twt_sc = float(2 * z_scatterer / v_ice)

_clssa_kw = dict(f_min_GHz=0.3, f_max_GHz=4.5, win_ns=2.0)

_SMALL_SEP_THRESH   = 0.25       # lambda — scenarios below this get zoom + clip
_ZOOM_HALF_WIDTH_NS = 1.5        # ns around twt_sc to zoom in on
_CLIP_SMALL         = np.pi / 6  # rad (~30°) colorbar clip for small-sep cases
_CLIP_FULL          = np.pi      # rad (180°) colorbar clip otherwise

for i in range(1, len(scenarios)):
    base_img = gazdag[0]
    mon_img  = gazdag[i]
    name     = str(scenarios[i])
    true_dx  = float(separation_lambda[i]) * lam
    x_sc0    = float(x_s1[0])

    ix0 = int(np.argmin(np.abs(x_traces - x_sc0)))
    base_trace = base_img[:, ix0].astype(float)
    mon_trace  = mon_img[:,  ix0].astype(float)

    A_b, ph_b, _, freqs_c, t_ax, _ = clssa_phase_decomposition(base_trace, dt_twt, **_clssa_kw)
    A_m, ph_m, _, freqs_c, _,    _ = clssa_phase_decomposition(mon_trace,  dt_twt, **_clssa_kw)

    # Cross-phase ΔΦ(f,τ) = φ_base − φ_mon, wrapped to (−π, +π]
    phi_xs = np.angle(np.exp(1j * np.deg2rad(ph_b - ph_m)))   # (n_f, n_t) radians

    is_small_sep = float(separation_lambda[i]) < _SMALL_SEP_THRESH
    clip = _CLIP_SMALL if is_small_sep else _CLIP_FULL

    ext = [freqs_c[0], freqs_c[-1], t_ax[-1], t_ax[0]]

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(phi_xs.T, aspect='auto', origin='upper', extent=ext,
                   cmap='hsv', vmin=-clip, vmax=clip, interpolation='nearest')
    ax.axhline(twt_sc, color='white', ls='--', lw=1.1, alpha=0.85,
               label=f'τ_sc ({twt_sc:.3f} ns)')

    _zoom_note = ''
    if is_small_sep:
        ax.set_ylim(twt_sc + _ZOOM_HALF_WIDTH_NS, twt_sc - _ZOOM_HALF_WIDTH_NS)
        _zoom_note = f'  [zoomed ±{_ZOOM_HALF_WIDTH_NS:.1f} ns, clip=±{np.degrees(clip):.0f}°]'

    ax.set_xlabel('Frequency [GHz]', fontsize=10)
    ax.set_ylabel('TWT [ns]', fontsize=10)
    ax.set_title(
        f'Cross-Phase Spectrogram (CLSSA)  ΔΦ(τ,f)  —  Gazdag  |  {name}\n'
        f'x = {x_sc0:.3f} m   (Δx = {true_dx*1e3:.1f} mm = {true_dx/lam:.4f}λ){_zoom_note}',
        fontsize=9.5, fontweight='bold'
    )
    ax.legend(fontsize=8)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='ΔΦ  [rad]')
    if is_small_sep:
        cbar.set_ticks([-clip, -clip / 2, 0, clip / 2, clip])
        cbar.set_ticklabels([f'{np.degrees(-clip):.0f}°', f'{np.degrees(-clip / 2):.0f}°',
                              '0', f'{np.degrees(clip / 2):.0f}°', f'{np.degrees(clip):.0f}°'])
    else:
        cbar.set_ticks([-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
        cbar.set_ticklabels(['-π', '-π/2', '0', 'π/2', 'π'])

    plt.tight_layout()
    plt.show()
