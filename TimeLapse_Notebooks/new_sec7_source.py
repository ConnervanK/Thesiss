# ─── Section 7: CLSSA Spectral Analysis ──────────────────────────────────────
# 4 rows × 5 columns per scenario.
# Rows 0–1: baseline antenna position (x_sc0)
# Rows 2–3: monitor  antenna position (x_sc_mon)
# Within each pair:
#   Even row (0, 2): baseline survey — BEFORE scatterer moves
#   Odd  row (1, 3): monitor  survey — AFTER  scatterer moves
#
# Col 0: wiggle trace (primary survey solid, other survey faint grey)
# Col 1: CLSSA amplitude  A(f,τ)      [dB]
# Col 2: CLSSA phase      φ(f,τ)      [°, raw amplitude-masked]
# Col 3: CLSSA phase gather S'(θ,τ)  — carrier-detrended
# Col 4 (even row): ΔΦ(f,τ) = φ_mon − φ_base  [2-D heatmap, joint-amp mask]
# Col 4 (odd  row): ΔΦ(f) at τ_sc  +  theoretical −2πfΔt  [1-D slice]

import sys as _sys
_sys.path.insert(0, r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks')
from phase_decomposition import clssa_phase_decomposition

dt_twt = float(2 * dz_mig / v_ice)
twt_ns = (2 * z_img / v_ice).astype(float)
twt_sc = float(2 * z_scatterer / v_ice)
y_lo   = float(twt_ns[-1])
y_hi   = float(twt_ns[0])

# Dominant-phase log: (scenario, row label, theta_dom [deg]) — filled below,
# printed once after all scenarios are processed so values can be compared.
_dominant_theta_log  = []
_ROW_LABELS_SHORT    = ['Baseline-Before', 'Baseline-After',
                         'Monitor-Before',  'Monitor-After']

for i in range(1, len(scenarios)):
    base_img = gazdag[0]
    mon_img  = gazdag[i]
    name     = str(scenarios[i])
    true_dx  = float(separation_lambda[i]) * lam
    x_sc0    = float(x_s1[0])
    x_sc_mon = float(x_s1[i])

    ix0 = int(np.argmin(np.abs(x_traces - x_sc0)))
    ixm = int(np.argmin(np.abs(x_traces - x_sc_mon)))

    base_col0 = base_img[:, ix0].astype(float)
    mon_col0  = mon_img[:,  ix0].astype(float)
    base_colm = base_img[:, ixm].astype(float)
    mon_colm  = mon_img[:,  ixm].astype(float)

    # ── CLSSA for all four traces ─────────────────────────────────────────
    _kw = dict(f_min_GHz=0.3, f_max_GHz=4.5, win_ns=2.0)
    A_b0, ph_b0, pg_b0, freqs_c, t_ax, theta_ax = clssa_phase_decomposition(base_col0, dt_twt, **_kw)
    A_m0, ph_m0, pg_m0, freqs_c, _,    theta_ax = clssa_phase_decomposition(mon_col0,  dt_twt, **_kw)
    A_bm, ph_bm, pg_bm, freqs_c, _,    theta_ax = clssa_phase_decomposition(base_colm, dt_twt, **_kw)
    A_mm, ph_mm, pg_mm, freqs_c, _,    theta_ax = clssa_phase_decomposition(mon_colm,  dt_twt, **_kw)

    it_sc = int(np.argmin(np.abs(t_ax - twt_sc)))

    # ── 1-D ΔΦ slices at τ_sc ────────────────────────────────────────────
    def _dphi_1d(ph_mon, ph_base, it):
        return np.rad2deg(np.angle(np.exp(1j * np.deg2rad(
            ph_mon[:, it] - ph_base[:, it]))))

    dphi_1d_0 = _dphi_1d(ph_m0, ph_b0, it_sc)
    dphi_1d_m = _dphi_1d(ph_mm, ph_bm, it_sc)

    # ── 2-D ΔΦ(f,τ) heatmaps, masked where joint amplitude < 5 % ─────────
    def _dphi_2d(ph_mon, ph_base, A_mon, A_base, thr=0.05):
        raw = np.rad2deg(np.angle(np.exp(1j * np.deg2rad(ph_mon - ph_base))))
        A_j = np.minimum(A_mon, A_base)
        return np.where(A_j / (A_j.max() + 1e-30) < thr, np.nan, raw)

    dphi_2d_0 = _dphi_2d(ph_m0, ph_b0, A_m0, A_b0)
    dphi_2d_m = _dphi_2d(ph_mm, ph_bm, A_mm, A_bm)

    # ── Theoretical Δt [ns] per antenna position ──────────────────────────
    _hyp   = float(np.sqrt(true_dx**2 + z_scatterer**2))
    dt_th0 = 2.0 * (_hyp - z_scatterer) / v_ice   # x_sc0:    scatterer moved away
    dt_thm = 2.0 * (z_scatterer - _hyp) / v_ice   # x_sc_mon: scatterer moved to

    # ── Carrier-phase detrending for the phase gather ─────────────────────
    _d_theta_pg = float(theta_ax[1] - theta_ax[0])
    _n_theta_pg = len(theta_ax)
    _shift_bins = (np.round(360.0 * f_c * t_ax / _d_theta_pg).astype(int)
                   % _n_theta_pg)

    def _detrend_gather(pg):
        return np.stack(
            [np.roll(pg[:, _it], -_shift_bins[_it]) for _it in range(pg.shape[1])],
            axis=1
        )

    ext_tf = [freqs_c[0], freqs_c[-1], y_lo, y_hi]
    ext_pg = [theta_ax[0], theta_ax[-1], y_lo, y_hi]

    # ── Row specifications ────────────────────────────────────────────────
    # Each tuple:
    #   (title, trace_primary, trace_ref, A, ph, pg,
    #    col4_type, dt_th, dphi_1d, dphi_2d, A_b_sc, A_m_sc)
    _ROW_SPECS = [
        (
            f'Baseline pos  x={x_sc0:.3f} m  —  BEFORE  (baseline survey)',
            base_col0, mon_col0,
            A_b0, ph_b0, pg_b0,
            '2d', dt_th0, dphi_1d_0, dphi_2d_0,
            A_b0[:, it_sc], A_m0[:, it_sc],
        ),
        (
            f'Baseline pos  x={x_sc0:.3f} m  —  AFTER   (monitor survey)',
            mon_col0, base_col0,
            A_m0, ph_m0, pg_m0,
            '1d', dt_th0, dphi_1d_0, dphi_2d_0,
            A_b0[:, it_sc], A_m0[:, it_sc],
        ),
        (
            f'Monitor pos   x={x_sc_mon:.3f} m  —  BEFORE  (baseline survey)',
            base_colm, mon_colm,
            A_bm, ph_bm, pg_bm,
            '2d', dt_thm, dphi_1d_m, dphi_2d_m,
            A_bm[:, it_sc], A_mm[:, it_sc],
        ),
        (
            f'Monitor pos   x={x_sc_mon:.3f} m  —  AFTER   (monitor survey)',
            mon_colm, base_colm,
            A_mm, ph_mm, pg_mm,
            '1d', dt_thm, dphi_1d_m, dphi_2d_m,
            A_bm[:, it_sc], A_mm[:, it_sc],
        ),
    ]

    # ── Figure ─────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(4, 5, figsize=(28, 20))
    fig.suptitle(
        f'Spectral Line (CLSSA)  —  Gazdag  |  {name}'
        f'   (Δx = {true_dx*1e3:.1f} mm = {true_dx/lam:.4f}λ)',
        fontsize=12, fontweight='bold',
    )

    # Faint background: blue tint for baseline-pos rows, orange tint for monitor-pos rows
    for _ri in range(4):
        _bg = '#eef2ff' if _ri < 2 else '#fff3ee'
        for _ci in range(5):
            axes[_ri, _ci].set_facecolor(_bg)

    for row_idx, (
        row_title, trace_p, trace_r,
        A_cl, ph_cl, pg_cl, col4_type, dt_th,
        dphi_1d, dphi_2d, A_b_sc, A_m_sc,
    ) in enumerate(_ROW_SPECS):

        _clr_p = 'steelblue' if row_idx % 2 == 0 else 'tomato'

        # ── Col 0: wiggle trace ─────────────────────────────────────────────
        ax = axes[row_idx, 0]
        ax.plot(trace_r, twt_ns, color='grey',  lw=0.8, alpha=0.28, label='other survey')
        ax.plot(trace_p, twt_ns, color=_clr_p,  lw=1.5, label='this survey')
        ax.fill_betweenx(twt_ns, trace_p, 0, where=(trace_p > 0),
                         color=_clr_p, alpha=0.18)
        ax.axhline(twt_sc, color='grey', ls='--', lw=0.9, alpha=0.7,
                   label=f'z_sc ({twt_sc:.3f} ns)')
        ax.axvline(0, color='k', lw=0.5, alpha=0.4)
        ax.set_ylim(y_lo, y_hi)
        ax.set_ylabel('TWT [ns]', fontsize=9)
        ax.set_xlabel('Amplitude [a.u.]', fontsize=9)
        ax.set_title(row_title, fontsize=8.5)
        ax.legend(fontsize=7.5, loc='lower right')
        ax.grid(True, alpha=0.25)

        # ── Col 1: CLSSA amplitude spectrum ────────────────────────────────
        ax = axes[row_idx, 1]
        A_dB = np.clip(20 * np.log10(A_cl / (A_cl.max() + 1e-30) + 1e-30), -40, 0)
        im1 = ax.imshow(A_dB.T, aspect='auto', origin='upper', extent=ext_tf,
                        cmap='jet', vmin=-40, vmax=0, interpolation='bilinear')
        ax.axhline(twt_sc, color='white', ls='--', lw=0.9, alpha=0.8,
                   label=f'z_sc ({twt_sc:.3f} ns)')
        ax.set_xlabel('Frequency [GHz]', fontsize=9)
        ax.set_ylabel('TWT [ns]', fontsize=9)
        ax.set_title('CLSSA Amplitude [dB]', fontsize=9)
        ax.legend(fontsize=7.5)
        fig.colorbar(im1, ax=ax, fraction=0.046, pad=0.04, label='[dB]')

        # ── Col 2: CLSSA phase spectrum φ(f,τ) — raw, amplitude-masked ──────────────
        ax = axes[row_idx, 2]
        ph_msk = np.where(A_cl / (A_cl.max() + 1e-30) > 0.05, ph_cl, np.nan)
        im2 = ax.imshow(ph_msk.T, aspect='auto', origin='upper', extent=ext_tf,
                        cmap='hsv', vmin=-180, vmax=180, interpolation='bilinear')
        ax.axhline(twt_sc, color='white', ls='--', lw=0.9, alpha=0.8,
                   label=f'z_sc ({twt_sc:.3f} ns)')
        ax.set_xlabel('Frequency [GHz]', fontsize=9)
        ax.set_ylabel('TWT [ns]', fontsize=9)
        ax.set_title(r'CLSSA Phase  $\phi(f,\tau)$  [°]', fontsize=9)
        ax.legend(fontsize=7.5)
        fig.colorbar(im2, ax=ax, fraction=0.046, pad=0.04, label='[°]')

        # ── Col 3: CLSSA phase gather — carrier-detrended ──────────────────
        ax = axes[row_idx, 3]
        pg_dt   = _detrend_gather(pg_cl)
        pg_clip = float(np.percentile(np.abs(pg_dt), 99)) or 1.0
        im3 = ax.imshow(pg_dt.T, aspect='auto', origin='upper', extent=ext_pg,
                        cmap='RdBu', vmin=-pg_clip, vmax=pg_clip,
                        interpolation='bilinear')
        ax.axhline(twt_sc, color='grey', ls='--', lw=0.9, alpha=0.7)
        for ph_line in (-90, 90):
            ax.axvline(ph_line, color='lime', lw=1.0, ls='--', alpha=0.85)
        ax.axvline(0, color='white', lw=0.6, ls=':', alpha=0.5)

        # Dominant phase: energy-weighted (integrated over all TWT) peak of
        # the detrended gather, so it is not tied to a single TWT sample.
        _theta_energy = np.nansum(np.abs(pg_dt), axis=1)   # (n_theta,)
        _idx_dom      = int(np.argmax(_theta_energy))
        theta_dom     = float(theta_ax[_idx_dom])
        ax.axvline(theta_dom, color='black', lw=1.4, ls='-', alpha=0.9, zorder=5)
        ax.text(0.02, 0.04, f'θ_dom = {theta_dom:.1f}°', transform=ax.transAxes,
                fontsize=8, color='black', ha='left', va='bottom',
                bbox=dict(boxstyle='round', fc='white', alpha=0.75, ec='none'))
        _dominant_theta_log.append((name, _ROW_LABELS_SHORT[row_idx], theta_dom))

        ax.set_xticks([-180, -90, 0, 90, 180])
        ax.set_xlabel('Phase [°]', fontsize=9)
        ax.set_ylabel('TWT [ns]', fontsize=9)
        ax.set_title("Phase Gather S'(θ,t) — detrended", fontsize=9)
        fig.colorbar(im3, ax=ax, fraction=0.046, pad=0.04, label='Amplitude')

        # ── Col 4 ──────────────────────────────────────────────────────────
        ax = axes[row_idx, 4]

        if col4_type == '2d':
            # 2-D ΔΦ(f,τ) heatmap: φ_mon − φ_base, masked by joint amplitude
            im4 = ax.imshow(dphi_2d.T, aspect='auto', origin='upper',
                            extent=ext_tf, cmap='hsv',
                            vmin=-180, vmax=180, interpolation='bilinear')
            ax.axhline(twt_sc, color='white', ls='--', lw=0.9, alpha=0.8,
                       label=f'z_sc ({twt_sc:.3f} ns)')
            ax.set_xlabel('Frequency [GHz]', fontsize=9)
            ax.set_ylabel('TWT [ns]', fontsize=9)
            ax.set_title(
                r'$\Delta\Phi(f,\tau) = \phi_\mathrm{mon} - \phi_\mathrm{base}$  [°]',
                fontsize=9,
            )
            ax.legend(fontsize=7.5)
            fig.colorbar(im4, ax=ax, fraction=0.046, pad=0.04, label='[°]')

        else:
            # Group delay: Δt(f) = -(1/2π) ∂ΔΦ/∂f
            # Unwrap the masked ΔΦ(f) (deg -> rad) to remove ±180° wraps, then
            # differentiate wrt frequency. This recovers the frequency-dependent
            # arrival-time shift directly, which should converge to the constant
            # ray-theory Δt at high frequency and curve away from it at low
            # frequency (geometric dispersion of the Gazdag migration operator).
            A_joint  = np.minimum(A_b_sc, A_m_sc)
            thr_amp  = 0.05 * A_joint.max() if A_joint.max() > 0 else np.inf
            dphi_msk = np.where(A_joint >= thr_amp, dphi_1d, np.nan)

            _valid = ~np.isnan(dphi_msk)
            dt_f = np.full_like(freqs_c, np.nan)
            if np.count_nonzero(_valid) > 2:
                f_v   = freqs_c[_valid]
                phi_v = np.unwrap(np.deg2rad(dphi_msk[_valid]))
                dt_f[_valid] = -1.0 / (2 * np.pi) * np.gradient(phi_v, f_v)

            ax_amp = ax.twinx()
            ax_amp.fill_between(freqs_c, A_b_sc, alpha=0.14, color='steelblue',
                                label='A_base(f)')
            ax_amp.fill_between(freqs_c, A_m_sc, alpha=0.14, color='tomato',
                                label='A_mon(f)')
            ax_amp.set_ylabel('Amplitude [a.u.]', fontsize=8, color='grey')
            ax_amp.tick_params(axis='y', labelcolor='grey', labelsize=7)
            ax_amp.set_ylim(bottom=0)
            ax_amp.legend(fontsize=7, loc='upper left')

            ax.plot(freqs_c, dt_f, color='darkorange', lw=2.2, zorder=3,
                    label=r'$\Delta t(f)=-\frac{1}{2\pi}\partial\Delta\Phi/\partial f$')
            ax.axhline(dt_th, color='limegreen', lw=1.8, ls='--', zorder=4,
                       label=f'ray theory  (Δt={dt_th:.6f} ns)')
            ax.axhline(0, color='k', lw=0.7, alpha=0.5, ls=':')

            _valid_dt = dt_f[~np.isnan(dt_f)]
            _dmax = float(np.max(np.abs(_valid_dt))) if len(_valid_dt) > 0 else abs(dt_th)
            _dmax = max(_dmax, abs(dt_th), 1e-6)
            _pad  = _dmax * 0.3
            ax.set_ylim(-(_dmax + _pad), _dmax + _pad)

            ax.set_xlabel('Frequency [GHz]', fontsize=9)
            ax.set_ylabel('Δt(f)  [ns]', fontsize=9)
            ax.set_title(
                f'Group Delay  Δt(f)   (z_sc={z_scatterer*1e3:.0f} mm)',
                fontsize=9,
            )
            ax.legend(fontsize=7.5, loc='upper right')
            ax.grid(True, alpha=0.25, zorder=0)

    plt.tight_layout()
    plt.show()

# ── Cross-scenario comparison of the dominant phase-gather angle ────────────
print()
print("Dominant phase angle theta_dom [deg]  (energy-weighted peak of the "
      "detrended phase gather S'(theta,t), integrated over all TWT):")
print(f'{"Scenario":<10}{"Row":<18}{"theta_dom [deg]":>16}')
for _name, _row_lbl, _th in _dominant_theta_log:
    print(f'{_name:<10}{_row_lbl:<18}{_th:>16.2f}')
