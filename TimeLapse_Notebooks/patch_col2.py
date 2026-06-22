import json, sys

NB_PATH = r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\TimeLapse_Processing.ipynb'
CELL_ID = '2790c368'

sys.stdout.reconfigure(encoding='utf-8')

# ── exact old → new replacements ─────────────────────────────────────────────

OLD_HEADER = "#   Col 2: CLSSA phase spectrum φ(f,τ) — carrier-detrended [Frequency \xd7 TWT]\n"
NEW_HEADER = "#   Col 2: CLSSA phase spectrum φ(f,τ) [Frequency \xd7 TWT, colour = phase °]\n"

OLD_HELPERS = (
    "    # ── Carrier-phase detrending helpers ──────────────────────────────────\n"
    "    # The CLSSA phase φ(f, τ) accumulates the propagation phase 2πf·τ, which\n"
    "    # causes diagonal bands in both the phase spectrum and the phase gather.\n"
    "    # Subtracting 360·f·τ (degrees) from each (f,τ) pixel removes this and\n"
    "    # reveals only the reflector / wavelet phase φ_0(f).\n"
    "    #\n"
    "    # For the phase gather (marginalized over f), the dominant carrier is f_c,\n"
    "    # so we circularly shift each TWT column by −360·f_c·τ / Δθ bins.\n"
    "    _carrier_2d  = 360.0 * np.outer(freqs_c, t_ax)   # (n_f, n_t) degrees\n"
    "\n"
    "    _d_theta_pg  = float(theta_ax[1] - theta_ax[0])  # degrees per phase-gather bin\n"
    "    _n_theta_pg  = len(theta_ax)\n"
    "    _shift_bins  = (np.round(360.0 * f_c * t_ax / _d_theta_pg).astype(int)\n"
    "                    % _n_theta_pg)   # (n_t,)\n"
    "\n"
    "    def _detrend_phase(ph, A):\n"
    "        # subtract 360*f*t, wrap to +-180 deg, mask noise\n"
    "        ph_dt = np.rad2deg(np.angle(np.exp(1j * np.deg2rad(ph - _carrier_2d))))\n"
    "        return np.where(A / (A.max() + 1e-30) > 0.05, ph_dt, np.nan)\n"
    "\n"
    "    def _detrend_gather(pg):\n"
    "        # circularly shift each TWT column by -360*f_c*t bins\n"
    "        return np.stack(\n"
    "            [np.roll(pg[:, _it], -_shift_bins[_it]) for _it in range(pg.shape[1])],\n"
    "            axis=1\n"
    "        )\n"
)

NEW_HELPERS = (
    "    # ── Carrier-phase detrending for the phase gather ─────────────────────\n"
    "    # The phase gather accumulates propagation phase 2πf_c·τ, creating diagonal\n"
    "    # bands. Circularly shifting each TWT column by −360·f_c·τ/Δθ bins removes it.\n"
    "    _d_theta_pg  = float(theta_ax[1] - theta_ax[0])  # degrees per bin\n"
    "    _n_theta_pg  = len(theta_ax)\n"
    "    _shift_bins  = (np.round(360.0 * f_c * t_ax / _d_theta_pg).astype(int)\n"
    "                    % _n_theta_pg)   # (n_t,)\n"
    "\n"
    "    def _detrend_gather(pg):\n"
    "        # circularly shift each TWT column by -360*f_c*t bins\n"
    "        return np.stack(\n"
    "            [np.roll(pg[:, _it], -_shift_bins[_it]) for _it in range(pg.shape[1])],\n"
    "            axis=1\n"
    "        )\n"
)

OLD_COL2 = (
    "        # ── Col 2: CLSSA phase spectrum φ(f,τ)  — carrier-detrended ────────\n"
    "        ax = axes[row_idx, 2]\n"
    "        ph_msk = _detrend_phase(ph_cl, A_cl)   # subtract 360·f·τ, mask noise\n"
    "        im2 = ax.imshow(ph_msk.T, aspect='auto', origin='upper', extent=ext_tf,\n"
    "                        cmap='hsv', vmin=-180, vmax=180, interpolation='bilinear')\n"
    "        ax.axhline(twt_sc, color='white', ls='--', lw=0.9, alpha=0.8,\n"
    "                   label=f'z_sc ({twt_sc:.3f} ns)')\n"
    "        ax.set_xlabel('Frequency [GHz]', fontsize=9)\n"
    "        ax.set_ylabel('TWT [ns]', fontsize=9)\n"
    "        ax.set_title(r'CLSSA Phase $\\phi(f,\\tau)$ — detrended [°]', fontsize=9)\n"
    "        ax.legend(fontsize=7.5)\n"
    "        fig.colorbar(im2, ax=ax, fraction=0.046, pad=0.04, label='[°]')\n"
)

NEW_COL2 = (
    "        # ── Col 2: CLSSA phase spectrum φ(f,τ) ──────────────────────────────\n"
    "        ax = axes[row_idx, 2]\n"
    "        ph_msk = np.where(A_cl / (A_cl.max() + 1e-30) > 0.05, ph_cl, np.nan)\n"
    "        im2 = ax.imshow(ph_msk.T, aspect='auto', origin='upper', extent=ext_tf,\n"
    "                        cmap='hsv', vmin=-180, vmax=180, interpolation='bilinear')\n"
    "        ax.axhline(twt_sc, color='white', ls='--', lw=0.9, alpha=0.8,\n"
    "                   label=f'z_sc ({twt_sc:.3f} ns)')\n"
    "        ax.set_xlabel('Frequency [GHz]', fontsize=9)\n"
    "        ax.set_ylabel('TWT [ns]', fontsize=9)\n"
    "        ax.set_title(r'CLSSA Phase  $\\phi(f,\\tau)$  [°]', fontsize=9)\n"
    "        ax.legend(fontsize=7.5)\n"
    "        fig.colorbar(im2, ax=ax, fraction=0.046, pad=0.04, label='[°]')\n"
)

REPLACEMENTS = [
    (OLD_HEADER,   NEW_HEADER),
    (OLD_HELPERS,  NEW_HELPERS),
    (OLD_COL2,     NEW_COL2),
]

# ── load notebook ─────────────────────────────────────────────────────────────
with open(NB_PATH, encoding='utf-8') as f:
    nb = json.load(f)

patched = False
for c in nb['cells']:
    if c.get('id') == CELL_ID:
        src = c['source']
        if isinstance(src, list):
            src = ''.join(src)          # work on a single string
        for old, new in REPLACEMENTS:
            if old not in src:
                print(f'WARNING: pattern not found (first 80 chars): {old[:80]!r}')
            src = src.replace(old, new)
        c['source'] = src               # store as single string (valid .ipynb)
        c['outputs'] = []
        c['execution_count'] = None
        patched = True
        break

if not patched:
    print('ERROR: cell not found', file=sys.stderr)
    sys.exit(1)

with open(NB_PATH, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print('OK — col2 reverted to raw phase; _carrier_2d/_detrend_phase removed')
