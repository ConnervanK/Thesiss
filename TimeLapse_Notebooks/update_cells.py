"""Update B_Scan_Migration_Playground.ipynb cells rs-load and lsm-debug."""
import json

NEW_RS_LOAD = """\
# Resolution Study --- Data Loading + Migration (Kirchhoff, Gazdag, Backprop, LSM)
import io, contextlib
import pylops
from scipy.signal import decimate as _scipy_decimate, fftconvolve as _fftconvolve

wavelength_ice = v_ice / f_c_GHz   # approx 0.1126 m  (1.5 GHz in ice)

STUDY_ROOT = BASE / 'resolution_study'
BG_PATH    = BASE / 'B_Scan_Migration_data' / 'bg_bscan_merged.out'

LAMBDA_FACTORS = [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625]

def _study_dir(factor):
    label = f'{factor:g}'.replace('.', 'p')
    return STUDY_ROOT / f'study_{label}lambda'


# ── Pure DAS Kirchhoff PyLops operator (no wavelet convolution) ────────────────
class KirchhoffDAS_Op(pylops.LinearOperator):
    \"\"\"Pure delay-and-sum Kirchhoff operator: no wavelet, no obliquity.

    Adjoint (_rmatvec) exactly mirrors kirchhoff_poststack:
        model[ix, iz] = sum_tr  data_mf[tr, r(ix,iz,tr)/v_mig/dt]

    Use with matched-filtered data whose peaks sit at TWTT = 2r/v_ice from t=0.

    dims  = (n_x_img, n_z)   model on (x_img, z_img) grid
    dimsd = (n_tr, nt)       data on x_traces time grid
    \"\"\"

    def __init__(self, x_img, z_img, x_traces, nt, dt, v_mig,
                 aperture_deg=90.0, dtype='float64'):
        nx   = len(x_img)
        nz   = len(z_img)
        n_tr = len(x_traces)
        ni   = nx * nz

        xi = x_img[:, None, None]
        zi = z_img[None, :, None]
        xt = x_traces[None, None, :]
        dx = xi - xt
        r  = np.sqrt(dx**2 + zi**2)
        t_frac = r / v_mig / dt            # (nx, nz, n_tr)

        itrav = t_frac.astype(np.int32)
        travd = (t_frac - itrav).astype(np.float32)
        valid = (itrav >= 0) & (itrav < nt - 1)

        if aperture_deg < 90.0:
            safe_z = np.where(np.abs(zi) > 1e-9, np.abs(zi), 1e-9)
            angle_deg = np.rad2deg(np.abs(np.arctan2(np.abs(dx), safe_z)))
            valid &= (angle_deg <= aperture_deg)

        itrav = itrav.reshape(ni, n_tr)
        travd = travd.reshape(ni, n_tr)
        valid = valid.reshape(ni, n_tr)

        ii_v, i_v = np.where(valid)
        self._ii_v = ii_v
        self._i_v  = i_v
        self._i0_v = itrav[ii_v, i_v].astype(np.int32)
        self._td_v = travd[ii_v, i_v].astype(np.float64)

        self.ni     = ni
        self.n_tr   = n_tr
        self.nt     = nt
        self._itrav = itrav
        self._travd = travd.astype(np.float64)
        self._valid = valid
        self._tr    = np.arange(n_tr, dtype=np.intp)

        super().__init__(dtype=np.dtype(dtype), dims=(nx, nz),
                         dimsd=(n_tr, nt), name='KirchhoffDAS')

    def _matvec(self, x):
        m = x.ravel()
        d = np.zeros((self.n_tr, self.nt), dtype=self.dtype)
        m_v = m[self._ii_v]
        np.add.at(d, (self._i_v, self._i0_v),     m_v * (1.0 - self._td_v))
        np.add.at(d, (self._i_v, self._i0_v + 1), m_v * self._td_v)
        return d.ravel()

    def _rmatvec(self, x):
        d   = x.reshape(self.n_tr, self.nt)
        i0  = np.clip(self._itrav, 0, self.nt - 2)
        td  = self._travd
        samp = ((1.0 - td) * d[self._tr[np.newaxis, :], i0    ]
                +       td  * d[self._tr[np.newaxis, :], i0 + 1])
        return (samp * self._valid).sum(axis=1)


# ── Matched filter (absorbs t0, corrects PEC polarity) ────────────────────────
def matched_filter_bscan(bscan, dt_ns, t0_ns, f_c_GHz):
    \"\"\"Cross-correlate each trace with the gprMax source Ricker; negate for PEC.

    gprMax background-subtracted data:  d[tr,t] ~= -Ricker(t - TWTT - t0_ns)
    (PEC scatterer inverts polarity; Ricker peaks at t0_ns in gprMax time).

    After matched filter + negation:  mf[tr,t] ~= autocorr(Ricker)(t - TWTT)
    -> positive peaks at TWTT from t=0 (t0_ns absorbed, PEC polarity corrected).

    Returns mf_data with same shape as bscan.
    \"\"\"
    n_tr, nt = bscan.shape
    t_grid   = np.arange(nt) * dt_ns
    # Reference wavelet: gprMax Ricker on the data time grid (peak at t0_ns)
    t_c = t_grid - t0_ns
    fc  = f_c_GHz
    h   = (1.0 - 2*(np.pi*fc*t_c)**2) * np.exp(-(np.pi*fc*t_c)**2)
    h  /= np.linalg.norm(h)

    h_rev   = h[::-1].copy()                 # time-reverse (= h for symmetric Ricker)
    mf_data = np.empty_like(bscan, dtype=np.float64)
    for tr in range(n_tr):
        conv_full   = _fftconvolve(bscan[tr].astype(np.float64), h_rev, mode='full')
        # causal portion starts at index nt-1; negate to correct PEC polarity
        mf_data[tr] = -conv_full[nt - 1: 2*nt - 1]
    return mf_data


# ── LSM: matched-filter -> decimate -> KirchhoffDAS -> FISTA ──────────────────
def lsm_poststack(bscan, dt_ns, x_traces, v_mig, x_img, z_img, t0_ns,
                  f_c_GHz=1.5, aperture_deg=40.0, n_iter=200, eps_frac=0.08):
    \"\"\"LSM with matched-filter preprocessing and pure DAS Kirchhoff operator.

    v_mig : use v_ice/2  (same half-velocity convention as kirchhoff_poststack).
    The matched filter absorbs t0_ns and corrects PEC polarity so arrivals become
    positive autocorrelation peaks at TWTT = 2r/v_ice from t=0.
    The KirchhoffDAS_Op adjoint then equals kirchhoff_poststack exactly.
    \"\"\"
    # Step 1: matched filter
    mf_data = matched_filter_bscan(bscan, dt_ns, t0_ns, f_c_GHz)

    # Step 2: decimate to target dt ~ 0.15 ns (tractable FISTA)
    dec = max(1, int(0.15 / dt_ns))
    if dec > 1:
        mf_dec = _scipy_decimate(mf_data, dec, axis=1, zero_phase=True).astype(np.float64)
        dt_d   = dt_ns * dec
    else:
        mf_dec = mf_data
        dt_d   = dt_ns
    nt_d = mf_dec.shape[1]
    print(f'    MF+decimate {dec}x  dt={dt_d:.4f} ns  n_t={nt_d}')

    # Step 3: build pure DAS operator
    K = KirchhoffDAS_Op(x_img, z_img, x_traces, nt_d, dt_d, v_mig,
                        aperture_deg=aperture_deg)

    # Step 4: FISTA L1
    d_flat   = mf_dec.ravel()
    adj_flat = K.H @ d_flat
    adj_peak = float(np.max(np.abs(adj_flat)))
    eps      = eps_frac * adj_peak
    print(f'    adj_peak={adj_peak:.3e}  eps={eps:.3e}')
    m_flat   = pylops.optimization.sparsity.fista(
        K, d_flat, eps=eps, niter=n_iter, show=False
    )[0]

    img  = m_flat.reshape(len(x_img), len(z_img)).T   # (nz, nx_img)
    peak = np.max(np.abs(img))
    return img / (peak + 1e-30)


# ── Load background and run all separations ───────────────────────────────────
bscan_bg_rs, _dt_bg, _ = load_gprmax_out(BG_PATH)
print(f'Background B-scan: {bscan_bg_rs.shape}  dt={_dt_bg:.6f} ns')

rs_results = []

for _factor in LAMBDA_FACTORS:
    _delta = _factor * wavelength_ice
    _xs1   = round((2.0 - _delta / 2) / dx_grid) * dx_grid
    _xs2   = round((2.0 + _delta / 2) / dx_grid) * dx_grid

    _psf_path = _study_dir(_factor) / 'psf_bscan_merged.out'
    print(f'\\n{_factor:7g}lambda  Dx={_delta*1e3:.2f} mm  '
          f'xs1={_xs1:.4f} m  xs2={_xs2:.4f} m')

    if not _psf_path.exists():
        print(f'  WARNING: {_psf_path} not found --- skipping')
        continue

    _full, _dt, _tn = load_gprmax_out(_psf_path)
    _sc = _full - bscan_bg_rs

    _ikh  = kirchhoff_poststack(_sc, _dt, x_traces, v_mig, x_img, z_img, t0_ns)
    _igz  = gazdag_poststack(   _sc, _dt, x_traces, v_mig, x_img, z_img, t0_ns)
    _ibp  = backprop_poststack( _sc, _dt, x_traces, v_ice,  x_img, z_img, t0_ns)
    print('  LSM (MF + DAS + FISTA):')
    _ilsm = lsm_poststack(_sc, _dt, x_traces, v_mig, x_img, z_img, t0_ns)

    rs_results.append(dict(
        factor=_factor, delta_m=_delta,
        xs1=_xs1, xs2=_xs2,
        bscan_sc=_sc, dt_ns=_dt, time_ns=_tn,
        img_kh=_ikh, img_gz=_igz, img_bp=_ibp, img_lsm=_ilsm,
    ))
    print(f'  KH={np.max(np.abs(_ikh)):.3f}  GZ={np.max(np.abs(_igz)):.3f}  '
          f'BP={np.max(_ibp):.3f}  LSM={np.max(np.abs(_ilsm)):.3f}  done')

print(f'\\nLoaded {len(rs_results)} separations.')
"""

NEW_LSM_DEBUG = """\
# ── LSM Debug Cell ─── tune settings before running full resolution study ──────────────
# Row 0: adjoint alignment
#   [0,0] Kirchhoff DAS  (reference)
#   [0,1] KirchhoffDAS_Op adjoint on MATCHED-FILTERED + decimated data  <- should match [0,0]
#   [0,2] Matched-filter trace (central) vs raw trace  (time-domain sanity check)
# Row 1: FISTA with two eps values and LSM from rs_results
# ─────────────────────────────────────────────────────────────────────────────────────────
TEST_FACTOR  = 2.0
ANGLE_AP     = 40.0
EPS_FRAC_A   = 0.04
EPS_FRAC_B   = 0.08
N_ITER       = 200
# ─────────────────────────────────────────────────────────────────────────────────────────

import pylops
from scipy.signal import decimate as _scipy_decimate, fftconvolve as _fftconvolve

_tr_r = next((r for r in rs_results if r['factor'] == TEST_FACTOR), rs_results[0])
_sc, _dt = _tr_r['bscan_sc'], _tr_r['dt_ns']
_xs1, _xs2 = _tr_r['xs1'], _tr_r['xs2']
_dec = max(1, int(0.15 / _dt))
print(f"factor={_tr_r['factor']}lam  dt={_dt:.6f} ns  dec={_dec}x  "
      f"Dx={_tr_r['delta_m']*1e3:.1f} mm  xs1={_xs1:.4f}  xs2={_xs2:.4f}")

# ── Matched filter ────────────────────────────────────────────────────────────
print('Matched filter...')
_mf = matched_filter_bscan(_sc, _dt, t0_ns, f_c_GHz)

# Decimated MF data
_mf_dec = _scipy_decimate(_mf, _dec, axis=1, zero_phase=True).astype(np.float64)
_dt_dec = _dt * _dec
_time_dec = np.arange(_mf_dec.shape[1]) * _dt_dec
print(f'Decimated MF: n_t={_mf_dec.shape[1]}  dt={_dt_dec:.4f} ns')

# ── Build DAS operator ────────────────────────────────────────────────────────
_K = KirchhoffDAS_Op(x_img, z_img, x_traces, _mf_dec.shape[1], _dt_dec, v_mig,
                     aperture_deg=ANGLE_AP)

# DAS adjoint on MF data  (should match Kirchhoff DAS)
_d_flat  = _mf_dec.ravel()
_adj_flat = _K.H @ _d_flat
_adj_peak = float(np.max(np.abs(_adj_flat)))
_img_adj_mf = (_adj_flat.reshape(len(x_img), len(z_img)).T
               / (_adj_peak + 1e-30))
print(f'DAS adj_peak (MF, decimated) = {_adj_peak:.4e}')

# ── FISTA ─────────────────────────────────────────────────────────────────────
def _fista(K, d, eps_frac, adj_peak, n_iter):
    eps = eps_frac * adj_peak
    print(f'  eps={eps:.3e}  ({eps_frac} * {adj_peak:.3e})')
    m   = pylops.optimization.sparsity.fista(K, d, eps=eps, niter=n_iter, show=False)[0]
    return m.reshape(len(x_img), len(z_img)).T, eps

print(f'FISTA A  eps_frac={EPS_FRAC_A}...')
_img_fa, _eps_a = _fista(_K, _d_flat, EPS_FRAC_A, _adj_peak, N_ITER)
print(f'FISTA B  eps_frac={EPS_FRAC_B}...')
_img_fb, _eps_b = _fista(_K, _d_flat, EPS_FRAC_B, _adj_peak, N_ITER)

# ── Plot ──────────────────────────────────────────────────────────────────────
_ext  = [x_img[0], x_img[-1], z_img[-1], z_img[0]]
_xlim = (min(_xs1, _xs2) - 0.25, max(_xs1, _xs2) + 0.25)
_ylim = (z_scat + 0.20, z_scat - 0.20)
_ctr  = np.argmin(np.abs(x_traces - 2.0))   # central trace index

def _mark(ax):
    ax.axvline(_xs1, color='lime', ls='--', lw=0.9, alpha=0.8)
    ax.axvline(_xs2, color='lime', ls='--', lw=0.9, alpha=0.8)
    ax.plot(_xs1, z_scat, 'c*', ms=9, zorder=5)
    ax.plot(_xs2, z_scat, 'c*', ms=9, zorder=5)
    ax.set_xlim(*_xlim)
    ax.set_ylim(*_ylim)
    ax.set_xlabel('x [m]', fontsize=8)
    ax.set_ylabel('Depth [m]', fontsize=8)

fig, axs = plt.subplots(2, 3, figsize=(17, 9))

# Row 0: adjoint comparison
for ax, img, title in [
    (axs[0,0], _tr_r['img_kh'],
     'Kirchhoff DAS  (reference)'),
    (axs[0,1], _img_adj_mf,
     f'DAS adj on MF data  ({_dec}x dec  dt={_dt_dec:.3f} ns)'),
]:
    _p = np.percentile(np.abs(img), 99)
    ax.imshow(img, aspect='auto', cmap='seismic', vmin=-_p, vmax=_p,
              extent=_ext, origin='upper')
    ax.set_title(title, fontsize=8.5)
    _mark(ax)

# [0,2]: time-domain trace comparison (raw vs MF at central trace)
_t_raw = np.arange(_sc.shape[1]) * _dt
_raw_tr  = _sc[_ctr]
_mf_tr   = _mf[_ctr]
_mf_dec_tr = _mf_dec[_ctr]
_t_dec_tr  = _time_dec
ax02 = axs[0,2]
_norm = np.max(np.abs(_raw_tr)) + 1e-30
ax02.plot(_t_raw, _raw_tr / _norm, 'b', lw=0.7, alpha=0.7, label='raw (norm)')
ax02.plot(_t_raw, _mf_tr / (np.max(np.abs(_mf_tr)) + 1e-30),
          'r', lw=0.7, alpha=0.7, label='MF (norm)')
# TWTT expected for scatterer at depth z_scat along the central trace
_x_ctr = x_traces[_ctr]
_r_sc1 = np.sqrt((_x_ctr - _xs1)**2 + z_scat**2)
_r_sc2 = np.sqrt((_x_ctr - _xs2)**2 + z_scat**2)
_twtt1 = _r_sc1 / v_mig
_twtt2 = _r_sc2 / v_mig
ax02.axvline(_twtt1, color='lime', ls='--', lw=1.2, label=f'TWTT sc1={_twtt1:.2f} ns')
ax02.axvline(_twtt2, color='cyan', ls='--', lw=1.2, label=f'TWTT sc2={_twtt2:.2f} ns')
ax02.set_xlim(0, min(15, _t_raw[-1]))
ax02.set_xlabel('t [ns]', fontsize=8)
ax02.set_title(f'Central trace (x={_x_ctr:.2f} m)  raw vs MF', fontsize=8.5)
ax02.legend(fontsize=7, loc='upper right')
ax02.grid(True, alpha=0.3)

# Row 1: FISTA envelope
for ax, img, title in [
    (axs[1,0], _img_fa,
     f'FISTA  eps_frac={EPS_FRAC_A}  eps={_eps_a:.3e}  niter={N_ITER}'),
    (axs[1,1], _img_fb,
     f'FISTA  eps_frac={EPS_FRAC_B}  eps={_eps_b:.3e}  niter={N_ITER}'),
    (axs[1,2], _tr_r['img_lsm'],
     'LSM in rs_results  (re-run rs-load to update)'),
]:
    _env = hilbert_env(img)
    ax.imshow(_env, aspect='auto', cmap='hot_r', vmin=0, vmax=1,
              extent=_ext, origin='upper')
    ax.set_title(title, fontsize=8.5)
    _mark(ax)

fig.suptitle(
    f'LSM debug  |  factor={_tr_r["factor"]}lam  Dx={_tr_r["delta_m"]*1e3:.1f} mm  |  '
    f'aperture={ANGLE_AP}deg  |  dec={_dec}x  |  adj_peak={_adj_peak:.3e}',
    fontsize=10)
plt.tight_layout()
plt.show()
"""

with open('B_Scan_Migration_Playground.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

cells = nb['cells']

idx_rs  = next(i for i, c in enumerate(cells) if c.get('id') == 'rs-load')
idx_dbg = next(i for i, c in enumerate(cells) if c.get('id') == 'lsm-debug')

cells[idx_rs]['source']  = NEW_RS_LOAD
cells[idx_dbg]['source'] = NEW_LSM_DEBUG

with open('B_Scan_Migration_Playground.ipynb', 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f'Updated cell {idx_rs}  (rs-load)  and  cell {idx_dbg}  (lsm-debug)')
print('Done.')
