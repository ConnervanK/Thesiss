"""
Export one field-data (method, profile pair) for the Geophysical Letter figure notebook.

Runs ONLY the compute steps of TimeLapse_Notebooks/FieldData_Final_Pipeline.ipynb that
the cell-27 WLS workflow figure (build_wls_workflow_figure) needs, for just the two
profiles involved (+ the pre-injection reference), and saves a compact .npz that
New_Figures/Geophysical_Letter.ipynb plots from in the `gprmax` env:

    cell 3  -- survey constants                     cell 9  -- Gazdag migration (- reference)
    cell 5  -- raw MALA profiles, DEPTH_AXIS          cell 11 -- apply_post_imaging (amp gate)
    cell 8  -- 5-step conditioning chain             cell 15 -- difference = mon - base
    cells 18/21 -- saved napari k-space / B-scan picks (kspace_picks.pkl, bscan_picks.pkl)
    cell 25 -- _fit_inputs(..., 'kspace'): WLS fit on the manual k-space mask, wls_pow=3
    cell 16 -- _physical_dz sign correction, _kspace_display_crop (kx >= 0 half-plane)

Every profile is conditioned independently against the aligned reference, so
processing only (0, run_a, run_b) reproduces the pipeline's arrays for that pair.

Must run in the `gdp` conda env (gdp.data_io / gdp.preprocessing); no plotting happens
here, so the gdp env's broken matplotlib text rendering doesn't matter:

    C:/Users/Administrator/miniconda3/envs/gdp/python.exe New_Figures/export_fielddata_pair.py
"""
import os
import pickle
import sys
from pathlib import Path

os.environ.setdefault('MPLBACKEND', 'Agg')
os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

import numpy as np
from scipy.ndimage import gaussian_filter

REPO_ROOT = Path(__file__).resolve().parents[1]
NB_DIR = REPO_ROOT / 'TimeLapse_Notebooks'
sys.path.insert(0, str(NB_DIR))

from gdp.data_io import load_mala                                  # noqa: E402
from gdp.preprocessing.filtering import filter_data, remove_mean   # noqa: E402
from gdp.preprocessing.gain import linear_gain                     # noqa: E402
from gdp.preprocessing.image_processing import remove_svd          # noqa: E402
from gdp.preprocessing.trace_ops import align_traces               # noqa: E402
from helper_functions.migration import gazdag_migration, wls_phase_plane_fit  # noqa: E402

# ── What to export ─────────────────────────────────────────────────────────────────────
METHOD = 'gazdag'
RUN_A, RUN_B = 2, 3          # STAGE_PAIRS[0], the cell-27 example pair
KDISP_FAC = 3.0              # k-space crop |k| < KDISP_FAC * kz_cent (figure zooms further)
OUT = REPO_ROOT / 'New_Figures' / 'data' / f'fielddata_{METHOD}_{RUN_A}_{RUN_B}.npz'

# ── Cell 3: survey constants ───────────────────────────────────────────────────────────
DATA = NB_DIR / 'fielddata' / 'raw_data' / '060616'
OUT_DIR = NB_DIR / 'fielddata' / 'output'
_prof_name = lambda n: 'prof' if n == 0 else f'prof{n}'   # noqa: E731
V_PROP = 0.10        # m/ns
F0_MIG = 0.10        # GHz
DL = 0.05            # m, trace spacing
MAX_DEPTH = 85.0     # m
RAD_CUT = 300        # time samples kept for imaging

# ── Cell 16: per-method fit constants (Gazdag / Kirchhoff-BP branch) ───────────────────
KZ_CENT = 2 * np.pi * F0_MIG / V_PROP      # rad/m (full-velocity domain)
DZ_SIGN = -1                               # _physical_dz: + = downward for gazdag
WLS_POW = 3                                # WLS_FIT_DEFAULTS_BY_METHOD['gazdag']['wls_pow']


# ── Cell 5: raw profiles ───────────────────────────────────────────────────────────────
def load_raw_profile(run):
    return load_mala(str(DATA / _prof_name(run)), return_object=False)


reference_raw, ref_info = load_raw_profile(0)
SF = ref_info['frequency (GHz)']
N_TRACES = reference_raw.shape[1]
N_SAMPLES = reference_raw.shape[0]
T_AXIS = np.arange(1, N_SAMPLES + 1) / SF
DEPTH_AXIS = np.linspace(MAX_DEPTH, MAX_DEPTH - N_TRACES * DL, N_TRACES)
raw_profiles = {run: load_raw_profile(run)[0] for run in (RUN_A, RUN_B)}


# ── Cell 8: conditioning chain ─────────────────────────────────────────────────────────
def preprocess_profile(raw, align_to, *, align_reference=False, normalize_align=True):
    n = min(raw.shape[1], N_TRACES)
    x = filter_data(raw[:, :n], fq=(0.02, 0.20), sfreq=SF, btype='bandpass')
    x, _ = remove_mean(x, 299, 517)
    x, _, _ = align_traces(x, align_to[:, :n], upsample=5,
                           normalize=normalize_align, align_reference=align_reference)
    x, _ = remove_svd(x, low_s=0, high_s=1)
    x, _ = linear_gain(x, T_AXIS)
    return x[:RAD_CUT, :].T


_ref_bp = filter_data(reference_raw, fq=(0.02, 0.20), sfreq=SF, btype='bandpass')
_ref_dc, _ = remove_mean(_ref_bp, 299, 517)
ref_aligned, _, _ = align_traces(_ref_dc, _ref_dc, upsample=5, normalize=False, align_reference=True)
reference_bscan = preprocess_profile(reference_raw, ref_aligned, align_reference=True,
                                     normalize_align=False)
conditioned = {run: preprocess_profile(raw, ref_aligned) for run, raw in raw_profiles.items()}

# ── Cell 9: Gazdag migration minus the reference migration ─────────────────────────────
T_MIG = T_AXIS[:RAD_CUT]
X_IMG = np.linspace(0, V_PROP * T_MIG[-1] / 2, RAD_CUT)
DX_IMG = X_IMG[1] - X_IMG[0]


def migrate_gazdag(bscan, n_traces):
    x_bh = np.arange(n_traces) * DL
    return gazdag_migration(bscan.T, x_bh, T_MIG, X_IMG, V_PROP).T   # (n_depth, n_radial)


ref_gazdag = migrate_gazdag(reference_bscan, N_TRACES)
migrated = {}
for run, bscan in conditioned.items():
    n = bscan.shape[0]
    migrated[run] = migrate_gazdag(bscan, n) - ref_gazdag[:n, :]


# ── Cell 11: post-imaging amplitude gate (Gazdag: no radial taper) ─────────────────────
def apply_post_imaging(image, dz, dx, amp_gate_sigma_m=0.15, amp_gate_percentile=80):
    sigma_px = (amp_gate_sigma_m / dz, amp_gate_sigma_m / dx)
    envelope = gaussian_filter(np.abs(image), sigma_px)
    gate = np.clip(envelope / np.percentile(envelope, amp_gate_percentile), 0, 1)
    return image * gate


base = apply_post_imaging(migrated[RUN_A], dz=DL, dx=DX_IMG)
mon = apply_post_imaging(migrated[RUN_B], dz=DL, dx=DX_IMG)
diff = mon - base                                                     # cell 15

# ── Cells 18/21: saved napari picks for this (method, pair) ────────────────────────────
key = (METHOD, RUN_A, RUN_B)
print('loading k-space picks (large file) ...', flush=True)
with open(OUT_DIR / 'kspace_picks.pkl', 'rb') as f:
    kspace_mask = np.asarray(pickle.load(f)[key], dtype=bool)
with open(OUT_DIR / 'bscan_picks.pkl', 'rb') as f:
    bscan_mask = np.asarray(pickle.load(f)[key]['full_mask'], dtype=bool)

# ── Cell 25: WLS fit on the manual k-space mask (same padded grid it was painted on) ───
dz_raw, dx_est, phi_0, n_mask, diag = wls_phase_plane_fit(
    base, mon, DL, DX_IMG, KZ_CENT, taper='tukey', mask=kspace_mask, wls_pow=WLS_POW,
    return_diagnostics=True)
dz_phys = DZ_SIGN * dz_raw                                            # cell 16 _physical_dz

# ── Cell 16: _kspace_display_crop (kx >= 0 half-plane: the carrier/mirror axis is x) ───
kz_disp, kx_disp = np.fft.fftshift(diag['kz_ax']), np.fft.fftshift(diag['kx_ax'])
ix = np.where((kx_disp >= 0) & (kx_disp < KDISP_FAC * KZ_CENT))[0]
iz = np.where(np.abs(kz_disp) < KDISP_FAC * KZ_CENT)[0]
iz0, iz1, ix0, ix1 = iz.min(), iz.max() + 1, ix.min(), ix.max() + 1
crop = lambda a: np.fft.fftshift(a)[iz0:iz1, ix0:ix1]                 # noqa: E731

OUT.parent.mkdir(parents=True, exist_ok=True)
np.savez_compressed(
    OUT, base=base, mon=mon, diff=diff, depth_axis=DEPTH_AXIS[:base.shape[0]], x_img=X_IMG,
    bscan_mask=bscan_mask, phi=crop(diag['phi']), w=crop(diag['w']),
    kspace_mask=crop(kspace_mask), kz_ax=kz_disp[iz0:iz1], kx_ax=kx_disp[ix0:ix1],
    kz_cent=KZ_CENT, lam=V_PROP / F0_MIG, dz=DL, dx=DX_IMG,
    dz_raw=dz_raw, dz_phys=dz_phys, dx_est=dx_est, phi_0=phi_0, n_mask=n_mask,
    wls_pow=WLS_POW, method=METHOD, run_a=RUN_A, run_b=RUN_B)
print(f'{METHOD} {RUN_A}->{RUN_B}: dx_est = {dx_est:+.4f} m, dz (physical, +down) = '
      f'{dz_phys:+.4f} m, n_mask = {n_mask}  ->  {OUT}')
