"""
Figure 1 (Geophysical Letter) -- the cross-spectrum phase method on synthetic data.

Two panels, one figure per movement type (lateral / vertical / diagonal):

    (a) Time-lapse difference image, monitor - baseline, in physical units [m].
    (b) Cross-spectrum phase arg(XS), XS = B . conj(M)  (thesis eq. cross-spectrum-def),
        with k-axes normalised by the central wavenumber k_zc = 2*pi/lambda. Bins whose
        energy |XS|^2 is below ENERGY_THR * max(|XS|^2) are drawn semi-transparent and
        the above-threshold lobes are outlined.

Data and processing mirror Hypothesis_1.ipynb exactly:
    * migrated images are read from the same cached npz archives as load_study()
      (cell 2), with the baseline = scenario 0 and the monitor = the chosen shift;
    * both images are cropped laterally to +-2.5 lambda around the scatterer apex,
      found with the same envelope logic as estimate_displacement() (cell 4);
    * XS and the WLS displacement estimate come from
      helper_functions.WLS.estimate_shift_2d, the same call the thesis results use.

Note on the k-space window: in a zero-offset migrated depth image the wavelet is
compressed by the two-way travel path, so the energy lobes peak at |k_z| ~ 1.8 k_zc
(not 1.0 k_zc). K_LIM is therefore set to 3.5 so the full lobes are visible.

Usage (run with the gprMax conda env, which has numpy/scipy/matplotlib/seaborn):

    python fig1_cross_spectrum_phase.py                        # all three, lambda/4, Gazdag
    python fig1_cross_spectrum_phase.py --movement vertical --shift 0.125 --method Kirchhoff

Output: New_Figures/figures/fig1_<movement>_<method>_<shift>.png  (dpi=300)
"""

import argparse
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import hilbert

# ── Repository paths ───────────────────────────────────────────────────────────────────
# This script lives in <repo>/New_Figures/; the cached synthetic data and the
# helper_functions package both live in <repo>/TimeLapse_Notebooks/.
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
NB_DIR = REPO_ROOT / 'TimeLapse_Notebooks'
OUT_DIR = SCRIPT_DIR / 'figures'

if str(NB_DIR) not in sys.path:
    sys.path.insert(0, str(NB_DIR))
from helper_functions.WLS import estimate_shift_2d  # noqa: E402

# ── Physical constants (identical to Hypothesis_1.ipynb, cell 2) ───────────────────────
EPS_R = 3.15                                # relative permittivity of pure ice
V_ICE = 0.299792458 / np.sqrt(EPS_R)       # wave speed in ice [m/ns]
F_C_GHZ = 1.5                              # Ricker centre frequency [GHz]
LAM = V_ICE / F_C_GHZ                      # central wavelength in ice [m] (~0.1126 m)
KZ_C = 2 * np.pi / LAM                     # central wavenumber [rad/m]
FDTD_CELL_M = 0.001                        # gprMax grid cell; true shifts snap to it

# ── Figure / processing parameters ─────────────────────────────────────────────────────
CROP_HW_LAM = 2.5      # lateral half-width of the analysis window [lambda] (as in H1)
SEARCH_LAM = 3.0       # half-width of the apex search window [lambda] (as in H1)
ENERGY_THR = 0.10      # opaque where |XS|^2 >= ENERGY_THR * max(|XS|^2)
ALPHA_LOW = 0.25       # alpha of the below-threshold bins
K_LIM = 3.5            # half-width of the plotted k-window [units of k_zc]
VIEW_Z_HW_LAM = CROP_HW_LAM  # panel (a) display only: depth half-width around the
                             # difference peak [lambda] (= lateral half-width -> square
                             # panel); the XS is still computed on the full depth column
FIG_SIZE = (7.2, 3.4)  # inches; ~ full text width of a two-column letter
DPI = 300

METHOD_KEYS = {'Kirchhoff': 'kirchhoff', 'Gazdag': 'gazdag', 'Back-prop': 'backprop'}

# Where each movement type's cache lives, and which npz keys carry its shift ladder
# (lambda units). Same files/keys as STUDIES in Hypothesis_1.ipynb. `primary` is the
# axis the --shift argument refers to.
STUDIES = {
    'lateral': dict(
        root=NB_DIR / 'timelapse_study', migrated='migrated_results.npz',
        diff='difference_migrated_results.npz',
        shift_x='separation_lambda', shift_z=None, primary='x'),
    'vertical': dict(
        root=NB_DIR / 'vertical_timelapse_study', migrated='vertical_migrated_results.npz',
        diff='vertical_difference_migrated_results.npz',
        shift_x=None, shift_z='shift_lambda', primary='z'),
    'diagonal': dict(
        root=NB_DIR / 'diagonal_timelapse_study', migrated='diagonal_migrated_results.npz',
        diff='diagonal_difference_migrated_results.npz',
        shift_x='shift_lambda_x', shift_z='shift_lambda_y', primary='x'),
}


def fdtd_true(nominal_m):
    """Round a nominal displacement [m] to the nearest 1 mm FDTD cell -- the shift
    gprMax actually realised (same as fdtd_true() in Hypothesis_1.ipynb)."""
    return np.round(np.asarray(nominal_m, dtype=float) / FDTD_CELL_M) * FDTD_CELL_M


def lambda_label(frac):
    """0.25 -> r'\\lambda/4', 1.0 -> r'\\lambda', 2.0 -> r'2\\lambda' (mathtext)."""
    f = Fraction(float(frac)).limit_denominator(64)
    if f == 0:
        return '0'
    if f.denominator == 1:
        return r'\lambda' if f.numerator == 1 else rf'{f.numerator}\lambda'
    num = '' if f.numerator == 1 else str(f.numerator)
    return rf'{num}\lambda/{f.denominator}'


def load_pair(movement, method, shift_lam):
    """Load the baseline and one monitor migrated image plus their grid axes.

    Returns a dict with base, mon (n_z, n_x), x [m], z [m], and the true dx/dz [m]
    and lambda fractions of the selected scenario.
    """
    cfg = STUDIES[movement]
    mig = np.load(cfg['root'] / cfg['migrated'])
    # diagonal_migrated_results.npz has no grid axes -> fall back to the diff archive
    diff = np.load(cfg['root'] / cfg['diff'])
    x = mig['x_traces'] if 'x_traces' in mig.files else diff['x_traces']
    z = mig['z_img'] if 'z_img' in mig.files else diff['z_img']

    n = len(mig['scenarios'])
    sx = mig[cfg['shift_x']] if cfg['shift_x'] else np.zeros(n)
    sz = mig[cfg['shift_z']] if cfg['shift_z'] else np.zeros(n)

    # Pick the scenario whose shift along the primary axis matches --shift
    ladder = sx if cfg['primary'] == 'x' else sz
    hits = np.flatnonzero(np.isclose(ladder, shift_lam))
    hits = hits[hits > 0]                      # index 0 is the baseline itself
    if hits.size == 0:
        raise ValueError(f'{movement}: no scenario with a {cfg["primary"]}-shift of '
                         f'{shift_lam} lambda; available: {sorted(set(ladder[1:]))}')
    i = int(hits[0])

    stack = mig[METHOD_KEYS[method]]
    return dict(base=np.nan_to_num(stack[0]), mon=np.nan_to_num(stack[i]), x=x, z=z,
                frac_x=float(sx[i]), frac_z=float(sz[i]),
                true_dx=float(fdtd_true(sx[i] * LAM)), true_dz=float(fdtd_true(sz[i] * LAM)))


def apex_crop_window(base, mon, x):
    """Lateral index range [ix_lo, ix_hi) of the +-CROP_HW_LAM*lambda analysis window.

    Port of estimate_displacement() in Hypothesis_1.ipynb: find the trace where the
    two envelopes differ most, then the baseline envelope peak within +-SEARCH_LAM
    lambda of it (the baseline apex), and centre the window there.
    """
    env_base = np.abs(hilbert(base, axis=0))
    env_mon = np.abs(hilbert(mon, axis=0))
    _, ix_d = np.unravel_index(np.argmax(np.abs(env_mon - env_base)), env_base.shape)

    lo_s, hi_s = np.searchsorted(x, [x[ix_d] - SEARCH_LAM * LAM, x[ix_d] + SEARCH_LAM * LAM])
    _, ix_loc = np.unravel_index(np.argmax(env_base[:, lo_s:hi_s]), env_base[:, lo_s:hi_s].shape)
    x_apex = x[lo_s + ix_loc]

    ix_lo, ix_hi = np.searchsorted(x, [x_apex - CROP_HW_LAM * LAM, x_apex + CROP_HW_LAM * LAM])
    return int(ix_lo), int(ix_hi)


def edge_extent(ax_h, ax_v):
    """imshow extent [left, right, bottom, top] whose pixel *edges* (not centres)
    bracket the sample axes, so every pixel sits exactly on its coordinate."""
    hh = 0.5 * (ax_h[1] - ax_h[0])
    hv = 0.5 * (ax_v[1] - ax_v[0])
    return [ax_h[0] - hh, ax_h[-1] + hh, ax_v[0] - hv, ax_v[-1] + hv]


def make_figure(movement, method, shift_lam):
    d = load_pair(movement, method, shift_lam)

    # ── 1. Crop both images to the analysis window (what the WLS fit actually sees) ──
    ix_lo, ix_hi = apex_crop_window(d['base'], d['mon'], d['x'])
    x_c, z = d['x'][ix_lo:ix_hi], d['z']
    base_c = d['base'][:, ix_lo:ix_hi]
    mon_c = d['mon'][:, ix_lo:ix_hi]
    dx_g, dz_g = float(x_c[1] - x_c[0]), float(z[1] - z[0])

    # ── 2. Panel (a) data: difference image, normalised by the baseline peak so the
    #       colour bar reads as a fraction of the reflector's own amplitude ──────────
    diff_img = (mon_c - base_c) / np.abs(base_c).max()
    vmax_diff = np.abs(diff_img).max()

    # ── 3. Panel (b) data: cross-spectrum XS = FFT2(base) . conj(FFT2(mon)) from the
    #       thesis helper (Tukey-tapered), plus its WLS displacement estimate ─────────
    dz_est, dx_est, _, XS, kz_ax, kx_ax = estimate_shift_2d(base_c, mon_c, dz_g, dx_g, KZ_C)

    # fftshift -> monotonic k-axes with DC in the middle; normalise by k_zc
    XS = np.fft.fftshift(XS)
    kz_n = np.fft.fftshift(kz_ax) / KZ_C
    kx_n = np.fft.fftshift(kx_ax) / KZ_C

    phase = np.angle(XS)                      # arg(XS) in (-pi, pi]
    energy = np.abs(XS) ** 2
    energy_rel = energy / energy.max()        # 0..1

    # 2-D alpha array: opaque where energy >= threshold, faded elsewhere
    alpha = np.where(energy_rel >= ENERGY_THR, 1.0, ALPHA_LOW)

    # ── 4. Figure ─────────────────────────────────────────────────────────────────────
    sns.set_theme(style='ticks', context='paper', font_scale=1.0)
    plt.rcParams.update({'savefig.dpi': DPI, 'image.interpolation': 'nearest'})

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=FIG_SIZE, layout='constrained')

    # (a) Difference image: depth increases downwards -> origin='upper' with the
    #     extent's bottom edge at the deepest sample.
    ext_a = edge_extent(x_c, z)
    ext_a = [ext_a[0], ext_a[1], ext_a[3], ext_a[2]]  # [xmin, xmax, zmax, zmin]
    im_a = ax_a.imshow(diff_img, extent=ext_a, origin='upper', aspect='equal',
                       cmap='RdBu_r', vmin=-vmax_diff, vmax=vmax_diff)
    # Zoom the view (not the data) onto the depth of the strongest time-lapse change
    # (window slides, rather than shrinks, where it would run off the image)
    z_pk = z[np.argmax(np.abs(hilbert(diff_img, axis=0)).max(axis=1))]
    z_hw = VIEW_Z_HW_LAM * LAM
    z_top = np.clip(z_pk - z_hw, ext_a[3], max(ext_a[2] - 2 * z_hw, ext_a[3]))
    ax_a.set_ylim(min(z_top + 2 * z_hw, ext_a[2]), z_top)
    ax_a.set_xlabel('Distance (m)')
    ax_a.set_ylabel('Depth (m)')
    ax_a.set_title('(a) Difference: monitor $-$ baseline', loc='left')
    cb_a = fig.colorbar(im_a, cax=ax_a.inset_axes([1.04, 0, 0.05, 1]))
    cb_a.set_label(r'$\Delta$amplitude / max$|$baseline$|$')

    # (b) Cross-spectrum phase: k_z increases upwards -> origin='lower'.
    im_b = ax_b.imshow(phase, extent=edge_extent(kx_n, kz_n), origin='lower',
                       aspect='equal', cmap='twilight', vmin=-np.pi, vmax=np.pi,
                       alpha=alpha)
    # Outline the above-threshold lobes (contour uses pixel-centre coordinates)
    ax_b.contour(kx_n, kz_n, energy_rel, levels=[ENERGY_THR], colors='k', linewidths=0.6)
    ax_b.set_xlim(-K_LIM, K_LIM)
    ax_b.set_ylim(-K_LIM, K_LIM)
    ax_b.set_xlabel(r'$k_x \,/\, k_{z,c}$')
    ax_b.set_ylabel(r'$k_z \,/\, k_{z,c}$')
    ax_b.set_title(r'(b) Cross-spectrum phase $\arg(XS)$', loc='left')
    cb_b = fig.colorbar(im_b, cax=ax_b.inset_axes([1.04, 0, 0.05, 1]),
                        ticks=[-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi])
    cb_b.ax.set_yticklabels([r'$-\pi$', r'$-\pi/2$', '0', r'$\pi/2$', r'$\pi$'])
    cb_b.set_label('Phase (rad)')
    cb_b.solids.set_alpha(1.0)  # colour bar shows the opaque colours

    # True vs WLS-estimated displacement, in the (empty) lower-left corner of (b)
    def _mm(v):
        return round(v * 1e3, 1) + 0.0        # + 0.0 turns -0.0 into 0.0

    shift_txt = (rf'$\Delta x$: ${lambda_label(d["frac_x"])}$ = {d["true_dx"] * 1e3:.0f} mm, '
                 rf'est. {_mm(dx_est):.1f} mm' + '\n' +
                 rf'$\Delta z$: ${lambda_label(d["frac_z"])}$ = {d["true_dz"] * 1e3:.0f} mm, '
                 rf'est. {_mm(dz_est):.1f} mm')
    ax_b.text(0.03, 0.03, shift_txt, transform=ax_b.transAxes, fontsize=6.5,
              va='bottom', ha='left',
              bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='0.7', lw=0.5))

    # ── 5. Save ───────────────────────────────────────────────────────────────────────
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    shift_tag = str(Fraction(shift_lam).limit_denominator(64)).replace('/', 'over')
    out = OUT_DIR / f'fig1_{movement}_{METHOD_KEYS[method]}_{shift_tag}lambda.png'
    fig.savefig(out, dpi=DPI, bbox_inches='tight')
    plt.close(fig)

    print(f'{movement:8s} {method:9s} shift={shift_lam:g} lambda | '
          f'dx true {d["true_dx"] * 1e3:6.1f} est {dx_est * 1e3:6.2f} mm | '
          f'dz true {d["true_dz"] * 1e3:6.1f} est {dz_est * 1e3:6.2f} mm -> {out.name}')
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--movement', choices=[*STUDIES, 'all'], default='all')
    p.add_argument('--method', choices=list(METHOD_KEYS), default='Gazdag')
    p.add_argument('--shift', type=float, default=0.25,
                   help='shift along the primary axis in lambda (x for lateral/diagonal, '
                        'z for vertical); e.g. 0.5, 0.25, 0.125, 0.0625, 0.03125')
    args = p.parse_args()

    movements = list(STUDIES) if args.movement == 'all' else [args.movement]
    for mv in movements:
        make_figure(mv, args.method, args.shift)


if __name__ == '__main__':
    main()
