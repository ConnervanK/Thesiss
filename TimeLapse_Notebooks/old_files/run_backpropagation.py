#!/usr/bin/env python3
"""
run_backpropagation.py
======================
FDTD physical time-reversal migration using gprMax.

For each study in the resolution ladder (2λ … 0.0625λ), this script:
  1. Loads the scattered B-scan  (full − background)
  2. Time-reverses each trace and applies a cosine fade-out taper
  3. Writes all 269 reversed traces to a single excitation.txt file
     consumed by gprMax via  #excitation_file
  4. Generates a backprop.in file with:
       - ε_r = 12.6  (= 4 × 3.15 → v_bp = v_ice/2, exploding-reflector equivalence)
       - 269 simultaneous Hertzian dipoles (one per original Rx position)
       - no PEC scatterers — empty homogeneous ice_bp medium
       - 30 wavefield snapshots over the last 3 ns (12 ns → 14.99 ns)
  5. Runs gprMax once per study (no -n flag; all sources fire simultaneously)

Physics note — why ε_r × 4?
  v = c / √ε_r.  To halve v: ε_r_new = 4 ε_r = 12.6.
  One-way travel at v/2 to z = 0.5 m ≈ 5.92 ns.
  Two-way travel at v_ice   to z = 0.5 m ≈ 5.92 ns.  ✓
  Energy therefore collapses to the scatterer positions at t ≈ T = 15 ns.

Directory layout produced
--------------------------
  resolution_study/
    study_*/
      backprop/
        excitation.txt       ← multi-column waveform file (269 reversed traces)
        backprop.in          ← gprMax input
        backprop_log.txt     ← gprMax stdout / stderr
        snap_000.h5 … snap_029.h5   ← wavefield snapshots written by gprMax

Usage
-----
  python run_backpropagation.py

Re-run safe: studies whose last snapshot (snap_029.h5) already exists are skipped.
"""

import copy
import logging
import os
import re
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

# ─── Physical constants ───────────────────────────────────────────────────────
C        = 0.299792458          # speed of light [m/ns]
EPS_R    = 3.15                 # relative permittivity of ice (forward sim)
EPS_R_BP = 4.0 * EPS_R         # = 12.6  (halves the wave velocity)
V_ICE    = C / EPS_R    ** 0.5  # ≈ 0.16891 m/ns
V_BP     = C / EPS_R_BP ** 0.5  # ≈ 0.08446 m/ns  (= v_ice / 2)

# ─── gprMax grid ─────────────────────────────────────────────────────────────
DX = 7.454e-4    # [m] grid cell size — 0.7454 mm, 30 cells/min-wavelength in ice
DY = DX
DZ = DX

# ─── Simulation domain ───────────────────────────────────────────────────────
DOMAIN_X    = 4.0        # [m]
DOMAIN_Y    = DY         # 1 cell thick (quasi-2D TMy)
DOMAIN_Z    = 0.8        # [m]
TIME_WINDOW = 15.0e-9    # [s]  = 15 ns

# ─── Zero-offset B-scan survey geometry ──────────────────────────────────────
N_CELLS_STEP  = 15
TRACE_SPACING = N_CELLS_STEP * DX     # 0.0111810 m  (15 cells per step)
N_TRACES      = 269

#  Tx is at x = snap(0.5) = 0.5001634 m  →  Rx is one cell further east
TX_X0 = round(0.5 / DX) * DX         # 0.5001634 m
RX_X0 = TX_X0 + DX                   # 0.5009088 m  (virtual source positions)
Y_SRC = DY / 2.0                     # 0.0003727 m
Z_SRC = DX                           # 0.0007454 m  (1 cell below surface)

# ─── Time-reversal tapering ───────────────────────────────────────────────────
TAPER_FRAC = 0.05    # cosine fade-out applied to the last 5 % of reversed traces

# ─── Snapshot configuration ───────────────────────────────────────────────────
N_SNAPS      = 30
SNAP_START_S = 12.0e-9     # [s]  start capturing at 12 ns (last 3 ns of 15 ns window)
SNAP_END_S   = 14.99e-9    # [s]  stop just before T to avoid off-by-one at final step
SNAP_MULT    = 4           # snapshot spatial coarsening factor (4× → 16× smaller files)
SNAP_DX      = SNAP_MULT * DX   # 0.0029816 m

# ─── File system paths ────────────────────────────────────────────────────────
BASE       = Path(r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks')
STUDY_ROOT = BASE / 'resolution_study'
BG_OUT     = BASE / 'B_Scan_Migration_data' / 'bg_bscan_merged.out'

GPRMAX_PYTHON = r'C:\Users\Administrator\miniconda3\envs\gprMax\python.exe'
MSVC_BIN = (
    r'C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools'
    r'\VC\Tools\MSVC\14.44.35207\bin\Hostx64\x64'
)

STUDY_DIRS = [
    'study_2lambda',
    'study_1lambda',
    'study_0p5lambda',
    'study_0p25lambda',
    'study_0p125lambda',
    'study_0p0625lambda',
]

# ─── Logging ──────────────────────────────────────────────────────────────────
STUDY_ROOT.mkdir(parents=True, exist_ok=True)
LOG_PATH = STUDY_ROOT / 'backprop_log.txt'
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOG_PATH, mode='a', encoding='utf-8'),
        logging.StreamHandler(sys.stdout),
    ],
)
log = logging.getLogger(__name__)


# ─── Helpers ──────────────────────────────────────────────────────────────────

def make_env():
    """Return os.environ copy with MSVC bin prepended so nvcc can find cl.exe."""
    env = copy.copy(os.environ)
    if MSVC_BIN not in env.get('PATH', ''):
        env['PATH'] = MSVC_BIN + ';' + env.get('PATH', '')
    env['PYTHONIOENCODING'] = 'utf-8'
    return env


def load_bscan(path: Path, component: str = 'Ey'):
    """
    Load a gprMax merged HDF5 B-scan.

    Returns
    -------
    bscan : ndarray, shape (N_traces, N_t)
    dt_s  : float  — time step in seconds
    """
    with h5py.File(path, 'r') as f:
        dt_s = float(f.attrs['dt'])
        rxs  = sorted(f['rxs'].keys(), key=lambda s: int(s.replace('rx', '')))
        bscan = np.stack([f['rxs'][r][component][:] for r in rxs], axis=0)
    return bscan, dt_s


def time_reverse_and_taper(bscan: np.ndarray, taper_frac: float = TAPER_FRAC):
    """
    Time-reverse and cosine-taper a B-scan.

    The cosine fade-out is applied to the LAST taper_frac of the reversed
    array.  In the original recording that corresponds to the early-time
    direct wave, which — when reversed — would fire as a giant simultaneous
    burst at t ≈ T and overwhelm the focused reflection energy.

    Parameters
    ----------
    bscan : ndarray, shape (N_traces, N_t)

    Returns
    -------
    rev : ndarray, shape (N_traces, N_t)  — time-reversed + tapered
    """
    rev = np.flip(bscan, axis=1).copy()       # reverse along time axis

    n_t   = rev.shape[1]
    n_tap = max(1, int(taper_frac * n_t))     # ≈ 426 samples ≈ 0.75 ns
    k     = np.arange(n_tap, dtype=float)
    # Hann-style fade: 1.0 at k=0, ~0 at k=n_tap−1
    taper = 0.5 * (1.0 + np.cos(np.pi * k / n_tap))
    rev[:, -n_tap:] *= taper[np.newaxis, :]   # apply in-place to end of each trace

    return rev


def write_excitation_file(bscan_rev: np.ndarray, dt_s: float, backprop_dir: Path) -> Path:
    """
    Write excitation.txt for gprMax's #excitation_file command.

    File format (space-separated):
      Row 0 :  time  waveform_000  waveform_001  …  waveform_268
      Rows 1+ : <time_s>  <amp_000>  <amp_001>  …  <amp_268>

    gprMax detects the leading 'time' token, uses column 0 as the time axis
    and builds a scipy.interpolate.interp1d for each subsequent column.

    Why N_EXTRA zero-padded samples?
      dt_bp (recomputed by gprMax via CFL) can differ from dt_fwd (read from
      the .out file) by a tiny floating-point amount.  If dt_bp > dt_fwd the
      simulation's last time point falls just above the excitation file's time
      range and interp1d raises ValueError (bounds_error=True by default).
      Extending by N_EXTRA * dt_s > TIME_WINDOW guarantees the file range
      always covers the full simulation window with no out-of-bounds access.

    Returns the path to the written file.
    """
    exc_path = backprop_dir / 'excitation.txt'
    n_traces, n_t = bscan_rev.shape

    N_EXTRA = 20   # zero-padded samples appended beyond the last real sample

    t_col  = np.arange(n_t + N_EXTRA, dtype=np.float64) * dt_s    # (N_t + N_EXTRA,)
    pad    = np.zeros((N_EXTRA, n_traces), dtype=bscan_rev.dtype)
    data   = np.hstack([
        t_col[:, np.newaxis],
        np.vstack([bscan_rev.T, pad]),
    ])   # shape: (N_t + N_EXTRA, 1 + N_TRACES)

    header = 'time ' + ' '.join(f'waveform_{i:03d}' for i in range(n_traces))

    log.info(
        f'  Writing excitation file: {exc_path}  '
        f'({n_t + N_EXTRA} rows × {n_traces + 1} cols)'
    )
    np.savetxt(exc_path, data, fmt='%.8e', header=header, comments='')
    return exc_path


def parse_scatterer_positions(in_path: Path):
    """
    Extract (xs1, xs2) from the title line of a gprMax .in file.
    Falls back to parsing the first two  #box … pec  lines.

    Returns
    -------
    xs1, xs2 : float  — x-coordinates [m] of the two PEC scatterers
    """
    content = in_path.read_text(encoding='utf-8')

    # Primary: title line contains  xs1=…m  xs2=…m
    m = re.search(r'xs1=([\d.]+)m\s+xs2=([\d.]+)m', content)
    if m:
        return float(m.group(1)), float(m.group(2))

    # Fallback: first two  #box: X … pec  lines give the west edge of each scatterer
    boxes = re.findall(
        r'#box:\s+([\d.]+)\s+0\.0\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+pec',
        content,
    )
    if len(boxes) >= 2:
        return float(boxes[0]), float(boxes[1])

    raise ValueError(f'Cannot parse scatterer positions from {in_path}')


def write_backprop_in(
    backprop_dir: Path,
    label: str,
    xs1: float,
    xs2: float,
) -> Path:
    """
    Write the gprMax backpropagation .in file.

    Key differences from the forward B-scan .in:
    - ε_r = 12.6 (4× forward) → wave speed halved → exploding-reflector equivalence
    - No PEC scatterers (empty homogeneous ice_bp medium)
    - No #rx, no #src_steps/#rx_steps (single simultaneous run)
    - 269 #hertzian_dipole sources at the original Rx positions, all driven by
      their respective reversed trace via #excitation_file
    - 30 #snapshot commands spanning [12 ns, 14.99 ns]
    """
    in_path = backprop_dir / 'backprop.in'

    snap_times = np.linspace(SNAP_START_S, SNAP_END_S, N_SNAPS)  # [s]

    lines = [
        f'#title: FDTD time-reversal backprop  {label}'
        f'  xs1={xs1:.6f}m  xs2={xs2:.6f}m  eps_r_bp={EPS_R_BP}',
        '',
        f'#domain: {DOMAIN_X} {DOMAIN_Y:.7f} {DOMAIN_Z}',
        f'#dx_dy_dz: {DX} {DY:.7f} {DZ}',
        f'#time_window: {TIME_WINDOW}',
        '',
        f'#material: {EPS_R_BP} 0.0 1.0 0.0 ice_bp',
        f'#box: 0.0 0.0 0.0 {DOMAIN_X} {DOMAIN_Y:.7f} {DOMAIN_Z} ice_bp',
        '',
        '#excitation_file: excitation.txt',
        '',
    ]

    for i in range(N_TRACES):
        x = RX_X0 + i * TRACE_SPACING
        lines.append(
            f'#hertzian_dipole: y {x:.7f} {Y_SRC:.7f} {Z_SRC:.7f} waveform_{i:03d}'
        )

    lines.append('')

    for k, t_s in enumerate(snap_times):
        lines.append(
            f'#snapshot: '
            f'0.0 0.0 0.0 '
            f'{DOMAIN_X} {DOMAIN_Y:.7f} {DOMAIN_Z} '
            f'{SNAP_DX:.7f} {DOMAIN_Y:.7f} {SNAP_DX:.7f} '
            f'{t_s:.6e} snap_{k:03d}'
        )

    in_path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    log.info(f'  Wrote .in file: {in_path}')
    return in_path


def run_gprmax_backprop(in_path: Path, label: str, log_path: Path) -> bool:
    """
    Run a single gprMax backpropagation simulation.

    No -n flag: all 269 dipoles fire simultaneously in one run.
    stdout/stderr are saved to log_path for post-run inspection.

    Returns True on success, False on failure/timeout.
    """
    log.info(f'[{label}] Starting gprMax backprop  →  {in_path}')
    try:
        result = subprocess.run(
            [GPRMAX_PYTHON, '-m', 'gprMax', str(in_path), '-gpu'],
            capture_output=True,
            encoding='utf-8',   # gprMax progress bars are UTF-8; Windows default cp1252 would crash
            errors='replace',   # replace any un-decodable bytes with ? rather than crashing
            env=make_env(),
            timeout=18000,
        )
    except subprocess.TimeoutExpired:
        log.error(f'[{label}] TIMEOUT after 18000 s')
        return False
    except Exception as exc:
        log.error(f'[{label}] Launch error: {exc}')
        return False

    # Save full gprMax output for inspection
    with open(log_path, 'w', encoding='utf-8') as fh:
        fh.write('=== STDOUT ===\n')
        fh.write(result.stdout or '(no output captured)')
        fh.write('\n=== STDERR ===\n')
        fh.write(result.stderr or '(no output captured)')

    if result.returncode != 0:
        log.error(
            f'[{label}] gprMax exited with code {result.returncode}. '
            f'See {log_path}'
        )
        log.error((result.stderr or '')[-2000:])
        return False

    log.info(f'[{label}] gprMax finished successfully')
    return True


# ─── Main pipeline ────────────────────────────────────────────────────────────

def main():
    log.info('=' * 72)
    log.info('run_backpropagation.py  —  FDTD time-reversal migration pipeline')
    log.info(f'  ε_r_bp      = {EPS_R_BP}  (= 4 × {EPS_R})')
    log.info(f'  v_bp        = {V_BP:.5f} m/ns  (= v_ice / 2)')
    log.info(f'  Snapshots   = {N_SNAPS} frames  over  '
             f'{SNAP_START_S*1e9:.0f} ns → {SNAP_END_S*1e9:.2f} ns')
    log.info('=' * 72)

    # Validate background file once before processing any study
    if not BG_OUT.exists():
        log.error(f'Background file not found: {BG_OUT}')
        log.error('Run the forward simulation / playground notebook first.')
        sys.exit(1)

    log.info(f'Loading background B-scan: {BG_OUT}')
    bscan_bg, _ = load_bscan(BG_OUT)
    log.info(f'  Background shape: {bscan_bg.shape}')

    for dir_name in STUDY_DIRS:
        study_dir    = STUDY_ROOT / dir_name
        psf_out      = study_dir  / 'psf_bscan_merged.out'
        psf_in       = study_dir  / 'psf_bscan.in'
        backprop_dir = study_dir  / 'backprop'
        skip_sentinel = backprop_dir / 'snap_029.h5'

        log.info('')
        log.info(f'─── {dir_name} {"─"*(50-len(dir_name))}')

        if not psf_out.exists():
            log.warning(f'  Forward B-scan missing: {psf_out} — skipping')
            continue

        if skip_sentinel.exists():
            log.info(f'  Last snapshot exists ({skip_sentinel.name}) — skipping')
            continue

        # ── Step 1: Load and scatter-field-subtract B-scan ───────────────────
        log.info(f'  Loading: {psf_out}')
        bscan_full, dt_s = load_bscan(psf_out)
        log.info(f'  B-scan shape: {bscan_full.shape}  dt = {dt_s:.4e} s')

        bscan_sc = bscan_full - bscan_bg   # scattered field only

        # ── Step 2: Time-reverse + taper ─────────────────────────────────────
        bscan_rev = time_reverse_and_taper(bscan_sc)
        log.info(f'  Time-reversed and tapered  (last {TAPER_FRAC*100:.0f}% cosine fade)')

        # ── Step 3: Write excitation file ─────────────────────────────────────
        backprop_dir.mkdir(exist_ok=True)
        write_excitation_file(bscan_rev, dt_s, backprop_dir)

        # ── Step 4: Parse scatterer positions and write .in file ──────────────
        xs1, xs2 = parse_scatterer_positions(psf_in)
        log.info(f'  Scatterers: xs1 = {xs1:.6f} m,  xs2 = {xs2:.6f} m')
        in_path = write_backprop_in(backprop_dir, dir_name, xs1, xs2)

        # ── Step 5: Run gprMax ────────────────────────────────────────────────
        gprmax_log = backprop_dir / 'gprmax_run.log'
        ok = run_gprmax_backprop(in_path, dir_name, gprmax_log)

        if not ok:
            log.error(f'  [{dir_name}] Backprop run FAILED — continuing with next study')

    log.info('')
    log.info('All studies processed.')
    log.info(
        'Visualise results by running the backpropagation cells in '
        'B_Scan_Migration_Playground.ipynb'
    )


if __name__ == '__main__':
    main()
