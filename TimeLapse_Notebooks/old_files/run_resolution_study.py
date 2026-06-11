#!/usr/bin/env python3
"""
run_resolution_study.py
=======================
Automated Rayleigh resolution sensitivity study for two PEC point scatterers
in a homogeneous ice medium, using gprMax zero-offset B-scan simulations.

Resolution ladder
-----------------
  2λ  →  1λ  →  0.5λ  →  0.25λ  →  0.125λ  →  0.0625λ   (λ = v_ice / f_c)

Directory layout produced
-------------------------
  TimeLapse_Notebooks/
    resolution_study/
      simulation_log.txt          ← timestamped log of every action
      background/
          bg_bscan.in
          bg_bscan_merged.out     ← generated once; reused by all jobs
      study_2lambda/
          psf_bscan.in
          psf_bscan_merged.out
      study_1lambda/   ...
      study_0p5lambda/ ...
      study_0p25lambda/...
      study_0p125lambda/...
      study_0p0625lambda/...

Usage
-----
  python run_resolution_study.py

The script is safe to interrupt and re-run: jobs whose merged .out file
already exists are skipped automatically.
"""

import copy
import os
import sys
import time
import logging
import subprocess
from pathlib import Path
from datetime import datetime

# ─── Physical constants ───────────────────────────────────────────────────────
C        = 0.299792458          # speed of light [m/ns]
EPS_R    = 3.15                 # relative permittivity of ice
V_ICE    = C / EPS_R ** 0.5    # ≈ 0.16891 m/ns
F_C_GHz  = 1.5                 # centre frequency [GHz]
LAMBDA   = V_ICE / F_C_GHz     # wavelength in ice ≈ 0.1126 m

# ─── gprMax grid ─────────────────────────────────────────────────────────────
DX = 7.454e-4    # grid cell size [m] = 0.7454 mm  (30 cells/min-wavelength in ice)
DY = DX          # 2-D TMy simulation: one cell thick in y

# ─── Simulation domain ───────────────────────────────────────────────────────
DOMAIN_X    = 4.0       # [m]
DOMAIN_Z    = 0.8       # [m]
TIME_WINDOW = 15e-9     # [s]

# ─── Zero-offset B-scan survey ───────────────────────────────────────────────
X_START       = 0.5     # first Tx/Rx midpoint [m]
X_END         = 3.5     # last  Tx/Rx midpoint [m]
N_CELLS_STEP  = 15      # grid cells between successive traces
TRACE_SPACING = N_CELLS_STEP * DX                                    # ≈ 0.01118 m
N_TRACES      = int(round((X_END - X_START) / TRACE_SPACING)) + 1   # 269

# ─── Scatterer geometry ──────────────────────────────────────────────────────
Z_SCAT   = 0.5              # scatterer depth [m] — both on the same horizontal plane
X_CENTER = DOMAIN_X / 2.0  # 2.0 m — domain centre = survey centre

# ─── Resolution ladder ───────────────────────────────────────────────────────
LAMBDA_FACTORS = [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625]

# ─── File system paths ───────────────────────────────────────────────────────
BASE         = Path(r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks')
STUDY_ROOT   = BASE / 'resolution_study'
BG_DIR       = STUDY_ROOT / 'background'

# Reuse the playground background if it already exists (identical geometry)
PLAYGROUND_BG = BASE / 'B_Scan_Migration_data' / 'bg_bscan_merged.out'

GPRMAX_PYTHON = r'C:\Users\Administrator\miniconda3\envs\gprMax\python.exe'
MERGE_SCRIPT  = str(BASE / 'merge_bscan.py')
MSVC_BIN      = (
    r'C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools'
    r'\VC\Tools\MSVC\14.44.35207\bin\Hostx64\x64'
)

# ─── Logging ─────────────────────────────────────────────────────────────────
STUDY_ROOT.mkdir(parents=True, exist_ok=True)
BG_DIR.mkdir(parents=True, exist_ok=True)

LOG_PATH = STUDY_ROOT / 'simulation_log.txt'
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


# ─── Helpers ─────────────────────────────────────────────────────────────────

def snap(val):
    """Snap a coordinate to the nearest gprMax grid-cell boundary."""
    return round(val / DX) * DX


def make_env():
    """Return os.environ with MSVC bin prepended so nvcc can find cl.exe."""
    env = copy.copy(os.environ)
    if MSVC_BIN not in env.get('PATH', ''):
        env['PATH'] = MSVC_BIN + ';' + env.get('PATH', '')
    env['PYTHONIOENCODING'] = 'utf-8'
    return env


def write_in_file(path, xs1=None, xs2=None):
    """
    Write a gprMax zero-offset B-scan .in file.

    xs1, xs2 : snapped x-coordinates [m] of the two PEC scatterers.
               Pass None for both to write a background (empty ice) file.
    """
    tx_x0 = snap(X_START)          # Tx starts at first survey position
    rx_x0 = tx_x0 + DX             # Rx one cell to the right (≈ zero offset)
    src_z = snap(DX)               # source / receiver one cell below surface

    has_scatterers = xs1 is not None

    if has_scatterers:
        delta_mm = abs(xs2 - xs1) * 1e3
        title = f'PSF B-scan  delta={delta_mm:.3f}mm  xs1={xs1:.6f}m  xs2={xs2:.6f}m'
    else:
        title = 'Background B-scan  (empty ice, no scatterers)'

    lines = [
        f'#title: {title}',
        '',
        f'#domain: {DOMAIN_X} {DY:.7f} {DOMAIN_Z}',
        f'#dx_dy_dz: {DX} {DY:.7f} {DX}',
        f'#time_window: {TIME_WINDOW}',
        '',
        '#material: 3.15 0.0 1.0 0.0 ice',
        f'#box: 0.0 0.0 0.0 {DOMAIN_X} {DY:.7f} {DOMAIN_Z} ice',
        '',
    ]

    if has_scatterers:
        zs = snap(Z_SCAT)
        for xs in [xs1, xs2]:
            lines += [
                f'#box: {xs:.7f} 0.0 {zs:.7f} {xs + DX:.7f} {DY:.7f} {zs + DX:.7f} pec',
                '',
            ]

    lines += [
        f'#waveform: ricker 1.0 {F_C_GHz * 1e9:.1f} myricker',
        f'#hertzian_dipole: y {tx_x0:.7f} {DY / 2:.7f} {src_z:.7f} myricker',
        f'#rx: {rx_x0:.7f} {DY / 2:.7f} {src_z:.7f}',
        f'#src_steps: {TRACE_SPACING:.7f} 0.0 0.0',
        f'#rx_steps:  {TRACE_SPACING:.7f} 0.0 0.0',
    ]

    path.write_text('\n'.join(lines), encoding='utf-8')


def run_gprmax_and_merge(in_path, merged_path, label):
    """
    Run gprMax for N_TRACES traces (GPU), then merge per-trace output files.

    Returns True on success, False on any failure.
    Skips silently if merged_path already exists.
    """
    if merged_path.exists():
        log.info(f'[{label}] SKIP  — {merged_path.name} already exists '
                 f'({merged_path.stat().st_size / 1e6:.1f} MB)')
        return True

    env = make_env()
    t_wall = time.perf_counter()

    # ── gprMax simulation ─────────────────────────────────────────────────
    log.info(f'[{label}] gprMax START  ({N_TRACES} traces)  →  {in_path}')
    try:
        result = subprocess.run(
            [GPRMAX_PYTHON, '-m', 'gprMax', str(in_path),
             '-n', str(N_TRACES), '-gpu'],
            capture_output=True,
            text=True,
            env=env,
            timeout=18000,   # 5-hour hard cap per job
        )
    except subprocess.TimeoutExpired:
        log.error(f'[{label}] TIMEOUT after 5 h — job abandoned')
        return False
    except Exception as exc:
        log.error(f'[{label}] Could not launch gprMax: {exc}')
        return False

    elapsed = time.perf_counter() - t_wall
    if result.returncode != 0:
        log.error(f'[{label}] gprMax returned code {result.returncode} '
                  f'after {elapsed / 60:.1f} min')
        log.error(f'[{label}] STDERR tail:\n{result.stderr[-1000:]}')
        return False

    log.info(f'[{label}] gprMax OK  ({elapsed / 60:.1f} min)')

    # ── Merge per-trace .out files ────────────────────────────────────────
    stem = str(in_path.parent / in_path.stem)
    log.info(f'[{label}] Merging {N_TRACES} output files …')
    try:
        merge = subprocess.run(
            [GPRMAX_PYTHON, MERGE_SCRIPT, stem, str(N_TRACES), '--remove-files'],
            capture_output=True,
            text=True,
            env=env,
            timeout=600,
        )
    except Exception as exc:
        log.error(f'[{label}] Merge launch failed: {exc}')
        return False

    if merge.returncode != 0:
        log.error(f'[{label}] Merge error (code {merge.returncode}): '
                  f'{merge.stderr[-500:]}')
        return False

    if not merged_path.exists():
        log.error(f'[{label}] Merged file not found after merge: {merged_path}')
        return False

    log.info(f'[{label}] Merged → {merged_path.name}  '
             f'({merged_path.stat().st_size / 1e6:.1f} MB)')
    return True


# ─── Build job queue ──────────────────────────────────────────────────────────

jobs = []
for factor in LAMBDA_FACTORS:
    delta        = factor * LAMBDA               # target separation [m]
    xs1          = snap(X_CENTER - delta / 2)    # left scatterer (grid-snapped)
    xs2          = snap(X_CENTER + delta / 2)    # right scatterer (grid-snapped)
    delta_actual = xs2 - xs1                     # realised separation after snapping

    # Directory name: 2.0 → "study_2lambda", 0.5 → "study_0p5lambda"
    factor_str = f'{factor:g}'.replace('.', 'p')
    dir_name   = f'study_{factor_str}lambda'

    jobs.append({
        'factor':        factor,
        'dir_name':      dir_name,
        'label':         f'{factor:g}λ',
        'delta_target':  delta,
        'delta_actual':  delta_actual,
        'xs1':           xs1,
        'xs2':           xs2,
        'zs':            snap(Z_SCAT),
        'job_dir':       STUDY_ROOT / dir_name,
    })


# ─── Header / geometry table ──────────────────────────────────────────────────

log.info('')
log.info('=' * 72)
log.info('RESOLUTION STUDY  —  START')
log.info(f'  Script     : {Path(__file__).resolve()}')
log.info(f'  Study root : {STUDY_ROOT}')
log.info(f'  Log file   : {LOG_PATH}')
log.info(f'  Launched   : {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
log.info('')
log.info('Medium & survey parameters')
log.info(f'  ε_r         = {EPS_R}')
log.info(f'  v_ice       = {V_ICE:.5f} m/ns')
log.info(f'  f_c         = {F_C_GHz} GHz')
log.info(f'  λ           = {LAMBDA * 1e3:.3f} mm')
log.info(f'  Domain      = {DOMAIN_X} m × {DOMAIN_Z} m,  t_window = {TIME_WINDOW * 1e9:.1f} ns')
log.info(f'  Survey      = {X_START} – {X_END} m,  {N_TRACES} traces,  '
         f'step = {TRACE_SPACING * 1e3:.4f} mm')
log.info(f'  z_scat      = {Z_SCAT} m')
log.info(f'  x_center    = {X_CENTER} m  (domain centre = survey centre)')
log.info('')
log.info('GEOMETRY VERIFICATION TABLE')
log.info(f'  {"Label":<10}  {"Δ target [mm]":>14}  {"Δ actual [mm]":>14}  '
         f'{"Δ cells":>8}  {"x1 [m]":>10}  {"x2 [m]":>10}  {"z [m]":>7}')
log.info('  ' + '-' * 68)
for j in jobs:
    n_cells = round(j['delta_actual'] / DX)
    log.info(
        f'  {j["label"]:<10}  {j["delta_target"] * 1e3:>14.3f}  '
        f'{j["delta_actual"] * 1e3:>14.3f}  '
        f'{n_cells:>8d}  '
        f'{j["xs1"]:>10.6f}  {j["xs2"]:>10.6f}  {j["zs"]:>7.4f}'
    )
log.info('=' * 72)
log.info('')


# ─── Background scan (shared by all jobs) ─────────────────────────────────────

bg_merged = BG_DIR / 'bg_bscan_merged.out'

if PLAYGROUND_BG.exists():
    # The playground already ran an identical background — use a symlink or note
    log.info(f'BACKGROUND  Reusing playground background:')
    log.info(f'  {PLAYGROUND_BG}  ({PLAYGROUND_BG.stat().st_size / 1e6:.1f} MB)')
    bg_merged = PLAYGROUND_BG   # point to the existing file directly
else:
    log.info('BACKGROUND  Playground background not found — generating now.')
    bg_in = BG_DIR / 'bg_bscan.in'
    write_in_file(bg_in, xs1=None, xs2=None)
    log.info(f'  .in written → {bg_in}')
    ok = run_gprmax_and_merge(bg_in, bg_merged, 'BACKGROUND')
    if not ok:
        log.critical('Background simulation FAILED — cannot continue without it.')
        sys.exit(1)

log.info('')


# ─── Sequential job runner ────────────────────────────────────────────────────

results = {}   # label → 'OK' | 'SKIP' | 'FAIL' | error string

for idx, j in enumerate(jobs):
    label    = j['label']
    job_dir  = j['job_dir']
    in_path  = job_dir / 'psf_bscan.in'
    merged   = job_dir / 'psf_bscan_merged.out'

    log.info('─' * 72)
    log.info(f'JOB {idx + 1}/{len(jobs)}  START  {label}')
    log.info(f'  Target Δ  = {j["delta_target"] * 1e3:.3f} mm  '
             f'({j["factor"]:g} × λ = {j["factor"]:g} × {LAMBDA * 1e3:.3f} mm)')
    log.info(f'  Actual Δ  = {j["delta_actual"] * 1e3:.3f} mm  '
             f'({round(j["delta_actual"] / DX)} grid cells)')
    log.info(f'  (x1, z)   = ({j["xs1"]:.6f} m, {j["zs"]:.4f} m)')
    log.info(f'  (x2, z)   = ({j["xs2"]:.6f} m, {j["zs"]:.4f} m)')
    log.info(f'  Directory : {job_dir}')

    t_job = time.perf_counter()

    try:
        job_dir.mkdir(parents=True, exist_ok=True)

        # Write gprMax .in file
        write_in_file(in_path, xs1=j['xs1'], xs2=j['xs2'])
        log.info(f'  .in written → {in_path.name}')

        # Run gprMax + merge
        ok = run_gprmax_and_merge(in_path, merged, label)

        elapsed = time.perf_counter() - t_job
        if ok:
            size_mb = merged.stat().st_size / 1e6 if merged.exists() else 0
            log.info(f'JOB {idx + 1}/{len(jobs)}  OK    {label}  '
                     f'({elapsed / 60:.1f} min, {size_mb:.1f} MB)')
            results[label] = 'OK' if elapsed > 1 else 'SKIP (already existed)'
        else:
            log.error(f'JOB {idx + 1}/{len(jobs)}  FAIL  {label}  '
                      f'({elapsed / 60:.1f} min)')
            results[label] = 'FAIL'

    except Exception as exc:
        elapsed = time.perf_counter() - t_job
        log.error(f'JOB {idx + 1}/{len(jobs)}  ERROR {label}: {exc}  '
                  f'({elapsed / 60:.1f} min)', exc_info=True)
        results[label] = f'ERROR: {exc}'

    # Stabilisation pause between jobs (skip after last job)
    if idx < len(jobs) - 1:
        log.info('  Stabilisation pause 5 s …')
        time.sleep(5)


# ─── Final summary ────────────────────────────────────────────────────────────

log.info('')
log.info('=' * 72)
log.info('STUDY COMPLETE')
log.info(f'  Finished : {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
log.info(f'  Log file : {LOG_PATH}')
log.info('')
log.info('  Job results:')
all_ok = True
for lbl, status in results.items():
    icon = '✓' if 'OK' in status or 'SKIP' in status else '✗'
    log.info(f'    {icon}  {lbl:<12}  {status}')
    if 'OK' not in status and 'SKIP' not in status:
        all_ok = False

log.info('')
log.info('  Output files:')
log.info(f'    Background : {bg_merged}')
for j in jobs:
    merged = j['job_dir'] / 'psf_bscan_merged.out'
    exists = '✓' if merged.exists() else '✗ MISSING'
    size   = f'{merged.stat().st_size / 1e6:.1f} MB' if merged.exists() else ''
    log.info(f'    {exists}  {j["label"]:<12}  {merged}  {size}')

log.info('')
if all_ok:
    log.info('  All jobs completed successfully.')
else:
    log.info('  WARNING: one or more jobs failed — check the log above.')
log.info('=' * 72)
