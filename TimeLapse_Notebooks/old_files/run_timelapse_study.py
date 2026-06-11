#!/usr/bin/env python3
"""
run_timelapse_study.py
======================
Automated Phase 2 Time-Lapse Sensitivity Study for two PEC point scatterers
in a homogeneous ice medium using gprMax zero-offset B-scan simulations.

Experimental Design
-------------------
  - Left Scatterer (Internal Control): Remains completely frozen at a baseline position.
  - Right Scatterer (Moving Target): Starts at a baseline position far away, then 
    shifts leftward by fractions of the subsurface wavelength toward the control.

Time-Lapse Shift Ladder
-----------------------
  2λ  →  1λ  →  0.5λ  →  0.25λ  →  0.125λ  →  0.0625λ   (λ = v_ice / f_c)

Directory layout produced
-------------------------
  TimeLapse_Notebooks/
    timelapse_study/                  ← Specific Phase 2 destination folder
      simulation_log.txt              ← Timestamped log of execution
      baseline/
          base_bscan.in
          base_bscan_merged.out       ← Core baseline data matrix
      shift_2lambda/
          shift_bscan.in
          shift_bscan_merged.out
          timelapse_data.npz          ← Contains raw baseline and timelapse arrays
      shift_1lambda/   ...
      shift_0p5lambda/ ...
      shift_0p25lambda/...
      shift_0p125lambda/...
      shift_0p0625lambda/...

Usage
-----
  python run_timelapse_study.py
"""

import copy
import os
import sys
import time
import logging
import subprocess
import numpy as np
import h5py
from pathlib import Path
from datetime import datetime

# ─── Physical constants ───────────────────────────────────────────────────────
C        = 0.299792458          # speed of light [m/ns]
EPS_R    = 3.15                 # relative permittivity of ice
V_ICE    = C / EPS_R ** 0.5    # ≈ 0.16891 m/ns
F_C_GHz  = 1.5                 # centre frequency [GHz]
LAMBDA   = V_ICE / F_C_GHz     # wavelength in ice ≈ 0.1126 m

# ─── gprMax grid ─────────────────────────────────────────────────────────────
DX = 7.454e-4    # grid cell size [m] = 0.7454 mm (30 cells/min-wavelength in ice)
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

# ─── Scatterer geometry (Phase 2 Time-Lapse Specific) ────────────────────────
Z_SCAT           = 0.5   # scatterer depth [m]
X_LEFT_FIXED     = 1.0   # Left scatterer position (STATIONARY CONTROL) [m]
X_RIGHT_BASELINE = 2.5   # Right scatterer baseline starting position [m]

# ─── Time-Lapse Shift Ladder ──────────────────────────────────────────────────
LAMBDA_FACTORS = [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625]  # Multiples of wavelength for right scatterer shifts

# ─── File system paths ───────────────────────────────────────────────────────
BASE          = Path(r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks')
STUDY_ROOT    = BASE / 'timelapse_study'    # Dedicated proper folder for Phase 2
BASELINE_DIR  = STUDY_ROOT / 'baseline'

GPRMAX_PYTHON = r'C:\Users\Administrator\miniconda3\envs\gprMax\python.exe'
MERGE_SCRIPT  = str(BASE / 'merge_bscan.py')
MSVC_BIN      = (
    r'C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools'
    r'\VC\Tools\MSVC\14.44.35207\bin\Hostx64\x64'
)

# ─── Logging Setup ───────────────────────────────────────────────────────────
STUDY_ROOT.mkdir(parents=True, exist_ok=True)
BASELINE_DIR.mkdir(parents=True, exist_ok=True)

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


def write_in_file(path, xl, xr, label_info=""):
    """Write a gprMax zero-offset B-scan input (.in) file."""
    tx_x0 = snap(X_START)
    rx_x0 = tx_x0 + DX
    src_z = snap(DX)
    zs    = snap(Z_SCAT)

    title = f'Phase 2 Time-Lapse Setup [{label_info}] - XL_Fixed={xl:.6f}m, XR={xr:.6f}m'

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
        f'#box: {xl:.7f} 0.0 {zs:.7f} {xl + DX:.7f} {DY:.7f} {zs + DX:.7f} pec',
        '',
        f'#box: {xr:.7f} 0.0 {zs:.7f} {xr + DX:.7f} {DY:.7f} {zs + DX:.7f} pec',
        '',
        f'#waveform: ricker 1.0 {F_C_GHz * 1e9:.1f} myricker',
        f'#hertzian_dipole: y {tx_x0:.7f} {DY / 2:.7f} {src_z:.7f} myricker',
        f'#rx: {rx_x0:.7f} {DY / 2:.7f} {src_z:.7f}',
        f'#src_steps: {TRACE_SPACING:.7f} 0.0 0.0',
        f'#rx_steps:  {TRACE_SPACING:.7f} 0.0 0.0',
    ]

    path.write_text('\n'.join(lines), encoding='utf-8')


def run_gprmax_and_merge(in_path, merged_path, label):
    """Run gprMax forward loop, then combine per-trace outputs."""
    if merged_path.exists():
        log.info(f'[{label}] SKIP  — {merged_path.name} already exists '
                 f'({merged_path.stat().st_size / 1e6:.1f} MB)')
        return True

    env = make_env()
    t_wall = time.perf_counter()

    log.info(f'[{label}] gprMax START  ({N_TRACES} traces)  →  {in_path}')
    try:
        result = subprocess.run(
            [GPRMAX_PYTHON, '-m', 'gprMax', str(in_path), '-n', str(N_TRACES), '-gpu'],
            capture_output=True, text=True, env=env, timeout=18000
        )
    except subprocess.TimeoutExpired:
        log.error(f'[{label}] TIMEOUT after 5 h — job abandoned')
        return False
    except Exception as exc:
        log.error(f'[{label}] Could not launch gprMax: {exc}')
        return False

    elapsed = time.perf_counter() - t_wall
    if result.returncode != 0:
        log.error(f'[{label}] gprMax returned code {result.returncode} after {elapsed / 60:.1f} min')
        log.error(f'[{label}] STDERR tail:\n{result.stderr[-1000:]}')
        return False

    log.info(f'[{label}] gprMax OK  ({elapsed / 60:.1f} min)')

    stem = str(in_path.parent / in_path.stem)
    log.info(f'[{label}] Merging {N_TRACES} files …')
    try:
        merge = subprocess.run(
            [GPRMAX_PYTHON, MERGE_SCRIPT, stem, str(N_TRACES), '--remove-files'],
            capture_output=True, text=True, env=env, timeout=600
        )
    except Exception as exc:
        log.error(f'[{label}] Merge launch failed: {exc}')
        return False

    if merge.returncode != 0:
        log.error(f'[{label}] Merge error (code {merge.returncode}): {merge.stderr[-500:]}')
        return False

    log.info(f'[{label}] Merged → {merged_path.name} ({merged_path.stat().st_size / 1e6:.1f} MB)')
    return True


def extract_bscan_data(merged_hdf5_path):
    """Extract the B-scan trace matrix from merged gprMax HDF5 output files."""
    with h5py.File(merged_hdf5_path, 'r') as f:
        nrx = int(f.attrs['nrx'])
        dt = float(f.attrs['dt'])
        n_samples = int(f.attrs['Iterations'])
        traces = [f[f'rxs/rx{i}/Ey'][:] for i in range(1, nrx + 1)]
        t = np.arange(n_samples) * dt
    data = np.column_stack(traces)   # shape: (n_samples, n_traces)
    return data, t


# ─── Build time-lapse shift queue ───────────────────────────────────────────

xl_snapped   = snap(X_LEFT_FIXED)
xr_base_snap = snap(X_RIGHT_BASELINE)

jobs = []
for factor in LAMBDA_FACTORS:
    shift_val       = factor * LAMBDA
    xr_shifted      = snap(X_RIGHT_BASELINE - shift_val)
    actual_shift    = xr_base_snap - xr_shifted

    factor_str = f'{factor:g}'.replace('.', 'p')
    dir_name   = f'shift_{factor_str}lambda'

    jobs.append({
        'factor':       factor,
        'dir_name':     dir_name,
        'label':        f'{factor:g}λ_Shift',
        'shift_target': shift_val,
        'shift_actual': actual_shift,
        'xl':           xl_snapped,
        'xr':           xr_shifted,
        'job_dir':      STUDY_ROOT / dir_name
    })


# ─── Screen Header and Geometry Table ────────────────────────────────────────

log.info('')
log.info('=' * 72)
log.info('PHASE 2 TIME-LAPSE STUDY  —  START')
log.info(f'  Output root : {STUDY_ROOT}')
log.info(f'  Wavelength  : {LAMBDA * 1e3:.3f} mm')
log.info(f'  Fixed Control (Left)  x: {xl_snapped:.4f} m, z: {snap(Z_SCAT):.4f} m')
log.info(f'  Baseline Target(Right) x: {xr_base_snap:.4f} m, z: {snap(Z_SCAT):.4f} m')
log.info('')
log.info('TIME-LAPSE DISPLACEMENT GEOMETRY TABLE')
log.info(f'  {"Label":<12}  {"Shift Target[mm]":>16}  {"Shift Actual[mm]":>16}  {"Target X2 [m]":>14}')
log.info('  ' + '-' * 64)
for j in jobs:
    log.info(f'  {j["label"]:<12}  {j["shift_target"] * 1e3:>16.3f}  {j["shift_actual"] * 1e3:>16.3f}  {j["xr"]:>14.6f}')
log.info('=' * 72)
log.info('')


# ─── STEP 1: Execute Core Baseline Model ────────────────────────────────────

base_in     = BASELINE_DIR / 'base_bscan.in'
base_merged = BASELINE_DIR / 'base_bscan_merged.out'

log.info('Executing Primary baseline simulation...')
write_in_file(base_in, xl=xl_snapped, xr=xr_base_snap, label_info="BASELINE")
ok = run_gprmax_and_merge(base_in, base_merged, 'BASELINE_RUN')
if not ok:
    log.critical('Baseline configuration simulation FAILED. Terminating process.')
    sys.exit(1)

# Extract baseline data matrix for storage
base_data, time_vector = extract_bscan_data(base_merged)
log.info(f'Successfully loaded Baseline matrix shape: {base_data.shape}')
log.info('')


# ─── STEP 2: Loop Through Time-Lapse Shifts ─────────────────────────────────

results = {}

for idx, j in enumerate(jobs):
    label    = j['label']
    job_dir  = j['job_dir']
    in_path  = job_dir / 'shift_bscan.in'
    merged   = job_dir / 'shift_bscan_merged.out'
    npz_path = job_dir / 'timelapse_data.npz'

    log.info('─' * 72)
    log.info(f'TIME-LAPSE SHIFT JOB {idx + 1}/{len(jobs)}: {label}')
    log.info(f'  Target Displacement = {j["shift_target"] * 1e3:.3f} mm')
    log.info(f'  Actual Displacement = {j["shift_actual"] * 1e3:.3f} mm')
    log.info(f'  New Right Scatterer Coordinate = {j["xr"]:.6f} m')

    t_job = time.perf_counter()

    try:
        job_dir.mkdir(parents=True, exist_ok=True)
        
        # Write input file
        write_in_file(in_path, xl=j['xl'], xr=j['xr'], label_info=label)
        
        # Execute forward model
        ok = run_gprmax_and_merge(in_path, merged, label)
        
        if ok:
            log.info(f'[{label}] Extracting time-lapse dataset...')
            tl_data, _ = extract_bscan_data(merged)
            
            # Calculate physical antenna positions along tracking profile
            antenna_x = np.linspace(X_START, X_END, N_TRACES)
            
            # Save raw datasets side-by-side (NO subtraction performed here)
            np.savez_compressed(
                npz_path,
                baseline=base_data,
                timelapse=tl_data,
                time=time_vector,
                antenna_x=antenna_x,
                xl_fixed=j['xl'],
                xr_moved=j['xr']
            )
            log.info(f'[{label}] Exported raw arrays to compressed archive → {npz_path.name}')
            
            elapsed = time.perf_counter() - t_job
            results[label] = 'OK'
        else:
            results[label] = 'FAIL'

    except Exception as exc:
        log.error(f'[{label}] Error encountered during process loop: {exc}', exc_info=True)
        results[label] = f'ERROR: {exc}'

    if idx < len(jobs) - 1:
        time.sleep(2)

# ─── Summary Report ──────────────────────────────────────────────────────────
log.info('')
log.info('=' * 72)
log.info('TIME-LAPSE RESOLUTION RUN COMPLETE')
log.info(f'Master Folder: {STUDY_ROOT}')
for lbl, status in results.items():
    icon = '✓' if status == 'OK' else '✗'
    log.info(f'  {icon}  {lbl:<15} : {status}')
log.info('=' * 72)