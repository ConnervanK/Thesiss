import os
import glob
import shutil
import numpy as np
import stat
import time
import sys
from fractions import Fraction

def _rmtree_onerror(func, path, exc_info):
    try:
        os.chmod(path, stat.S_IWRITE)
        func(path)
    except Exception:
        pass

def safe_remove_dir(path, retries=3, delay_s=0.4):
    if not os.path.isdir(path):
        return True
    for attempt in range(1, retries + 1):
        try:
            shutil.rmtree(path, onerror=_rmtree_onerror)
            return True
        except PermissionError as e:
            if attempt < retries:
                time.sleep(delay_s)
                continue
            return False
        except Exception:
            return False
    return False

def clear_results(directory="."):
    for ext in ["*.vti", "*.in", "*.out", "*.png"]:
        for file in glob.glob(os.path.join(directory, ext)):
            try:
                os.remove(file)
            except Exception:
                pass
    for folder in glob.glob(os.path.join(directory, "*_snaps*")):
        if os.path.isdir(folder):
            try:
                safe_remove_dir(folder)
            except Exception:
                pass


def create_input_file(
        f_central=1.5 * 1e9,
        c=3 * 1e8,
        permittivity_ice=3.15,
        permittivity_air=1,
        permittivity_fracture=80,
        conductivity_ice=1e-6,
        conductivity_air=0,
        conductivity_fracture=0.01,
        fracture_depth=0.6,
        depth_below_fracture=0.1,
        air_thickness=0.1,
        rx_count=50,
        rx_spread=2.0,
        block_size=None,
        block_size_fraction=(1/4),
        time_lapse_shift=0.0
    ):
    print(f"Generating gprMax input file (mode: time-lapse)...")
    
    c_ice = c / np.sqrt(permittivity_ice)
    wavelength_ice = c_ice / f_central

    permittivity_fracture_plus = 1.0
    permittivity_fracture_minus = 80

    conductivity_fracture_plus = 0
    conductivity_fracture_minus = 0.01

    if block_size is None:
        block_size = block_size_fraction * wavelength_ice

    frac = block_size / wavelength_ice
    fr = Fraction(frac).limit_denominator(32)
    if fr.numerator == 0:
        block_label = "0Lambda"
        block_label_print = "0Lambda"
    else:
        block_label = f"{fr.numerator}_{fr.denominator}Lambda"
        block_label_print = f"{fr.numerator}/{fr.denominator}Lambda"

    wavelength_plus = (c / np.sqrt(permittivity_fracture_plus)) / f_central
    wavelength_minus = (c / np.sqrt(permittivity_fracture_minus)) / f_central
    min_wavelength = min(wavelength_plus, wavelength_minus)
    
    thickness_fracture = min_wavelength / 5

    dx_min = min_wavelength / 30
    dx_dy_dz = dx_min

    # Fixed geometry: domain always 4 m wide, source at x=2 m
    domain_width = 4.0
    domain_height = air_thickness + fracture_depth + thickness_fracture + depth_below_fracture

    domain_width = np.ceil(domain_width / dx_dy_dz) * dx_dy_dz
    domain_height = np.ceil(domain_height / dx_dy_dz) * dx_dy_dz

    TW = (5 / f_central) + (3 * (domain_height - air_thickness) / c_ice)
    n_snapshots = 36
    snapshot_time = 0.25e-9

    fracture_top = air_thickness + fracture_depth
    fracture_bottom = fracture_top + thickness_fracture

    y_air_bottom = domain_height - air_thickness
    y_frac_top_gpr = domain_height - fracture_top
    y_frac_bottom_gpr = domain_height - fracture_bottom
    y_frac_min = min(y_frac_top_gpr, y_frac_bottom_gpr)
    y_frac_max = max(y_frac_top_gpr, y_frac_bottom_gpr)

    n_blocks = int(np.round(domain_width / block_size))
    if n_blocks < 1: n_blocks = 1
    if n_blocks % 2 != 0: n_blocks += 1
    block_width = domain_width / n_blocks

    # Fixed transmitter position at center
    tx_start_x = 2.0
    rx_y = y_air_bottom
    
    rx_min_x = tx_start_x - (rx_spread / 2.0)
    rx_max_x = tx_start_x + (rx_spread / 2.0)
    rx_spacing = rx_spread / max(1, rx_count - 1)

    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_file_base = os.path.join(base_dir, 'configs', 'horizontal_scattering_baseline.in')

    with open(out_file_base, 'w') as f:
        f.write('#title: Horizontal Scattering Baseline')
        f.write('\n#domain: {} {} {}'.format(domain_width, domain_height, dx_dy_dz))
        f.write('\n#dx_dy_dz: {} {} {}'.format(dx_dy_dz, dx_dy_dz, dx_dy_dz))
        f.write('\n#time_window: {}'.format(TW))

        f.write('\n#material: {} {} 1 0 ice'.format(permittivity_ice, conductivity_ice))
        f.write('\n#material: {} {} 1 0 air'.format(permittivity_air, conductivity_air))
        f.write('\n#material: {} {} 1 0 fracture_plus'.format(permittivity_fracture_plus, conductivity_fracture_plus))
        f.write('\n#material: {} {} 1 0 fracture_minus'.format(permittivity_fracture_minus, conductivity_fracture_minus))

        f.write('\n#waveform: ricker 1 {} my_ricker'.format(f_central))
        f.write('\n#hertzian_dipole: z {} {} 0 my_ricker'.format(tx_start_x, rx_y))
        
        for i in range(rx_count):
            f.write('\n#rx: {} {} 0'.format(rx_min_x + i * rx_spacing, rx_y))

        f.write('\n#box: 0 {} 0 {} {} {} air'.format(y_air_bottom, domain_width, domain_height, dx_dy_dz))
        f.write('\n#box: 0 0 0 {} {} {} ice'.format(domain_width, y_air_bottom, dx_dy_dz))

        for i in range(n_blocks):
            x0 = i * block_width
            x1 = (i + 1) * block_width if i < n_blocks - 1 else domain_width
            x1 = min(x1, domain_width)
            mat = 'fracture_plus' if i % 2 == 0 else 'fracture_minus'
            f.write('\n#box: {} {} 0 {} {} {} {}'.format(x0, y_frac_min, x1, y_frac_max, dx_dy_dz, mat))

        f.write('\n#python:')
        f.write('\nfor i in range(1, {} + 1):'.format(n_snapshots))
        f.write(
            "\n    print('#snapshot: 0 0 0 {} {} {} {} {} {} {{}} snapshot_mid_x_{{}}'.format(i*{}, i))".format(
                domain_width, domain_height, dx_dy_dz, dx_dy_dz, dx_dy_dz, dx_dy_dz, snapshot_time
            )
        )
        f.write('\n#end_python:')

    print(f"Created: {out_file_base}")

    out_file_tl = os.path.join(base_dir, 'configs', 'horizontal_scattering_timelapse.in')
    with open(out_file_tl, 'w') as f:
        f.write('#title: Horizontal Scattering Time-Lapse')
        f.write('\n#domain: {} {} {}'.format(domain_width, domain_height, dx_dy_dz))
        f.write('\n#dx_dy_dz: {} {} {}'.format(dx_dy_dz, dx_dy_dz, dx_dy_dz))
        f.write('\n#time_window: {}'.format(TW))

        f.write('\n#material: {} {} 1 0 ice'.format(permittivity_ice, conductivity_ice))
        f.write('\n#material: {} {} 1 0 air'.format(permittivity_air, conductivity_air))
        f.write('\n#material: {} {} 1 0 fracture_plus'.format(permittivity_fracture_plus, conductivity_fracture_plus))
        f.write('\n#material: {} {} 1 0 fracture_minus'.format(permittivity_fracture_minus, conductivity_fracture_minus))

        f.write('\n#waveform: ricker 1 {} my_ricker'.format(f_central))
        f.write('\n#hertzian_dipole: z {} {} 0 my_ricker'.format(tx_start_x, rx_y))
        
        for i in range(rx_count):
            f.write('\n#rx: {} {} 0'.format(rx_min_x + i * rx_spacing, rx_y))

        f.write('\n#box: 0 {} 0 {} {} {} air'.format(y_air_bottom, domain_width, domain_height, dx_dy_dz))
        f.write('\n#box: 0 0 0 {} {} {} ice'.format(domain_width, y_air_bottom, dx_dy_dz))

        for i in range(n_blocks):
            x0 = i * block_width + time_lapse_shift
            x1 = (i + 1) * block_width + time_lapse_shift if i < n_blocks - 1 else domain_width + time_lapse_shift
            x0 = min(max(0, x0), domain_width)
            x1 = min(max(0, x1), domain_width)
            if x0 >= x1: continue
            mat = 'fracture_plus' if i % 2 == 0 else 'fracture_minus'
            f.write('\n#box: {} {} 0 {} {} {} {}'.format(x0, y_frac_min, x1, y_frac_max, dx_dy_dz, mat))

        f.write('\n#python:')
        f.write('\nfor i in range(1, {} + 1):'.format(n_snapshots))
        f.write(
            "\n    print('#snapshot: 0 0 0 {} {} {} {} {} {} {{}} snapshot_mid_x_{{}}'.format(i*{}, i))".format(
                domain_width, domain_height, dx_dy_dz, dx_dy_dz, dx_dy_dz, dx_dy_dz, snapshot_time
            )
        )
        f.write('\n#end_python:')

    print(f"Created: {out_file_tl}")

    output_context = {
        'domain_width': domain_width,
        'domain_height': domain_height,
        'air_thickness': air_thickness,
        'fracture_top': fracture_top,
        'fracture_bottom': fracture_bottom,
        'snapshot_time': snapshot_time,
        'dx_dy_dz': dx_dy_dz,
        'n_blocks': n_blocks,
        'block_width': block_width,
        'tx_start_x': tx_start_x,
        'rx_start_x': rx_min_x,
        'total_rx': rx_count,
        'actual_traces': 1,
        'block_label': block_label,
        'block_label_print': block_label_print
    }

    return output_context


def run_gprmax(traces_count=1):
    import subprocess
    import sys
    print("Running gprMax with GPU support via PowerShell...")
    import os
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    py_exe = sys.executable
    base_in = os.path.join(base_dir, "configs", "horizontal_scattering_baseline.in")
    tl_in = os.path.join(base_dir, "configs", "horizontal_scattering_timelapse.in")

    for snaps_dir in ["horizontal_scattering_baseline_snaps", "horizontal_scattering_timelapse_snaps"]:
        path = os.path.join(base_dir, "data", "outputs", snaps_dir)
        if os.path.exists(path): safe_remove_dir(path)
        path = os.path.join(base_dir, "configs", snaps_dir)
        if os.path.exists(path): safe_remove_dir(path)

    vcvars_script = r"C:\Progra~2\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
    ps_command = f"""
    $envDump = cmd /c '"{vcvars_script}" && set'
    $envDump | ForEach-Object {{ $p = $_ -split '=',2; if($p.Length -eq 2){{ Set-Item -Path env:$($p[0]) -Value $p[1] }} }}
    
    & "{py_exe}" -m gprMax "{base_in}" -n 1 -gpu
    if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
    
    & "{py_exe}" -m gprMax "{tl_in}" -n 1 -gpu
    if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
    """
    subprocess.run(["powershell", "-Command", ps_command], check=True)
    print("gprMax execution completed.")

def run_gprmax_to_subdir(traces_count=1, output_subdir=None):
    run_gprmax(traces_count=traces_count)
    if output_subdir is None:
        return
    import os
    import shutil
    import glob
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dest_dir = os.path.join(base_dir, "data", "outputs", output_subdir)
    os.makedirs(dest_dir, exist_ok=True)
    candidates = [
        os.path.join(base_dir, "configs", "horizontal_scattering_baseline.out"),
        os.path.join(base_dir, "configs", "horizontal_scattering_timelapse.out"),
    ]
    for src in candidates:
        if os.path.exists(src):
            try:
                shutil.move(src, os.path.join(dest_dir, os.path.basename(src)))
            except Exception:
                pass
    snapshot_roots = [
        os.path.join(base_dir, "configs", "horizontal_scattering_baseline_snaps"),
        os.path.join(base_dir, "configs", "horizontal_scattering_timelapse_snaps"),
    ]
    for root in snapshot_roots:
        for src in glob.glob(root + "*"):
            if os.path.isdir(src):
                try:
                    shutil.move(src, os.path.join(dest_dir, os.path.basename(src)))
                except Exception:
                    pass


def create_bscan_input_file(
        f_central=1.5 * 1e9,
        c=3 * 1e8,
        permittivity_ice=3.15,
        permittivity_air=1,
        permittivity_fracture=80,
        conductivity_ice=1e-6,
        conductivity_air=0,
        conductivity_fracture=0.01,
        fracture_depth=0.6,
        depth_below_fracture=0.1,
        air_thickness=0.1,
        rx_count=50,
        rx_spread=2.0,
        tx_rx_offset=0.1,
        block_size=None,
        block_size_fraction=(1/4),
        time_lapse_shift=0.0
    ):
    print(f"Generating gprMax B-scan input files...")

    c_ice = c / np.sqrt(permittivity_ice)
    wavelength_ice = c_ice / f_central

    permittivity_fracture_plus = 1.0
    permittivity_fracture_minus = 80
    conductivity_fracture_plus = 0
    conductivity_fracture_minus = 0.01

    if block_size is None:
        block_size = block_size_fraction * wavelength_ice

    frac = block_size / wavelength_ice
    fr = Fraction(frac).limit_denominator(32)
    if fr.numerator == 0:
        block_label = "0Lambda"
        block_label_print = "0Lambda"
    else:
        block_label = f"{fr.numerator}_{fr.denominator}Lambda"
        block_label_print = f"{fr.numerator}/{fr.denominator}Lambda"

    wavelength_plus = (c / np.sqrt(permittivity_fracture_plus)) / f_central
    wavelength_minus = (c / np.sqrt(permittivity_fracture_minus)) / f_central
    min_wavelength = min(wavelength_plus, wavelength_minus)
    thickness_fracture = min_wavelength / 5

    dx_dy_dz = min_wavelength / 30

    domain_width = 4.0
    domain_height = air_thickness + fracture_depth + thickness_fracture + depth_below_fracture
    domain_width = np.ceil(domain_width / dx_dy_dz) * dx_dy_dz
    domain_height = np.ceil(domain_height / dx_dy_dz) * dx_dy_dz

    TW = (5 / f_central) + (3 * (domain_height - air_thickness) / c_ice)

    fracture_top = air_thickness + fracture_depth
    fracture_bottom = fracture_top + thickness_fracture

    y_air_bottom = domain_height - air_thickness
    y_frac_top_gpr = domain_height - fracture_top
    y_frac_bottom_gpr = domain_height - fracture_bottom
    y_frac_min = min(y_frac_top_gpr, y_frac_bottom_gpr)
    y_frac_max = max(y_frac_top_gpr, y_frac_bottom_gpr)

    n_blocks = int(np.round(domain_width / block_size))
    if n_blocks < 1: n_blocks = 1
    if n_blocks % 2 != 0: n_blocks += 1
    block_width = domain_width / n_blocks

    # Midpoints of Tx-Rx pair sweep the same positions as the original receivers
    midpoint_start_x = 2.0 - (rx_spread / 2.0)
    midpoint_spacing = rx_spread / max(1, rx_count - 1)
    tx_x_start = midpoint_start_x - tx_rx_offset / 2.0
    rx_x_start = midpoint_start_x + tx_rx_offset / 2.0
    rx_y = y_air_bottom

    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    def _write_bscan_config(f, title, block_shift):
        f.write(f'#title: {title}')
        f.write('\n#domain: {} {} {}'.format(domain_width, domain_height, dx_dy_dz))
        f.write('\n#dx_dy_dz: {} {} {}'.format(dx_dy_dz, dx_dy_dz, dx_dy_dz))
        f.write('\n#time_window: {}'.format(TW))
        f.write('\n#material: {} {} 1 0 ice'.format(permittivity_ice, conductivity_ice))
        f.write('\n#material: {} {} 1 0 air'.format(permittivity_air, conductivity_air))
        f.write('\n#material: {} {} 1 0 fracture_plus'.format(permittivity_fracture_plus, conductivity_fracture_plus))
        f.write('\n#material: {} {} 1 0 fracture_minus'.format(permittivity_fracture_minus, conductivity_fracture_minus))
        f.write('\n#waveform: ricker 1 {} my_ricker'.format(f_central))
        f.write('\n#hertzian_dipole: z {} {} 0 my_ricker'.format(tx_x_start, rx_y))
        f.write('\n#rx: {} {} 0'.format(rx_x_start, rx_y))
        f.write('\n#src_steps: {} 0 0'.format(midpoint_spacing))
        f.write('\n#rx_steps: {} 0 0'.format(midpoint_spacing))
        f.write('\n#box: 0 {} 0 {} {} {} air'.format(y_air_bottom, domain_width, domain_height, dx_dy_dz))
        f.write('\n#box: 0 0 0 {} {} {} ice'.format(domain_width, y_air_bottom, dx_dy_dz))
        for i in range(n_blocks):
            x0 = i * block_width + block_shift
            x1 = (i + 1) * block_width + block_shift if i < n_blocks - 1 else domain_width + block_shift
            x0 = min(max(0, x0), domain_width)
            x1 = min(max(0, x1), domain_width)
            if x0 >= x1: continue
            mat = 'fracture_plus' if i % 2 == 0 else 'fracture_minus'
            f.write('\n#box: {} {} 0 {} {} {} {}'.format(x0, y_frac_min, x1, y_frac_max, dx_dy_dz, mat))

    out_file_base = os.path.join(base_dir, 'configs', 'horizontal_scattering_baseline_bscan.in')
    with open(out_file_base, 'w') as f:
        _write_bscan_config(f, 'Horizontal Scattering B-scan Baseline', block_shift=0.0)
    print(f"Created: {out_file_base}")

    out_file_tl = os.path.join(base_dir, 'configs', 'horizontal_scattering_timelapse_bscan.in')
    with open(out_file_tl, 'w') as f:
        _write_bscan_config(f, 'Horizontal Scattering B-scan Time-Lapse', block_shift=time_lapse_shift)
    print(f"Created: {out_file_tl}")

    return {
        'domain_width': domain_width,
        'domain_height': domain_height,
        'air_thickness': air_thickness,
        'fracture_top': fracture_top,
        'fracture_bottom': fracture_bottom,
        'dx_dy_dz': dx_dy_dz,
        'n_blocks': n_blocks,
        'block_width': block_width,
        'midpoint_start_x': midpoint_start_x,
        'midpoint_spacing': midpoint_spacing,
        'tx_rx_offset': tx_rx_offset,
        'n_traces': rx_count,
        'block_label': block_label,
        'block_label_print': block_label_print,
    }


def run_gprmax_bscan(n_traces=50):
    import subprocess
    import sys
    print("Running gprMax B-scan with GPU support via PowerShell...")
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    py_exe = sys.executable
    base_in = os.path.join(base_dir, "configs", "horizontal_scattering_baseline_bscan.in")
    tl_in = os.path.join(base_dir, "configs", "horizontal_scattering_timelapse_bscan.in")

    vcvars_script = r"C:\Progra~2\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
    ps_command = f"""
    $envDump = cmd /c '"{vcvars_script}" && set'
    $envDump | ForEach-Object {{ $p = $_ -split '=',2; if($p.Length -eq 2){{ Set-Item -Path env:$($p[0]) -Value $p[1] }} }}

    & "{py_exe}" -m gprMax "{base_in}" -n {n_traces} -gpu
    if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}

    & "{py_exe}" -m gprMax "{tl_in}" -n {n_traces} -gpu
    if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
    """
    subprocess.run(["powershell", "-Command", ps_command], check=True)
    print("gprMax B-scan execution completed.")


def merge_bscan_outputs(configs_dir, base_names):
    """
    Merge per-trace gprMax B-scan output files (basename1.out, basename2.out, ...)
    into a single combined .out file and remove the originals.
    """
    import h5py
    import glob
    import re

    for base_name in base_names:
        pattern = os.path.join(configs_dir, f"{base_name}[0-9]*.out")
        individual_files = glob.glob(pattern)

        if not individual_files:
            print(f"No files found matching {pattern}")
            continue

        def _trace_num(path):
            m = re.search(r'(\d+)\.out$', path)
            return int(m.group(1)) if m else 0

        individual_files.sort(key=_trace_num)
        print(f"Merging {len(individual_files)} traces for {base_name}...")

        merged_path = os.path.join(configs_dir, f"{base_name}.out")
        with h5py.File(merged_path, 'w') as fout:
            with h5py.File(individual_files[0], 'r') as f0:
                for k, v in f0.attrs.items():
                    fout.attrs[k] = v

            rxs_grp = fout.create_group('rxs')
            for i, fpath in enumerate(individual_files):
                with h5py.File(fpath, 'r') as fi:
                    src_rxs = fi['rxs'] if 'rxs' in fi else fi
                    src_rx = (src_rxs['rx1'] if 'rx1' in src_rxs
                              else next(iter(src_rxs.values())))
                    dst_rx = rxs_grp.create_group(f'rx{i + 1}')
                    for field in src_rx.keys():
                        dst_rx.create_dataset(field, data=src_rx[field][:])
                    for k, v in src_rx.attrs.items():
                        dst_rx.attrs[k] = v

        for fpath in individual_files:
            try:
                os.remove(fpath)
            except Exception:
                pass

        print(f"Merged: {merged_path}")


def run_gprmax_bscan_to_subdir(n_traces=50, output_subdir=None):
    run_gprmax_bscan(n_traces=n_traces)
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    configs_dir = os.path.join(base_dir, "configs")

    merge_bscan_outputs(configs_dir, [
        "horizontal_scattering_baseline_bscan",
        "horizontal_scattering_timelapse_bscan",
    ])

    if output_subdir is None:
        return
    dest_dir = os.path.join(base_dir, "data", "outputs", output_subdir)
    os.makedirs(dest_dir, exist_ok=True)
    for fname in [
        "horizontal_scattering_baseline_bscan.out",
        "horizontal_scattering_timelapse_bscan.out",
    ]:
        src = os.path.join(configs_dir, fname)
        if os.path.exists(src):
            try:
                shutil.move(src, os.path.join(dest_dir, fname))
            except Exception:
                pass