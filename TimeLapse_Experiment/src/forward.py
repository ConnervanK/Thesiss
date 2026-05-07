import os
import glob
import shutil
import numpy as np
import subprocess
import stat
import time
import sys
from fractions import Fraction


def _rmtree_onerror(func, path, exc_info):
    """Best-effort handler for read-only files during shutil.rmtree on Windows."""
    try:
        os.chmod(path, stat.S_IWRITE)
        func(path)
    except Exception:
        pass


def safe_remove_dir(path, retries=3, delay_s=0.4):
    """Remove a directory tree robustly on Windows without hard-failing the run."""
    if not os.path.isdir(path):
        return True

    for attempt in range(1, retries + 1):
        try:
            shutil.rmtree(path, onerror=_rmtree_onerror)
            return True
        except PermissionError as e:
            if attempt < retries:
                print(f"Directory busy/locked, retrying ({attempt}/{retries}) for {path}: {e}")
                time.sleep(delay_s)
                continue
            print(f"Warning: could not remove locked directory {path}: {e}")
            return False
        except Exception as e:
            print(f"Warning: could not remove directory {path}: {e}")
            return False

    return False

def clear_results(directory="."):
    """
    Clears all result files (.vti, .in, .out) and output folders from the directory,
    but keeps the python scripts.
    """
    print(f"Clearing results in {os.path.abspath(directory)}...")
    for ext in ["*.vti", "*.in", "*.out", "*.png"]:
        for file in glob.glob(os.path.join(directory, ext)):
            try:
                os.remove(file)
                print(f"Removed file: {file}")
            except Exception as e:
                print(f"Error removing {file}: {e}")
                
    # Optionally also remove generated snapshot directories
    for folder in glob.glob(os.path.join(directory, "*_snaps")):
        if os.path.isdir(folder):
            try:
                removed = safe_remove_dir(folder)
                if not removed:
                    continue
                print(f"Removed directory: {folder}")
            except Exception as e:
                print(f"Error removing directory {folder}: {e}")

def create_input_file(
        f_central=1.5 * 1e9,
        c=3 * 1e8,
        permittivity_ice=3.15,
        permittivity_air=1,
        permittivity_fracture=80,
        conductivity_ice=1e-6,
        conductivity_air=0,
        conductivity_fracture=0.01,
        fracture_depth=0.3,
        depth_below_fracture=0.1,
        air_thickness=0.05,
        mode='static',          # 'static' for one measurement with array across domain, 'bscan' for moving Tx-Rx
        rx_per_block=1,         # Used only if mode == 'static'
        rx_count=1,             # Number of receivers used if mode == 'bscan'
        rx_spacing=0.05,        # Receiver spacing if mode == 'bscan' and rx_count > 1
        bscan_traces=50,        # Number of traces for 'bscan'
        bscan_step_x=0.05,      # Movement step for 'bscan'
        block_size=None,        # Absolute block size in metres. If None, uses block_size_fraction * wavelength_ice
        block_size_fraction=(1/4),
        time_lapse_shift=0.0    # Lateral shift (m) for the time-lapse model (must be < wavelength_ice)
    ):
    print(f"Generating gprMax input file (mode: {mode})...")
    
    c_ice   = c / np.sqrt(permittivity_ice)
    c_water = c / np.sqrt(permittivity_fracture)
    wavelength_ice = c_ice / f_central

    

    permittivity_fracture_plus  = 1.0 # air permittivity_fracture * 0.9
    permittivity_fracture_minus = 80 # water permittivity_fracture * 1.0

    conductivity_fracture_plus  = 0 # air conductivity_fracture * 1.1
    conductivity_fracture_minus = 0.01 # water conductivity_fracture * 0.9

    permittivity_fracture_average = (permittivity_fracture_plus + permittivity_fracture_minus) / 2
    conductivity_fracture_average = (conductivity_fracture_plus + conductivity_fracture_minus) / 2

    # Determine block size (meters). Prefer explicit `block_size` if provided, else use fraction.
    if block_size is None:
        block_size = block_size_fraction * wavelength_ice

    # Create a human- and filesystem-friendly label for this block size using a rational fraction of wavelength
    frac = block_size / wavelength_ice
    fr = Fraction(frac).limit_denominator(32)
    if fr.numerator == 0:
        block_label = f"0Lambda"
        block_label_print = "0Lambda"
    else:
        block_label = f"{fr.numerator}_{fr.denominator}Lambda"
        block_label_print = f"{fr.numerator}/{fr.denominator}Lambda"

    wavelength_plus = (c / np.sqrt(permittivity_fracture_plus)) / f_central
    wavelength_minus = (c / np.sqrt(permittivity_fracture_minus)) / f_central
    min_wavelength = min(wavelength_plus, wavelength_minus)

    
    thickness_fracture = min(wavelength_plus, wavelength_minus) / 5

    dx_min = min_wavelength / 30
    dx_dy_dz = dx_min
    source_receiver_steps = wavelength_ice / 10

    domain_width   = 12 * wavelength_ice

    domain_height = air_thickness + fracture_depth + thickness_fracture + depth_below_fracture

    TW = (5 / f_central) + (3 * (domain_height - air_thickness) / c_ice)
    n_snapshots = 36
    snapshot_time = 0.25e-9  # quarter of a nanosecond

    fracture_top = air_thickness + fracture_depth
    fracture_bottom = fracture_top + thickness_fracture

    y_air_bottom = domain_height - air_thickness
    y_frac_top_gpr = domain_height - fracture_top
    y_frac_bottom_gpr = domain_height - fracture_bottom
    y_frac_min = min(y_frac_top_gpr, y_frac_bottom_gpr)
    y_frac_max = max(y_frac_top_gpr, y_frac_bottom_gpr)

    wavelength_fracture = min_wavelength
    x_first_measurement = 1/2 * wavelength_fracture
    
    n_blocks = int(np.round(domain_width / block_size))
    if n_blocks < 1: n_blocks = 1
    if n_blocks % 2 != 0: n_blocks += 1
    block_width = domain_width / n_blocks

    target_model_run = 56 # Assuming from user's -restart 55
    tx_rx_offset = source_receiver_steps
    rx_y = y_air_bottom
    
    if mode == 'bscan':
        tx_start_x = domain_width / 2 #x_first_measurement
        rx_start_x = tx_start_x + tx_rx_offset
        actual_traces = bscan_traces
        
        # Check domain limits and fix traces if needed
        max_dist = rx_start_x + (rx_count - 1) * rx_spacing + (bscan_traces - 1) * bscan_step_x
        if max_dist > domain_width - dx_dy_dz:
            allowed_traces = int((domain_width - dx_dy_dz - (rx_start_x + (rx_count - 1) * rx_spacing)) / bscan_step_x) + 1
            print(f"Warning: B-scan would exceed domain width. Reducing traces from {bscan_traces} to {allowed_traces}.")
            actual_traces = max(1, allowed_traces)

        total_rx = rx_count
    else:
        # Puts the source in the middle to measure the block diffractions.
        tx_start_x = domain_width / 2
        rx_start_x = tx_start_x + tx_rx_offset
        actual_traces = 1
        total_rx = n_blocks * rx_per_block

    # 1. Generate Alternating Input File
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out_file_alt = os.path.join(base_dir, 'configs', 'horizontal_scattering_0p5lambda.in')

    with open(out_file_alt, 'w') as f:
        f.write('#title: Horizontal Scattering Experiment')
        f.write('\n#domain: {} {} {}'.format(domain_width, domain_height, dx_dy_dz))
        f.write('\n#dx_dy_dz: {} {} {}'.format(dx_dy_dz, dx_dy_dz, dx_dy_dz))
        f.write('\n#time_window: {}'.format(TW))

        f.write('\n#material: {} {} 1 0 ice'.format(permittivity_ice, conductivity_ice))
        f.write('\n#material: {} {} 1 0 air'.format(permittivity_air, conductivity_air))
        f.write('\n#material: {} {} 1 0 fracture_plus'.format(permittivity_fracture_plus, conductivity_fracture_plus))
        f.write('\n#material: {} {} 1 0 fracture_minus'.format(permittivity_fracture_minus, conductivity_fracture_minus))

        f.write('\n#waveform: ricker 1 {} my_ricker'.format(f_central))
        f.write('\n#hertzian_dipole: z {} {} 0 my_ricker'.format(tx_start_x, rx_y))
        
        if mode == 'static':
            rx_dx = block_width / rx_per_block
            for i in range(n_blocks * rx_per_block):
                rx_x_init = (i + 0.5) * rx_dx
                if rx_x_init > dx_dy_dz and rx_x_init < domain_width - dx_dy_dz:
                    f.write('\n#rx: {} {} 0'.format(rx_x_init, rx_y))
        elif mode == 'bscan':
            f.write('\n#src_steps: {} 0 0'.format(bscan_step_x))
            f.write('\n#rx_steps: {} 0 0'.format(bscan_step_x))
            for i in range(rx_count):
                rx_x_init = rx_start_x + i * rx_spacing
                f.write('\n#rx: {} {} 0'.format(rx_x_init, rx_y))

        f.write('\n#box: 0 {} 0 {} {} {} air'.format(y_air_bottom, domain_width, domain_height, dx_dy_dz))
        f.write('\n#box: 0 0 0 {} {} {} ice'.format(domain_width, y_air_bottom, dx_dy_dz))

        for i in range(n_blocks):
            x0 = i * block_width
            x1 = (i + 1) * block_width if i < n_blocks - 1 else domain_width
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

        # Comment out geometry_view as an int8 overflow occurs in gprMax 
        # when adding more than exactly 127 receivers to the VTK file export.
        # f.write('\n#geometry_view: 0 0 0 {} {} {} {} {} {} horizontal_scattering_0p5lambda n'.format(domain_width, domain_height, dx_dy_dz, dx_dy_dz, dx_dy_dz, dx_dy_dz))

    print(f"Created: {out_file_alt}  (block size: {block_size:.4e} m -> {block_label_print})")

    # 2. Generate Homogeneous Input File
    out_file_homo = os.path.join(base_dir, 'configs', 'horizontal_scattering_homogeneous.in')

    with open(out_file_homo, 'w') as f:
        f.write('#title: Horizontal Scattering Homogeneous')
        f.write('\n#domain: {} {} {}'.format(domain_width, domain_height, dx_dy_dz))
        f.write('\n#dx_dy_dz: {} {} {}'.format(dx_dy_dz, dx_dy_dz, dx_dy_dz))
        f.write('\n#time_window: {}'.format(TW))

        f.write('\n#material: {} {} 1 0 ice'.format(permittivity_ice, conductivity_ice))
        f.write('\n#material: {} {} 1 0 air'.format(permittivity_air, conductivity_air))
        f.write('\n#material: {} {} 1 0 fracture_homo'.format(permittivity_fracture_average, conductivity_fracture_average))

        f.write('\n#waveform: ricker 1 {} my_ricker'.format(f_central))
        f.write('\n#hertzian_dipole: z {} {} 0 my_ricker'.format(tx_start_x, rx_y))
        
        if mode == 'static':
            rx_dx = block_width / rx_per_block
            for i in range(n_blocks * rx_per_block):
                rx_x_init = (i + 0.5) * rx_dx
                if rx_x_init > dx_dy_dz and rx_x_init < domain_width - dx_dy_dz: 
                    f.write('\n#rx: {} {} 0'.format(rx_x_init, rx_y))
        elif mode == 'bscan':
            f.write('\n#src_steps: {} 0 0'.format(bscan_step_x))
            f.write('\n#rx_steps: {} 0 0'.format(bscan_step_x))
            for i in range(rx_count):
                rx_x_init = rx_start_x + i * rx_spacing
                f.write('\n#rx: {} {} 0'.format(rx_x_init, rx_y))

        f.write('\n#box: 0 {} 0 {} {} {} air'.format(y_air_bottom, domain_width, domain_height, dx_dy_dz))
        f.write('\n#box: 0 0 0 {} {} {} ice'.format(domain_width, y_air_bottom, dx_dy_dz))

        f.write('\n#box: 0 {} 0 {} {} {} fracture_homo'.format(y_frac_min, domain_width, y_frac_max, dx_dy_dz))

        f.write('\n#python:')
        f.write('\nfor i in range(1, {} + 1):'.format(n_snapshots))
        f.write(
            "\n    print('#snapshot: 0 0 0 {} {} {} {} {} {} {{}} snapshot_mid_x_{{}}'.format(i*{}, i))".format(
                domain_width, domain_height, dx_dy_dz, dx_dy_dz, dx_dy_dz, dx_dy_dz, snapshot_time
            )
        )
        f.write('\n#end_python:')

        # Comment out geometry_view to prevent int8 receiver overflow in gprMax
        # f.write('\n#geometry_view: 0 0 0 {} {} {} {} {} {} horizontal_scattering_homogeneous n'.format(domain_width, domain_height, dx_dy_dz, dx_dy_dz, dx_dy_dz, dx_dy_dz))

    print(f"Created: {out_file_homo}")

    # 3. Optionally generate a time-lapse shifted alternating file when a small shift is provided
    timelapse_file = None
    if time_lapse_shift and abs(time_lapse_shift) < wavelength_ice:
        timelapse_file = os.path.join(base_dir, 'configs', 'horizontal_scattering_0p5lambda_timelapse.in')
        with open(timelapse_file, 'w') as f:
            f.write('#title: Horizontal Scattering Experiment (Time-Lapse Shifted)')
            f.write('\n#domain: {} {} {}'.format(domain_width, domain_height, dx_dy_dz))
            f.write('\n#dx_dy_dz: {} {} {}'.format(dx_dy_dz, dx_dy_dz, dx_dy_dz))
            f.write('\n#time_window: {}'.format(TW))

            f.write('\n#material: {} {} 1 0 ice'.format(permittivity_ice, conductivity_ice))
            f.write('\n#material: {} {} 1 0 air'.format(permittivity_air, conductivity_air))
            f.write('\n#material: {} {} 1 0 fracture_plus'.format(permittivity_fracture_plus, conductivity_fracture_plus))
            f.write('\n#material: {} {} 1 0 fracture_minus'.format(permittivity_fracture_minus, conductivity_fracture_minus))

            f.write('\n#waveform: ricker 1 {} my_ricker'.format(f_central))
            f.write('\n#hertzian_dipole: z {} {} 0 my_ricker'.format(tx_start_x, rx_y))

            if mode == 'static':
                rx_dx = block_width / rx_per_block
                for i in range(n_blocks * rx_per_block):
                    rx_x_init = (i + 0.5) * rx_dx
                    if rx_x_init > dx_dy_dz and rx_x_init < domain_width - dx_dy_dz:
                        f.write('\n#rx: {} {} 0'.format(rx_x_init, rx_y))
            elif mode == 'bscan':
                f.write('\n#src_steps: {} 0 0'.format(bscan_step_x))
                f.write('\n#rx_steps: {} 0 0'.format(bscan_step_x))
                for i in range(rx_count):
                    rx_x_init = rx_start_x + i * rx_spacing
                    f.write('\n#rx: {} {} 0'.format(rx_x_init, rx_y))

            f.write('\n#box: 0 {} 0 {} {} {} air'.format(y_air_bottom, domain_width, domain_height, dx_dy_dz))
            f.write('\n#box: 0 0 0 {} {} {} ice'.format(domain_width, y_air_bottom, dx_dy_dz))

            # Shift the blocks by the requested (small) amount
            for i in range(n_blocks):
                x0 = i * block_width + time_lapse_shift
                x1 = (i + 1) * block_width + time_lapse_shift if i < n_blocks - 1 else domain_width + time_lapse_shift
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

        print(f"Created time-lapse config: {timelapse_file}")

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
        'rx_start_x': rx_start_x,
        'total_rx': total_rx,
        'actual_traces': actual_traces,
        'mode': mode,
        'rx_per_block': rx_per_block
    }

    # Include block label metadata for downstream routing of outputs
    output_context['block_label'] = block_label
    output_context['block_label_print'] = block_label_print
    return output_context

def run_gprmax(traces_count=1):
    print("Running gprMax with GPU support via PowerShell...")
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    py_exe = sys.executable
    alt_in = os.path.join(base_dir, "configs", "horizontal_scattering_0p5lambda.in")
    homo_in = os.path.join(base_dir, "configs", "horizontal_scattering_homogeneous.in")

    snapshots_dir_alt = os.path.join(base_dir, "data", "outputs", "horizontal_scattering_0p5lambda_snaps")
    if os.path.exists(snapshots_dir_alt):
        print(f"Removing old snapshots directory: {snapshots_dir_alt}")
        safe_remove_dir(snapshots_dir_alt)

    snapshots_dir_homo = "horizontal_scattering_homogeneous_snaps"
    if os.path.exists(snapshots_dir_homo):
        print(f"Removing old snapshots directory: {snapshots_dir_homo}")
        safe_remove_dir(snapshots_dir_homo)

    ps_command = f'''
    $envDump = cmd /c '"C:\\Progra~2\\Microsoft Visual Studio\\2022\\BuildTools\\VC\\Auxiliary\\Build\\vcvars64.bat" && set'
    $envDump | ForEach-Object {{ $p = $_ -split '=',2; if($p.Length -eq 2){{ Set-Item -Path env:$($p[0]) -Value $p[1] }} }}
    
    & "{py_exe}" -m gprMax "{alt_in}" -n {traces_count} -gpu
    if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
    
    & "{py_exe}" -m gprMax "{homo_in}" -n {traces_count} -gpu
    if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
    '''
    
    subprocess.run(["powershell", "-Command", ps_command], check=True)
    print("gprMax execution completed.")
    
    if traces_count > 1:
        print("Merging multi-trace output files into a single B-scan...")
        from tools.outputfiles_merge import merge_files
        
        # Merge alternating file
        alt_base = os.path.join(base_dir, "configs", "horizontal_scattering_0p5lambda")
        merge_files(alt_base, removefiles=True)
        merged_alt = os.path.join(base_dir, "configs", "horizontal_scattering_0p5lambda_merged.out")
        target_alt = os.path.join(base_dir, "data", "outputs", "horizontal_scattering_0p5lambda.out")
        if os.path.exists(merged_alt):
            shutil.move(merged_alt, target_alt)
        elif os.path.exists(os.path.join(base_dir, "configs", "horizontal_scattering_0p5lambda.out")):
            shutil.move(os.path.join(base_dir, "configs", "horizontal_scattering_0p5lambda.out"), target_alt)
            
        # Merge homogeneous file
        homo_base = os.path.join(base_dir, "configs", "horizontal_scattering_homogeneous")
        merge_files(homo_base, removefiles=True)
        merged_homo = os.path.join(base_dir, "configs", "horizontal_scattering_homogeneous_merged.out")
        target_homo = os.path.join(base_dir, "data", "outputs", "horizontal_scattering_homogeneous.out")
        if os.path.exists(merged_homo):
            shutil.move(merged_homo, target_homo)
        elif os.path.exists(os.path.join(base_dir, "configs", "horizontal_scattering_homogeneous.out")):
            shutil.move(os.path.join(base_dir, "configs", "horizontal_scattering_homogeneous.out"), target_homo)

        # Move snapshot folders to data/outputs
        snaps_alt_src = os.path.join(base_dir, "configs", "horizontal_scattering_0p5lambda_snaps")
        snaps_alt_dst = os.path.join(base_dir, "data", "outputs", "horizontal_scattering_0p5lambda_snaps")
        if os.path.exists(snaps_alt_src):
            shutil.move(snaps_alt_src, snaps_alt_dst)

        snaps_homo_src = os.path.join(base_dir, "configs", "horizontal_scattering_homogeneous_snaps")
        snaps_homo_dst = os.path.join(base_dir, "data", "outputs", "horizontal_scattering_homogeneous_snaps")
        if os.path.exists(snaps_homo_src):
            shutil.move(snaps_homo_src, snaps_homo_dst)
            
        print("Merging complete.")


def run_gprmax_to_subdir(traces_count=1, output_subdir=None):
    """
    Run gprMax and move produced .out and snapshot folders
    into `data/outputs/<output_subdir>` if provided.
    """
    run_gprmax(traces_count=traces_count)

    if output_subdir is None:
        return

    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dest_dir = os.path.join(base_dir, "data", "outputs", output_subdir)
    os.makedirs(dest_dir, exist_ok=True)

    # Move known .out candidates from configs into the destination
    candidates = [
        os.path.join(base_dir, "configs", "horizontal_scattering_0p5lambda_merged.out"),
        os.path.join(base_dir, "configs", "horizontal_scattering_0p5lambda.out"),
        os.path.join(base_dir, "configs", "horizontal_scattering_homogeneous_merged.out"),
        os.path.join(base_dir, "configs", "horizontal_scattering_homogeneous.out"),
    ]
    for src in candidates:
        if os.path.exists(src):
            try:
                shutil.move(src, os.path.join(dest_dir, os.path.basename(src)))
            except Exception:
                pass

    # Move any snapshot folders created by gprMax into the labeled subdir.
    # gprMax writes numbered siblings like `horizontal_scattering_0p5lambda_snaps1`,
    # `horizontal_scattering_0p5lambda_snaps2`, etc. so we move the entire family.
    snapshot_roots = [
        os.path.join(base_dir, "configs", "horizontal_scattering_0p5lambda_snaps"),
        os.path.join(base_dir, "configs", "horizontal_scattering_homogeneous_snaps"),
        os.path.join(base_dir, "data", "outputs", "horizontal_scattering_0p5lambda_snaps"),
        os.path.join(base_dir, "data", "outputs", "horizontal_scattering_homogeneous_snaps"),
    ]
    for root in snapshot_roots:
        for src in glob.glob(root + "*"):
            if os.path.isdir(src):
                try:
                    shutil.move(src, os.path.join(dest_dir, os.path.basename(src)))
                except Exception:
                    pass

    print(f"Moved gprMax outputs into: {dest_dir}")
