import os
import glob
import shutil
import numpy as np
import subprocess

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
                shutil.rmtree(folder)
                print(f"Removed directory: {folder}")
            except Exception as e:
                print(f"Error removing directory {folder}: {e}")

def create_input_file(
        f_central=1.5 * 1e9,
        c=3 * 1e8,
        permittivity_ice=6,
        permittivity_air=1,
        permittivity_fracture=6,
        conductivity_ice=1e-6,
        conductivity_air=0,
        conductivity_fracture=0.001,
        fracture_depth=0.3,
        depth_below_fracture=0.1,
        air_thickness=0.05,
        rx_per_block=4
    ):
    print("Generating gprMax input file...")
    
    c_ice   = c / np.sqrt(permittivity_ice)
    c_water = c / np.sqrt(permittivity_fracture)
    c_fracture = c / np.sqrt(permittivity_fracture)
    wavelength_fracture = c_fracture / f_central

    thickness_fracture = wavelength_fracture / 10

    permittivity_fracture_plus  = permittivity_fracture * 1.0
    permittivity_fracture_minus = permittivity_fracture * 1.0
    conductivity_fracture_plus  = conductivity_fracture * 1.1
    conductivity_fracture_minus = conductivity_fracture * 0.9

    block_size = (1/4) * wavelength_fracture

    dx_min = wavelength_fracture / 20
    dx_dy_dz = dx_min
    source_receiver_steps = wavelength_fracture / 10

    domain_width   = 12 * wavelength_fracture

    domain_height = air_thickness + fracture_depth + thickness_fracture + depth_below_fracture

    TW = (5 / f_central) + (2 * (domain_height - air_thickness) / c_ice)
    snapshot_time = TW / 10
    n_snapshots = 10

    fracture_top = air_thickness + fracture_depth
    fracture_bottom = fracture_top + thickness_fracture

    y_air_bottom = domain_height - air_thickness
    y_frac_top_gpr = domain_height - fracture_top
    y_frac_bottom_gpr = domain_height - fracture_bottom
    y_frac_min = min(y_frac_top_gpr, y_frac_bottom_gpr)
    y_frac_max = max(y_frac_top_gpr, y_frac_bottom_gpr)

    x_first_measurement = 1/2 * wavelength_fracture
    # Pre-calculate positions for step 55 (54 steps offset)
    tx_start_x = x_first_measurement + 54 * source_receiver_steps
    tx_rx_offset = source_receiver_steps
    rx_start_x = tx_start_x + tx_rx_offset
    rx_y = y_air_bottom
    
    n_blocks = int(np.round(domain_width / block_size))
    if n_blocks < 1: n_blocks = 1
    if n_blocks % 2 != 0: n_blocks += 1
    block_width = domain_width / n_blocks

    target_model_run = 56 # Assuming from user's -restart 55

    # 1. Generate Alternating Input File
    out_file_alt = 'horizontal_scattering_0p5lambda.in'

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
        f.write('\n#rx: {} {} 0'.format(rx_start_x, rx_y))
        rx_dx = block_width / rx_per_block
        for i in range(n_blocks * rx_per_block):
            rx_x_init = (i + 0.5) * rx_dx
            # Ensure receiver starts within domain
            if rx_x_init > dx_dy_dz and rx_x_init < domain_width - dx_dy_dz:
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

    print(f"Created: {out_file_alt}")

    # 2. Generate Homogeneous Input File
    out_file_homo = 'horizontal_scattering_homogeneous.in'

    with open(out_file_homo, 'w') as f:
        f.write('#title: Horizontal Scattering Homogeneous')
        f.write('\n#domain: {} {} {}'.format(domain_width, domain_height, dx_dy_dz))
        f.write('\n#dx_dy_dz: {} {} {}'.format(dx_dy_dz, dx_dy_dz, dx_dy_dz))
        f.write('\n#time_window: {}'.format(TW))

        f.write('\n#material: {} {} 1 0 ice'.format(permittivity_ice, conductivity_ice))
        f.write('\n#material: {} {} 1 0 air'.format(permittivity_air, conductivity_air))
        f.write('\n#material: {} {} 1 0 fracture_homo'.format(permittivity_fracture, conductivity_fracture))

        f.write('\n#waveform: ricker 1 {} my_ricker'.format(f_central))
        f.write('\n#hertzian_dipole: z {} {} 0 my_ricker'.format(tx_start_x, rx_y))
        f.write('\n#rx: {} {} 0'.format(rx_start_x, rx_y))
        rx_dx = block_width / rx_per_block
        for i in range(n_blocks * rx_per_block):
            rx_x_init = (i + 0.5) * rx_dx
            if rx_x_init > dx_dy_dz and rx_x_init < domain_width - dx_dy_dz:
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

    return domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time, dx_dy_dz, n_blocks, block_width, target_model_run * source_receiver_steps, rx_per_block

def run_gprmax():
    print("Running gprMax with GPU support via PowerShell...")
    snapshots_dir_alt = "horizontal_scattering_0p5lambda_snaps"
    if os.path.exists(snapshots_dir_alt):
        print(f"Removing old snapshots directory: {snapshots_dir_alt}")
        shutil.rmtree(snapshots_dir_alt)

    snapshots_dir_homo = "horizontal_scattering_homogeneous_snaps"
    if os.path.exists(snapshots_dir_homo):
        print(f"Removing old snapshots directory: {snapshots_dir_homo}")
        shutil.rmtree(snapshots_dir_homo)

    ps_command = r'''
    $envDump = cmd /c '"C:\Progra~2\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" && set'
    $envDump | ForEach-Object { $p = $_ -split '=',2; if($p.Length -eq 2){ Set-Item -Path env:$($p[0]) -Value $p[1] } }
    
    python -m gprMax "C:\Users\Administrator\Thesis\Diffraction_Experiment\horizontal_scattering_0p5lambda.in" -n 1 -gpu
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    
    python -m gprMax "C:\Users\Administrator\Thesis\Diffraction_Experiment\horizontal_scattering_homogeneous.in" -n 1 -gpu
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    '''
    
    subprocess.run(["powershell", "-Command", ps_command], check=True)
    print("gprMax execution completed.")
