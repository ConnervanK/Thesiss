import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pyvista as pv
import subprocess
import h5py



def create_input_file():
    print("Generating gprMax input file...")
    # Parameters from notebook
    f_central  = 1.5 * 1e9 # 1.5 GHz
    c = 3 * 1e8 # speed of light in m/s

    permittivity_ice      = 6
    permittivity_air      = 1
    permittivity_fracture = 6

    conductivity_ice      = 1e-6
    conductivity_air      = 0
    conductivity_fracture = 0.001 #0.001

    c_ice   = c / np.sqrt(permittivity_ice)
    c_water = c / np.sqrt(permittivity_fracture)
    c_fracture = c / np.sqrt(permittivity_fracture)
    wavelength_fracture = c_fracture / f_central

    thickness_fracture = wavelength_fracture / 10

    permittivity_fracture_plus  = permittivity_fracture * 1.0
    permittivity_fracture_minus = permittivity_fracture * 1.0
    conductivity_fracture_plus  = conductivity_fracture * 1.1
    conductivity_fracture_minus = conductivity_fracture * 0.9

    block_size = (1/2) * wavelength_fracture

    dx_min = wavelength_fracture / 20
    dx_dy_dz = dx_min
    source_receiver_steps = wavelength_fracture / 10

    domain_width   = 12 * wavelength_fracture
    fracture_depth = 0.3
    depth_below_fracture = 0.1
    air_thickness  = 0.05

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
        for i in range(n_blocks):
            rx_x_init = (i + 0.5) * block_width
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

        f.write('\n#geometry_view: 0 0 0 {} {} {} {} {} {} horizontal_scattering_0p5lambda n'.format(domain_width, domain_height, dx_dy_dz, dx_dy_dz, dx_dy_dz, dx_dy_dz))

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
        for i in range(n_blocks):
            rx_x_init = (i + 0.5) * block_width
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

        f.write('\n#geometry_view: 0 0 0 {} {} {} {} {} {} horizontal_scattering_homogeneous n'.format(domain_width, domain_height, dx_dy_dz, dx_dy_dz, dx_dy_dz, dx_dy_dz))

    print(f"Created: {out_file_homo}")

    return domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time, dx_dy_dz, n_blocks, block_width, target_model_run * source_receiver_steps

import os
import shutil

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
    
    python -m gprMax "C:\Users\Administrator\Thesis\True_and_Baseline_Models\Horizontal_scattering_experiment\small_conductivity\horizontal_scattering_0p5lambda.in" -n 1 -gpu
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    
    python -m gprMax "C:\Users\Administrator\Thesis\True_and_Baseline_Models\Horizontal_scattering_experiment\small_conductivity\horizontal_scattering_homogeneous.in" -n 1 -gpu
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    '''
    
    subprocess.run(["powershell", "-Command", ps_command], check=True)
    print("gprMax execution completed.")

def get_snapshots_data(snapshot_folder, snapshot_prefix, snapshot_indices):
    loaded_snapshots = []
    for snap_num in snapshot_indices:
        snapshot_path = os.path.join(snapshot_folder, f'{snapshot_prefix}{snap_num}.vti')
        if not os.path.exists(snapshot_path):
            print(f"Warning: {snapshot_path} not found. Skipping plot.")
            return []

        mesh = pv.read(snapshot_path)
        dims = mesh.dimensions

        field_name = 'E-field'
        data = mesh[field_name]
        if data.shape[1] == 3:
            Ez_data = data[:, 2]
        else:
            Ez_data = data.flatten()

        data_2d = Ez_data.reshape(dims[1] - 1, dims[0] - 1)
        loaded_snapshots.append((snap_num, data_2d))
    return loaded_snapshots

def do_plot(loaded_snapshots, domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time, save_filename, title):
    if not loaded_snapshots:
        return

    length = domain_width
    depth  = domain_height
    upper_layer_thickness = air_thickness
    fracture_y0 = domain_height - fracture_bottom
    fracture_y1 = domain_height - fracture_top

    early_data = [d for s, d in loaded_snapshots if s <= 4]
    mid_data = [d for s, d in loaded_snapshots if 5 <= s <= 6]
    late_data = [d for s, d in loaded_snapshots if s > 6]

    early_abs_max = max(np.abs(d).max() for d in early_data) if early_data else 1.0
    mid_pct_ref = max(np.percentile(np.abs(d), 99.2) for d in mid_data) if mid_data else early_abs_max
    late_pct_ref = max(np.percentile(np.abs(d), 99.7) for d in late_data) if late_data else early_abs_max

    mid_abs_max = min(mid_pct_ref, early_abs_max * 0.22)
    mid_abs_max = max(mid_abs_max, early_abs_max * 0.035)

    late_abs_max = min(late_pct_ref, early_abs_max * 0.32)
    late_abs_max = max(late_abs_max, early_abs_max * 0.06)

    n_plots = len(loaded_snapshots)
    ncols = 2
    nrows = int(np.ceil(n_plots / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 3.3 * nrows), sharex=True, sharey=True)
    axes = np.array(axes).reshape(-1)

    for i, (snap_num, data_2d) in enumerate(loaded_snapshots):
        ax = axes[i]

        if snap_num <= 4:
            vlim = early_abs_max
            scale_tag = 'regular scale'
        elif snap_num <= 6:
            vlim = mid_abs_max
            scale_tag = 'reflection-enhanced'
        else:
            vlim = late_abs_max
            scale_tag = 'reflection-enhanced'

        im = ax.imshow(
            data_2d,
            aspect='auto',
            cmap='seismic',
            origin='lower',
            vmin=-vlim,
            vmax=vlim,
            extent=[0, length, 0, depth],
        )

        time_ns = snap_num * snapshot_time * 1e9
        ax.set_title(f'Snapshot {snap_num} (t = {time_ns:.2f} ns) | {scale_tag}', fontsize=10, weight='bold')

        ax.add_patch(Rectangle((0, depth - upper_layer_thickness), length, upper_layer_thickness, facecolor='#cfe8ff', edgecolor='blue', alpha=0.2))
        ax.add_patch(Rectangle((0, fracture_y0), length, fracture_y1 - fracture_y0, facecolor='#ffe6e6', edgecolor='red', linewidth=1.0, alpha=0.2, label='Fracture' if i == 0 else ''))

        ax.set_xlabel('Length (m)', fontsize=10)
        ax.set_ylabel('Depth (m)', fontsize=10)

        if i == 0:
            ax.legend(loc='upper right', fontsize=9)

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
        cbar.set_label('Ez (V/m)', fontsize=8)
        cbar.ax.tick_params(labelsize=8)

    for j in range(n_plots, len(axes)):
        axes[j].axis('off')

    fig.suptitle(title, fontsize=14, weight='bold', y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    plt.savefig(save_filename, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Plot saved to: {save_filename}")

def plot_time_traces(out_alt, out_homo, n_blocks, block_width, rx_offset, save_filename, title):
    fig, ax = plt.subplots(figsize=(10, 8))
    try:
        with h5py.File(out_alt, 'r') as fa, h5py.File(out_homo, 'r') as fh:
            iterations = fa.attrs['Iterations']
            dt = fa.attrs['dt']
            time = np.arange(iterations) * dt * 1e9 # ns
            
            diff_traces = []
            for i in range(n_blocks):
                # The first receiver was rx1. The next blocks start at rx2.
                rx_name = f'rx{i+2}'
                try:
                    ez_alt = np.array(fa['rxs'][rx_name]['Ez'])
                    ez_homo = np.array(fh['rxs'][rx_name]['Ez'])
                    diff_traces.append(ez_alt - ez_homo)
                except KeyError:
                    pass
            
            if not diff_traces:
                return

            diff_traces = np.array(diff_traces)
            
            max_val = np.max(np.abs(diff_traces))
            if max_val == 0: max_val = 1
            scale = (block_width * 0.8) / max_val  # Amplify amplitudes
            
            for i, trace in enumerate(diff_traces):
                # We don't have the 54 * step offset anymore
                x_base = (i + 0.5) * block_width
                scaled_trace = x_base + trace * scale
                
                ax.plot(scaled_trace, time, 'k-', linewidth=0.8)
                ax.fill_betweenx(time, x_base, scaled_trace, where=(scaled_trace > x_base), facecolor='k', alpha=0.5)
                ax.axvline(x_base, color='k', linestyle=':', linewidth=0.5, alpha=0.3)
                
            ax.set_ylim(time[-1], time[0])
            ax.set_xlim(0, n_blocks * block_width)
            ax.set_xlabel('Length (m)', fontsize=12)
            ax.set_ylabel('Time (ns)', fontsize=12)
            ax.set_title(title, fontsize=14, weight='bold')
            
        plt.tight_layout()
        plt.savefig(save_filename, bbox_inches='tight', dpi=300)
    except Exception as e:
        print(f"Error plotting traces: {e}")
    finally:
        plt.close(fig)
        print(f"Plot saved to: {save_filename}")

def plot_snapshots(domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time, dx_dy_dz, n_blocks, block_width, rx_offset):
    print("Plotting snapshots...")
    snapshot_prefix = 'snapshot_mid_x_'
    snapshot_indices = list(range(1, 11))

    # 1. Alternating snapshots
    folder_alt = r'horizontal_scattering_0p5lambda_snaps'
    alt_data = get_snapshots_data(folder_alt, snapshot_prefix, snapshot_indices)
    do_plot(
        alt_data, domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time,
        'gpr_snapshots_result_alt.png', 'GPR Forward Modeling Snapshots (Alternating Block Fracture)'
    )

    # 2. Homogeneous snapshots
    folder_homo = r'horizontal_scattering_homogeneous_snaps'
    homo_data = get_snapshots_data(folder_homo, snapshot_prefix, snapshot_indices)
    do_plot(
        homo_data, domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time,
        'gpr_snapshots_result_homo.png', 'GPR Forward Modeling Snapshots (Homogeneous Fracture)'
    )

    # 3. Difference snapshots (Alternating - Homogeneous)
    if alt_data and homo_data:
        diff_data = []
        for (snap_num, data_alt), (_, data_homo) in zip(alt_data, homo_data):
            diff_data.append((snap_num, data_alt - data_homo))
        do_plot(
            diff_data, domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time,
            'gpr_snapshots_result_diff.png', 'GPR Forward Modeling Snapshots (Diff: Alternating - Homogeneous)'
        )
        
    # 4. Difference Time Traces (Alternating - Homogeneous) from .out files
    plot_time_traces(
        'horizontal_scattering_0p5lambda.out',
        'horizontal_scattering_homogeneous.out',
        n_blocks, block_width, rx_offset,
        'gpr_snapshots_result_diff_traces.png', 'GPR Difference Time Traces (Alternating - Homogeneous)'
    )

if __name__ == "__main__":
    import sys
    # Change to correct directory to match notebook behavior
    target_dir = r"C:\Users\Administrator\Thesis\True_and_Baseline_Models\Horizontal_scattering_experiment\small_conductivity"
    if os.path.exists(target_dir):
        os.chdir(target_dir)
        
    width, height, air_thick, f_top, f_bottom, snap_time, dx_dy_dz, n_blocks, block_width, rx_offset = create_input_file()
    run_gprmax()
    plot_snapshots(width, height, air_thick, f_top, f_bottom, snap_time, dx_dy_dz, n_blocks, block_width, rx_offset)
