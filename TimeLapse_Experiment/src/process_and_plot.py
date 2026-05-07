import os
import numpy as np
from utils.processing import GPRModelData
from utils.plotting import plot_wiggle_traces, plot_csd, plot_xwt_phase_arrows
from utils.plotting import plot_aligned_traces, get_snapshots_data, do_plot

def main():
    target_dir = r"c:\Users\Administrator\OneDrive\Thesis\TimeLapse_Experiment"
    os.chdir(target_dir)

    block_label = "1_4Lambda"
    
    baseline_file = os.path.join(target_dir, "data", "outputs", block_label, "horizontal_scattering_baseline.out")
    timelapse_file = os.path.join(target_dir, "data", "outputs", block_label, "horizontal_scattering_timelapse.out")

    print(f"Loading data from {block_label}...")
    
    # Re-compute EXACT parameters from forward.py to ensure plots align centrally:
    f_central = 1.5e9
    c = 3e8
    permittivity_ice = 3.15
    c_ice = c / np.sqrt(permittivity_ice)
    wavelength_ice = c_ice / f_central
    rx_spread = 2.0
    domain_width_approx = max(12 * wavelength_ice, rx_spread + 2.0)
    
    air_thickness = 0.1
    fracture_depth = 0.6
    permittivity_fracture_plus = 1.0
    permittivity_fracture_minus = 80
    wavelength_plus = (c / np.sqrt(permittivity_fracture_plus)) / f_central
    wavelength_minus = (c / np.sqrt(permittivity_fracture_minus)) / f_central
    min_wavelength = min(wavelength_plus, wavelength_minus)
    thickness_fracture = min_wavelength / 5
    depth_below_fracture = 0.1
    domain_height_approx = air_thickness + fracture_depth + thickness_fracture + depth_below_fracture
    
    dx_dy_dz = min_wavelength / 30
    domain_width = np.ceil(domain_width_approx / dx_dy_dz) * dx_dy_dz
    domain_height = np.ceil(domain_height_approx / dx_dy_dz) * dx_dy_dz

    # Use the correctly snapped bounds matching forward.py
    tx_start_x = domain_width / 2
    rx_count = 50
    rx_start_x = tx_start_x - (rx_spread / 2.0)
    rx_spacing = rx_spread / max(1, rx_count - 1)
    
    gpr_model = GPRModelData(
        baseline_file=baseline_file,
        timelapse_file=timelapse_file,
        tx_start_x=tx_start_x,
        rx_start_x=rx_start_x,
        rx_spacing=rx_spacing,
        rx_count=rx_count
    )
    
    gpr_model.load_data()

    out_plot_dir = os.path.join("data", "outputs", block_label)
    
    # 1. Visualize Traces (Baseline, Timelapse, and Difference)
    diff_traces = gpr_model.subtract()
    
    plot_wiggle_traces(gpr_model, gpr_model.baseline_traces, 
                       "Baseline Traces", os.path.join(out_plot_dir, "baseline_traces.png"),
                       domain_width=domain_width)
    plot_wiggle_traces(gpr_model, gpr_model.timelapse_traces, 
                       "Time-Lapse Traces", os.path.join(out_plot_dir, "timelapse_traces.png"),
                       domain_width=domain_width)
    plot_wiggle_traces(gpr_model, diff_traces, 
                       "Difference Traces", os.path.join(out_plot_dir, "difference_traces.png"),
                       domain_width=domain_width)

    # 2. Trace Alignment
    v_ice = c_ice
    ref_idx = rx_count // 2
    
    print("Aligning traces...")
    best_velocity, _ = gpr_model.estimate_best_alignment_velocity(
        diff_traces,
        depth=fracture_depth,
        velocity_init=v_ice,
        tx_x=tx_start_x,
        reference_idx=ref_idx,
        window_width=3e-9,
        search_factors=(0.85, 1.15),
        n_trials=31,
    )

    aligned_diff = gpr_model.align_all_traces(
        diff_traces, 
        depth=fracture_depth, 
        velocity=best_velocity,
        tx_x=tx_start_x, 
        reference_idx=ref_idx
    )
    
    plot_wiggle_traces(gpr_model, aligned_diff, 
                       "Aligned Difference Traces", os.path.join(out_plot_dir, "aligned_diff_traces.png"),
                       domain_width=domain_width)

    # 3. Cross Spectral Density
    print("Computing Cross Spectral Density...")
    t1_idx = ref_idx - 5
    t2_idx = ref_idx + 5
    
    if len(aligned_diff) > max(t1_idx, t2_idx):
        fs = 1.0 / gpr_model.dt
        plot_csd(
            aligned_diff[t1_idx], 
            aligned_diff[t2_idx],
            fs,
            "Cross Spectral Density (Symmetric Traces)",
            os.path.join(out_plot_dir, "csd_symmetric.png"),
            source_freq_hz=1.5e9
        )
    
    # 4. XWT Phase
    print("Computing Cross Wavelet Transform Phase...")
    if len(aligned_diff) > max(t1_idx, t2_idx):
        freqs_xwt, power_xwt, phase_xwt = gpr_model.compute_xwt(
            aligned_diff[t1_idx],
            aligned_diff[t2_idx]
        )
        plot_xwt_phase_arrows(
            gpr_model.time, 
            freqs_xwt, 
            power_xwt, 
            phase_xwt, 
            "XWT Phase Analysis", 
            os.path.join(out_plot_dir, "xwt_phase.png")
        )

    # 5. Plot Snapshots
    print("Plotting Snapshots...")
    snapshot_time = 0.25e-9
    snapshot_prefix = 'snapshot_mid_x_'
    snapshot_indices = list(range(1, 37))
    n_blocks = int(np.round(domain_width / (0.25 * wavelength_ice)))
    if n_blocks < 1: n_blocks = 1
    if n_blocks % 2 != 0: n_blocks += 1
    block_width = domain_width / n_blocks

    fracture_top = air_thickness + fracture_depth
    fracture_bottom = fracture_top + thickness_fracture

    baseline_snaps_dir = os.path.join(out_plot_dir, "horizontal_scattering_baseline_snaps")
    baseline_snaps_data = get_snapshots_data(baseline_snaps_dir, snapshot_prefix, snapshot_indices)
    if baseline_snaps_data:
        do_plot(
            baseline_snaps_data, domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time,
            os.path.join(out_plot_dir, "baseline_snapshots.png"), "Baseline GPR Snapshots",
            n_blocks=n_blocks, block_width=block_width
        )
    
    timelapse_snaps_dir = os.path.join(out_plot_dir, "horizontal_scattering_timelapse_snaps")
    timelapse_snaps_data = get_snapshots_data(timelapse_snaps_dir, snapshot_prefix, snapshot_indices)
    if timelapse_snaps_data:
        do_plot(
            timelapse_snaps_data, domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time,
            os.path.join(out_plot_dir, "timelapse_snapshots.png"), "Time-Lapse GPR Snapshots",
            n_blocks=n_blocks, block_width=block_width
        )
        
        # Difference Snapshots
        diff_snaps_data = []
        for i in range(len(baseline_snaps_data)):
            snap_num = baseline_snaps_data[i][0]
            b_data = baseline_snaps_data[i][1]
            t_data = timelapse_snaps_data[i][1]
            diff_data = t_data - b_data
            diff_snaps_data.append((snap_num, diff_data))
        
        do_plot(
            diff_snaps_data, domain_width, domain_height, air_thickness, fracture_top, fracture_bottom, snapshot_time,
            os.path.join(out_plot_dir, "difference_snapshots.png"), "Difference GPR Snapshots",
            n_blocks=n_blocks, block_width=block_width, is_diff=True
        )

    print("Process and Plot complete.")

if __name__ == "__main__":
    main()
