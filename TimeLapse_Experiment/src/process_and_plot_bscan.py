import os
import numpy as np
from utils.processing_bscan import GPRBScanData
from utils.plotting import plot_bscan_section, plot_migrated_image, plot_csd, plot_xwt_phase_arrows


def main():
    target_dir = r"c:\Users\Administrator\OneDrive\Thesis\TimeLapse_Experiment"
    os.chdir(target_dir)

    block_label = "1_8Lambda"

    baseline_file = os.path.join(target_dir, "data", "outputs", block_label,
                                 "horizontal_scattering_baseline_bscan.out")
    timelapse_file = os.path.join(target_dir, "data", "outputs", block_label,
                                  "horizontal_scattering_timelapse_bscan.out")

    # Physical constants — must match run_experiments_bscan.py
    f_central = 1.5e9
    c = 3e8
    permittivity_ice = 3.15
    c_ice = c / np.sqrt(permittivity_ice)

    air_thickness = 0.1
    fracture_depth = 0.6
    permittivity_fracture_minus = 80
    wavelength_minus = (c / np.sqrt(permittivity_fracture_minus)) / f_central
    min_wavelength = min((c / np.sqrt(1.0)) / f_central, wavelength_minus)
    thickness_fracture = min_wavelength / 5
    depth_below_fracture = 0.1

    # B-scan geometry — must match run_experiments_bscan.py
    rx_count = 50
    rx_spread = 2.0
    tx_rx_offset = 0.1
    midpoint_start_x = 2.0 - rx_spread / 2.0      # = 1.0
    midpoint_spacing = rx_spread / max(1, rx_count - 1)

    print(f"Loading B-scan data from {block_label}...")
    bscan = GPRBScanData(
        baseline_file=baseline_file,
        timelapse_file=timelapse_file,
        midpoint_start_x=midpoint_start_x,
        midpoint_spacing=midpoint_spacing,
        n_traces=rx_count,
        tx_rx_offset=tx_rx_offset,
    )
    bscan.load_data()

    # Sync count with what was actually loaded
    rx_count = bscan.n_traces
    midpoint_spacing = rx_spread / max(1, rx_count - 1)
    bscan.midpoint_spacing = midpoint_spacing

    out_dir = os.path.join("data", "outputs", block_label)
    midpoints = bscan.get_midpoint_x_array()

    # 1. B-scan sections
    diff_traces = bscan.subtract()

    plot_bscan_section(
        bscan.baseline_traces, bscan.time, midpoints,
        "Baseline B-scan", os.path.join(out_dir, "bscan_baseline.png"),
        fracture_depth=fracture_depth, velocity=c_ice, tx_rx_offset=tx_rx_offset,
    )
    plot_bscan_section(
        bscan.timelapse_traces, bscan.time, midpoints,
        "Time-Lapse B-scan", os.path.join(out_dir, "bscan_timelapse.png"),
        fracture_depth=fracture_depth, velocity=c_ice, tx_rx_offset=tx_rx_offset,
    )
    plot_bscan_section(
        diff_traces, bscan.time, midpoints,
        "Difference B-scan", os.path.join(out_dir, "bscan_difference.png"),
        fracture_depth=fracture_depth, velocity=c_ice, tx_rx_offset=tx_rx_offset,
    )

    # 2. Kirchhoff migration of baseline and timelapse B-scans
    print("Running Kirchhoff migration on baseline and timelapse B-scans...")
    max_migration_depth = fracture_depth + thickness_fracture + depth_below_fracture
    migrated_baseline, z_array = bscan.kirchhoff_migration_bscan(
        bscan.baseline_traces,
        velocity=c_ice,
        max_depth=max_migration_depth,
        dz=0.005,
    )
    plot_migrated_image(
        migrated_baseline, z_array,
        rx_x_array=midpoints,
        tx_x=midpoints[rx_count // 2],
        title="Kirchhoff Migration — Baseline B-scan",
        save_filename=os.path.join(out_dir, "bscan_kirchhoff_baseline.png"),
        fracture_depth=fracture_depth,
    )

    migrated_timelapse, _ = bscan.kirchhoff_migration_bscan(
        bscan.timelapse_traces,
        velocity=c_ice,
        max_depth=max_migration_depth,
        dz=0.005,
    )
    plot_migrated_image(
        migrated_timelapse, z_array,
        rx_x_array=midpoints,
        tx_x=midpoints[rx_count // 2],
        title="Kirchhoff Migration — Time-Lapse B-scan",
        save_filename=os.path.join(out_dir, "bscan_kirchhoff_timelapse.png"),
        fracture_depth=fracture_depth,
    )

    migrated_diff = migrated_baseline - migrated_timelapse
    plot_migrated_image(
        migrated_diff, z_array,
        rx_x_array=midpoints,
        tx_x=midpoints[rx_count // 2],
        title="Kirchhoff Migration — Difference (Baseline − Time-Lapse)",
        save_filename=os.path.join(out_dir, "bscan_kirchhoff_difference.png"),
        fracture_depth=fracture_depth,
    )

    # 4. CSD between two symmetric traces around the profile centre
    print("Computing Cross Spectral Density...")
    ref_idx = rx_count // 2
    t1_idx = max(0, ref_idx - 5)
    t2_idx = min(rx_count - 1, ref_idx + 5)
    if len(diff_traces) > max(t1_idx, t2_idx):
        assert bscan.dt is not None
        fs = 1.0 / bscan.dt
        plot_csd(
            diff_traces[t1_idx],
            diff_traces[t2_idx],
            fs,
            "Cross Spectral Density — Symmetric B-scan Traces",
            os.path.join(out_dir, "bscan_csd.png"),
            source_freq_hz=f_central,
        )

    # 5. XWT phase between the same two traces
    print("Computing Cross Wavelet Transform Phase...")
    if len(diff_traces) > max(t1_idx, t2_idx):
        freqs_xwt, power_xwt, phase_xwt = bscan.compute_xwt(
            diff_traces[t1_idx], diff_traces[t2_idx]
        )
        assert bscan.time is not None
        plot_xwt_phase_arrows(
            bscan.time, freqs_xwt, power_xwt, phase_xwt,
            "XWT Phase — B-scan",
            os.path.join(out_dir, "bscan_xwt_phase.png"),
        )

    print("B-scan process and plot complete.")


if __name__ == "__main__":
    main()
