import os
import numpy as np

# Import modules from our new files
from plotting import plot_snapshots
from processing_OO import GPRModelData
from plotting_OO import plot_wiggle_traces, plot_fk_image, plot_cwt_image, plot_migrated_image

def main():
    # Setup directory
    target_dir = r"C:\Users\Administrator\Thesis\Diffraction_Experiment"
    if os.path.exists(target_dir):
        os.chdir(target_dir)
        
    # Define model parameters to give you control
    model_params = {
        'f_central': 1.5 * 1e9,              # 1.5 GHz
        'c': 3 * 1e8,                        # speed of light in m/s
        'permittivity_ice': 6,
        'permittivity_air': 1,
        'permittivity_fracture': 6,
        'conductivity_ice': 1e-6,
        'conductivity_air': 0,
        'conductivity_fracture': 0.001,
        'fracture_depth': 0.3,
        'depth_below_fracture': 0.1,
        'air_thickness': 0.05,
        'rx_per_block': 4                    # Measure 4 times per spatial block
    }

    # ==========================
    # Object-Oriented Processing
    # ==========================
    # 1. Instantiate the Model Object WITH parameters
    gpr_model = GPRModelData(
        out_alt_file='horizontal_scattering_0p5lambda.out', 
        out_homo_file='horizontal_scattering_homogeneous.out',
        model_params=model_params
    )

    # 2. Run the simulation through the object
    gpr_model.run_simulation()

    # (Optional) Plot the full domain snapshots using the attributes stored in gpr_model
    plot_snapshots(
        gpr_model.width, gpr_model.height, gpr_model.air_thick, gpr_model.f_top, 
        gpr_model.f_bottom, gpr_model.snap_time, gpr_model.dx_dy_dz, 
        gpr_model.n_blocks, gpr_model.block_width, gpr_model.rx_offset
    )

    if len(gpr_model.time) == 0:
        print("Traces could not be loaded. Ensure the simulation generated .out files.")
        return

    # 3. Raw Difference Traces Plot
    plot_wiggle_traces(gpr_model, gpr_model.diff_traces,
                       "Difference Traces (Raw)", "plot_diff_raw.png")

    # 4. Apply SVD Filter
    n_comp_mute = 1 # Number of horizontal modes to drop
    svd_filtered_traces = gpr_model.apply_svd_filter(n_comp_mute)
    plot_wiggle_traces(gpr_model, svd_filtered_traces,
                       f"Difference Traces (SVD Muted = {n_comp_mute})", "plot_diff_svd.png")

    # 5. Apply Cross-Correlation
    cc_traces = gpr_model.cross_correlate()
    plot_wiggle_traces(gpr_model, cc_traces,
                       "Cross-Correlation (Avg Homo vs Diff)", "plot_diff_cc.png")

    # 6. Transform to F-K Domain
    # We apply this to the svd_filtered_traces to see the diffraction energy
    k_array, fk_freqs, fk_mag = gpr_model.fk_transform(svd_filtered_traces)
    if k_array is not None:
        plot_fk_image(k_array, fk_freqs, fk_mag,
                      "F-K Transform of Diffractions", "plot_diff_fk.png",
                      block_width=gpr_model.block_width)  # Pass block_width for k-lines!

    # 7. Continuous Wavelet Transform (CWT)
    cwt_freqs, cwt_image = gpr_model.compute_cwt_image()
    if cwt_freqs is not None:
        plot_cwt_image(gpr_model.time, cwt_freqs, cwt_image,
                       "Wavelet Transform of Global Avg Difference Trace", "plot_diff_cwt.png")

    # 8. Prestack Kirchhoff Depth Migration + Envelope
    # Geometry Setup
    v_ice = model_params['c'] / np.sqrt(model_params['permittivity_ice'])
    wavelength_fracture = (model_params['c'] / np.sqrt(model_params['permittivity_fracture'])) / model_params['f_central']
    source_receiver_steps = wavelength_fracture / 10
    x_first_measurement = 1/2 * wavelength_fracture
    tx_start_x = x_first_measurement + 54 * source_receiver_steps
    max_depth = 0.5 

    # E. Multi-Band Migration
    # Low-freq for mean fracture (using full alt field)
    alt_low_traces = gpr_model.bandpass_filter(gpr_model.alt_traces, lowcut=1e8, highcut=1e9)
    migrated_low, depths_low = gpr_model.migrate(
        alt_low_traces, tx_x=tx_start_x, velocity=v_ice, max_depth=max_depth
    )
    plot_migrated_image(gpr_model, gpr_model.apply_envelope(migrated_low), depths_low,
                        "Low-Freq Migration (Mean Fracture)", "plot_migrated_low.png")

    # High-freq for internal heterogeneity (using diffracted field)
    diff_high_traces = gpr_model.bandpass_filter(svd_filtered_traces, lowcut=1e9, highcut=4e9)
    migrated_high, depths_high = gpr_model.migrate(
        diff_high_traces, tx_x=tx_start_x, velocity=v_ice, max_depth=max_depth
    )
    plot_migrated_image(gpr_model, gpr_model.apply_envelope(migrated_high), depths_high,
                        "High-Freq Migration (Diffracted Interior)", "plot_migrated_high.png")

    # Full band diffraction migration
    migrated_img, depths = gpr_model.migrate(
        svd_filtered_traces, tx_x=tx_start_x, velocity=v_ice, max_depth=max_depth
    )
    
    # Envelope gives positive analytic magnitude
    migrated_envelope = gpr_model.apply_envelope(migrated_img)
    plot_migrated_image(gpr_model, migrated_envelope, depths,
                        "Depth Migrated Diffractions (Envelope)", "plot_diff_migrated.png")
                        
    # 9. Evaluate Lateral Tuning Theory Metrics & Spectral Inversion
    wavelength_bg = v_ice / model_params['f_central']
    eta = gpr_model.block_width / wavelength_bg
    diff_energy = gpr_model.calculate_diffraction_energy(svd_filtered_traces)
    
    # Inversion: Estimate block width (d) from peak k_x
    k_peak, peak_mag = gpr_model.extract_fk_peaks(k_array, fk_freqs, fk_mag, f_min=1e9, f_max=2.5e9)
    # The alternating blocks (+A, -A, +A, -A) form a spatial period of 2*d.
    # Therefore, the dominant spatial wavenumber is K = 2*pi / (2*d) = pi / d.
    d_est = (np.pi / k_peak) if k_peak else 0.0

    print("\n" + "="*50)
    print("LATERAL TUNING METRICS & SPECTRAL INVERSION")
    print("="*50)
    print(f"Background Wavelength (lambda) : {wavelength_bg:.4f} m")
    print(f"Internal Block Width (True d)  : {gpr_model.block_width:.4f} m")
    print(f"Lateral Tuning Parameter (eta) : {eta:.4f} (d / lambda)")
    print(f"Total Diffraction Energy       : {diff_energy:.5e}")
    if k_peak:
        print(f"Extracted Peak F-K (K_0)       : {k_peak:.4f} rad/m")
        print(f"Inverted Block Width (Est. d)  : {d_est:.4f} m")
        print(f"Inversion Error                : {abs(d_est - gpr_model.block_width)/gpr_model.block_width * 100:.1f}%")
    else:
        print("No valid spatial wavenumber peaks found.")
    print("="*50 + "\n")

if __name__ == "__main__":
    main()
