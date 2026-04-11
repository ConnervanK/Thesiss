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
        
        # --- Mode selection: 'static' or 'bscan' ---
        'mode': 'static',                    # Switch to 'bscan' to run moving Tx-Rx array 
        'rx_per_block': 1,                   # Used only if mode == 'static'
        
        # --- B-scan parameters (used if mode == 'bscan') ---
        'rx_count': 1,                       # Number of receivers in moving array
        'rx_spacing': 0.02,                  # Metres between receivers in the array
        'bscan_traces': 20,                  # Number of traces for the B-scan
        'bscan_step_x': 0.02                 # Movement step size in metres
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

    # A. Diffraction Energy Calibration
    # Sub-wavelength scattering energy scales ~ (d/\lambda)^a.
    # Using existing empirical constant for this dataset setup.
    energy_cal_factor = 280.0  # Assumed calibration constant for demonstration
    d_est_energy = wavelength_bg * np.sqrt(diff_energy / energy_cal_factor)
    
    # B. High-Resolution Spatial Spectral Estimation (MUSIC)
    # Resolves spatial frequencies beyond the Rayleigh limit of standard FFTs.
    k_music, music_spectrum = gpr_model.compute_music_spectrum(svd_filtered_traces, num_sources=3)
    if k_music is not None:
        music_peak_idx = np.argmax(music_spectrum)
        k_music_peak = k_music[music_peak_idx]
        d_est_music = np.pi / k_music_peak
    else:
        k_music_peak, d_est_music = 0.0, 0.0

    # C. Spectral Centroid / Frequency Shift Analysis
    # Compare the centroid frequency of the difference field to the background frequency.
    # Sub-wavelength features act as high-pass or resonant filters (e.g., Rayleigh scattering ~ f^4)
    centroid_homo = gpr_model.compute_spectral_centroid(gpr_model.homo_traces)
    centroid_diff = gpr_model.compute_spectral_centroid(svd_filtered_traces)
    
    # Simple empirical shift mapping
    centroid_shift_ratio = centroid_diff / (centroid_homo + 1e-9)
    # Target central freq is around 1.5e9. Assume scaling relationship. Placeholder conversion.
    # Actual relations depend on forward modeling calibration.
    d_est_centroid = wavelength_bg * (1.0 / centroid_shift_ratio) if centroid_diff > 0 else 0.0

    # D. Amplitude Versus Offset (AVO)
    # Extract peak reflection amplitudes across the receiver array
    rx_x, avo_amps, offsets = gpr_model.extract_avo(svd_filtered_traces, tx_x=tx_start_x)
    if rx_x is not None and len(offsets) > 1:
        # Calculate AVO gradient (simplified linear fit of amplitude vs offset)
        # Small 'd' scatters more uniformly (isotropic, low gradient), large 'd' reflects specularly (sharp decay).
        p = np.polyfit(offsets, avo_amps, 1)
        avo_gradient = p[0]
        # Empirical conversion placeholder
        d_est_avo = wavelength_bg * np.exp(avo_gradient * 50) 
    else:
        avo_gradient = 0.0
        d_est_avo = 0.0

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
        print(f"F-K Inversion Error            : {abs(d_est - gpr_model.block_width)/gpr_model.block_width * 100:.1f}%")
        
    print("-"*50)
    print(f"Energy-Calibrated Block Width  : {d_est_energy:.4f} m")
    print(f"Energy Inversion Error         : {abs(d_est_energy - gpr_model.block_width)/gpr_model.block_width * 100:.1f}%")
    
    if k_music_peak > 0.0:
        print("-"*50)
        print(f"MUSIC Spectrum Peak (K_0)      : {k_music_peak:.4f} rad/m")
        print(f"MUSIC Inverted Block Width     : {d_est_music:.4f} m")
        print(f"MUSIC Inversion Error          : {abs(d_est_music - gpr_model.block_width)/gpr_model.block_width * 100:.1f}%")
    else:
        print("No valid spatial wavenumber peaks found.")
        
    print("-"*50)
    print(f"Homogeneous Spectral Centroid  : {centroid_homo/1e9:.4f} GHz")
    print(f"Diffraction Spectral Centroid  : {centroid_diff/1e9:.4f} GHz")
    print(f"Centroid-Inverted Block Width  : {d_est_centroid:.4f} m")
    print(f"Centroid Inversion Error       : {abs(d_est_centroid - gpr_model.block_width)/gpr_model.block_width * 100:.1f}%")

    print("-"*50)
    print(f"AVO Gradient (Amp/m)           : {avo_gradient:.5e}")
    print(f"AVO-Inverted Block Width       : {d_est_avo:.4f} m")
    print(f"AVO Inversion Error            : {abs(d_est_avo - gpr_model.block_width)/gpr_model.block_width * 100:.1f}%")
        
    print("="*50 + "\n")

if __name__ == "__main__":
    main()
