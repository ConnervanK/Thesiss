import os
import numpy as np

# Import modules from our new files
from utils.plotting import plot_snapshots
from utils.processing import GPRModelData
from utils.plotting import plot_wiggle_traces, plot_fk_image, plot_cwt_image, plot_migrated_image, plot_cwt_cross_sections, plot_xwt_phase_arrows, plot_aligned_traces

def main():
    # Setup directory 
    target_dir = r"C:\Users\Administrator\OneDrive\Thesis\Diffraction_Experiment"
    os.makedirs(os.path.join(target_dir, "data", "outputs"), exist_ok=True)
    os.chdir(target_dir)
        
    # Define model parameters to give you control
    model_params = {
        'f_central': 1.5 * 1e9,              # 1.5 GHz
        'c': 3 * 1e8,                        # speed of light in m/s
        'permittivity_ice': 3.15,
        'permittivity_air': 1,
        'permittivity_fracture': 10, #80,
        'conductivity_ice': 1e-6,
        'conductivity_air': 0,
        'conductivity_fracture': 1e-3, #1,
        'fracture_depth': 0.6,
        'depth_below_fracture': 0.1,
        'air_thickness': 0.1,
        
        # --- Mode selection: 'static' or 'bscan' ---
        'mode': 'static',                    # Switch to 'bscan' to run moving Tx-Rx array 
        'rx_per_block': 2,                   # Used only if mode == 'static'
        
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
        out_alt_file=r'configs\horizontal_scattering_0p5lambda.out', 
        out_homo_file=r'configs\horizontal_scattering_homogeneous.out',
        model_params=model_params
    )

    # 2. Run the simulation through the object
    # Force rerun so that modifications in forward.py take effect
    gpr_model.run_simulation(force_rerun=False)

    # (Optional) Plot the full domain snapshots using the attributes stored in gpr_model
    plot_snapshots(
        gpr_model.width, gpr_model.height, gpr_model.air_thick, gpr_model.f_top, 
        gpr_model.f_bottom, gpr_model.snap_time, gpr_model.dx_dy_dz, 
        gpr_model.n_blocks, gpr_model.block_width, gpr_model.rx_offset,
        rx_per_block=gpr_model.rx_per_block
    )

    if len(gpr_model.time) == 0:
        print("Traces could not be loaded. Ensure the simulation generated .out files.")
        return

    # 3. Raw Difference Traces Plot
    plot_wiggle_traces(gpr_model, gpr_model.diff_traces,
                       "Difference Traces (Raw Amplitude)", os.path.join("data", "outputs", "plot_diff_raw.png"))

    # 3.1 Normal Moveout (NMO) Correction
    # Get velocity of the background medium (ice)
    v_ice = model_params['c'] / np.sqrt(model_params['permittivity_ice'])
    # Assume source is at center of domain since tx_start_x is not dynamically fetched yet
    tx_start_x = getattr(gpr_model, 'tx_start_x', gpr_model.width / 2.0)
    
    # nmo_traces = gpr_model.apply_nmo(gpr_model.diff_traces, velocity=v_ice * 1.1, tx_x=tx_start_x, max_stretch=0.5)
    # plot_wiggle_traces(gpr_model, nmo_traces,
    #                    "Difference Traces (NMO Corrected)", os.path.join("data", "outputs", "plot_diff_nmo.png"))

    # phase_diff_traces = gpr_model.compute_instantaneous_phase(gpr_model.diff_traces)
    # plot_wiggle_traces(gpr_model, phase_diff_traces,
    #                    "Difference Traces (Instantaneous Phase)", os.path.join("data", "outputs", "plot_diff_phase.png"))

    # phase_alt_traces = gpr_model.compute_instantaneous_phase(gpr_model.alt_traces)
    # plot_wiggle_traces(gpr_model, phase_alt_traces,
    #                    "Alternating Traces (Instantaneous Phase)", os.path.join("data", "outputs", "plot_alt_phase.png"))

    # # 4. Apply SVD Filter
    # n_comp_mute = 1 # Number of horizontal modes to drop
    # svd_filtered_traces = gpr_model.apply_svd_filter(n_comp_mute)
    # plot_wiggle_traces(gpr_model, svd_filtered_traces,
    #                    f"Difference Traces (SVD Muted = {n_comp_mute})", os.path.join("data", "outputs", "plot_diff_svd.png"))

    # # 5. Apply Cross-Correlation
    # cc_traces = gpr_model.cross_correlate()
    # plot_wiggle_traces(gpr_model, cc_traces,
    #                    "Cross-Correlation (Avg Homo vs Diff)", os.path.join("data", "outputs", "plot_diff_cc.png"))

    auto_cc_traces = gpr_model.auto_correlate_diff()
    plot_wiggle_traces(gpr_model, auto_cc_traces,
                       "Auto-Correlation of Difference Traces", os.path.join("data", "outputs", "plot_diff_autocc.png"))

    # # 6. Transform to F-K Domain
    # # We apply this to the svd_filtered_traces to see the diffraction energy
    # k_array, fk_freqs, fk_mag = gpr_model.fk_transform(svd_filtered_traces)
    # if k_array is not None:
    #     plot_fk_image(k_array, fk_freqs, fk_mag,
    #                   "F-K Transform of Diffractions", os.path.join("data", "outputs", "plot_diff_fk.png"),
    #                   block_width=gpr_model.block_width)  # Pass block_width for k-lines!

    # 7. Continuous Wavelet Transform (CWT)
    cwt_freqs, cwt_diff_image = gpr_model.compute_cwt_image(traces=gpr_model.diff_traces)
    if cwt_freqs is not None:
        plot_cwt_image(gpr_model.time, cwt_freqs, cwt_diff_image,
                       "Wavelet Transform of Global Avg Difference Trace", os.path.join("data", "outputs", "plot_diff_cwt.png"))

    cwt_freqs_alt, cwt_alt_image = gpr_model.compute_cwt_image(traces=gpr_model.alt_traces)
    if cwt_freqs_alt is not None:
        plot_cwt_image(gpr_model.time, cwt_freqs_alt, cwt_alt_image,
                       "Wavelet Transform of Global Avg Alternating Trace", os.path.join("data", "outputs", "plot_alt_cwt.png"))

    cwt_freqs_diff_phase, cwt_diff_phase_image = gpr_model.compute_cwt_image(traces=gpr_model.diff_traces, return_phase=True)
    if cwt_freqs_diff_phase is not None:
        plot_cwt_image(gpr_model.time, cwt_freqs_diff_phase, cwt_diff_phase_image,
                       "Wavelet Transform Phase of Global Avg Difference Trace", os.path.join("data", "outputs", "plot_diff_cwt_phase.png"), cmap='hsv')
                       
        # Visualize explicit wave cycles bounded between -1 and 1 without phase jumps
        # Scaled by normalized CWT magnitude to hide phase where signal energy is zero
        mag_norm_diff = cwt_diff_image / (np.max(cwt_diff_image) + 1e-12)
        wave_cycles_diff = np.cos(cwt_diff_phase_image) * mag_norm_diff
        
        plot_cwt_image(gpr_model.time, cwt_freqs_diff_phase, wave_cycles_diff,
                       "Normalized Wave Cycles (Diff)", os.path.join("data", "outputs", "plot_diff_cwt_phase_unwrapped.png"), cmap='PuOr', vmin=-1, vmax=1)
        
        # Cross sections at 1.5 GHz (Central) and 2.5 GHz (Higher offset)
        plot_cwt_cross_sections(gpr_model.time, cwt_freqs_diff_phase, wave_cycles_diff,
                                f1=1.5e9, f2=2.5e9,
                                title=("Wave Cycle Evolution over Time (1.5GHz vs 2.5GHz)"),
                                save_filename=os.path.join("data", "outputs", "plot_diff_phase_cross_section.png"))
                       
        # Compute and plot the time derivative of the phase
        unwrapped_phase = np.unwrap(cwt_diff_phase_image, axis=1)
        dt = gpr_model.time[1] - gpr_model.time[0]
        phase_derivative = np.gradient(unwrapped_phase, dt, axis=1)
        
        plot_cwt_image(gpr_model.time, cwt_freqs_diff_phase, phase_derivative,
                       "Time Derivative of WT Phase (Instantaneous Frequency)", os.path.join("data", "outputs", "plot_diff_cwt_phase_derivative.png"))

        # Compute and plot the frequency derivative of the phase
        unwrapped_phase_freq = np.unwrap(cwt_diff_phase_image, axis=0)
        phase_derivative_freq = np.gradient(unwrapped_phase_freq, cwt_freqs_diff_phase, axis=0)
        
        plot_cwt_image(gpr_model.time, cwt_freqs_diff_phase, phase_derivative_freq,
                       "Frequency Derivative of WT Phase", os.path.join("data", "outputs", "plot_diff_cwt_phase_derivative_freq.png"))

        # Scaled Frequency Derivative of the phase (hides noise where signal is zero)
        phase_derivative_freq_scaled = phase_derivative_freq * mag_norm_diff
        plot_cwt_image(gpr_model.time, cwt_freqs_diff_phase, phase_derivative_freq_scaled,
                       "Scaled Frequency Derivative of WT Phase", os.path.join("data", "outputs", "plot_diff_cwt_phase_derivative_freq_scaled.png"))

    cwt_freqs_alt_phase, cwt_alt_phase_image = gpr_model.compute_cwt_image(traces=gpr_model.alt_traces, return_phase=True)
    if cwt_freqs_alt_phase is not None:
        plot_cwt_image(gpr_model.time, cwt_freqs_alt_phase, cwt_alt_phase_image,
                       "Wavelet Transform Phase of Global Avg Alternating Trace", os.path.join("data", "outputs", "plot_alt_cwt_phase.png"), cmap='hsv')
        
        # Visualize explicit wave cycles for alt trace
        mag_norm_alt = cwt_alt_image / (np.max(cwt_alt_image) + 1e-12)
        wave_cycles_alt = np.cos(cwt_alt_phase_image) * mag_norm_alt
        
        plot_cwt_image(gpr_model.time, cwt_freqs_alt_phase, wave_cycles_alt,
                       "Normalized Wave Cycles (Alt)", os.path.join("data", "outputs", "plot_alt_cwt_phase_unwrapped.png"), cmap='PuOr', vmin=-1, vmax=1)
        
        plot_cwt_cross_sections(gpr_model.time, cwt_freqs_alt_phase, wave_cycles_alt,
                                f1=1.5e9, f2=2.5e9,
                                title=("Wave Cycle Evolution over Time (Alt, 1.5GHz vs 2.5GHz)"),
                                save_filename=os.path.join("data", "outputs", "plot_alt_phase_cross_section.png"))

    # 7.1 Phase-Weighted Stacking (PWS) CWT
    pws_diff = gpr_model.phase_weighted_stack(traces=gpr_model.diff_traces, power=2)
    cwt_freqs_pws, cwt_pws_image = gpr_model.compute_cwt_image(traces=np.array([pws_diff]))
    if cwt_freqs_pws is not None:
        plot_cwt_image(gpr_model.time, cwt_freqs_pws, cwt_pws_image,
                       "Wavelet Transform of PWS Difference Stack (Power=2)", os.path.join("data", "outputs", "plot_diff_cwt_pws.png"))

    cwt_freqs_pws_phase, cwt_pws_phase_image = gpr_model.compute_cwt_image(traces=np.array([pws_diff]), return_phase=True)
    if cwt_freqs_pws_phase is not None:
        plot_cwt_image(gpr_model.time, cwt_freqs_pws_phase, cwt_pws_phase_image,
                       "Wavelet Transform Phase of PWS Difference Stack", os.path.join("data", "outputs", "plot_diff_cwt_pws_phase.png"))

    # # 7.2 3D CWT Example
    # # Resolves CWT for each individual trace to avoid average-loss. Array shape: (num_traces, num_freqs, num_time)
    # cwt_freqs_3d, cwt_3d_diff = gpr_model.compute_cwt_3d(traces=gpr_model.diff_traces)
    # if cwt_freqs_3d is not None and len(cwt_3d_diff) > 0:
    #     # Plotting the CWT for the very first trace as an example
    #     plot_cwt_image(gpr_model.time, cwt_freqs_3d, cwt_3d_diff[0],
    #                    "Wavelet Transform of Difference Trace 1", os.path.join("data", "outputs", "plot_diff_cwt_trace_1.png"))

    # # 7.2 3D CWT Example
    # # Resolves CWT for each individual trace to avoid average-loss. Array shape: (num_traces, num_freqs, num_time)
    # cwt_freqs_3d, cwt_3d_diff = gpr_model.compute_cwt_3d(traces=gpr_model.diff_traces)
    
    # if cwt_freqs_3d is not None and len(cwt_3d_diff) > 0:
    #     # Create a subfolder to store all individual trace CWT images
    #     cwt_output_dir = "cwt_3d_traces"
    #     os.makedirs(cwt_output_dir, exist_ok=True)
        
    #     # Loop through all available traces
    #     for i in range(len(cwt_3d_diff)):
    #         trace_num = i + 1
    #         output_filepath = os.path.join(cwt_output_dir, os.path.join("data", "outputs", f"plot_diff_cwt_trace_{trace_num}.png"))
            
    #         # Plot and save each trace's CWT into the subfolder
    #         plot_cwt_image(gpr_model.time, cwt_freqs_3d, cwt_3d_diff[i],
    #                        f"Wavelet Transform of Difference Trace {trace_num}", 
    #                        output_filepath)
            
    #     print(f"[{len(cwt_3d_diff)}] 3D CWT trace plots saved in the '{cwt_output_dir}' folder.")

    # 8. Shift-and-Correlate XWT Strategy
    print("Computing Geometric Time Shift and XWT...")
    
    # Define custom trace pairs to evaluate (e.g. adjacent and non-adjacent)
    num_diff_traces = len(gpr_model.diff_traces)
    trace_pairs_to_test = [(10, 11)] # Default adjacent example
    if num_diff_traces > 5:
        trace_pairs_to_test.append((50, 51)) # Another adjacent example
        trace_pairs_to_test.append((50, 60)) # Example of non-adjacent (gap 4)
    if num_diff_traces > 10:
        trace_pairs_to_test.append((50, 80)) # Example of further non-adjacent (gap 9)
        
    freqs_xwt, power_xwt, phase_xwt, shifted_windows, aligned_traces_full = gpr_model.run_shift_and_correlate(
        gpr_model.diff_traces, 
        depth=model_params['fracture_depth'], 
        velocity=v_ice, 
        tx_x=tx_start_x, 
        window_width=3e-9,
        trace_pairs=trace_pairs_to_test
    )
    
    if freqs_xwt is not None and len(power_xwt) > 0:
        for idx, (t1, t2) in enumerate(trace_pairs_to_test):
            # Power describes magnitude similarity, Phase arrows denote angular lead/lag
            plot_xwt_phase_arrows(gpr_model.time * 1e9, freqs_xwt, power_xwt[idx], phase_xwt[idx],
                           f"Cross-Wavelet Transform (XWT) - Traces {t1+1} & {t2+1}", os.path.join("data", "outputs", f"plot_diff_xwt_t{t1+1}_t{t2+1}.png"))
                           
            plot_aligned_traces(gpr_model.time, shifted_windows[idx][0], shifted_windows[idx][1],
                                f"Windowed - Traces {t1+1} & {t2+1}", os.path.join("data", "outputs", f"plot_diff_xwt_windowed_t{t1+1}_t{t2+1}.png"))
                                
            plot_aligned_traces(gpr_model.time, aligned_traces_full[idx][0], aligned_traces_full[idx][1],
                                f"Full Aligned - Traces {t1+1} & {t2+1}", os.path.join("data", "outputs", f"plot_diff_xwt_full_aligned_t{t1+1}_t{t2+1}.png"))

    # 8.5 Align all traces globally to visualize the flattened hyperbola reflection
    print("Aligning all traces geometrically to visualize the flattened hyperbola...")
    ref_idx = len(gpr_model.diff_traces) // 2
    best_velocity, best_score = gpr_model.estimate_best_alignment_velocity(
        gpr_model.diff_traces,
        depth=model_params['fracture_depth'],
        velocity_init=v_ice,
        tx_x=tx_start_x,
        reference_idx=ref_idx,
        window_width=3e-9,
        search_factors=(0.85, 1.15),
        n_trials=31,
    )
    print(f"Best alignment velocity: {best_velocity:.3e} m/s (coherence score: {best_score:.4f})")

    aligned_all = gpr_model.align_all_traces(
        gpr_model.diff_traces, 
        depth=model_params['fracture_depth'], 
        velocity=best_velocity,
        tx_x=tx_start_x, 
        reference_idx=ref_idx
    )
    plot_wiggle_traces(
        gpr_model, 
        aligned_all, 
        "All Geometrically Aligned Traces (Flat Hyperbola)", 
        os.path.join("data", "outputs", "plot_diff_aligned_bscan.png")
    )

    # 8.6 Fine-tune the alignment using empirical cross-correlation (perfectly flatten it)
    print("Running data-driven residual alignment to perfectly flatten the reflection...")
    # Identify the arrival time of the reference trace to window around the target reflection
    ref_offset = np.abs(gpr_model.get_rx_x_array()[ref_idx] - tx_start_x)
    ref_arrival_time = np.sqrt(ref_offset**2 + (4 * (model_params['fracture_depth']**2))) / best_velocity

    aligned_all_fine = gpr_model.residual_align_traces(
        aligned_all, 
        reference_idx=ref_idx,
        window_center=ref_arrival_time, 
        window_width=3e-9
    )
    plot_wiggle_traces(
        gpr_model, 
        aligned_all_fine, 
        "Fine-Tuned Geometrical & Residual Aligned Traces", 
        os.path.join("data", "outputs", "plot_diff_aligned_bscan_fine.png")
    )
                           
    # # 9. Prestack Kirchhoff Depth Migration + Envelope
    # # Geometry Setup
    # v_ice = model_params['c'] / np.sqrt(model_params['permittivity_ice'])
    # wavelength_fracture = (model_params['c'] / np.sqrt(model_params['permittivity_fracture'])) / model_params['f_central']
    # source_receiver_steps = wavelength_fracture / 10
    # x_first_measurement = 1/2 * wavelength_fracture
    # tx_start_x = x_first_measurement + 54 * source_receiver_steps
    # max_depth = 0.5 

    # # E. Multi-Band Migration
    # # Low-freq for mean fracture (using full alt field)
    # alt_low_traces = gpr_model.bandpass_filter(gpr_model.alt_traces, lowcut=1e8, highcut=1e9)
    # migrated_low, depths_low = gpr_model.migrate(
    #     alt_low_traces, tx_x=tx_start_x, velocity=v_ice, max_depth=max_depth
    # )
    # plot_migrated_image(gpr_model, gpr_model.apply_envelope(migrated_low), depths_low,
    #                     "Low-Freq Migration (Mean Fracture)", os.path.join("data", "outputs", "plot_migrated_low.png"))

    # # High-freq for internal heterogeneity (using diffracted field)
    # diff_high_traces = gpr_model.bandpass_filter(svd_filtered_traces, lowcut=1e9, highcut=4e9)
    # migrated_high, depths_high = gpr_model.migrate(
    #     diff_high_traces, tx_x=tx_start_x, velocity=v_ice, max_depth=max_depth
    # )
    # plot_migrated_image(gpr_model, gpr_model.apply_envelope(migrated_high), depths_high,
    #                     "High-Freq Migration (Diffracted Interior)", os.path.join("data", "outputs", "plot_migrated_high.png"))

    # # Full band diffraction migration
    # migrated_img, depths = gpr_model.migrate(
    #     svd_filtered_traces, tx_x=tx_start_x, velocity=v_ice, max_depth=max_depth
    # )
    
    # # Envelope gives positive analytic magnitude
    # migrated_envelope = gpr_model.apply_envelope(migrated_img)
    # plot_migrated_image(gpr_model, migrated_envelope, depths,
    #                     "Depth Migrated Diffractions (Envelope)", os.path.join("data", "outputs", "plot_diff_migrated.png"))
                        
    # # 9. Evaluate Lateral Tuning Theory Metrics & Spectral Inversion
    # wavelength_bg = v_ice / model_params['f_central']
    # eta = gpr_model.block_width / wavelength_bg
    # diff_energy = gpr_model.calculate_diffraction_energy(svd_filtered_traces)
    
    # # Inversion: Estimate block width (d) from peak k_x
    # k_peak, peak_mag = gpr_model.extract_fk_peaks(k_array, fk_freqs, fk_mag, f_min=1e9, f_max=2.5e9)
    # # The alternating blocks (+A, -A, +A, -A) form a spatial period of 2*d.
    # # Therefore, the dominant spatial wavenumber is K = 2*pi / (2*d) = pi / d.
    # d_est = (np.pi / k_peak) if k_peak else 0.0

    # # A. Diffraction Energy Calibration
    # # Sub-wavelength scattering energy scales ~ (d/\lambda)^a.
    # # Using existing empirical constant for this dataset setup.
    # energy_cal_factor = 280.0  # Assumed calibration constant for demonstration
    # d_est_energy = wavelength_bg * np.sqrt(diff_energy / energy_cal_factor)
    
    # # B. High-Resolution Spatial Spectral Estimation (MUSIC)
    # # Resolves spatial frequencies beyond the Rayleigh limit of standard FFTs.
    # k_music, music_spectrum = gpr_model.compute_music_spectrum(svd_filtered_traces, num_sources=3)
    # if k_music is not None:
    #     music_peak_idx = np.argmax(music_spectrum)
    #     k_music_peak = k_music[music_peak_idx]
    #     d_est_music = np.pi / k_music_peak
    # else:
    #     k_music_peak, d_est_music = 0.0, 0.0

    # # C. Spectral Centroid / Frequency Shift Analysis
    # # Compare the centroid frequency of the difference field to the background frequency.
    # # Sub-wavelength features act as high-pass or resonant filters (e.g., Rayleigh scattering ~ f^4)
    # centroid_homo = gpr_model.compute_spectral_centroid(gpr_model.homo_traces)
    # centroid_diff = gpr_model.compute_spectral_centroid(svd_filtered_traces)
    
    # # Simple empirical shift mapping
    # centroid_shift_ratio = centroid_diff / (centroid_homo + 1e-9)
    # # Target central freq is around 1.5e9. Assume scaling relationship. Placeholder conversion.
    # # Actual relations depend on forward modeling calibration.
    # d_est_centroid = wavelength_bg * (1.0 / centroid_shift_ratio) if centroid_diff > 0 else 0.0

    # # D. Amplitude Versus Offset (AVO)
    # # Extract peak reflection amplitudes across the receiver array
    # rx_x, avo_amps, offsets = gpr_model.extract_avo(svd_filtered_traces, tx_x=tx_start_x)
    # if rx_x is not None and len(offsets) > 1:
    #     # Calculate AVO gradient (simplified linear fit of amplitude vs offset)
    #     # Small 'd' scatters more uniformly (isotropic, low gradient), large 'd' reflects specularly (sharp decay).
    #     p = np.polyfit(offsets, avo_amps, 1)
    #     avo_gradient = p[0]
    #     # Empirical conversion placeholder
    #     d_est_avo = wavelength_bg * np.exp(avo_gradient * 50) 
    # else:
    #     avo_gradient = 0.0
    #     d_est_avo = 0.0

    # print("\n" + "="*50)
    # print("LATERAL TUNING METRICS & SPECTRAL INVERSION")
    # print("="*50)
    # print(f"Background Wavelength (lambda) : {wavelength_bg:.4f} m")
    # print(f"Internal Block Width (True d)  : {gpr_model.block_width:.4f} m")
    # print(f"Lateral Tuning Parameter (eta) : {eta:.4f} (d / lambda)")
    # print(f"Total Diffraction Energy       : {diff_energy:.5e}")
    # if k_peak:
    #     print(f"Extracted Peak F-K (K_0)       : {k_peak:.4f} rad/m")
    #     print(f"Inverted Block Width (Est. d)  : {d_est:.4f} m")
    #     print(f"F-K Inversion Error            : {abs(d_est - gpr_model.block_width)/gpr_model.block_width * 100:.1f}%")
        
    # print("-"*50)
    # print(f"Energy-Calibrated Block Width  : {d_est_energy:.4f} m")
    # print(f"Energy Inversion Error         : {abs(d_est_energy - gpr_model.block_width)/gpr_model.block_width * 100:.1f}%")
    
    # if k_music_peak > 0.0:
    #     print("-"*50)
    #     print(f"MUSIC Spectrum Peak (K_0)      : {k_music_peak:.4f} rad/m")
    #     print(f"MUSIC Inverted Block Width     : {d_est_music:.4f} m")
    #     print(f"MUSIC Inversion Error          : {abs(d_est_music - gpr_model.block_width)/gpr_model.block_width * 100:.1f}%")
    # else:
    #     print("No valid spatial wavenumber peaks found.")
        
    # print("-"*50)
    # print(f"Homogeneous Spectral Centroid  : {centroid_homo/1e9:.4f} GHz")
    # print(f"Diffraction Spectral Centroid  : {centroid_diff/1e9:.4f} GHz")
    # print(f"Centroid-Inverted Block Width  : {d_est_centroid:.4f} m")
    # print(f"Centroid Inversion Error       : {abs(d_est_centroid - gpr_model.block_width)/gpr_model.block_width * 100:.1f}%")

    # print("-"*50)
    # print(f"AVO Gradient (Amp/m)           : {avo_gradient:.5e}")
    # print(f"AVO-Inverted Block Width       : {d_est_avo:.4f} m")
    # print(f"AVO Inversion Error            : {abs(d_est_avo - gpr_model.block_width)/gpr_model.block_width * 100:.1f}%")
        
    # print("="*50 + "\n")

if __name__ == "__main__":
    main()
