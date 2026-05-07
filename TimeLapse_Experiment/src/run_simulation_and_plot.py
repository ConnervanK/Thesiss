import os
import numpy as np

# Import modules from our new files
from utils.plotting import plot_snapshots
from utils.processing import GPRModelData
from utils.plotting import plot_wiggle_traces, plot_fk_image, plot_cwt_image, plot_migrated_image, plot_cwt_cross_sections, plot_xwt_phase_arrows, plot_aligned_traces, plot_csd
import forward

def main():
    # Setup directory 
    target_dir = r"C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Experiment"
    os.makedirs(os.path.join(target_dir, "data", "outputs"), exist_ok=True)
    os.chdir(target_dir)
        
    # Define model parameters to give you control
    model_params = {
        'f_central': 1.5 * 1e9,              # 1.5 GHz
        'c': 3 * 1e8,                        # speed of light in m/s
        'permittivity_ice': 3.15,
        'permittivity_air': 1,
        'permittivity_fracture': 80, #80,
        'conductivity_ice': 1e-6,
        'conductivity_air': 0,
        'conductivity_fracture': 0.01, #1,
        'fracture_depth': 0.6,
        'depth_below_fracture': 0.1,
        'air_thickness': 0.1,
        
        # --- Mode selection: 'static' or 'bscan' ---
        'mode': 'static',                    # Switch to 'bscan' to run moving Tx-Rx array 
        'rx_per_block': 1,                   # Used only if mode == 'static'
        
        # --- B-scan parameters (used if mode == 'bscan') ---
        'rx_count': 1,                       # Number of receivers in moving array
        'rx_spacing': 0.02,                  # Metres between receivers in the array
        'bscan_traces': 20,                  # Number of traces for the B-scan
        'bscan_step_x': 0.02                 # Movement step size in metres
        ,
        # --- Block size control ---
        'block_size': None,                  # Absolute block size (m). If None uses fraction below.
        'block_size_fraction': 1/6,          # Fraction of ice wavelength used when block_size is None
        'time_lapse_shift': 0.5              # Lateral shift (m) for time-lapse variant; must be < wavelength_ice
    }

    # ==========================
    # Object-Oriented Processing
    # ==========================
    # 0. Generate input files and run forward model (optional). Compute block size.
    # Compute ice wavelength for block fraction calculation
    c = model_params['c']
    f_central = model_params['f_central']
    permittivity_ice = model_params['permittivity_ice']
    c_ice = c / np.sqrt(permittivity_ice)
    wavelength_ice = c_ice / f_central

    # Determine block_size in meters
    block_size = model_params.get('block_size', None)
    if block_size is None:
        block_size = model_params.get('block_size_fraction', 1/8) * wavelength_ice

    # Allow user to set a small time-lapse shift (meters)
    time_lapse_shift = model_params.get('time_lapse_shift', 0.1)

    # Call forward.create_input_file to write .in files (and optional time-lapse .in)
    call_params = model_params.copy()
    call_params['block_size'] = block_size
    call_params['time_lapse_shift'] = time_lapse_shift
    output_context = forward.create_input_file(**call_params)

    # Run gprMax and move outputs into a labeled subfolder for this block size
    block_label = output_context.get('block_label', None)
    if block_label is None:
        block_label = f"{int(round(block_size*1e6))}um"
    
    # For static mode: one measurement with fixed Tx and array of Rx -> traces_count=1
    # For bscan mode: moving Tx-Rx pair -> traces_count=bscan_traces
    if call_params.get('mode') == 'static':
        traces_count = 1
    else:
        traces_count = call_params.get('bscan_traces', 1)
    forward.run_gprmax_to_subdir(traces_count=traces_count, output_subdir=block_label)

    # 1. Instantiate the Model Object WITH parameters (pointing to output files in labeled subfolder)
    out_alt_path = os.path.join('data', 'outputs', block_label, 'horizontal_scattering_0p5lambda.out')
    out_homo_path = os.path.join('data', 'outputs', block_label, 'horizontal_scattering_homogeneous.out')

    gpr_model = GPRModelData(
        out_alt_file=out_alt_path,
        out_homo_file=out_homo_path,
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
        rx_per_block=gpr_model.rx_per_block,
        output_dir=os.path.join('data', 'outputs', block_label)
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
    
    auto_cc_traces = gpr_model.auto_correlate_diff()
    plot_wiggle_traces(gpr_model, auto_cc_traces,
                       "Auto-Correlation of Difference Traces", os.path.join("data", "outputs", "plot_diff_autocc.png"))

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
            # Only plot if this index was actually computed
            if idx >= len(power_xwt):
                print(f"Skipping traces {t1+1} & {t2+1}: not computed (only {len(power_xwt)} pairs available)")
                continue
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
    
    # 10. Cross Spectral Density
    print("Computing Cross Spectral Density on two selected traces...")
    t1_idx = 20
    t2_idx = 30
    if len(aligned_all_fine) > 11:
        fs_hz = 1.0 / gpr_model.dt
        plot_csd(
            aligned_all_fine[t1_idx], 
            aligned_all_fine[t2_idx], 
            fs_hz, 
            f"Cross Spectral Density: Trace {t1_idx+1} vs {t2_idx+1}", 
            os.path.join("data", "outputs", f"plot_diff_csd_t{t1_idx+1}_t{t2_idx+1}.png"),
            source_freq_hz=model_params['f_central']
        )
                           

if __name__ == "__main__":
    main()
