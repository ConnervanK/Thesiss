import numpy as np
from scipy.signal.windows import kaiser
from numpy.fft import ifft

def gpr_calc(data, n_freq=8192, freq_min=0.0, freq_max=14.0e9):
    freq = np.linspace(300e3, 14e9, n_freq)

    # Create a mask for the frequencies within the desired range
    ind_ = (freq >= freq_min) & (freq <= freq_max)

    # Filter the frequency array
    freq = freq[ind_]
    data = data[ind_]

    n_freq = len(data)
    ov = 40

    nf = ov * n_freq
    # gpr = np.zeros(nf)
    # time = np.zeros(nf)

    # Apply Kaiser window to the signal
    sig_f = data * kaiser(n_freq, 6)

    # Calculate frequency difference
    df = freq[1] - freq[0]

    # Zero-pad the signal
    sig_f2 = np.zeros(nf, dtype=complex)
    sig_f2[:n_freq] = sig_f
    sig_f = sig_f2

    # Perform inverse FFT
    h_t = ifft(sig_f)*len(sig_f)

    # Time step and time vector
    dt = 1 / (nf * df)
    time = np.arange(nf) * dt

    # Compute GPR values
    gpr = np.abs(2 * ov * h_t)

    return gpr, time, freq

def hrtr_method_fin(s11, freq, dsr, nt, time_min=0.0e-9, time_max=5.0e-9 + 0.01e-9, dt=0.01e-9):

    data_ds = s11[::dsr]
    freq_ds = freq[::dsr]
    freq_ds = freq_ds.reshape(-1, 1)

    time_tr = np.arange(time_min, time_max, dt)

    L = len(data_ds)
    M = L // 2
    N = L - M + 1

    R_ = np.zeros((N, N), dtype=complex)

    r = data_ds
    for j0 in range(1, M + 1):
        H_temp = r[j0 - 1:j0 + N - 1]
        temp = np.outer(H_temp, H_temp.conj())
        R_ += temp

    R_ /= M

    U, S, _ = np.linalg.svd(R_, full_matrices=False)

    En = U[:, nt + 1:]
    En_EnH = En @ En.conj().T
    p_music3 = np.zeros(len(time_tr), dtype=complex)


    for j0, t_sample in enumerate(time_tr):

        av = np.exp(-1j * 2.0 * np.pi * freq_ds[:N] * t_sample)

        p_music3[j0] = ((np.conj(av).T @ av) / (np.conj(av).T @ En_EnH @ av)).item()

    return p_music3, time_tr

if __name__ == '__main__':
    import os
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend to avoid Qt hanging issues
    import matplotlib.pyplot as plt
    from processing_OO import GPRModelData

    # Ensure we are in the correct directory
    target_dir = r"C:\Users\Administrator\Thesis\Diffraction_Experiment"
    if os.path.exists(target_dir):
        os.chdir(target_dir)

    model_params = {
        'f_central': 1.5 * 1e9,              
        'c': 3 * 1e8,                        
        'permittivity_ice': 6,
        'permittivity_air': 1,
        'permittivity_fracture': 6,
        'conductivity_ice': 1e-6,
        'conductivity_air': 0,
        'conductivity_fracture': 0.001,
        'fracture_depth': 0.3,
        'depth_below_fracture': 0.1,
        'air_thickness': 0.05,
        'mode': 'bscan',          
        'rx_count': 1,            
        'rx_spacing': 0.02,       
        'bscan_traces': 20,       
        'bscan_step_x': 0.02 
    }

    print("Initializing GPRModelData...")
    gpr_model = GPRModelData(
        out_alt_file='horizontal_scattering_0p5lambda.out', 
        out_homo_file='horizontal_scattering_homogeneous.out',
        model_params=model_params
    )

    print("Loading simulation data...")
    gpr_model.run_simulation()

    print("Applying HRTR...")
    diff_bscan = gpr_model.diff_traces
    dt = gpr_model.dt
    time_array = gpr_model.time * 1e-9 # Convert from ns back to seconds for HRTR
    
    hrtr_bscan = []
    
    for i in range(diff_bscan.shape[0]):
        trace = diff_bscan[i]
        
        # Compute Frequency Spectrum (S11) and Frequencies
        # The trace is real, compute rfft
        s11 = np.fft.fft(trace)
        freq = np.fft.fftfreq(len(trace), d=dt)
        
        # Use positive frequencies only and limit the bandwidth around the center frequency
        pos_mask = (freq >= 0.5e9) & (freq <= 3.0e9)
        s11 = s11[pos_mask]
        freq = freq[pos_mask]

        dsr = 1  # Downsampling rate 
        nt = 1   # Expected number of targets (noise space thresholding)
        
        p_music3, time_tr = hrtr_method_fin(
            s11, 
            freq, 
            dsr, 
            nt, 
            time_min=0.0e-9, 
            time_max=time_array[-1], 
            dt=dt
        )
        
        hrtr_bscan.append(np.abs(p_music3))

    hrtr_bscan = np.array(np.log(hrtr_bscan).T) # Transpose for plotting
    
    # # Normalize the amplitude
    # hrtr_bscan = hrtr_bscan / np.max(hrtr_bscan)

    print("Plotting results...")
    plt.figure(figsize=(10, 6))
    plt.imshow(hrtr_bscan, aspect='auto', extent=[0, diff_bscan.shape[0], time_array[-1]*1e9, 0], cmap='viridis')
    plt.colorbar(label='Log(HRTR Amplitude) [dB]')
    plt.title('HRTR Applied to B-scan')
    plt.xlabel('Trace Number')
    plt.ylabel('Time (ns)')
    plt.tight_layout()
    plt.savefig('hrtr_results.png', dpi=300)
    print("Plot saved to hrtr_results.png")
    
    # Using 'Agg' backend, plt.show() is omitted.

