import numpy as np
from scipy.signal import hilbert
import pywt

def apply_svd_filter(traces, n_components_to_mute=1):
    """
    Performs SVD on the time traces (receivers x time_steps).
    Mutes the first N singular components and reconstructs the image.
    This is often used to remove strong horizontal events (like direct waves).
    """
    if len(traces) == 0:
        return traces
        
    # Perform SVD
    # traces is expected to be a 2D numpy array (num_receivers x time_steps)
    U, S, Vh = np.linalg.svd(traces, full_matrices=False)
    
    # Mute the strongest singular values
    S_filtered = S.copy()
    num_to_mute = min(n_components_to_mute, len(S_filtered))
    if num_to_mute > 0:
        S_filtered[:num_to_mute] = 0.0
        
    # Reconstruct the traces
    reconstructed_traces = U @ np.diag(S_filtered) @ Vh
    return reconstructed_traces

def envelope(traces):
    """
    Applies Hilbert transform to calculate the instantaneous amplitude (envelope) of the traces.
    traces: 2D array (num_receivers, num_time_steps)
    Returns: 2D array of the same shape containing strictly positive envelope values.
    """
    if len(traces) == 0:
        return traces
    analytic_signal = hilbert(traces, axis=1) # Transform along the time axis
    amplitude_envelope = np.abs(analytic_signal)
    return amplitude_envelope

def kirchhoff_migration(traces, time_array, rx_x_array, tx_x, velocity, max_depth, dz=0.005):
    """
    Prestack Depth Migration for a common shot gather.
    traces: 2D array (num_receivers, num_time_steps)
    time_array: 1D array of time (seconds)
    rx_x_array: 1D array of receiver x-coordinates (meters)
    tx_x: float, transmitter x-coordinate (meters)
    velocity: propagation velocity (meters/second)
    max_depth: float, maximum depth to migrate (meters)
    dz: vertical step size for output (meters)
    Returns:
        migrated_image: 2D array (depth_steps, rx_x_steps)
        z_array: 1D array of depths
    """
    dt = time_array[1] - time_array[0]
    z_array = np.arange(0, max_depth + dz, dz)
    
    num_rx = len(rx_x_array)
    num_out_x = num_rx
    x_out_array = rx_x_array # We migrate to the same horizontal positions as the receivers
    
    migrated_image = np.zeros((len(z_array), num_out_x))
    
    # We do a simple diffraction summation (Kirchhoff)
    for iz, z in enumerate(z_array):
        if z == 0:
            continue
        for ix_out, x_out in enumerate(x_out_array):
            # For this image point (x_out, z), calculate distance from TX
            dist_tx = np.sqrt((x_out - tx_x)**2 + z**2)
            
            # Loop over all receivers (traces) and sum the corresponding travel time amplitude
            for irx, rx_x in enumerate(rx_x_array):
                dist_rx = np.sqrt((x_out - rx_x)**2 + z**2)
                
                # Total travel time from tx -> diffractor -> rx
                t_total_sec = (dist_tx + dist_rx) / velocity
                
                # time_array is in nanoseconds, so we convert t_total to nanoseconds
                t_total_ns = t_total_sec * 1e9
                
                # Find the index in the time trace
                t_idx = int(np.round((t_total_ns - time_array[0]) / dt))
                
                if 0 <= t_idx < len(time_array):
                    # Add the trace amplitude to the output image point
                    migrated_image[iz, ix_out] += traces[irx, t_idx]
                    
    return migrated_image, z_array

def calculate_cross_correlation(homo_traces, diff_traces):
    """
    Takes the average trace of the homogeneous results and cross correlates it 
    with the difference traces. Returns the correlation outputs.
    """
    if len(homo_traces) == 0 or len(diff_traces) == 0:
        return []
        
    # Get the average trace of the homogeneous results
    avg_homo_trace = np.mean(homo_traces, axis=0)
    
    # Cross correlate the average homogeneous trace with each difference trace
    cc_traces = []
    for diff_trace in diff_traces:
        # Cross correlate and keep the 'same' length array
        cc = np.correlate(diff_trace, avg_homo_trace, mode='same')
        cc_traces.append(cc)
        
    return np.array(cc_traces)

def compute_cwt_image(traces, dt, wavelet='cmor1.5-1.0', freqs=None):
    """
    Computes the Continuous Wavelet Transform (CWT) of the traces.
    To show a 2D image of Time vs Frequency, this usually averages the CWT 
    magnitude of all traces.
    
    Returns: freq_array, cwt_image (magnitude)
    """
    if len(traces) == 0:
        return None, None
        
    avg_trace = np.mean(traces, axis=0) # Take the global average trace to represent the survey
    
    if freqs is None:
        # Frequencies from 0.1 GHz to 5 GHz
        freqs = np.linspace(0.1e9, 3.5e9, 100) 
        
    # Calculate scales for CWT corresponding to target frequencies
    sampling_freq = 1.0 / dt
    scales = pywt.frequency2scale(wavelet, freqs / sampling_freq)
    
    # Compute CWT
    coefs, freqs_out = pywt.cwt(avg_trace, scales, wavelet, sampling_period=dt)
    
    # The output is complex for Morlet, so we take magnitude
    cwt_image = np.abs(coefs)
    
    return freqs, cwt_image

