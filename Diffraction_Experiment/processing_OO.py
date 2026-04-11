import os
import numpy as np
import h5py
from scipy.signal import hilbert, butter, filtfilt, find_peaks
import pywt
from forward import clear_results, create_input_file, run_gprmax

class GPRModelData:
    """
    Object-oriented wrapper for a GPR forward modeling result.
    It encapsulates the data loading, difference tracing, and all processing steps.
    """
    def __init__(self, out_alt_file, out_homo_file, model_params=None):
        self.out_alt_file = out_alt_file
        self.out_homo_file = out_homo_file
        self.model_params = model_params
        
        self.n_blocks = None
        self.block_width = None
        self.rx_offset = None
        self.rx_per_block = 1
        
        self.time = None
        self.dt = None
        self.diff_traces = None
        self.homo_traces = None
        self.alt_traces = None

    def run_simulation(self):
        """Runs the full GPRMax pipeline and loads the resulting data."""
        if self.model_params is None:
            raise ValueError("No model parameters specified.")
        
        if os.path.exists(self.out_alt_file) and os.path.exists(self.out_homo_file):
            print(f"Loading existing `{self.out_alt_file}` and `{self.out_homo_file}`, skipping simulation.")
            out_ctx = create_input_file(**self.model_params)
        else:
            clear_results()
            out_ctx = create_input_file(**self.model_params)
            run_gprmax(traces_count=out_ctx['actual_traces'])

        # Store attributes
        self.width = out_ctx['domain_width']
        self.height = out_ctx['domain_height']
        self.air_thick = out_ctx['air_thickness']
        self.f_top = out_ctx['fracture_top']
        self.f_bottom = out_ctx['fracture_bottom']
        self.snap_time = out_ctx['snapshot_time']
        self.dx_dy_dz = out_ctx['dx_dy_dz']
        self.n_blocks = out_ctx['n_blocks']
        self.block_width = out_ctx['block_width']
        self.rx_offset = out_ctx['rx_start_x'] - out_ctx['tx_start_x']
        self.rx_per_block = out_ctx.get('rx_per_block', 1)
        self.mode = out_ctx.get('mode', 'static')
        
        self.total_rx = out_ctx['total_rx']
        self.actual_traces = out_ctx['actual_traces']

        self._load_data()

    def _load_data(self):
        with h5py.File(self.out_alt_file, 'r') as fa, h5py.File(self.out_homo_file, 'r') as fh:
            iterations = fa.attrs['Iterations']
            self.dt = fa.attrs['dt']
            self.time = np.arange(iterations) * self.dt * 1e9 # in ns
            
            diff_traces = []
            homo_traces = []
            alt_traces = []
            
            # The receivers in the output are simply numbered rx1, rx2... up to total_rx
            for i in range(self.total_rx):
                rx_name = f'rx{i+1}'
                # Inside 'static' mode, we used to skip the rx1 because rx1 was positioned at source location. 
                # Our rewritten logic outputs the exact array receivers natively, so we just start from rx1!
                try:
                    ez_alt = np.array(fa['rxs'][rx_name]['Ez'])
                    ez_homo = np.array(fh['rxs'][rx_name]['Ez'])
                    
                    if self.actual_traces > 1:
                        # B-scan matrix: shape (time, traces). 
                        # We transpose so Shape: (traces, time), which handles as an ensemble of classical trace arrays.
                        ez_alt = ez_alt.T
                        ez_homo = ez_homo.T
                        
                        # In the single rx case (bscan with 1 source, 1 rx), 'diff_traces' 
                        # is logically 'ez_alt_transposed', behaving like a single pseudo-spatial array 
                        # just like 'static'.
                        if self.total_rx == 1:
                            alt_traces = ez_alt
                            homo_traces = ez_homo
                            diff_traces = ez_alt - ez_homo
                        else:
                            alt_traces.append(ez_alt)
                            homo_traces.append(ez_homo)
                            diff_traces.append(ez_alt - ez_homo)
                    else:
                        # shape is (time,), just append
                        alt_traces.append(ez_alt)
                        homo_traces.append(ez_homo)
                        diff_traces.append(ez_alt - ez_homo)
                except KeyError as e:
                    print(f"Warning: {rx_name} not found in output files.")
                    pass
                    
            self.diff_traces = np.array(diff_traces)
            self.homo_traces = np.array(homo_traces)
            self.alt_traces = np.array(alt_traces)
            
    def get_rx_x_array(self):
        """Returns the x-coordinates of the receiver array."""
        if hasattr(self, 'mode') and self.mode == 'bscan':
            # If a B-scan, the 'x_array' for plotting AVO/Migrations corresponds to the moving array steps.
            bscan_step_x = self.model_params.get('bscan_step_x', 0.05)
            # Adjust offset relative to start 
            start_x = self.model_params.get('tx_start_x', 0.0) 
            return np.array([start_x + i * bscan_step_x for i in range(self.actual_traces)])
        else:
            rx_dx = self.block_width / self.rx_per_block
            return np.array([(i + 0.5) * rx_dx for i in range(self.total_rx)])

    def apply_svd_filter(self, n_components_to_mute=1):
        """Mutes the first N singular components of the difference traces."""
        traces = self.diff_traces
        if len(traces) == 0: return traces
        U, S, Vh = np.linalg.svd(traces, full_matrices=False)
        S_filtered = S.copy()
        num_to_mute = min(n_components_to_mute, len(S_filtered))
        if num_to_mute > 0:
            S_filtered[:num_to_mute] = 0.0
        return U @ np.diag(S_filtered) @ Vh

    def apply_envelope(self, traces):
        """Applies Hilbert transform to calculate the positive amplitude envelope."""
        if len(traces) == 0: return traces
        return np.abs(hilbert(traces, axis=1))

    def calculate_diffraction_energy(self, traces):
        """Calculates the total energy (L2 norm squared) of the diffracted field traces."""
        return np.sum(np.square(traces))
        
    def bandpass_filter(self, traces, lowcut, highcut, order=4):
        """Applies a zero-phase Butterworth bandpass filter to the traces."""
        if len(traces) == 0: return traces
        nyq = 0.5 / self.dt
        low = lowcut / nyq
        high = highcut / nyq
        # Handle cases where highcut is above Nyquist safely
        if high >= 1.0: high = 0.99
        
        b, a = butter(order, [low, high], btype='band')
        filtered_traces = filtfilt(b, a, traces, axis=1)
        return filtered_traces

    def extract_fk_peaks(self, k, freqs, fk_mag, f_min=1e9, f_max=2e9):
        """
        Extracts the dominant lateral wavenumber (k_x) ridge in a specific frequency band.
        This provides an estimate of the grating spatial frequency K_0 = 2*pi/d.
        """
        if fk_mag is None: return None, None
        
        # Limit to target frequency band
        f_idx = np.where((freqs >= f_min) & (freqs <= f_max))[0]
        if len(f_idx) == 0: return None, None
        
        # Average fk magnitude over the frequency band to find spatial dominant peaks
        avg_fk_k = np.mean(fk_mag[:, f_idx], axis=1)
        
        # Find peaks in the spatial wavenumber domain
        peaks, _ = find_peaks(avg_fk_k, prominence=np.max(avg_fk_k)*0.1)
        
        if len(peaks) > 0:
            # Sort peaks by amplitude
            sorted_peaks = sorted(peaks, key=lambda p: avg_fk_k[p], reverse=True)
            # Find the peak corresponding to positive k (skip dc k=0)
            for p in sorted_peaks:
                if k[p] > 0.1: # Skip exactly zero
                    return k[p], avg_fk_k[p]
        return None, None
        
    def cross_correlate(self):
        """Cross correlates the average homogeneous trace with each difference trace."""
        if len(self.homo_traces) == 0 or len(self.diff_traces) == 0: return []
        avg_homo = np.mean(self.homo_traces, axis=0)
        cc_traces = []
        for diff_trace in self.diff_traces:
            cc = np.correlate(diff_trace, avg_homo, mode='same')
            cc_traces.append(cc)
        return np.array(cc_traces)

    def compute_cwt_image(self, wavelet='cmor1.5-1.0', freqs=None):
        """Computes Continuous Wavelet Transform (CWT) magnitude of the global average difference trace."""
        if len(self.diff_traces) == 0: return None, None
        avg_trace = np.mean(self.diff_traces, axis=0) 
        
        if freqs is None:
            freqs = np.linspace(0.1e9, 3.5e9, 100) 
            
        sampling_freq = 1.0 / self.dt
        scales = pywt.frequency2scale(wavelet, freqs / sampling_freq)
        coefs, _ = pywt.cwt(avg_trace, scales, wavelet, sampling_period=self.dt)
        return freqs, np.abs(coefs)

    def compute_spectral_centroid(self, traces):
        """
        Computes the average spectral centroid (frequency focus) of the given traces.
        """
        if len(traces) == 0: return 0.0
        
        # Compute FFT along the time axis (axis=1)
        spectra = np.abs(np.fft.rfft(traces, axis=1))
        freqs = np.fft.rfftfreq(traces.shape[1], d=self.dt)
        
        # Calculate centroid for each trace
        centroids = np.sum(freqs * spectra, axis=1) / (np.sum(spectra, axis=1) + 1e-12)
        
        # Return the mean centroid across all traces
        return np.mean(centroids)

    def extract_avo(self, traces, tx_x=None):
        """
        Extracts the maximum amplitude of the envelope for each trace, representing
        Amplitude vs Offset.
        Returns: (rx_x_array, amplitudes, offsets)
        """
        if len(traces) == 0: return None, None, None
        rx_x_array = self.get_rx_x_array()
        
        # Get envelope of traces
        env_traces = self.apply_envelope(traces)
        
        # Maximum amplitude per trace
        amplitudes = np.max(env_traces, axis=1)
        
        if tx_x is None:
            tx_x = getattr(self, 'width', 0.5) / 2.0
            
        # Ensure sizes match
        num_traces = len(traces)
        if len(rx_x_array) > num_traces:
            rx_x_array = rx_x_array[:num_traces]
            
        offsets = np.abs(rx_x_array - tx_x)
        return rx_x_array, amplitudes, offsets

    def fk_transform(self, traces):
        """
        Transforms traces (x, t) to Frequency-Wavenumber (F-K) domain.
        Returns: (wavenumbers, frequencies, FK_magnitude_image)
        """
        if len(traces) == 0: return None, None, None
        
        # traces shape: (num_rx, num_t)
        fk_data = np.fft.fft2(traces)
        fk_data = np.fft.fftshift(fk_data)
        
        # Frequencies (axis 1 of traces -> columns of FFT) in Hz
        freqs = np.fft.fftshift(np.fft.fftfreq(traces.shape[1], self.dt))
        # Wavenumbers (axis 0 of traces -> rows of FFT) in rad/m
        rx_dx = self.block_width / self.rx_per_block
        k = np.fft.fftshift(np.fft.fftfreq(traces.shape[0], rx_dx)) * 2 * np.pi
        
        fk_mag = np.abs(fk_data)
        return k, freqs, fk_mag

    def compute_music_spectrum(self, traces, num_sources=2):
        """
        Computes the high-resolution spatial spectrum using the MUSIC algorithm.
        This resolves sub-wavelength block spacing much better than standard F-K transform.
        """
        if len(traces) == 0: return None, None
        
        # Calculate spatial covariance matrix of the traces
        # traces shape: (num_rx, num_t)
        # R shape: (num_rx, num_rx)
        R = traces @ traces.T / traces.shape[1]
        
        # Eigenvalue decomposition
        eigenvalues, eigenvectors = np.linalg.eigh(R)
        
        # Sort eigenvalues and eigenvectors in descending order
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        
        # Extract noise subspace
        En = eigenvectors[:, num_sources:]
        
        # Wavenumbers to search (rad/m) - focus on sub-wavelength scale
        k_search = np.linspace(10, 500, 2000)
        
        rx_x_array = self.get_rx_x_array()
        music_spectrum = np.zeros_like(k_search)
        
        for i, k in enumerate(k_search):
            # Steering vector (assuming plane waves or far-field)
            a = np.exp(-1j * k * rx_x_array)
            # P_MUSIC = 1 / (a^H * En * En^H * a)
            den = (a.conj().T @ En) @ (En.conj().T @ a)
            music_spectrum[i] = 1.0 / np.abs(den + 1e-10)
            
        # Normalize spectrum
        music_spectrum = music_spectrum / np.max(music_spectrum)
        
        return k_search, music_spectrum

    def migrate(self, traces, tx_x, velocity, max_depth, dz=0.002):
        """
        Prestack Depth Migration for a common shot gather, collapsing hyperbolas using Kirchhoff summation.
        """
        rx_x_array = self.get_rx_x_array()
        dt_ns = self.time[1] - self.time[0] if len(self.time) > 1 else 1.0
        z_array = np.arange(0, max_depth + dz, dz)
        
        migrated_image = np.zeros((len(z_array), len(rx_x_array)))
        
        for iz, z in enumerate(z_array):
            if z == 0: continue
            for ix_out, x_out in enumerate(rx_x_array):
                dist_tx = np.sqrt((x_out - tx_x)**2 + z**2)
                for irx, rx_x in enumerate(rx_x_array):
                    if irx >= traces.shape[0]:
                        continue
                    
                    dist_rx = np.sqrt((x_out - rx_x)**2 + z**2)
                    
                    t_total_sec = (dist_tx + dist_rx) / velocity
                    t_total_ns = t_total_sec * 1e9
                    
                    t_idx = int(np.round((t_total_ns - self.time[0]) / dt_ns))
                    if 0 <= t_idx < len(self.time):
                        migrated_image[iz, ix_out] += traces[irx, t_idx]
                        
        return migrated_image, z_array