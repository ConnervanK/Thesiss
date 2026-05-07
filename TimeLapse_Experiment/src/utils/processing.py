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

    def run_simulation(self, force_rerun=False):
        """Runs the full GPRMax pipeline and loads the resulting data."""
        if self.model_params is None:
            raise ValueError("No model parameters specified.")
        
        if not force_rerun and os.path.exists(self.out_alt_file) and os.path.exists(self.out_homo_file):
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
            domain_width = self.width
            dx_dy_dz = self.dx_dy_dz
            valid_x = []
            for i in range(self.n_blocks * self.rx_per_block):
                rx_x_init = (i + 0.5) * rx_dx
                if rx_x_init > dx_dy_dz and rx_x_init < domain_width - dx_dy_dz:
                    valid_x.append(rx_x_init)
            
            # Match the length to what was actually loaded (in case partial loading happened)
            return np.array(valid_x)[:len(self.diff_traces)]

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

    def compute_instantaneous_phase(self, traces):
        """Applies Hilbert transform to calculate the instantaneous phase."""
        if len(traces) == 0: return traces
        return np.angle(hilbert(traces, axis=1))

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

    def auto_correlate_diff(self):
        """Autocorrelates each difference trace with itself."""
        if len(self.diff_traces) == 0: return []
        cc_traces = []
        for diff_trace in self.diff_traces:
            cc = np.correlate(diff_trace, diff_trace, mode='same')
            cc_traces.append(cc)
        return np.array(cc_traces)

    def compute_cwt_image(self, traces=None, wavelet='cmor1.5-1.0', freqs=None, return_phase=False):
        """Computes Continuous Wavelet Transform (CWT) magnitude or phase of the global average trace."""
        if traces is None:
            traces = self.diff_traces
        if len(traces) == 0: return None, None
        avg_trace = np.mean(traces, axis=0) 
        
        if freqs is None:
            freqs = np.linspace(0.1e9, 5.5e9, 100)
            
        sampling_freq = 1.0 / self.dt
        scales = pywt.frequency2scale(wavelet, freqs / sampling_freq)
        coefs, _ = pywt.cwt(avg_trace, scales, wavelet, sampling_period=self.dt)
        if return_phase:
            phase = np.angle(coefs)
            # Return raw unscaled phase so 'np.unwrap' can mathematically detect the exact 2*pi jumps
            return freqs, phase
        return freqs, np.abs(coefs)

    def compute_cwt_3d(self, traces=None, wavelet='cmor1.5-1.0', freqs=None, return_phase=False):
        """
        Computes the CWT for each individual trace, forming a 3D array: 
        (num_traces, num_freqs, num_time)
        """
        if traces is None:
            traces = self.diff_traces
        if len(traces) == 0: return None, None
        
        if freqs is None:
            freqs = np.linspace(0.1e9, 5.5e9, 100) 
            
        sampling_freq = 1.0 / self.dt
        scales = pywt.frequency2scale(wavelet, freqs / sampling_freq)
        
        cwt_3d = []
        for trace in traces:
            coefs, _ = pywt.cwt(trace, scales, wavelet, sampling_period=self.dt)
            if return_phase:
                phase = np.angle(coefs)
                cwt_3d.append(phase)
            else:
                cwt_3d.append(np.abs(coefs))
                
        return freqs, np.array(cwt_3d)

    def phase_weighted_stack(self, traces=None, power=2):
        """
        Performs Phase-Weighted Stacking (PWS) on an ensemble of traces.
        The linear stack is weighted by the instantaneous phase coherence.
        Returns a single 1D trace.
        """
        if traces is None:
            traces = self.diff_traces
        if len(traces) == 0: return None
        
        # 1. Compute linear stack (average trace)
        linear_stack = np.mean(traces, axis=0)
        
        # 2. Compute instantaneous phase for all traces
        analytic_signal = hilbert(traces, axis=1)
        instantaneous_phase = np.angle(analytic_signal)
        
        # 3. Compute phase coherence c(t) = |(1/N) * sum(e^{i*phi(t)})|
        phase_coherence = np.abs(np.mean(np.exp(1j * instantaneous_phase), axis=0))
        
        # 4. Weight the stack by coherence
        pws = linear_stack * (phase_coherence ** power)
        return pws

    def apply_nmo(self, traces, velocity, tx_x=None, max_stretch=1.0):
        """
        Applies Normal Moveout (NMO) correction to the given traces.
        Corrects for the travel-time delay due to the offset between Tx and Rx.
        """
        if traces is None or len(traces) == 0:
            return traces
            
        nmo_traces = np.zeros_like(traces)
        
        rx_x_array = self.get_rx_x_array()
        if len(rx_x_array) > len(traces):
            rx_x_array = rx_x_array[:len(traces)]
            
        if tx_x is None:
            # If not provided, assume Tx is in the center or calculate it based on out_ctx
            tx_x = getattr(self, 'tx_start_x', self.width / 2.0)
            
        # Physical offset for each trace
        offsets = np.abs(rx_x_array - tx_x)
        
        # Time axis
        t = self.time
        dt = t[1] - t[0] if len(t) > 1 else 1e-11
        
        for i in range(len(traces)):
            x = offsets[i]
            trace = traces[i]
            
            # The NMO equation: t^2 = t0^2 + (x/v)^2
            # Note: t is in nanoseconds (ns), so x/v (seconds) must be multiplied by 1e9 to match units!
            # To find the amplitude at t0, we evaluate the trace at t = sqrt(t0^2 + (x/v)^2)
            t_squared = t**2 + ((x / velocity) * 1e9)**2
            t_lookup = np.sqrt(t_squared)
            
            # Interpolate the trace at calculated times
            nmo_trace = np.interp(t_lookup, t, trace, left=0.0, right=0.0)
            
            # NMO stretch Mute: stretch = (t_lookup - t) / t
            t_safe = np.where(t == 0, 1e-12, t)
            stretch = (t_lookup - t) / t_safe
            nmo_trace[stretch > max_stretch] = 0.0
            
            nmo_traces[i] = nmo_trace
            
        return nmo_traces

    def frequency_domain_shift(self, trace, delta_t):
        """
        Shifts a 1D array `trace` in the time domain by `delta_t` (in seconds)
        using the Fourier Shift Theorem to allow for sub-sample shifts.
        """
        N = len(trace)
        F = np.fft.rfft(trace)
        freqs = np.fft.rfftfreq(N, d=self.dt)
        # Shift theorem: F(w) * exp(-j * 2 * pi * f * delta_t)
        F_shifted = F * np.exp(-1j * 2 * np.pi * freqs * delta_t)
        return np.fft.irfft(F_shifted, n=N)

    def extract_window(self, trace, center_time, window_width):
        """
        Extracts a window of a trace around a center time.
        Any points outside the window are set to zero, or we could taper it.
        """
        # smooth taper (Tukey window) or a simple rect window
        from scipy.signal.windows import tukey
        window = np.zeros_like(trace)
        start_idx = max(0, int((center_time - window_width/2) / self.dt))
        end_idx = min(len(trace), int((center_time + window_width/2) / self.dt))
        
        if end_idx > start_idx:
            # apply tukey alpha 0.2 to smooth edges
            taper = tukey(end_idx - start_idx, alpha=0.2)
            window[start_idx:end_idx] = trace[start_idx:end_idx] * taper
            
        return window

    def compute_xwt(self, trace1, trace2, wavelet='cmor1.5-1.0', freqs=None):
        """
        Computes the Cross-Wavelet Transform (XWT) between two traces.
        """
        if freqs is None:
            freqs = np.linspace(0.1e9, 5.5e9, 100)
            
        sampling_freq = 1.0 / self.dt
        scales = pywt.frequency2scale(wavelet, freqs / sampling_freq)
        
        coefs1, _ = pywt.cwt(trace1, scales, wavelet, sampling_period=self.dt)
        coefs2, _ = pywt.cwt(trace2, scales, wavelet, sampling_period=self.dt)
        
        # Cross wavelet transform
        xwt = coefs1 * np.conj(coefs2)
        power = np.abs(xwt)
        phase = np.angle(xwt)
        
        return freqs, power, phase

    def align_all_traces(self, traces, depth, velocity, tx_x=None, reference_idx=0):
        """
        Shifts all traces geometrically in the frequency domain so that the reflection
        from a horizontal fracture at given depth becomes perfectly horizontal, 
        aligned to the arrival time of the reference trace.
        """
        if len(traces) == 0: return traces
        if tx_x is None:
            tx_x = getattr(self, 'tx_start_x', self.width / 2.0)
            
        rx_array = self.get_rx_x_array()
        if len(rx_array) > len(traces):
            rx_array = rx_array[:len(traces)]
            
        # Physical offset for each trace
        offsets = np.abs(rx_array - tx_x)
        
        # Expected travel time for a horizontal reflector at given depth
        t_arrivals = np.sqrt(offsets**2 + (4 * (depth**2))) / velocity
        
        target_t = t_arrivals[reference_idx]
        
        aligned_traces = []
        for i in range(len(traces)):
            # delta_tau is how much later trace i arrives compared to reference
            delta_tau = t_arrivals[i] - target_t
            # Shift trace i backwards (negative delta_tau) so it aligns with target_t
            shifted = self.frequency_domain_shift(traces[i], -delta_tau)
            aligned_traces.append(shifted)
            
        return np.array(aligned_traces)

    def estimate_best_alignment_velocity(
        self,
        traces,
        depth,
        velocity_init,
        tx_x=None,
        reference_idx=0,
        window_width=3e-9,
        search_factors=(0.85, 1.15),
        n_trials=31,
    ):
        """
        Finds the best velocity for geometric flattening by maximizing
        trace-to-trace coherence in a reflection window after alignment.
        """
        if traces is None or len(traces) == 0:
            return velocity_init, 0.0
        if tx_x is None:
            tx_x = getattr(self, 'tx_start_x', self.width / 2.0)

        rx_array = self.get_rx_x_array()
        if len(rx_array) > len(traces):
            rx_array = rx_array[:len(traces)]

        lo, hi = search_factors
        vel_candidates = np.linspace(velocity_init * lo, velocity_init * hi, n_trials)

        best_velocity = velocity_init
        best_score = -np.inf

        ref_offset = np.abs(rx_array[reference_idx] - tx_x)

        for vel in vel_candidates:
            aligned = self.align_all_traces(
                traces,
                depth=depth,
                velocity=vel,
                tx_x=tx_x,
                reference_idx=reference_idx,
            )

            ref_arrival = np.sqrt(ref_offset**2 + (4 * (depth**2))) / vel
            windows = np.array([
                self.extract_window(tr, ref_arrival, window_width) for tr in aligned
            ])

            # Coherence score: ratio between stacked energy and average single-trace energy.
            stack = np.mean(windows, axis=0)
            coherent_energy = np.sum(stack**2)
            average_energy = np.mean(np.sum(windows**2, axis=1)) + 1e-12
            score = coherent_energy / average_energy

            if score > best_score:
                best_score = score
                best_velocity = vel

        return best_velocity, best_score

    def residual_align_traces(self, traces, reference_idx=0, window_center=6.0e-9, window_width=3.0e-9):
        """
        Fine-tunes the alignment of given traces using cross-correlation with a reference trace.
        Calculates and applies an empirical (data-driven) residual time shift to optimally flatten events.
        """
        if len(traces) == 0: return traces
        
        ref_trace = traces[reference_idx]
        
        # Extract window if specified
        if window_center is not None and window_width is not None:
            ref_window = self.extract_window(ref_trace, window_center, window_width)
        else:
            ref_window = ref_trace
            
        aligned = []
        for i in range(len(traces)):
            if i == reference_idx:
                aligned.append(traces[i])
                continue
                
            trace = traces[i]
            if window_center is not None and window_width is not None:
                tr_window = self.extract_window(trace, window_center, window_width)
            else:
                tr_window = trace
                
            # Cross correlate
            cc = np.correlate(tr_window, ref_window, mode='full')
            lags = np.arange(-len(tr_window)+1, len(tr_window))
            
            # Find the peak of the cross-correlation
            peak_idx = np.argmax(cc)
            lag_samples = lags[peak_idx]
            
            # Parabolic interpolation for sub-sample accuracy
            if 0 < peak_idx < len(cc) - 1:
                alpha, beta, gamma = cc[peak_idx-1], cc[peak_idx], cc[peak_idx+1]
                p = 0.5 * (alpha - gamma) / (alpha - 2*beta + gamma + 1e-12)
                sub_lag = lag_samples + p
            else:
                sub_lag = lag_samples
                
            # delta_tau is the delay of 'trace' relative to 'ref_trace' in seconds
            delta_tau = sub_lag * self.dt
            
            # Shift the trace backward (negative delta_tau) so it aligns with ref_trace
            shifted = self.frequency_domain_shift(traces[i], -delta_tau)
            aligned.append(shifted)
            
        return np.array(aligned)

    def run_shift_and_correlate(self, traces, depth, velocity, tx_x=None, window_width=2e-9, trace_pairs=None):
        """
        Calculates the Geometric Time Shift for a horizontal fracture at depth d.
        Applying the Shift-and-Correlate method, traces are shifted relative to specified pairs, 
        windowed, and their XWT is computed. Returns the power and phase of the XWT.
        """
        if len(traces) < 2:
            return None, None, None, None

        if tx_x is None:
            tx_x = getattr(self, 'tx_start_x', self.width / 2.0)
            
        rx_array = self.get_rx_x_array()
        offsets = np.abs(rx_array - tx_x)
        
        # t(x) = sqrt(x^2 + (2d)^2) / v
        t_arrivals = np.sqrt(offsets**2 + (4 * (depth**2))) / velocity
        
        xwt_power_list = []
        xwt_phase_list = []
        shifted_windows_list = []
        aligned_traces_full_list = []
        freqs_out = None
        
        if trace_pairs is None:
            trace_pairs = [(i, i+1) for i in range(len(traces) - 1)]
            
        for (i, j) in trace_pairs:
            if i >= len(traces) or j >= len(traces):
                continue
                
            t1 = t_arrivals[i]
            t2 = t_arrivals[j]
            delta_tau = t2 - t1
            
            trace1 = traces[i]
            trace2 = traces[j]
            
            # Sub-sample frequency domain shift (move trace2 backwards by delta_tau)
            # advancing trace 2 so that its arrival perfectly matches trace 1
            shifted_trace2 = self.frequency_domain_shift(trace2, -delta_tau)
            
            # Extract window around identical arrival time t1 for BOTH
            # After trace2 was shifted to align with t1, both have arrival aligned to t1
            window1 = self.extract_window(trace1, t1, window_width)
            window2 = self.extract_window(shifted_trace2, t1, window_width)
            
            # Compute Cross-Wavelet Transform (XWT) on the aligned windows
            freqs_out, power, phase = self.compute_xwt(window1, window2)
            
            xwt_power_list.append(power)
            xwt_phase_list.append(phase)
            shifted_windows_list.append((window1, window2))
            aligned_traces_full_list.append((trace1, shifted_trace2))
            
        return freqs_out, np.array(xwt_power_list), np.array(xwt_phase_list), shifted_windows_list, aligned_traces_full_list

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

