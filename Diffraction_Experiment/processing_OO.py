import numpy as np
import h5py
from scipy.signal import hilbert
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
        
        clear_results()
        
        (self.width, self.height, self.air_thick, self.f_top, self.f_bottom, 
         self.snap_time, self.dx_dy_dz, self.n_blocks, self.block_width, 
         self.rx_offset, self.rx_per_block) = create_input_file(**self.model_params)
         
        run_gprmax()
        self._load_data()

    def _load_data(self):
        with h5py.File(self.out_alt_file, 'r') as fa, h5py.File(self.out_homo_file, 'r') as fh:
            iterations = fa.attrs['Iterations']
            self.dt = fa.attrs['dt']
            self.time = np.arange(iterations) * self.dt * 1e9 # in ns
            
            diff_traces = []
            homo_traces = []
            alt_traces = []
            
            total_rx = self.n_blocks * self.rx_per_block
            for i in range(total_rx):
                # The first receiver is rx1 (usually next to the tx). The array receivers start at rx2.
                rx_name = f'rx{i+2}'
                try:
                    ez_alt = np.array(fa['rxs'][rx_name]['Ez'])
                    ez_homo = np.array(fh['rxs'][rx_name]['Ez'])
                    
                    alt_traces.append(ez_alt)
                    homo_traces.append(ez_homo)
                    diff_traces.append(ez_alt - ez_homo)
                except KeyError:
                    pass
                    
            self.diff_traces = np.array(diff_traces)
            self.homo_traces = np.array(homo_traces)
            self.alt_traces = np.array(alt_traces)
            
    def get_rx_x_array(self):
        """Returns the x-coordinates of the receiver array."""
        rx_dx = self.block_width / self.rx_per_block
        return np.array([(i + 0.5) * rx_dx for i in range(self.n_blocks * self.rx_per_block)])

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
        k = np.fft.fftshift(np.fft.fftfreq(traces.shape[0], self.block_width)) * 2 * np.pi
        
        fk_mag = np.abs(fk_data)
        return k, freqs, fk_mag

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
                    dist_rx = np.sqrt((x_out - rx_x)**2 + z**2)
                    
                    t_total_sec = (dist_tx + dist_rx) / velocity
                    t_total_ns = t_total_sec * 1e9
                    
                    t_idx = int(np.round((t_total_ns - self.time[0]) / dt_ns))
                    if 0 <= t_idx < len(self.time):
                        migrated_image[iz, ix_out] += traces[irx, t_idx]
                        
        return migrated_image, z_array