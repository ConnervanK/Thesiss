import numpy as np
import h5py
from scipy.signal import hilbert, butter, filtfilt
import pywt

from utils.processing import apply_svd_filter, envelope, compute_cwt_image


class GPRBScanData:
    """
    Wrapper for a GPR B-scan result (common-offset, moving Tx-Rx pair).

    Each trace i was recorded with:
      Tx at  midpoint_x[i] - tx_rx_offset/2
      Rx at  midpoint_x[i] + tx_rx_offset/2
    where midpoint_x spans the same positions as the original multi-receiver survey.
    """

    def __init__(self, baseline_file, timelapse_file,
                 midpoint_start_x=1.0, midpoint_spacing=None,
                 n_traces=50, tx_rx_offset=0.1):
        self.baseline_file = baseline_file
        self.timelapse_file = timelapse_file
        self.midpoint_start_x = midpoint_start_x
        self.midpoint_spacing = midpoint_spacing if midpoint_spacing is not None else 2.0 / max(1, n_traces - 1)
        self.n_traces = n_traces
        self.tx_rx_offset = tx_rx_offset

        self.time = None
        self.dt = None
        self.baseline_traces = None
        self.timelapse_traces = None
        self.diff_traces = None

    # ------------------------------------------------------------------
    # Geometry helpers
    # ------------------------------------------------------------------

    def get_midpoint_x_array(self):
        return np.linspace(
            self.midpoint_start_x,
            self.midpoint_start_x + self.midpoint_spacing * (self.n_traces - 1),
            self.n_traces,
        )

    def get_tx_x_array(self):
        return self.get_midpoint_x_array() - self.tx_rx_offset / 2.0

    def get_rx_x_array(self):
        return self.get_midpoint_x_array() + self.tx_rx_offset / 2.0

    # ------------------------------------------------------------------
    # Data I/O
    # ------------------------------------------------------------------

    def load_data(self):
        with h5py.File(self.baseline_file, 'r') as fa, h5py.File(self.timelapse_file, 'r') as fh:
            rx_group_base = fa['rxs'] if 'rxs' in fa else fa
            rx_group_tl = fh['rxs'] if 'rxs' in fh else fh

            def sorted_rx_keys(grp):
                keys = [k for k in grp.keys() if k.startswith('rx') and k != 'rxs']
                keys.sort(key=lambda k: int(k[2:]))
                return keys

            keys_base = sorted_rx_keys(rx_group_base)
            keys_tl = sorted_rx_keys(rx_group_tl)

            self.dt = fa.attrs['dt']
            iterations = fa.attrs['Iterations']
            self.time = np.arange(iterations) * self.dt

            base_traces, tl_traces = [], []
            for kb, kt in zip(keys_base, keys_tl):
                ez_b = (np.array(rx_group_base[kb]['Ez'])
                        if isinstance(rx_group_base[kb], h5py.Group)
                        else np.array(rx_group_base[kb]))
                ez_t = (np.array(rx_group_tl[kt]['Ez'])
                        if isinstance(rx_group_tl[kt], h5py.Group)
                        else np.array(rx_group_tl[kt]))
                base_traces.append(ez_b)
                tl_traces.append(ez_t)

        self.baseline_traces = np.array(base_traces)
        self.timelapse_traces = np.array(tl_traces)

        if self.baseline_traces.ndim == 3:
            self.baseline_traces = self.baseline_traces[:, 0, :]
            self.timelapse_traces = self.timelapse_traces[:, 0, :]

        self.n_traces = len(self.baseline_traces)
        self.diff_traces = self.subtract()

    def subtract(self):
        return self.baseline_traces - self.timelapse_traces

    # ------------------------------------------------------------------
    # Signal processing
    # ------------------------------------------------------------------

    def apply_svd_filter(self, traces, n_components_to_mute=1):
        return apply_svd_filter(traces, n_components_to_mute)

    def apply_envelope(self, traces):
        return np.abs(hilbert(traces, axis=1))

    def bandpass_filter(self, traces, lowcut, highcut, order=4):
        nyq = 0.5 / self.dt
        low, high = lowcut / nyq, highcut / nyq
        if high >= 1.0:
            high = 0.99
        b, a = butter(order, [low, high], btype='band')
        return filtfilt(b, a, traces, axis=1)

    def frequency_domain_shift(self, trace, delta_t):
        N = len(trace)
        F = np.fft.rfft(trace)
        freqs = np.fft.rfftfreq(N, d=self.dt)
        return np.fft.irfft(F * np.exp(-1j * 2 * np.pi * freqs * delta_t), n=N)

    def compute_cwt_image(self, traces=None, wavelet='cmor1.5-1.0', freqs=None):
        if traces is None:
            traces = self.diff_traces
        return compute_cwt_image(traces, self.dt, wavelet=wavelet, freqs=freqs)

    def compute_xwt(self, trace1, trace2, wavelet='cmor1.5-1.0', freqs=None):
        if freqs is None:
            freqs = np.linspace(0.1e9, 5.5e9, 100)
        sampling_freq = 1.0 / self.dt
        scales = pywt.frequency2scale(wavelet, freqs / sampling_freq)
        coefs1, _ = pywt.cwt(trace1, scales, wavelet, sampling_period=self.dt)
        coefs2, _ = pywt.cwt(trace2, scales, wavelet, sampling_period=self.dt)
        xwt = coefs1 * np.conj(coefs2)
        return freqs, np.abs(xwt), np.angle(xwt)

    # ------------------------------------------------------------------
    # Migration
    # ------------------------------------------------------------------

    def kirchhoff_migration_bscan(self, traces, velocity, max_depth, dz=0.005,
                                   max_angle_deg=65.0, apply_obliquity=True):
        """
        Common-offset Kirchhoff depth migration.

        For each image point (x_out, z) the contribution of trace i is looked up at
        t = (dist(x_out, tx_i) + dist(x_out, rx_i)) / velocity,
        where tx_i and rx_i are the positions of the moving Tx/Rx pair for trace i.

        self.time is in seconds; velocity in m/s.
        max_angle_deg: aperture mute — skip contributions beyond this angle from vertical.
        apply_obliquity: weight each contribution by the average cosine of incidence angles.
        """
        assert self.time is not None, "load_data() must be called first"
        time_ns = self.time * 1e9
        dt_ns = time_ns[1] - time_ns[0]

        midpoints = self.get_midpoint_x_array()
        tx_arr = self.get_tx_x_array()
        rx_arr = self.get_rx_x_array()

        z_array = np.arange(0, max_depth + dz, dz)
        migrated = np.zeros((len(z_array), len(midpoints)))
        cos_min = np.cos(np.deg2rad(max_angle_deg))

        for iz, z in enumerate(z_array):
            if z == 0:
                continue
            for ix_out, x_out in enumerate(midpoints):
                for i in range(len(traces)):
                    dist_tx = np.sqrt((x_out - tx_arr[i]) ** 2 + z ** 2)
                    dist_rx = np.sqrt((x_out - rx_arr[i]) ** 2 + z ** 2)
                    cos_tx = z / dist_tx
                    cos_rx = z / dist_rx

                    if cos_tx < cos_min or cos_rx < cos_min:
                        continue

                    t_ns = (dist_tx + dist_rx) / velocity * 1e9
                    t_idx = int(np.round((t_ns - time_ns[0]) / dt_ns))
                    if 0 <= t_idx < traces.shape[1]:
                        w = (cos_tx + cos_rx) * 0.5 if apply_obliquity else 1.0
                        migrated[iz, ix_out] += w * traces[i, t_idx]

        return migrated, z_array
