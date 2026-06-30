import numpy as np


def _as_2d(data: np.ndarray) -> np.ndarray:
	array = np.asarray(data)
	if array.ndim == 1:
		return array[:, None]
	return array


def filter_data(data, fq, sfreq, btype='bandpass'):
	array = _as_2d(np.asarray(data, dtype=float))
	n_samples = array.shape[0]
	freqs = np.fft.rfftfreq(n_samples, d=1.0 / float(sfreq))
	spectrum = np.fft.rfft(array, axis=0)

	if btype == 'bandpass':
		f_low, f_high = fq
		keep = (freqs >= f_low) & (freqs <= f_high)
	elif btype == 'lowpass':
		f_high = float(fq)
		keep = freqs <= f_high
	elif btype == 'highpass':
		f_low = float(fq)
		keep = freqs >= f_low
	else:
		raise ValueError(f"Unsupported filter type: {btype}")

	spectrum[~keep, ...] = 0.0
	filtered = np.fft.irfft(spectrum, n=n_samples, axis=0)
	return filtered.reshape(np.asarray(data, dtype=float).shape)


def remove_mean(data, start, stop):
	array = np.asarray(data, dtype=float)
	window = array[start:stop]
	mean = np.mean(window, axis=0, keepdims=True)
	return array - mean, mean.squeeze()


def linear_gain(data, t):
	array = np.asarray(data, dtype=float)
	time = np.asarray(t, dtype=float)
	if array.ndim == 1:
		gain = time
	else:
		gain = time[:, None]
	return array * gain, gain


def remove_svd(data, low_s=0, high_s=1):
	array = np.asarray(data, dtype=float)
	u, singular_values, vt = np.linalg.svd(array, full_matrices=False)
	filtered_singular_values = singular_values.copy()
	filtered_singular_values[low_s:high_s] = 0.0
	reconstructed = (u * filtered_singular_values) @ vt
	return reconstructed, singular_values


def _shift_trace(trace, shift):
	trace = np.asarray(trace, dtype=float)
	indices = np.arange(trace.shape[0], dtype=float)
	return np.interp(indices - shift, indices, trace, left=0.0, right=0.0)


def align_traces(data, reference, upsample=5, normalize=False, align_reference=False):
	array = np.asarray(data, dtype=float)
	ref = np.asarray(reference, dtype=float)

	if array.ndim == 1:
		array = array[:, None]
	if ref.ndim == 1:
		ref = ref[:, None]

	n_samples, n_traces = array.shape
	aligned = np.zeros_like(array)
	shifts = np.zeros(n_traces, dtype=float)
	correlations = np.zeros(n_traces, dtype=float)

	for idx in range(n_traces):
		trace = array[:, idx]
		ref_idx = idx if ref.shape[1] > 1 else 0
		ref_trace = ref[:, min(ref_idx, ref.shape[1] - 1)]

		trace_work = trace - np.mean(trace)
		ref_work = ref_trace - np.mean(ref_trace)

		if normalize:
			trace_std = np.std(trace_work)
			ref_std = np.std(ref_work)
			if trace_std > 0:
				trace_work = trace_work / trace_std
			if ref_std > 0:
				ref_work = ref_work / ref_std

		correlation = np.correlate(trace_work, ref_work, mode='full')
		lag_samples = np.arange(-n_samples + 1, n_samples, dtype=float)
		peak = int(np.argmax(correlation))
		shift = lag_samples[peak]

		if 0 < peak < correlation.size - 1:
			y0, y1, y2 = correlation[peak - 1], correlation[peak], correlation[peak + 1]
			denom = y0 - 2.0 * y1 + y2
			if denom != 0.0:
				shift += 0.5 * (y0 - y2) / denom

		shifts[idx] = shift
		correlations[idx] = correlation[peak]
		aligned[:, idx] = _shift_trace(trace, shift)

	if data.ndim == 1:
		aligned = aligned[:, 0]

	if align_reference:
		return aligned, shifts, correlations

	return aligned, shifts, correlations
