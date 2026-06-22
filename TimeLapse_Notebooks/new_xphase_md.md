## Section 8 — Cross-Phase Spectrogram (TWT vs. Frequency)

For the central trace (x = x_sc0, the scatterer's baseline position — the fixed
reference column shared by every monitor scenario), use CLSSA (the same narrowband
decomposition as Section 7) to get the baseline and monitor phase spectra
$\phi_\text{base}(f,\tau)$ and $\phi_\text{mon}(f,\tau)$, and form the cross-phase

$$\Delta\Phi(\tau, f) = \angle\Big[\exp\big(i\,(\phi_\text{base} - \phi_\text{mon})\big)\Big]
= \angle\,XS(\tau,f), \qquad XS = S_\text{base}\cdot S_\text{mon}^{*}$$

without ever forming the complex STFT explicitly — CLSSA samples frequency directly, so
unlike a generic STFT it isn't limited to a handful of coarse FFT bins regardless of
window length.

Plotting $\Delta\Phi(\tau,f)$ as a 2-D heatmap (frequency on the x-axis, TWT on the
y-axis, phase in $[-\pi, +\pi]$ on the colorbar) tests whether the measured phase shift
is tied to the actual target rather than to background noise:

- **Outside** the target's reflection time, the background should look like chaotic,
  pixelated salt-and-pepper noise — randomly wrapped phase from incoherent background
  scattering.
- **Exactly at** the TWT of the scatterer, a clean, coherent window of structured phase
  should appear across the antenna's frequency band.
- If the target **moved**, that window shows a smooth vertical colour gradient
  (a fringe pattern) — the phase ramps with frequency, consistent with the $-2\pi f\Delta t$
  Fourier-shift theorem.
- If a diffuse change occurred instead (e.g. fluid filling a pore space), the window would
  instead show a solid, near-uniform block of a single colour.

For separations **below ¼λ** the phase change at the scatterer is very subtle, so those
scenarios additionally zoom the TWT axis in on a window around the scatterer depth and
clip the colorbar to a much smaller range, so the subtle contrast is still visible.
