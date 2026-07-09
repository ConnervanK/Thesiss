#import "../template.typ": *

= Hypothesis 1.5 (Optional): Local Phase-Gradient Methods <sec:hyp1-h15>

The phase-plane fit of @sec:hyp1-phaseplane is a _global_ operation: it acts
on the whole 2D spectrum at once. This appendix explores a _local_ alternative
--- the spatial instantaneous phase, extracted trace-by-trace with a
Riesz/Hilbert transform, and its time-frequency relatives --- to test
Hypothesis 1.5: whether the local phase gradients $partial phi \/ partial x$,
$partial phi \/ partial y$ offer an equivalent, simpler route to the same
displacement estimate.

#draftnote[this appendix is exploratory: the underlying
Riesz-transform/monogenic-signal machinery is implemented in
`CWT_playground.ipynb`, but it has not yet been validated quantitatively
against the global WLS fit of @sec:hyp1-phaseplane on the same datasets ---
treat the results below as a qualitative complement to, not a replacement for,
Hypothesis 1's primary evidence.]

== Instantaneous Phase Imaging --- Lateral Displacement <sec:tlp-instphase-horizontal>

The 2D Riesz transform of @eq:kx-inst and @eq:dphi-lateral is applied
directly to the Gazdag-migrated baseline and monitor images, producing an
instantaneous-phase difference map $#Dphi (z,x)$ and, at the scatterer depth,
a 1D phase cross-section whose slope at the zero-crossing encodes $#Dx$. This
analysis is repeated at all seven lateral scales from $2 lambda$ to
$1 \/ 32 lambda$; @fig:tlp-instphase-horiz-example shows the representative
$1\/4 lambda$ case ($#Dx = 28.0 "mm"$), and
@fig:tlp-instphase-horiz-summary summarises the fitted zero-crossing slope
across all seven scales. The full seven-scale sweep is given in
@app:sec-instphase-horizontal.

#figure(
  subfigs(cols: 2,
    img("TLP_025_Instantaneous_Phase__Gazdag____Baseline_vs_¼λ___Δx__280_mm.png"),
    img("TLP_026_Phase_cross-section__z__676_mm____Gazdag____¼λ__Δx__280_mm.png"),
    img("TLP_027_Phase_cross-section__zoomed__½λ_around_scatterers.png"),
  ),
  caption: [Instantaneous-phase analysis for the representative
    $#Dx = 1\/4 lambda = 28.0 "mm"$ lateral displacement, Gazdag-migrated:
    (a) the 2D wrapped phase-difference map; (b) the 1D phase cross-section
    at the scatterer depth $z = 676 "mm"$; (c) the same cross-section zoomed
    to $plus.minus 1\/2 lambda$ around the scatterer.],
) <fig:tlp-instphase-horiz-example>

#figure(
  img("TLP_037_Gazdag__Δφ_zero-crossing_slope_vs_scatterer_separation.png", width: 70%),
  caption: [Fitted zero-crossing slope of the lateral instantaneous
    phase-difference cross-section, as a function of true scatterer
    displacement, across all seven scales from $2 lambda$ to $1\/32 lambda$
    (Gazdag-migrated).],
) <fig:tlp-instphase-horiz-summary>

== Instantaneous Phase Imaging --- Vertical Displacement <sec:tlp-instphase-vertical>

The same instantaneous-phase analysis is repeated for the vertical
displacement dataset, with the cross-section now taken along $z$ at the
fixed lateral position $x = 2.0 "m"$ of the scatterer, across six scales
from $1 lambda$ to $1 \/ 32 lambda$. @fig:tlp-instphase-vert-example shows
the representative $1\/4 lambda$ case ($#Dz = 28.0 "mm"$), and
@fig:tlp-instphase-vert-summary summarises the mean phase difference across
all six scales. The full sweep is given in @app:sec-instphase-vertical.

#figure(
  subfigs(cols: 2,
    img("TLP_044_VerticalTimeLapse__Instantaneous_Phase__Gazdag____Baseline_v.png"),
    img("TLP_045_Phase_cross-section__x__2000_mm____Gazdag____¼λ__Δz__280_mm.png"),
    img("TLP_046_Phase_cross-section__zoomed__½λ_around_scatterer_depths.png"),
  ),
  caption: [Instantaneous-phase analysis for the representative
    $#Dz = 1\/4 lambda = 28.0 "mm"$ vertical displacement, Gazdag-migrated:
    (a) the 2D wrapped phase-difference map; (b) the 1D phase cross-section
    at $x = 2000 "mm"$; (c) the same cross-section zoomed to
    $plus.minus 1\/2 lambda$ around the scatterer depths.],
) <fig:tlp-instphase-vert-example>

#figure(
  img("TLP_056_VerticalTimeLapse__Gazdag__Mean_Δφ_vs_scatterer_vertical_shi.png", width: 70%),
  caption: [Mean instantaneous phase difference $#Dphi$ as a function of true
    vertical scatterer shift, across all six scales from $1 lambda$ to
    $1\/32 lambda$ (Gazdag-migrated). Unlike the lateral case
    (@fig:tlp-instphase-horiz-summary), this is a plateau level rather than a
    cross-section slope, consistent with the asymmetry derived in
    @sec:th-duality.],
) <fig:tlp-instphase-vert-summary>

== Spectral-Line Phase Analysis <sec:tlp-spectral-line>

The localised-Fourier-shift spectral-line equation @eq:local-phase-line is
verified directly on individual traces using the CLSSA decomposition of
`phase_decomposition.py`, applied both to the Gazdag-migrated image and to
the raw, unmigrated B-scan, at the column position containing each scenario's
scatterer. @fig:tlp-spectral-line-example compares the representative
$1\/4 lambda$ case in both domains; the full seven-scale sweep for both is
given in @app:sec-spectral-line.

#figure(
  subfigs(cols: 1,
    img("TLP_060_Spectral_Line_CLSSA____Gazdag____¼λ___Δx__280_mm__02500λ.png"),
    img("TLP_067_Spectral_Line_CLSSA____Raw_Unmigrated____¼λ___Δx__280_mm__02.png"),
  ),
  caption: [CLSSA spectral-line decomposition (amplitude spectrum, phase
    spectrum, and phase gather, @sec:th-local) at the representative
    $#Dx = 1\/4 lambda = 28.0 "mm"$ scale: (a) Gazdag-migrated trace; (b) raw,
    unmigrated trace at the same scenario.],
) <fig:tlp-spectral-line-example>

#draftnote[state whether the migrated and raw spectral lines in
@fig:tlp-spectral-line-example give consistent $#Dt$/slope estimates once
converted with @eq:dt-dz, which would confirm that migration is not required
for the temporal spectral-line method to work, only for the 2D wavenumber
plane fit of @sec:th-wls.]

== Cross-Phase Spectrograms <sec:tlp-spectrogram>

@fig:tlp-spectrogram-example shows the cross-phase spectrogram
$#DPhi (tau, f)$ of @sec:th-local for two contrasting scales: a large, easily
visible $1\/4 lambda$ displacement and the smallest, $1\/32 lambda$, where the
coherent window around the scatterer's two-way time becomes much harder to
distinguish from the incoherent background. The full seven-scale sweep is
given in @app:sec-spectrogram.

#figure(
  subfigs(cols: 2,
    img("TLP_074_Cross-Phase_Spectrogram_CLSSA__ΔΦτf____Gazdag____¼λ.png"),
    img("TLP_077_Cross-Phase_Spectrogram_CLSSA__ΔΦτf____Gazdag____¹₃₂λ.png"),
  ),
  caption: [Cross-phase spectrogram $#DPhi (tau, f)$ (two-way time against
    frequency) at two contrasting lateral displacement scales,
    Gazdag-migrated: (a) $#Dx = 1\/4 lambda$; (b) $#Dx = 1\/32 lambda$.],
) <fig:tlp-spectrogram-example>

== Localised Fourier Shift via Short-Time Fourier Transform <sec:tlp-stft>

Finally, the three-panel localised-STFT decomposition of @eq:local-shift and
@eq:local-phase-line --- amplitude spectrum $A(tau,f)$, phase spectrum
$#DPhi (tau,f)$, and phase-angle domain $A(tau, theta)$ --- is applied with a
Gaussian analysis window. @fig:tlp-stft-example shows two of the seven outputs
produced by this sweep; the full set is given in @app:sec-stft.

#figure(
  subfigs(cols: 1,
    img("TLP_078_Localized_Fourier_Shift_STFT_Gaussian_win12_smp____Gazdag.png"),
    img("TLP_081_Localized_Fourier_Shift_STFT_Gaussian_win12_smp____Gazdag.png"),
  ),
  caption: [Localised-STFT phase decomposition (Gaussian window,
    Gazdag-migrated) for two of the seven lateral displacement scales:
    (a) scale 1 of 7; (b) scale 4 of 7.
    #draftnote[TLP\_078--084 all share an identical auto-generated filename
    stem and do not encode which $#Dx$ each panel corresponds to --- confirm
    the scale for each of the seven figures against the
    `TimeLapse_Processing.ipynb` cell order before finalising these
    captions.]],
) <fig:tlp-stft-example>
