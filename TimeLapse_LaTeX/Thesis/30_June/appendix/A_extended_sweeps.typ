#import "../template.typ": *

= Extended Wavelength-Sweep Figures <app:extended-sweeps>

@sec:hyp1-phaseplane showed one representative scale from each of five
repeated parameter sweeps produced by `TimeLapse_Processing.ipynb`, to keep
the main narrative focused. This appendix gives the complete sweep for each
of those five analyses, across every wavelength scale that was simulated.
Figures here intentionally repeat the scale already highlighted in
@sec:hyp1-phaseplane, so that each full sweep is shown as a single,
self-contained reference figure.

== Full Sweep: Lateral Instantaneous-Phase Analysis <app:sec-instphase-horizontal>

@fig:app-instphase-horiz-full, @fig:app-instphase-horiz-cs-full, and
@fig:app-instphase-horiz-zoom-full give the complete seven-scale lateral sweep
underlying @sec:tlp-instphase-horizontal: the 2D instantaneous
phase-difference map, the phase cross-section at the scatterer depth
$z = 676 "mm"$, and the same cross-section zoomed to
$plus.minus 1\/2 lambda$, for every scale from $2 lambda$ to $1\/32 lambda$.

#figure(
  grid(columns: (1fr, 1fr, 1fr, 1fr), gutter: 0.5em,
    img("TLP_016_Instantaneous_Phase__Gazdag____Baseline_vs_2λ___Δx__2240_mm.png"),
    img("TLP_019_Instantaneous_Phase__Gazdag____Baseline_vs_1λ___Δx__1120_mm.png"),
    img("TLP_022_Instantaneous_Phase__Gazdag____Baseline_vs_½λ___Δx__560_mm.png"),
    img("TLP_025_Instantaneous_Phase__Gazdag____Baseline_vs_¼λ___Δx__280_mm.png"),
    img("TLP_028_Instantaneous_Phase__Gazdag____Baseline_vs_⅛λ___Δx__140_mm.png"),
    img("TLP_031_Instantaneous_Phase__Gazdag____Baseline_vs_¹₁₆λ___Δx__70_mm.png"),
    img("TLP_034_Instantaneous_Phase__Gazdag____Baseline_vs_¹₃₂λ___Δx__35_mm.png"),
  ),
  caption: [Full seven-scale sweep: 2D instantaneous phase-difference map,
    Gazdag-migrated, lateral displacement from $2 lambda$ (top-left) down to
    $1\/32 lambda$ (bottom row): (a) $2 lambda$; (b) $1 lambda$;
    (c) $1\/2 lambda$; (d) $1\/4 lambda$; (e) $1\/8 lambda$;
    (f) $1\/16 lambda$; (g) $1\/32 lambda$.],
) <fig:app-instphase-horiz-full>

#figure(
  grid(columns: (1fr, 1fr, 1fr, 1fr), gutter: 0.5em,
    img("TLP_017_Phase_cross-section__z__676_mm____Gazdag____2λ__Δx__2240_mm.png"),
    img("TLP_020_Phase_cross-section__z__676_mm____Gazdag____1λ__Δx__1120_mm.png"),
    img("TLP_023_Phase_cross-section__z__676_mm____Gazdag____½λ__Δx__560_mm.png"),
    img("TLP_026_Phase_cross-section__z__676_mm____Gazdag____¼λ__Δx__280_mm.png"),
    img("TLP_029_Phase_cross-section__z__676_mm____Gazdag____⅛λ__Δx__140_mm.png"),
    img("TLP_032_Phase_cross-section__z__676_mm____Gazdag____¹₁₆λ__Δx__70_mm.png"),
    img("TLP_035_Phase_cross-section__z__676_mm____Gazdag____¹₃₂λ__Δx__35_mm.png"),
  ),
  caption: [Full seven-scale sweep: phase cross-section at $z = 676 "mm"$,
    Gazdag-migrated, lateral displacement from $2 lambda$ down to
    $1\/32 lambda$: (a) $2 lambda$; (b) $1 lambda$; (c) $1\/2 lambda$;
    (d) $1\/4 lambda$; (e) $1\/8 lambda$; (f) $1\/16 lambda$;
    (g) $1\/32 lambda$.],
) <fig:app-instphase-horiz-cs-full>

#figure(
  grid(columns: (1fr, 1fr, 1fr, 1fr), gutter: 0.5em,
    img("TLP_018_Phase_cross-section__zoomed__½λ_around_scatterers.png"),
    img("TLP_021_Phase_cross-section__zoomed__½λ_around_scatterers.png"),
    img("TLP_024_Phase_cross-section__zoomed__½λ_around_scatterers.png"),
    img("TLP_027_Phase_cross-section__zoomed__½λ_around_scatterers.png"),
    img("TLP_030_Phase_cross-section__zoomed__½λ_around_scatterers.png"),
    img("TLP_033_Phase_cross-section__zoomed__½λ_around_scatterers.png"),
    img("TLP_036_Phase_cross-section__zoomed__½λ_around_scatterers.png"),
  ),
  caption: [Full seven-scale sweep: phase cross-section zoomed to
    $plus.minus 1\/2 lambda$ around the scatterers, Gazdag-migrated, lateral
    displacement from $2 lambda$ down to $1\/32 lambda$: (a) $2 lambda$;
    (b) $1 lambda$; (c) $1\/2 lambda$; (d) $1\/4 lambda$; (e) $1\/8 lambda$;
    (f) $1\/16 lambda$; (g) $1\/32 lambda$.],
) <fig:app-instphase-horiz-zoom-full>

== Full Sweep: Vertical Instantaneous-Phase Analysis <app:sec-instphase-vertical>

@fig:app-instphase-vert-full, @fig:app-instphase-vert-cs-full, and
@fig:app-instphase-vert-zoom-full give the complete six-scale vertical sweep
underlying @sec:tlp-instphase-vertical, from $1 lambda$ to $1\/32 lambda$.

#figure(
  grid(columns: (1fr, 1fr, 1fr), gutter: 0.5em,
    img("TLP_038_VerticalTimeLapse__Instantaneous_Phase__Gazdag____Baseline_v.png"),
    img("TLP_041_VerticalTimeLapse__Instantaneous_Phase__Gazdag____Baseline_v.png"),
    img("TLP_044_VerticalTimeLapse__Instantaneous_Phase__Gazdag____Baseline_v.png"),
    img("TLP_047_VerticalTimeLapse__Instantaneous_Phase__Gazdag____Baseline_v.png"),
    img("TLP_050_VerticalTimeLapse__Instantaneous_Phase__Gazdag____Baseline_v.png"),
    img("TLP_053_VerticalTimeLapse__Instantaneous_Phase__Gazdag____Baseline_v.png"),
  ),
  caption: [Full six-scale sweep: 2D instantaneous phase-difference map,
    Gazdag-migrated, vertical displacement from $1 lambda$ down to
    $1\/32 lambda$: (a) $1 lambda$; (b) $1\/2 lambda$; (c) $1\/4 lambda$;
    (d) $1\/8 lambda$; (e) $1\/16 lambda$; (f) $1\/32 lambda$.],
) <fig:app-instphase-vert-full>

#figure(
  grid(columns: (1fr, 1fr, 1fr), gutter: 0.5em,
    img("TLP_039_Phase_cross-section__x__2000_mm____Gazdag____1λ__Δz__1120_mm.png"),
    img("TLP_042_Phase_cross-section__x__2000_mm____Gazdag____½λ__Δz__560_mm.png"),
    img("TLP_045_Phase_cross-section__x__2000_mm____Gazdag____¼λ__Δz__280_mm.png"),
    img("TLP_048_Phase_cross-section__x__2000_mm____Gazdag____⅛λ__Δz__140_mm.png"),
    img("TLP_051_Phase_cross-section__x__2000_mm____Gazdag____¹₁₆λ__Δz__70_mm.png"),
    img("TLP_054_Phase_cross-section__x__2000_mm____Gazdag____¹₃₂λ__Δz__35_mm.png"),
  ),
  caption: [Full six-scale sweep: phase cross-section at $x = 2000 "mm"$,
    Gazdag-migrated, vertical displacement from $1 lambda$ down to
    $1\/32 lambda$: (a) $1 lambda$; (b) $1\/2 lambda$; (c) $1\/4 lambda$;
    (d) $1\/8 lambda$; (e) $1\/16 lambda$; (f) $1\/32 lambda$.],
) <fig:app-instphase-vert-cs-full>

#figure(
  grid(columns: (1fr, 1fr, 1fr), gutter: 0.5em,
    img("TLP_040_Phase_cross-section__zoomed__½λ_around_scatterer_depths.png"),
    img("TLP_043_Phase_cross-section__zoomed__½λ_around_scatterer_depths.png"),
    img("TLP_046_Phase_cross-section__zoomed__½λ_around_scatterer_depths.png"),
    img("TLP_049_Phase_cross-section__zoomed__½λ_around_scatterer_depths.png"),
    img("TLP_052_Phase_cross-section__zoomed__½λ_around_scatterer_depths.png"),
    img("TLP_055_Phase_cross-section__zoomed__½λ_around_scatterer_depths.png"),
  ),
  caption: [Full six-scale sweep: phase cross-section zoomed to
    $plus.minus 1\/2 lambda$ around the scatterer depths, Gazdag-migrated,
    vertical displacement from $1 lambda$ down to $1\/32 lambda$:
    (a) $1 lambda$; (b) $1\/2 lambda$; (c) $1\/4 lambda$; (d) $1\/8 lambda$;
    (e) $1\/16 lambda$; (f) $1\/32 lambda$.],
) <fig:app-instphase-vert-zoom-full>

== Full Sweep: Spectral-Line Phase Analysis <app:sec-spectral-line>

@fig:app-spectral-migrated-full and @fig:app-spectral-raw-full give the
complete seven-scale CLSSA spectral-line sweep underlying
@sec:tlp-spectral-line, for the Gazdag-migrated and raw unmigrated traces
respectively.

#figure(
  grid(columns: (1fr, 1fr, 1fr, 1fr), gutter: 0.5em,
    img("TLP_057_Spectral_Line_CLSSA____Gazdag____2λ___Δx__2240_mm__20000λ.png"),
    img("TLP_058_Spectral_Line_CLSSA____Gazdag____1λ___Δx__1120_mm__10000λ.png"),
    img("TLP_059_Spectral_Line_CLSSA____Gazdag____½λ___Δx__560_mm__05000λ.png"),
    img("TLP_060_Spectral_Line_CLSSA____Gazdag____¼λ___Δx__280_mm__02500λ.png"),
    img("TLP_061_Spectral_Line_CLSSA____Gazdag____⅛λ___Δx__140_mm__01250λ.png"),
    img("TLP_062_Spectral_Line_CLSSA____Gazdag____¹₁₆λ___Δx__70_mm__00625λ.png"),
    img("TLP_063_Spectral_Line_CLSSA____Gazdag____¹₃₂λ___Δx__35_mm__00312λ.png"),
  ),
  caption: [Full seven-scale sweep: CLSSA spectral-line decomposition,
    Gazdag-migrated trace, lateral displacement from $2 lambda$ down to
    $1\/32 lambda$: (a) $2 lambda$; (b) $1 lambda$; (c) $1\/2 lambda$;
    (d) $1\/4 lambda$; (e) $1\/8 lambda$; (f) $1\/16 lambda$;
    (g) $1\/32 lambda$.],
) <fig:app-spectral-migrated-full>

#figure(
  grid(columns: (1fr, 1fr, 1fr, 1fr), gutter: 0.5em,
    img("TLP_064_Spectral_Line_CLSSA____Raw_Unmigrated____2λ___Δx__2240_mm__2.png"),
    img("TLP_065_Spectral_Line_CLSSA____Raw_Unmigrated____1λ___Δx__1120_mm__1.png"),
    img("TLP_066_Spectral_Line_CLSSA____Raw_Unmigrated____½λ___Δx__560_mm__05.png"),
    img("TLP_067_Spectral_Line_CLSSA____Raw_Unmigrated____¼λ___Δx__280_mm__02.png"),
    img("TLP_068_Spectral_Line_CLSSA____Raw_Unmigrated____⅛λ___Δx__140_mm__01.png"),
    img("TLP_069_Spectral_Line_CLSSA____Raw_Unmigrated____¹₁₆λ___Δx__70_mm__0.png"),
    img("TLP_070_Spectral_Line_CLSSA____Raw_Unmigrated____¹₃₂λ___Δx__35_mm__0.png"),
  ),
  caption: [Full seven-scale sweep: CLSSA spectral-line decomposition, raw
    unmigrated trace, lateral displacement from $2 lambda$ down to
    $1\/32 lambda$: (a) $2 lambda$; (b) $1 lambda$; (c) $1\/2 lambda$;
    (d) $1\/4 lambda$; (e) $1\/8 lambda$; (f) $1\/16 lambda$;
    (g) $1\/32 lambda$.],
) <fig:app-spectral-raw-full>

== Full Sweep: Cross-Phase Spectrograms <app:sec-spectrogram>

@fig:app-spectrogram-full gives the complete seven-scale cross-phase
spectrogram sweep underlying @sec:tlp-spectrogram.

#figure(
  grid(columns: (1fr, 1fr, 1fr, 1fr), gutter: 0.5em,
    img("TLP_071_Cross-Phase_Spectrogram_CLSSA__ΔΦτf____Gazdag____2λ.png"),
    img("TLP_072_Cross-Phase_Spectrogram_CLSSA__ΔΦτf____Gazdag____1λ.png"),
    img("TLP_073_Cross-Phase_Spectrogram_CLSSA__ΔΦτf____Gazdag____½λ.png"),
    img("TLP_074_Cross-Phase_Spectrogram_CLSSA__ΔΦτf____Gazdag____¼λ.png"),
    img("TLP_075_Cross-Phase_Spectrogram_CLSSA__ΔΦτf____Gazdag____⅛λ.png"),
    img("TLP_076_Cross-Phase_Spectrogram_CLSSA__ΔΦτf____Gazdag____¹₁₆λ.png"),
    img("TLP_077_Cross-Phase_Spectrogram_CLSSA__ΔΦτf____Gazdag____¹₃₂λ.png"),
  ),
  caption: [Full seven-scale sweep: cross-phase spectrogram $#DPhi (tau, f)$,
    Gazdag-migrated, lateral displacement from $2 lambda$ down to
    $1\/32 lambda$: (a) $2 lambda$; (b) $1 lambda$; (c) $1\/2 lambda$;
    (d) $1\/4 lambda$; (e) $1\/8 lambda$; (f) $1\/16 lambda$;
    (g) $1\/32 lambda$.],
) <fig:app-spectrogram-full>

== Full Sweep: Localised Fourier Shift via STFT <app:sec-stft>

@fig:app-stft-full gives the complete seven-output sweep underlying
@sec:tlp-stft.

#draftnote[TLP\_078--084 share an identical auto-generated filename stem and
do not encode which displacement scale each panel corresponds to. The
ordering below is assumed to follow the same $2 lambda arrow.r 1\/32 lambda$
sequence as every other sweep in this appendix, on the basis of the cell
order in `TimeLapse_Processing.ipynb`, but this must be confirmed against the
notebook before the per-panel labels below are trusted.]

#figure(
  grid(columns: (1fr, 1fr, 1fr, 1fr), gutter: 0.5em,
    img("TLP_078_Localized_Fourier_Shift_STFT_Gaussian_win12_smp____Gazdag.png"),
    img("TLP_079_Localized_Fourier_Shift_STFT_Gaussian_win12_smp____Gazdag.png"),
    img("TLP_080_Localized_Fourier_Shift_STFT_Gaussian_win12_smp____Gazdag.png"),
    img("TLP_081_Localized_Fourier_Shift_STFT_Gaussian_win12_smp____Gazdag.png"),
    img("TLP_082_Localized_Fourier_Shift_STFT_Gaussian_win12_smp____Gazdag.png"),
    img("TLP_083_Localized_Fourier_Shift_STFT_Gaussian_win12_smp____Gazdag.png"),
    img("TLP_084_Localized_Fourier_Shift_STFT_Gaussian_win12_smp____Gazdag.png"),
  ),
  caption: [Full seven-output sweep: localised-STFT phase decomposition
    (Gaussian window, Gazdag-migrated). Scale labels are provisional --- see
    draft note above: (a) $2 lambda$ (assumed); (b) $1 lambda$ (assumed);
    (c) $1\/2 lambda$ (assumed); (d) $1\/4 lambda$ (assumed);
    (e) $1\/8 lambda$ (assumed); (f) $1\/16 lambda$ (assumed);
    (g) $1\/32 lambda$ (assumed).],
) <fig:app-stft-full>
