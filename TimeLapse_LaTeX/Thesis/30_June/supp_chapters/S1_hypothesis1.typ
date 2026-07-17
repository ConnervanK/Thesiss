#import "../template.typ": *

= Hypothesis 1: Supplementary Figures <supp:hyp1>

This chapter collects the Baseline-versus-Monitor point-spread-function
(PSF) comparisons and phase-plane shift-estimation diagnostics for
Chapter 4, "Hypothesis 1: Can We Infer Subwavelength Movement from Phase
Changes?", of the main thesis. Both figure types recur for every scenario
set and migration method and are already summarised in the main text by
the Rayleigh-criterion-ratio and WLS-displacement-error tables, so they are
collected here rather than reproduced in the main chapter. Section numbers
below mirror Chapter 4's own section names one-for-one.

== Lateral Movement <supp:hyp1-lateral>

=== Amplitude Test <supp:hyp1-lat-amplitude>

#figure(
  img("H1_006_Lateral_--_Amplitude_PSF_zoom_Baseline_vs_Monitor_RAW_migrat.png", width: 80%),
  caption: [Zoomed Baseline-versus-Monitor point-spread-function comparison
    for the raw (non-difference) migrated lateral images: peak-normalised
    Baseline (blue) and Monitor (red) 1D slices, with the Baseline FWHM
    shaded and both peaks marked, for every scenario (rows) and migration
    method (columns). Referenced from the "Lateral Movement -- Amplitude
    Test" section of Chapter 4.],
) <fig:supp-h1-lat-amp-zoom>

=== Phase Test <supp:hyp1-lat-phase>

#figure(
  img("H1_007_Lateral_--_Phase-plane_shift_estimation_--_Kirchhoff____Base.png", width: 100%),
  caption: [Kirchhoff-migrated phase-plane shift estimation, lateral study:
    Baseline versus each of the seven displacement scenarios. Referenced
    from the "Lateral Movement -- Phase Test" section of Chapter 4.],
) <fig:supp-h1-lat-phase-kirchhoff>

#figure(
  img("H1_008_Lateral_--_Phase-plane_shift_estimation_--_Gazdag____Baselin.png", width: 100%),
  caption: [Gazdag-migrated phase-plane shift estimation, lateral study:
    Baseline versus each of the seven displacement scenarios. Referenced
    from the "Lateral Movement -- Phase Test" section of Chapter 4.],
) <fig:supp-h1-lat-phase-gazdag>

#figure(
  img("H1_009_Lateral_--_Phase-plane_shift_estimation_--_Back-prop____Base.png", width: 100%),
  caption: [Back-propagation-migrated phase-plane shift estimation, lateral
    study: Baseline versus each of the seven displacement scenarios.
    Referenced from the "Lateral Movement -- Phase Test" section of
    Chapter 4.],
) <fig:supp-h1-lat-phase-backprop>

== Vertical Movement <supp:hyp1-vertical>

=== Amplitude Test <supp:hyp1-vert-amplitude>

#figure(
  img("H1_015_Vertical_--_Amplitude_PSF_zoom_Baseline_vs_Monitor_RAW_migra.png", width: 80%),
  caption: [Zoomed Baseline-versus-Monitor point-spread-function comparison
    for the raw migrated vertical images, all scenarios and migration
    methods. Referenced from the "Vertical Movement -- Amplitude Test"
    section of Chapter 4.],
) <fig:supp-h1-vert-amp-zoom>

=== Phase Test <supp:hyp1-vert-phase>

#figure(
  img("H1_016_Vertical_--_Phase-plane_shift_estimation_--_Kirchhoff____Bas.png", width: 100%),
  caption: [Kirchhoff-migrated phase-plane shift estimation, vertical study:
    Baseline versus each of the six displacement scenarios. Referenced from
    the "Vertical Movement -- Phase Test" section of Chapter 4.],
) <fig:supp-h1-vert-phase-kirchhoff>

#figure(
  img("H1_017_Vertical_--_Phase-plane_shift_estimation_--_Gazdag____Baseli.png", width: 100%),
  caption: [Gazdag-migrated phase-plane shift estimation, vertical study:
    Baseline versus each of the six displacement scenarios. Referenced from
    the "Vertical Movement -- Phase Test" section of Chapter 4.],
) <fig:supp-h1-vert-phase-gazdag>

#figure(
  img("H1_018_Vertical_--_Phase-plane_shift_estimation_--_Back-prop____Bas.png", width: 100%),
  caption: [Back-propagation-migrated phase-plane shift estimation, vertical
    study: Baseline versus each of the six displacement scenarios.
    Referenced from the "Vertical Movement -- Phase Test" section of
    Chapter 4.],
) <fig:supp-h1-vert-phase-backprop>

== Diagonal Movement <supp:hyp1-diagonal>

=== Amplitude Test <supp:hyp1-diag-amplitude>

#figure(
  img("H1_024_Diagonal_--_Amplitude_PSF_zoom_Baseline_vs_Monitor_RAW_migra.png", width: 80%),
  caption: [Zoomed Baseline-versus-Monitor point-spread-function comparison
    for the raw migrated diagonal images, sampled along the motion direction,
    all scenarios and migration methods. Referenced from the "Diagonal
    Movement -- Amplitude Test" section of Chapter 4.],
) <fig:supp-h1-diag-amp-zoom>

=== Phase Test <supp:hyp1-diag-phase>

#figure(
  img("H1_025_Diagonal_--_Phase-plane_shift_estimation_--_Kirchhoff____Bas.png", width: 100%),
  caption: [Kirchhoff-migrated phase-plane shift estimation, diagonal study:
    Baseline versus each of the five displacement scenarios. Referenced
    from the "Diagonal Movement -- Phase Test" section of Chapter 4.],
) <fig:supp-h1-diag-phase-kirchhoff>

#figure(
  img("H1_026_Diagonal_--_Phase-plane_shift_estimation_--_Gazdag____Baseli.png", width: 100%),
  caption: [Gazdag-migrated phase-plane shift estimation, diagonal study:
    Baseline versus each of the five displacement scenarios. Referenced
    from the "Diagonal Movement -- Phase Test" section of Chapter 4.],
) <fig:supp-h1-diag-phase-gazdag>

#figure(
  img("H1_027_Diagonal_--_Phase-plane_shift_estimation_--_Back-prop____Bas.png", width: 100%),
  caption: [Back-propagation-migrated phase-plane shift estimation, diagonal
    study: Baseline versus each of the five displacement scenarios.
    Referenced from the "Diagonal Movement -- Phase Test" section of
    Chapter 4.],
) <fig:supp-h1-diag-phase-backprop>

== Fluid Flow <supp:hyp1-fluidflow>

=== Amplitude Test <supp:hyp1-ff-amplitude>

#figure(
  img("H1_033_FluidFlow_--_Amplitude_PSF_zoom_Baseline_vs_Monitor_RAW_migr.png", width: 80%),
  caption: [Zoomed Baseline-versus-Monitor point-spread-function comparison
    for the raw migrated fluid-flow images, all scenarios and migration
    methods. Referenced from the "Fluid Flow -- Amplitude Test" section of
    Chapter 4.],
) <fig:supp-h1-ff-amp-zoom>

=== Phase Test <supp:hyp1-ff-phase>

#figure(
  img("H1_034_FluidFlow_--_Phase-plane_shift_estimation_--_Kirchhoff____Ba.png", width: 100%),
  caption: [Kirchhoff-migrated phase-plane shift estimation, fluid-flow
    study: Baseline versus each of the seven displacement scenarios.
    Referenced from the "Fluid Flow -- Phase Test" section of Chapter 4.],
) <fig:supp-h1-ff-phase-kirchhoff>

#figure(
  img("H1_035_FluidFlow_--_Phase-plane_shift_estimation_--_Gazdag____Basel.png", width: 100%),
  caption: [Gazdag-migrated phase-plane shift estimation, fluid-flow study:
    Baseline versus each of the seven displacement scenarios. Referenced
    from the "Fluid Flow -- Phase Test" section of Chapter 4.],
) <fig:supp-h1-ff-phase-gazdag>

#figure(
  img("H1_036_FluidFlow_--_Phase-plane_shift_estimation_--_Back-prop____Ba.png", width: 100%),
  caption: [Back-propagation-migrated phase-plane shift estimation,
    fluid-flow study: Baseline versus each of the seven displacement
    scenarios. Referenced from the "Fluid Flow -- Phase Test" section of
    Chapter 4.],
) <fig:supp-h1-ff-phase-backprop>
