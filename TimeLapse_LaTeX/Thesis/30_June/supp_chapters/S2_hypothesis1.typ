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

=== Model Set Up, B-scans, and Migration Results <supp:hyp1-vert-setup>

#figure(
  subfigs(cols: 1,
    img("H1_010_Vertical_Movement_--_Model_Set_Up.png"),
    img("H1_011_Vertical_--_target_position_all_7_scenarios.png", width: 70%),
    table(
      columns: (auto, auto, auto),
      stroke: none,
      inset: (x: 0.8em, y: 0.3em),
      table.hline(stroke: 0.7pt),
      [*Scenario*], [*$#Dx$*], [*$#Dz$*],
      table.hline(stroke: 0.4pt),
      [Baseline], [$0$], [$0$],
      [$1 lambda$], [$0$], [$1 lambda$],
      [$1\/2 lambda$], [$0$], [$1\/2 lambda$],
      [$1\/4 lambda$], [$0$], [$1\/4 lambda$],
      [$1\/8 lambda$], [$0$], [$1\/8 lambda$],
      [$1\/16 lambda$], [$0$], [$1\/16 lambda$],
      [$1\/32 lambda$], [$0$], [$1\/32 lambda$],
      table.hline(stroke: 0.7pt),
    ),
  ),
  caption: [Forward-model setup for the vertical time-lapse study: (a) the
    gprMax domain and grid, fixed at $x = 2.0 "m"$; (b) the target depth for
    all seven scenarios (Baseline plus six displacements), colour-coded from
    $1 lambda$ to $1 \/ 32 lambda$; (c) the seven displacement scenarios: the
    scatterer moves downward in depth ($#Dz$) while its lateral position is
    fixed ($#Dx = 0$), swept from $1 lambda$ down to $1\/32 lambda$.
    Referenced from the "Vertical Movement" section of Chapter 4.],
) <fig:h1-vert-setup>

#figure(
  subfigs(cols: 1,
    img("H1_012_Vertical_Movement_--_Background-Subtracted_B-Scans.png"),
    img("H1_013_Vertical_--_TimeLapse_Migration_Comparison_Clean_--_Signed_A.png", width: 95%),
  ),
  caption: [(a) Background-subtracted B-scans for the vertical time-lapse
    study, all seven scenarios; (b) signed time-lapse-difference amplitude
    (monitor-minus-baseline) for all three migration algorithms, clean data.
    Referenced from the "Vertical Movement" section of Chapter 4.],
) <fig:h1-vert-summary-amp>

=== Amplitude Test <supp:hyp1-vert-amplitude>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$1 lambda$],     [6.096], [7.564], [6.528],
    [$1\/2 lambda$],  [3.180], [3.946], [3.406],
    [$1\/4 lambda$],  [1.590], [1.973], [1.703],
    [$1\/8 lambda$],  [0.795], [0.987], [0.851],
    [$1\/16 lambda$], [0.530], [0.658], [0.568],
    [$1\/32 lambda$], [0.265], [0.329], [0.284],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Vertical Rayleigh-criterion ratio (Baseline--Monitor peak
    separation / Baseline FWHM). Referenced from the "Vertical Movement --
    Amplitude Test" section of Chapter 4.],
  kind: table,
) <tab:h1-vert-amp>

#figure(
  img("H1_015_Vertical_--_Amplitude_PSF_zoom_Baseline_vs_Monitor_RAW_migra.png", width: 80%),
  caption: [Zoomed Baseline-versus-Monitor point-spread-function comparison
    for the raw migrated vertical images, all scenarios and migration
    methods. Referenced from the "Vertical Movement -- Amplitude Test"
    section of Chapter 4.],
) <fig:supp-h1-vert-amp-zoom>

=== Phase Test <supp:hyp1-vert-phase>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$1 lambda$],     [$-100.01 (-88.5%)$], [$-90.08 (-79.7%)$],  [$-86.93 (-76.9%)$],
    [$1\/2 lambda$],  [$-59.67 (-106.6%)$], [$-80.37 (-143.5%)$], [$-83.99 (-150.0%)$],
    [$1\/4 lambda$],  [$+0.03 (+0.1%)$],    [$-0.01 (-0.0%)$],    [$-0.00 (-0.0%)$],
    [$1\/8 lambda$],  [$+0.03 (+0.2%)$],    [$-0.01 (-0.0%)$],    [$-0.00 (-0.0%)$],
    [$1\/16 lambda$], [$+0.02 (+0.2%)$],    [$-0.00 (-0.0%)$],    [$-0.00 (-0.0%)$],
    [$1\/32 lambda$], [$+0.01 (+0.3%)$],    [$-0.00 (-0.0%)$],    [$-0.00 (-0.0%)$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Vertical phase-plane WLS displacement error, $#Dz$
    (estimated $-$ true), in millimetres, with the equivalent percentage of
    the true displacement in parentheses. Referenced from the "Vertical
    Movement -- Phase Test" section of Chapter 4.],
  kind: table,
) <tab:h1-vert-phase>

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

=== Model Set Up, B-scans, and Migration Results <supp:hyp1-diag-setup>

#figure(
  subfigs(cols: 1,
    img("H1_019_Diagonal_Movement_--_Model_Set_Up.png"),
    img("H1_020_Diagonal_--_target_position_all_6_scenarios.png", width: 70%),
    table(
      columns: (auto, auto, auto),
      stroke: none,
      inset: (x: 0.8em, y: 0.3em),
      table.hline(stroke: 0.7pt),
      [*Scenario*], [*$#Dx$ (right)*], [*$#Dz$ (down)*],
      table.hline(stroke: 0.4pt),
      [Baseline], [$0$], [$0$],
      [1], [$1 lambda$], [$1\/2 lambda$],
      [2], [$1\/2 lambda$], [$1\/4 lambda$],
      [3], [$1\/4 lambda$], [$1\/8 lambda$],
      [4], [$1\/8 lambda$], [$1\/16 lambda$],
      [5], [$1\/16 lambda$], [$1\/32 lambda$],
      table.hline(stroke: 0.7pt),
    ),
  ),
  caption: [Forward-model setup for the diagonal time-lapse study: (a) the
    gprMax domain and grid; (b) the target position for all six scenarios
    (Baseline plus five displacements) along the $2$:$1$ diagonal; (c) the
    five displacement scenarios: lateral ($#Dx$, rightward) and vertical
    ($#Dz$, downward) components, both expressed as fractions of the
    dominant wavelength $lambda$, in a fixed $2$:$1$ ratio. Referenced from
    the "Diagonal Movement" section of Chapter 4.],
) <fig:h1-diag-setup>

#figure(
  subfigs(cols: 1,
    img("H1_021_Diagonal_Movement_--_Background-Subtracted_B-Scans.png"),
    img("H1_022_Diagonal_--_TimeLapse_Migration_Comparison_Clean_--_Signed_A.png", width: 95%),
  ),
  caption: [(a) Background-subtracted B-scans for the diagonal time-lapse
    study, all six scenarios; (b) signed time-lapse-difference amplitude
    (monitor-minus-baseline) for all three migration algorithms, clean data,
    sampled along the diagonal motion direction. Referenced from the
    "Diagonal Movement" section of Chapter 4.],
) <fig:h1-diag-summary-amp>

=== Amplitude Test <supp:hyp1-diag-amplitude>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [1 ($1 lambda_x$, $1\/2 lambda_z$)],   [5.065], [4.169], [3.429],
    [2 ($1\/2 lambda_x$, $1\/4 lambda_z$)], [2.392], [2.025], [1.663],
    [3 ($1\/4 lambda_x$, $1\/8 lambda_z$)], [1.266], [0.953], [0.831],
    [4 ($1\/8 lambda_x$, $1\/16 lambda_z$)],[0.563], [0.596], [0.312],
    [5 ($1\/16 lambda_x$, $1\/32 lambda_z$)],[0.422],[0.238], [0.104],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Diagonal Rayleigh-criterion ratio (Baseline--Monitor peak
    separation / Baseline FWHM, sampled along the motion direction).
    Referenced from the "Diagonal Movement -- Amplitude Test" section of
    Chapter 4.],
  kind: table,
) <tab:h1-diag-amp>

#figure(
  img("H1_024_Diagonal_--_Amplitude_PSF_zoom_Baseline_vs_Monitor_RAW_migra.png", width: 80%),
  caption: [Zoomed Baseline-versus-Monitor point-spread-function comparison
    for the raw migrated diagonal images, sampled along the motion direction,
    all scenarios and migration methods. Referenced from the "Diagonal
    Movement -- Amplitude Test" section of Chapter 4.],
) <fig:supp-h1-diag-amp-zoom>

=== Phase Test <supp:hyp1-diag-phase>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [1], [$-59.49 (-106.2%)$], [$-56.09 (-100.2%)$], [$-54.30 (-97.0%)$],
    [2], [$-49.09 (-175.3%)$], [$-7.48 (-26.7%)$],   [$-12.54 (-44.8%)$],
    [3], [$-1.51 (-10.8%)$],   [$-0.01 (-0.0%)$],    [$-0.00 (-0.0%)$],
    [4], [$+0.01 (+0.2%)$],    [$-0.00 (-0.0%)$],    [$-0.00 (-0.0%)$],
    [5], [$+0.01 (+0.3%)$],    [$-0.00 (-0.0%)$],    [$-0.00 (-0.0%)$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Diagonal phase-plane WLS displacement error, $#Dz$
    (estimated $-$ true), in millimetres, with the equivalent percentage of
    the true displacement in parentheses. Referenced from the "Diagonal
    Movement -- Phase Test" section of Chapter 4.],
  kind: table,
) <tab:h1-diag-phase-dz>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [1], [$-121.17 (-107.2%)$], [$-184.50 (-163.3%)$], [$-146.94 (-130.0%)$],
    [2], [$-47.25 (-84.4%)$],   [$-70.89 (-126.6%)$],  [$-83.58 (-149.3%)$],
    [3], [$-1.32 (-4.7%)$],     [$-0.00 (-0.0%)$],     [$+0.02 (+0.1%)$],
    [4], [$-0.02 (-0.1%)$],     [$-0.00 (-0.0%)$],     [$+0.01 (+0.1%)$],
    [5], [$-0.01 (-0.1%)$],     [$-0.00 (-0.0%)$],     [$+0.01 (+0.1%)$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Diagonal phase-plane WLS displacement error, $#Dx$
    (estimated $-$ true), in millimetres, with the equivalent percentage of
    the true displacement in parentheses. Referenced from the "Diagonal
    Movement -- Phase Test" section of Chapter 4.],
  kind: table,
) <tab:h1-diag-phase-dx>

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

=== B-scans and Migration Results <supp:hyp1-ff-setup>

#figure(
  subfigs(cols: 1,
    img("H1_030_FluidFlow_Movement_--_Background-Subtracted_B-Scans.png"),
    img("H1_031_FluidFlow_--_TimeLapse_Migration_Comparison_Clean_--_Signed.png", width: 95%),
  ),
  caption: [(a) Background-subtracted B-scans for the fluid-flow time-lapse
    study, all eight scenarios; (b) signed time-lapse-difference amplitude
    (monitor-minus-baseline) for all three migration algorithms, clean data.
    Referenced from the "Fluid Flow" section of Chapter 4.],
) <fig:h1-ff-summary-amp>

=== Amplitude Test <supp:hyp1-ff-amplitude>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$2 lambda$],     [1.978], [0.781], [5.442],
    [$1 lambda$],     [1.032], [0.391], [2.721],
    [$1\/2 lambda$],  [0.516], [0.178], [1.237],
    [$1\/4 lambda$],  [0.258], [0.107], [0.495],
    [$1\/8 lambda$],  [0.172], [0.036], [0.247],
    [$1\/16 lambda$], [0.086], [0.0],   [0.0],
    [$1\/32 lambda$], [0.086], [0.0],   [0.0],
    table.hline(stroke: 0.7pt),
  ),
  caption: [FluidFlow Rayleigh-criterion ratio (Baseline--Monitor peak
    separation / Baseline FWHM). Referenced from the "Fluid Flow --
    Amplitude Test" section of Chapter 4.],
  kind: table,
) <tab:h1-ff-amp>

#figure(
  img("H1_033_FluidFlow_--_Amplitude_PSF_zoom_Baseline_vs_Monitor_RAW_migr.png", width: 80%),
  caption: [Zoomed Baseline-versus-Monitor point-spread-function comparison
    for the raw migrated fluid-flow images, all scenarios and migration
    methods. Referenced from the "Fluid Flow -- Amplitude Test" section of
    Chapter 4.],
) <fig:supp-h1-ff-amp-zoom>

=== Phase Test <supp:hyp1-ff-phase>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$2 lambda$],     [$-293.35 (-130.4%)$], [$-225.00 (-100.0%)$], [$-225.00 (-100.0%)$],
    [$1 lambda$],     [$-84.57 (-74.8%)$],   [$+0.15 (+0.1%)$],     [$+1.18 (+1.0%)$],
    [$1\/2 lambda$],  [$+0.47 (+0.8%)$],     [$+0.61 (+1.1%)$],     [$+0.90 (+1.6%)$],
    [$1\/4 lambda$],  [$+0.52 (+1.8%)$],     [$-0.07 (-0.2%)$],     [$+0.52 (+1.9%)$],
    [$1\/8 lambda$],  [$+0.29 (+2.0%)$],     [$-0.04 (-0.3%)$],     [$+0.07 (+0.5%)$],
    [$1\/16 lambda$], [$+0.34 (+4.8%)$],     [$-0.02 (-0.3%)$],     [$+0.04 (+0.6%)$],
    [$1\/32 lambda$], [$+0.15 (+3.8%)$],     [$-0.01 (-0.3%)$],     [$+0.03 (+0.7%)$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [FluidFlow phase-plane WLS front-displacement error, $#Dx$
    (estimated $-$ true, centroid-corrected $times 2$ for Kirchhoff/Gazdag),
    in millimetres, with the equivalent percentage of the true displacement
    in parentheses. Referenced from the "Fluid Flow -- Phase Test" section
    of Chapter 4.],
  kind: table,
) <tab:h1-ff-phase>

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
