#import "../template.typ": *

= Hypothesis 2: Supplementary Figures <supp:hyp2>

This chapter collects the Baseline-versus-Monitor point-spread-function
(PSF) comparisons and phase-plane shift-estimation diagnostics for
Chapter 5, "Hypothesis 2: Which Migration Technique Is Best Suited for
Noise-Robust Phase-Plane Tracking?", of the main thesis. As in
@supp:hyp1, both figure types recur for every scenario set and migration
method and are already summarised in the main text by the
Rayleigh-criterion-ratio and WLS-displacement-error tables, so they are
collected here rather than reproduced in the main chapter. Section numbers
below mirror Chapter 5's own section names one-for-one.

#block(breakable: false)[
== Noisy Lateral Movement <supp:hyp2-lateral>

=== Amplitude Test <supp:hyp2-lat-amplitude>

#figure(
  img("H2_010_Lateral_Noisy_--_Amplitude_PSF_zoom_Baseline_vs_Monitor_RAW.png", width: 65%),
  caption: [Zoomed Baseline-versus-Monitor point-spread-function comparison
    for the raw _noisy_ migrated lateral images, all scenarios and migration
    methods. Referenced from the "Noisy Lateral Movement -- Amplitude Test"
    section of Chapter 5.],
) <fig:supp-h2-lat-amp-zoom>
]

#block(breakable: false)[
=== Phase Test <supp:hyp2-lat-phase>

#figure(
  img("H2_011_Lateral_Noisy_--_Phase-plane_shift_estimation_--_Kirchhoff.png", width: 75%),
  caption: [Kirchhoff-migrated phase-plane shift estimation, _noisy_ lateral
    study: Baseline versus each of the seven displacement scenarios.
    Referenced from the "Noisy Lateral Movement -- Phase Test" section of
    Chapter 5.],
) <fig:supp-h2-lat-phase-kirchhoff>
]

#figure(
  img("H2_012_Lateral_Noisy_--_Phase-plane_shift_estimation_--_Gazdag____B.png", width: 75%),
  caption: [Gazdag-migrated phase-plane shift estimation, _noisy_ lateral
    study: Baseline versus each of the seven displacement scenarios.
    Referenced from the "Noisy Lateral Movement -- Phase Test" section of
    Chapter 5.],
) <fig:supp-h2-lat-phase-gazdag>

#figure(
  img("H2_013_Lateral_Noisy_--_Phase-plane_shift_estimation_--_Back-prop.png", width: 75%),
  caption: [Back-propagation-migrated (sign-bit) phase-plane shift
    estimation, _noisy_ lateral study: Baseline versus each of the seven
    displacement scenarios. Referenced from the "Noisy Lateral Movement --
    Phase Test" section of Chapter 5.],
) <fig:supp-h2-lat-phase-backprop>

#block(breakable: false)[
== Noisy Vertical Movement <supp:hyp2-vertical>

=== Migration Results <supp:hyp2-vert-migration>

#figure(
  img("H2_014_Vertical_--_TimeLapse_Migration_Comparison_Noisy_--_Signed_A.png", width: 69%),
  caption: [Signed time-lapse-difference amplitude (monitor-minus-baseline)
    for all three migration algorithms, vertical time-lapse study, _noisy_
    data. Referenced from the "Noisy Vertical Movement" section of
    Chapter 5.],
) <fig:h2-vert-summary-amp>
]

#block(breakable: false)[
=== Amplitude Test <supp:hyp2-vert-amplitude>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$1 lambda$],     [6.508], [5.695], [6.528],
    [$1\/2 lambda$],  [3.977], [3.480], [3.406],
    [$1\/4 lambda$],  [2.169], [1.582], [1.703],
    [$1\/8 lambda$],  [1.085], [0.633], [0.851],
    [$1\/16 lambda$], [1.085], [0.316], [0.568],
    [$1\/32 lambda$], [0.362], [0.0],   [0.284],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Vertical Rayleigh-criterion ratio, _noisy_ data (compare
    @tab:h1-vert-amp). Referenced from the "Noisy Vertical Movement --
    Amplitude Test" section of Chapter 5.],
  kind: table,
) <tab:h2-vert-amp>
]

#figure(
  img("H2_015_Vertical_Noisy_--_Amplitude_PSF_zoom_Baseline_vs_Monitor_RAW.png", width: 65%),
  caption: [Zoomed Baseline-versus-Monitor point-spread-function comparison
    for the raw _noisy_ migrated vertical images, all scenarios and
    migration methods. Referenced from the "Noisy Vertical Movement --
    Amplitude Test" section of Chapter 5.],
) <fig:supp-h2-vert-amp-zoom>

#block(breakable: false)[
=== Phase Test <supp:hyp2-vert-phase>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$1 lambda$],     [$-79.52 (-70.4%)$],  [$-105.61 (-93.5%)$],  [$-87.51 (-77.4%)$],
    [$1\/2 lambda$],  [$-63.82 (-114.0%)$], [$-76.54 (-136.7%)$],  [$-86.27 (-154.0%)$],
    [$1\/4 lambda$],  [$-0.82 (-2.9%)$],    [$-0.83 (-3.0%)$],     [$-0.01 (-0.0%)$],
    [$1\/8 lambda$],  [$-0.50 (-3.5%)$],    [$-2.93 (-20.9%)$],    [$+0.02 (+0.1%)$],
    [$1\/16 lambda$], [$-0.53 (-7.5%)$],    [$+6.95 (+99.3%)$],    [$+0.08 (+1.2%)$],
    [$1\/32 lambda$], [$-0.42 (-10.5%)$],   [$-1.31 (-32.7%)$],    [$+0.02 (+0.5%)$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Vertical phase-plane WLS displacement error, $#Dz$
    (estimated $-$ true), in millimetres, with the equivalent percentage of
    the true displacement in parentheses, _noisy_ data (compare
    @tab:h1-vert-phase). Referenced from the "Noisy Vertical Movement --
    Phase Test" section of Chapter 5.],
  kind: table,
) <tab:h2-vert-phase>
]

#figure(
  img("H2_016_Vertical_Noisy_--_Phase-plane_shift_estimation_--_Kirchhoff.png", width: 75%),
  caption: [Kirchhoff-migrated phase-plane shift estimation, _noisy_
    vertical study: Baseline versus each of the six displacement scenarios.
    Referenced from the "Noisy Vertical Movement -- Phase Test" section of
    Chapter 5.],
) <fig:supp-h2-vert-phase-kirchhoff>

#figure(
  img("H2_017_Vertical_Noisy_--_Phase-plane_shift_estimation_--_Gazdag.png", width: 75%),
  caption: [Gazdag-migrated phase-plane shift estimation, _noisy_ vertical
    study: Baseline versus each of the six displacement scenarios.
    Referenced from the "Noisy Vertical Movement -- Phase Test" section of
    Chapter 5.],
) <fig:supp-h2-vert-phase-gazdag>

#figure(
  img("H2_018_Vertical_Noisy_--_Phase-plane_shift_estimation_--_Back-prop.png", width: 75%),
  caption: [Back-propagation-migrated (sign-bit) phase-plane shift
    estimation, _noisy_ vertical study: Baseline versus each of the six
    displacement scenarios. Referenced from the "Noisy Vertical Movement --
    Phase Test" section of Chapter 5.],
) <fig:supp-h2-vert-phase-backprop>

#block(breakable: false)[
== Noisy Diagonal Movement <supp:hyp2-diagonal>

=== Migration Results <supp:hyp2-diag-migration>

#figure(
  img("H2_019_Diagonal_--_TimeLapse_Migration_Comparison_Noisy_--_Signed_A.png", width: 69%),
  caption: [Signed time-lapse-difference amplitude (monitor-minus-baseline)
    for all three migration algorithms, diagonal time-lapse study, _noisy_
    data. Referenced from the "Noisy Diagonal Movement" section of
    Chapter 5.],
) <fig:h2-diag-summary-amp>
]

#block(breakable: false)[
=== Amplitude Test <supp:hyp2-diag-amplitude>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [1 ($1 lambda_x$, $1\/2 lambda_z$)],   [5.295], [4.157], [3.430],
    [2 ($1\/2 lambda_x$, $1\/4 lambda_z$)], [2.492], [2.019], [1.663],
    [3 ($1\/4 lambda_x$, $1\/8 lambda_z$)], [1.090], [1.069], [0.831],
    [4 ($1\/8 lambda_x$, $1\/16 lambda_z$)],[0.623], [0.594], [0.312],
    [5 ($1\/16 lambda_x$, $1\/32 lambda_z$)],[0.156],[0.238], [0.104],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Diagonal Rayleigh-criterion ratio, _noisy_ data (compare
    @tab:h1-diag-amp). Referenced from the "Noisy Diagonal Movement --
    Amplitude Test" section of Chapter 5.],
  kind: table,
) <tab:h2-diag-amp>
]

#figure(
  img("H2_020_Diagonal_Noisy_--_Amplitude_PSF_zoom_Baseline_vs_Monitor_RAW.png", width: 65%),
  caption: [Zoomed Baseline-versus-Monitor point-spread-function comparison
    for the raw _noisy_ migrated diagonal images, all scenarios and
    migration methods. Referenced from the "Noisy Diagonal Movement --
    Amplitude Test" section of Chapter 5.],
) <fig:supp-h2-diag-amp-zoom>

#block(breakable: false)[
=== Phase Test <supp:hyp2-diag-phase>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [1], [$-56.97 (-101.7%)$], [$-58.82 (-105.0%)$], [$-55.46 (-99.0%)$],
    [2], [$-38.78 (-138.5%)$], [$-5.45 (-19.5%)$],   [$-12.38 (-44.2%)$],
    [3], [$-1.00 (-7.1%)$],   [$-0.27 (-1.9%)$],    [$-0.02 (-0.1%)$],
    [4], [$-0.23 (-3.3%)$],   [$-2.16 (-30.8%)$],   [$+0.01 (+0.2%)$],
    [5], [$-0.36 (-9.1%)$],   [$+6.77 (+169.3%)$],  [$+0.08 (+2.0%)$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Diagonal phase-plane WLS displacement error, $#Dz$
    (estimated $-$ true), in millimetres, with the equivalent percentage of
    the true displacement in parentheses, _noisy_ data (compare
    @tab:h1-diag-phase-dz). Referenced from the "Noisy Diagonal Movement --
    Phase Test" section of Chapter 5.],
  kind: table,
) <tab:h2-diag-phase-dz>
]

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [1], [$-114.95 (-101.7%)$], [$-117.97 (-104.4%)$], [$-150.68 (-133.3%)$],
    [2], [$-46.55 (-83.1%)$],   [$-44.21 (-79.0%)$],   [$-81.95 (-146.3%)$],
    [3], [$-1.05 (-3.8%)$],     [$-15.08 (-53.9%)$],   [$-0.02 (-0.1%)$],
    [4], [$+0.04 (+0.3%)$],     [$-41.00 (-292.9%)$],  [$+0.01 (+0.0%)$],
    [5], [$-0.05 (-0.8%)$],     [$-12.00 (-171.4%)$],  [$+0.13 (+1.9%)$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Diagonal phase-plane WLS displacement error, $#Dx$
    (estimated $-$ true), in millimetres, with the equivalent percentage of
    the true displacement in parentheses, _noisy_ data (compare
    @tab:h1-diag-phase-dx). Referenced from the "Noisy Diagonal Movement --
    Phase Test" section of Chapter 5.],
  kind: table,
) <tab:h2-diag-phase-dx>

#figure(
  img("H2_021_Diagonal_Noisy_--_Phase-plane_shift_estimation_--_Kirchhoff.png", width: 75%),
  caption: [Kirchhoff-migrated phase-plane shift estimation, _noisy_
    diagonal study: Baseline versus each of the five displacement scenarios.
    Referenced from the "Noisy Diagonal Movement -- Phase Test" section of
    Chapter 5.],
) <fig:supp-h2-diag-phase-kirchhoff>

#figure(
  img("H2_022_Diagonal_Noisy_--_Phase-plane_shift_estimation_--_Gazdag.png", width: 75%),
  caption: [Gazdag-migrated phase-plane shift estimation, _noisy_ diagonal
    study: Baseline versus each of the five displacement scenarios.
    Referenced from the "Noisy Diagonal Movement -- Phase Test" section of
    Chapter 5.],
) <fig:supp-h2-diag-phase-gazdag>

#figure(
  img("H2_023_Diagonal_Noisy_--_Phase-plane_shift_estimation_--_Back-prop.png", width: 75%),
  caption: [Back-propagation-migrated (sign-bit) phase-plane shift
    estimation, _noisy_ diagonal study: Baseline versus each of the five
    displacement scenarios. Referenced from the "Noisy Diagonal Movement --
    Phase Test" section of Chapter 5.],
) <fig:supp-h2-diag-phase-backprop>

#block(breakable: false)[
== Noisy Fluid Flow <supp:hyp2-fluidflow>

=== Migration Results <supp:hyp2-ff-migration>

#figure(
  img("H2_024_FluidFlow_--_TimeLapse_Migration_Comparison_Noisy_--_Signed.png", width: 69%),
  caption: [Signed time-lapse-difference amplitude (monitor-minus-baseline)
    for all three migration algorithms, fluid-flow time-lapse study, _noisy_
    data. Referenced from the "Noisy Fluid Flow" section of Chapter 5.],
) <fig:h2-ff-summary-amp>
]

#block(breakable: false)[
=== Amplitude Test <supp:hyp2-ff-amplitude>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$2 lambda$],     [4.853], [0.036], [5.450],
    [$1 lambda$],     [2.427], [0.214], [2.725],
    [$1\/2 lambda$],  [1.103], [0.605], [1.239],
    [$1\/4 lambda$],  [0.441], [0.107], [0.495],
    [$1\/8 lambda$],  [0.221], [0.071], [0.248],
    [$1\/16 lambda$], [0.0],   [0.0],   [0.0],
    [$1\/32 lambda$], [0.0],   [0.036], [0.0],
    table.hline(stroke: 0.7pt),
  ),
  caption: [FluidFlow Rayleigh-criterion ratio, _noisy_ data (compare
    @tab:h1-ff-amp). Referenced from the "Noisy Fluid Flow -- Amplitude
    Test" section of Chapter 5.],
  kind: table,
) <tab:h2-ff-amp>
]

#figure(
  img("H2_025_FluidFlow_Noisy_--_Amplitude_PSF_zoom_Baseline_vs_Monitor_RA.png", width: 65%),
  caption: [Zoomed Baseline-versus-Monitor point-spread-function comparison
    for the raw _noisy_ migrated fluid-flow images, all scenarios and
    migration methods. Referenced from the "Noisy Fluid Flow -- Amplitude
    Test" section of Chapter 5.],
) <fig:supp-h2-ff-amp-zoom>

#block(breakable: false)[
=== Phase Test <supp:hyp2-ff-phase>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$2 lambda$],     [$-211.14 (-93.8%)$], [$-225.00 (-100.0%)$], [$-225.00 (-100.0%)$],
    [$1 lambda$],     [$-105.94 (-93.8%)$], [$+7.52 (+6.7%)$],     [$+11.14 (+9.9%)$],
    [$1\/2 lambda$],  [$-57.02 (-101.8%)$], [$-0.33 (-0.6%)$],     [$-0.63 (-1.1%)$],
    [$1\/4 lambda$],  [$-0.72 (-2.6%)$],    [$+1.15 (+4.1%)$],     [$-0.28 (-1.0%)$],
    [$1\/8 lambda$],  [$-0.52 (-3.7%)$],    [$+1.50 (+10.7%)$],    [$+0.02 (+0.2%)$],
    [$1\/16 lambda$], [$-0.91 (-13.1%)$],   [$+0.50 (+7.2%)$],     [$-0.01 (-0.1%)$],
    [$1\/32 lambda$], [$-1.27 (-31.7%)$],   [$+0.96 (+24.0%)$],    [$-0.05 (-1.4%)$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [FluidFlow phase-plane WLS front-displacement error, $#Dx$
    (estimated $-$ true), in millimetres, with the equivalent percentage of
    the true displacement in parentheses, _noisy_ data (compare
    @tab:h1-ff-phase). Referenced from the "Noisy Fluid Flow -- Phase Test"
    section of Chapter 5.],
  kind: table,
) <tab:h2-ff-phase>
]

#figure(
  img("H2_026_FluidFlow_Noisy_--_Phase-plane_shift_estimation_--_Kirchhoff.png", width: 75%),
  caption: [Kirchhoff-migrated phase-plane shift estimation, _noisy_
    fluid-flow study: Baseline versus each of the seven displacement
    scenarios. Referenced from the "Noisy Fluid Flow -- Phase Test" section
    of Chapter 5.],
) <fig:supp-h2-ff-phase-kirchhoff>

#figure(
  img("H2_027_FluidFlow_Noisy_--_Phase-plane_shift_estimation_--_Gazdag.png", width: 75%),
  caption: [Gazdag-migrated phase-plane shift estimation, _noisy_ fluid-flow
    study: Baseline versus each of the seven displacement scenarios.
    Referenced from the "Noisy Fluid Flow -- Phase Test" section of
    Chapter 5.],
) <fig:supp-h2-ff-phase-gazdag>

#figure(
  img("H2_028_FluidFlow_Noisy_--_Phase-plane_shift_estimation_--_Back-prop.png", width: 75%),
  caption: [Back-propagation-migrated (sign-bit) phase-plane shift
    estimation, _noisy_ fluid-flow study: Baseline versus each of the seven
    displacement scenarios. Referenced from the "Noisy Fluid Flow -- Phase
    Test" section of Chapter 5.],
) <fig:supp-h2-ff-phase-backprop>

#block(breakable: false)[
== Noise Characterization <supp:hyp2-noise-stages>

#figure(
  img("H2_001_Noise_amplitude_distribution_by_processing_stage.png"),
  caption: [Noise amplitude distribution at every one of the eleven tracked
    processing stages, each with a Gaussian reference overlay. Referenced
    from the "Noise Creation" section of Chapter 5.],
) <fig:h2-noise-stages>
]
