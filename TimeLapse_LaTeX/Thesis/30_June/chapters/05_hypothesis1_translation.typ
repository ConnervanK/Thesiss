#import "../template.typ": *

= Hypothesis 1: Can We Infer Subwavelength Movement from Phase Changes? <ch:hyp1>

@sec:meth-resolution established that amplitude-based detection of a
_stationary_ point scatterer fails below roughly half a wavelength of
separation. This chapter scales that finding to the time-lapse setting: if a
single scatterer moves by a sub-wavelength amount between a baseline and a
monitor survey, can that motion be inferred? The answer is that amplitude
differencing fails below the same resolution floor --- but examining the
_phase change_ in the two-dimensional Fourier domain (@ch:theory) reveals
displacements as small as $1 \/ 32 lambda$. The displacement is swept from
$2 lambda$ down to $1 \/ 32 lambda$, tested independently in three
directions: lateral (@sec:hyp1-lateral), vertical (@sec:hyp1-vertical), and
diagonal (@sec:hyp1-diagonal). @sec:hyp1-phaseplane then applies the global
2D weighted least-squares (WLS) phase-plane fit, and @sec:hyp1-h15 explores
an optional, local alternative.

#para-head[Hypothesis 1.] Can multi-dimensional phase-plane regression infer
lateral, vertical, and diagonal subwavelength displacements from time-lapse
migrated GPR images, at scales where amplitude differencing has already
failed?

#para-head[Hypothesis 1.5 (Optional).] Tracking the local phase gradients
$partial phi \/ partial x$ and $partial phi \/ partial y$ gives an
equivalent, simpler alternative to the global 2D phase-plane fit for
inferring horizontal and vertical translation.

== Lateral Movement <sec:hyp1-lateral>

A single PEC cylinder ($r = 28 "mm"$, depth $0.676 "m"$) is held fixed in
the baseline survey and displaced laterally by each of eight scenarios
(@tab:tl-scenarios) in the monitor survey, from $2 lambda$ down to
$1 \/ 32 lambda$ (@fig:tl-setup), using the same gprMax domain and grid as
@sec:meth-resolution.

#figure(
  table(
    columns: (auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*$#Dx$*], [*$#Dz$*],
    table.hline(stroke: 0.4pt),
    [Baseline], [$0$], [$0$],
    [$2 lambda$], [$2 lambda$], [$0$],
    [$1 lambda$], [$1 lambda$], [$0$],
    [$1\/2 lambda$], [$1\/2 lambda$], [$0$],
    [$1\/4 lambda$], [$1\/4 lambda$], [$0$],
    [$1\/8 lambda$], [$1\/8 lambda$], [$0$],
    [$1\/16 lambda$], [$1\/16 lambda$], [$0$],
    [$1\/32 lambda$], [$1\/32 lambda$], [$0$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Lateral displacement scenarios: the scatterer moves laterally
    ($#Dx$) while depth is fixed ($#Dz = 0$), swept from $2 lambda$ down to
    $1\/32 lambda$.],
  kind: table,
) <tab:tl-scenarios>

#figure(
  subfigs(cols: 1,
    img("TL_001_TimeLapse_Study__Model_Geometry__domain_4010_m_Δx__1_mm_PML.png"),
    img("TL_002_Moving_scatterer_s1__all_8_scenarios__r__28_mm_depth__0676_m.png", width: 55%),
  ),
  caption: [Forward-model setup for the lateral time-lapse study: (a) the
    gprMax domain and grid; (b) the lateral displacement of scatterer `s1`
    across all eight scenarios, $r = 28 "mm"$, depth $0.676 "m"$.],
) <fig:tl-setup>

@fig:tl-bscans shows the background-baseline and time-lapsed raw B-scans,
the background-subtracted result, and (for later use in @ch:hyp2) the effect
of adding synthetic Laplace-distributed noise at $10%$ of the signal
standard deviation.

#figure(
  subfigs(cols: 1,
    img("TL_003_GPR_B-Scans__Background_Baseline_and_TimeLapsed_Models.png"),
    img("TL_004_GPR_B-Scans__Background_Subtracted.png"),
    img("TL_005_GPR_B-Scans__With_Synthetic_Laplace_Noise_10_of_signal_std.png"),
  ),
  caption: [Raw B-scans for the lateral time-lapse study: (a) background and
    time-lapsed models; (b) background-subtracted; (c) with synthetic
    Laplace-distributed noise at $10%$ of the signal standard deviation
    (used in @ch:hyp2).],
) <fig:tl-bscans>

The same tapering and $t_0$-shift conditioning as @sec:meth-resolution
(@sec:meth-conditioning) is applied before migration (@fig:tl-taper).

#figure(
  subfigs(cols: 1,
    img("TL_008_Effect_of_Tapering_and_t0_Shift__2λ_dataset_single_trace.png"),
    img("TL_009_B-scan_effect_of_tapering_and_t0_shift__2λ_dataset.png"),
  ),
  caption: [Effect of tapering and the $t_0$ shift on the $2 lambda$ lateral
    displacement dataset: (a) a single trace; (b) the full B-scan.],
) <fig:tl-taper>

Per-method migrated images and time-lapse differences for Kirchhoff, Gazdag,
and back-propagation migration are given in @app:hyp1-methods
(@fig:tl-kirchhoff, @fig:tl-gazdag, @fig:tl-backprop). @fig:tl-summary-amp
below instead overlays the signed difference from all three algorithms
directly.

=== Amplitude-Based Detectability Summary <sec:tl-detectability>

@fig:tl-summary-amp overlays the signed time-lapse-difference amplitude from
all three algorithms, and @fig:tl-summary-psf plots the normalised lateral
point-spread function of that difference at the true scatterer depth, as a
function of displacement.

#figure(
  img("TL_024_TimeLapse_Migration_Comparison__Signed_Amplitude____f_c15_GH.png"),
  caption: [Signed time-lapse-difference amplitude for all three migration
    algorithms at $f_c = 1.5 "GHz"$.],
) <fig:tl-summary-amp>

#figure(
  img("TL_025_Normalised_Lateral_PSF__TimeLapse_Difference_at_True_Scatter.png"),
  caption: [Normalised lateral PSF of the time-lapse-difference image at the
    true scatterer depth, swept across lateral displacements from $2 lambda$
    to $1 \/ 32 lambda$.],
) <fig:tl-summary-psf>

#draftnote[state the smallest displacement at which the difference image in
@fig:tl-summary-psf still shows a clear, unambiguous peak above the numerical
background, and compare it with the amplitude resolution floor found for
stationary scatterers in @sec:meth-resolution. The noisy version of this
sweep, and the noise-mitigation strategies tested against it, are presented
in @ch:hyp2 rather than here.]

== Vertical Movement <sec:hyp1-vertical>

This section mirrors @sec:hyp1-lateral, replacing lateral displacement of
the scatterer with _vertical_ (depth) displacement. The two directions are
not expected to behave identically: @sec:th-duality showed that a migrated
image separates lateral spatial frequencies like a prism but does not separate
vertical ones, so the amplitude-based detectability established here is an
important point of comparison for the phase-based vertical results of
@sec:hyp1-phaseplane.

A single PEC cylinder ($r = 28 "mm"$) sits at lateral position $x = 2.0 "m"$
and is displaced _downward_ from its baseline depth across seven scenarios
(@tab:vtl-scenarios), from $1 lambda$ down to $1 \/ 32 lambda$
(@fig:vtl-setup), using the same domain and grid as @sec:meth-resolution and
@sec:hyp1-lateral.

#figure(
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
  caption: [Vertical displacement scenarios: the scatterer moves downward in
    depth ($#Dz$) while its lateral position is fixed ($#Dx = 0$), swept from
    $1 lambda$ down to $1\/32 lambda$.],
  kind: table,
) <tab:vtl-scenarios>

#figure(
  subfigs(cols: 1,
    img("VTL_001_Vertical_TimeLapse_Study__Model_Geometry__domain_4010_m_Δx.png"),
    img("VTL_002_Moving_scatterer__all_7_scenarios__r__28_mm_x__20_mBaseline.png", width: 45%),
  ),
  caption: [Forward-model setup for the vertical time-lapse study: (a) the
    gprMax domain and grid; (b) the vertical displacement of the scatterer
    across all seven scenarios, $r = 28 "mm"$, fixed at $x = 2.0 "m"$.],
) <fig:vtl-setup>

@fig:vtl-bscans shows the background-baseline and vertically time-lapsed raw
B-scans and the background-subtracted result, with the expected arrival times
for each scenario marked, and @fig:vtl-taper the effect of the standard
tapering and $t_0$-shift conditioning on the $1 lambda$ dataset.

#figure(
  subfigs(cols: 1,
    img("VTL_003_GPR_B-Scans__Background_Baseline_and_Vertical_TimeLapsed_Mod.png"),
    img("VTL_004_GPR_B-Scans__Background_Subtracted__green_dashed__expected_a.png"),
  ),
  caption: [Raw B-scans for the vertical time-lapse study: (a) background and
    time-lapsed models; (b) background-subtracted, with expected arrival
    times marked (green, dashed).],
) <fig:vtl-bscans>

#figure(
  subfigs(cols: 1,
    img("VTL_007_Effect_of_Tapering_and_t0_Shift__1λ_dataset_single_trace.png"),
    img("VTL_008_B-scan_effect_of_tapering_and_t0_shift__1λ_dataset.png"),
  ),
  caption: [Effect of tapering and the $t_0$ shift on the $1 lambda$ vertical
    displacement dataset: (a) a single trace; (b) the full B-scan.],
) <fig:vtl-taper>

Per-method migrated images and time-lapse differences for Kirchhoff, Gazdag,
and back-propagation migration are given in @app:hyp1-methods
(@fig:vtl-kirchhoff, @fig:vtl-gazdag, @fig:vtl-backprop). @fig:vtl-summary-amp
below instead overlays the signed difference from all three algorithms
directly.

=== Vertical Detectability Summary <sec:vtl-detectability>

@fig:vtl-summary-amp overlays the signed time-lapse-difference amplitude from
all three algorithms, and @fig:vtl-summary-psf plots the normalised
_vertical_ point-spread function of that difference at $x = 2.0 "m"$, as a
function of vertical displacement.

#figure(
  img("VTL_023_Vertical_TimeLapse_Migration_Comparison__Signed_Amplitude.png"),
  caption: [Signed time-lapse-difference amplitude for all three migration
    algorithms.],
) <fig:vtl-summary-amp>

#figure(
  img("VTL_024_Normalised_Vertical_PSF__TimeLapse_Difference_at_x__20_m.png", width: 55%),
  caption: [Normalised vertical PSF of the time-lapse-difference image at
    $x = 2.0 "m"$, swept across vertical displacements from $1 lambda$ to
    $1 \/ 32 lambda$.],
) <fig:vtl-summary-psf>

#draftnote[compare the smallest reliably-detected vertical displacement in
@fig:vtl-summary-psf against the lateral result of @fig:tl-summary-psf
(@sec:tl-detectability) and discuss whether the lateral/vertical asymmetry
predicted by the instantaneous-phase theory of @sec:th-duality is also visible
at the amplitude level, or only emerges once the phase-plane estimator of
@sec:hyp1-phaseplane is applied.]

== Diagonal Movement <sec:hyp1-diagonal>

The third and final translation direction combines the previous two: the
scatterer moves simultaneously laterally and vertically, along a fixed $2:1$
diagonal ($#Dx = 2 thin #Dz$ in every scenario), so that all scenarios lie
on the same line through the baseline position. A single PEC cylinder at
$x = 2.0 "m"$, baseline depth $0.676 "m"$, is displaced diagonally across
five scenarios (@tab:dtl-scenarios), using the same domain and grid as
@sec:meth-resolution, @sec:hyp1-lateral, and @sec:hyp1-vertical.

#figure(
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
  caption: [Diagonal displacement scenarios: lateral ($#Dx$, rightward) and
    vertical ($#Dz$, downward) components, both expressed as fractions of the
    dominant wavelength $lambda$.],
  kind: table,
) <tab:dtl-scenarios>

#figure(
  subfigs(cols: 1,
    img("DTL_001_Diagonal_TimeLapse_Study__Model_Geometry__domain_4010_m_Δx.png"),
    img("DTL_002_Diagonal_scatterer_path__all_6_scenarios__r__28_mm.png", width: 55%),
  ),
  caption: [Forward-model setup for the diagonal time-lapse study: (a) the
    gprMax domain and grid; (b) the diagonal displacement path of the scatterer
    across all scenarios, $r = 28 "mm"$.],
) <fig:dtl-setup>

@fig:dtl-bscans shows the background-baseline and diagonally time-lapsed raw
B-scans, the background-subtracted result, and the effect of the synthetic
Laplace noise used later in @ch:hyp2; @fig:dtl-taper shows the standard
tapering and $t_0$-shift conditioning on the scenario-1 dataset.

#figure(
  subfigs(cols: 1,
    img("DTL_003_GPR_B-Scans__Background_Baseline_and_Diagonal_TimeLapsed_Mod.png"),
    img("DTL_004_GPR_B-Scans__Background_Subtracted__green_dashed__expected_a.png"),
    img("DTL_005_GPR_B-Scans__With_Synthetic_Laplace_Noise_10_of_signal_std.png"),
  ),
  caption: [Raw B-scans for the diagonal time-lapse study: (a) background and
    time-lapsed models; (b) background-subtracted, with expected arrival times
    marked; (c) with synthetic Laplace-distributed noise at $10%$ of the
    signal standard deviation (used in @ch:hyp2).],
) <fig:dtl-bscans>

#figure(
  subfigs(cols: 1,
    img("DTL_007_Effect_of_Tapering_and_t0_Shift__Scenario_1_dataset_single_t.png"),
    img("DTL_008_B-scan_effect_of_tapering_and_t0_shift__Scenario_1_dataset.png"),
  ),
  caption: [Effect of tapering and the $t_0$ shift on the scenario-1
    ($1 lambda, 1\/2 lambda$) diagonal displacement dataset: (a) a single
    trace; (b) the full B-scan.],
) <fig:dtl-taper>

Per-method migrated images and time-lapse differences for Kirchhoff, Gazdag,
and back-propagation migration are given in @app:hyp1-methods
(@fig:dtl-kirchhoff, @fig:dtl-gazdag, @fig:dtl-backprop). @fig:dtl-summary-amp
below instead overlays the signed difference from all three algorithms
directly.

=== Diagonal Detectability Summary <sec:dtl-detectability>

Because every scenario moves the scatterer along the same $2:1$ diagonal line
through the baseline, the detectability analysis samples each
time-lapse-difference image _along that line_ rather than along a single
Cartesian axis: @fig:dtl-summary-amp overlays the signed, diagonally-sampled
amplitude from all three algorithms, and @fig:dtl-summary-psf plots the
normalised diagonal point-spread function (signed amplitude and Hilbert
envelope) as a function of signed distance along the motion direction.

#figure(
  img("DTL_023_Diagonal_TimeLapse_Migration_Comparison__Signed_Amplitude.png"),
  caption: [Signed time-lapse-difference amplitude for all three migration
    algorithms, sampled along the diagonal motion direction.],
) <fig:dtl-summary-amp>

#figure(
  img("DTL_024_Normalised_Diagonal_PSF__TimeLapse_Difference_Along_Motion_D.png", width: 80%),
  caption: [Normalised diagonal PSF of the time-lapse-difference image along
    the motion direction, swept across the five diagonal scenarios of
    @tab:dtl-scenarios.],
) <fig:dtl-summary-psf>

#draftnote[state the smallest diagonal displacement at which
@fig:dtl-summary-psf still shows a clear, unambiguous peak above background,
and compare it with the lateral and vertical floors of @sec:tl-detectability
and @sec:vtl-detectability --- since the diagonal step combines a lateral and
a vertical component of different magnitude ($#Dx = 2 #Dz$), state whether
the effective detectability floor tracks the (tighter) lateral floor, the
(looser) vertical floor, or some combination of the two.]

== Global Phase-Plane Validation <sec:hyp1-phaseplane>

@sec:hyp1-lateral, @sec:hyp1-vertical, and @sec:hyp1-diagonal showed that
simple amplitude differencing of migrated images detects sub-wavelength
displacement only down to a finite floor, in all three directions. This
section applies the phase-plane shift estimator derived in
@sec:th-fourier-shift, @sec:th-wls, and @sec:th-material-change, and
implemented as described in @sec:meth-phaseplane, to the lateral and vertical
displacement datasets, and shows that it remains accurate at displacement
scales where the amplitude image is already featureless.

Because several of the analyses below are repeated identically across seven
(or six) wavelength scales, only one or two representative scales are shown
in this section; the complete sweep for every analysis is given in
@app:extended-sweeps, so that every one of the relevant figures produced by
`TimeLapse_Processing.ipynb` for these scenarios appears at least once in the
thesis.

=== Migrated Baseline Overview Across Scenarios <sec:tlp-overview>

@fig:tlp-overview shows the starting point for every analysis in this section:
the raw background-subtracted B-scans and the three migrated images
(Kirchhoff, Gazdag, back-propagation) for every scenario used below.

#figure(
  subfigs(cols: 2,
    img("TLP_001_Raw_Background-Subtracted_B-Scans__All_Scenarios.png"),
    img("TLP_002_Kirchhoff_Migration__All_Scenarios.png"),
    img("TLP_003_Gazdag_Migration__All_Scenarios.png"),
    img("TLP_004_Back-prop_Migration__All_Scenarios.png"),
  ),
  caption: [Overview of all scenarios re-used from @sec:hyp1-lateral in this
    section: (a) raw background-subtracted B-scans; (b--d) the same data
    migrated with Kirchhoff, Gazdag, and back-propagation migration.],
) <fig:tlp-overview>

=== Phase-Plane Fit: Lateral Shift Validation <sec:tlp-horizontal-validation>

For every lateral-displacement scenario of @sec:hyp1-lateral, the estimator
of @sec:meth-phaseplane is applied to the baseline/monitor pair, cropped to
$plus.minus 2.5 lambda$ around the scatterer apex (located from the Hilbert
envelope of the baseline image). @fig:tlp-horiz-validation shows the
recovered $(#Dz, #Dx)$ compared with the known ground truth for all three
migration methods.

#figure(
  subfigs(cols: 2,
    img("TLP_005_Phase-plane_shift_estimation__Kirchhoff____Baseline_vs_each.png"),
    img("TLP_006_Phase-plane_shift_estimation__Gazdag____Baseline_vs_each_sce.png"),
    img("TLP_007_Phase-plane_shift_estimation__Back-prop____Baseline_vs_each.png"),
  ),
  caption: [2D phase-plane shift estimation applied to the
    baseline-versus-each-scenario pairs of @sec:hyp1-lateral, for all three
    migration methods: (a) Kirchhoff; (b) Gazdag; (c) back-propagation.],
) <fig:tlp-horiz-validation>

#draftnote[read off @fig:tlp-horiz-validation the smallest lateral displacement
at which the estimated $#Dx$ still tracks the true value to within (state your
accepted tolerance, e.g. $plus.minus 5%$ or $plus.minus 1 "mm"$), for each
migration method, and contrast this explicitly with the amplitude floor found
in @sec:tl-detectability.]

=== Phase-Plane Fit: Vertical Shift Validation <sec:tlp-vertical>

@fig:tlp-vert-validation repeats the previous analysis for the
vertical-displacement dataset of @sec:hyp1-vertical, with the ROI crop and
cross-spectrum fit unchanged except that the displacement being recovered is
now $#Dz$ rather than $#Dx$.

#figure(
  subfigs(cols: 2,
    img("TLP_011_VerticalTimeLapse__Phase-plane_shift_estimation__Kirchhoff.png"),
    img("TLP_012_VerticalTimeLapse__Phase-plane_shift_estimation__Gazdag____B.png"),
    img("TLP_013_VerticalTimeLapse__Phase-plane_shift_estimation__Back-prop.png"),
  ),
  caption: [2D phase-plane shift estimation applied to the
    baseline-versus-each-scenario pairs of @sec:hyp1-vertical, for all three
    migration methods: (a) Kirchhoff; (b) Gazdag; (c) back-propagation.],
) <fig:tlp-vert-validation>

#draftnote[state whether the vertical phase-plane fit remains accurate to the
same smallest displacement found for the lateral case in
@sec:tlp-horizontal-validation, or whether the lateral/vertical asymmetry
predicted in @sec:th-duality (a constant vertical _plateau_ versus a lateral
_slope_) shows up as a difference in achievable precision between
@fig:tlp-horiz-validation and @fig:tlp-vert-validation.]

=== Phase-Plane Fit: Diagonal Shift Validation <sec:tlp-diagonal>

#draftnote[*Gap:* unlike the lateral and vertical cases above, the 2D WLS
phase-plane fit has not yet been run on the diagonal dataset of
@sec:hyp1-diagonal --- `Diagonal_TimeLapse_Playground.ipynb` currently only
produces the amplitude-based detectability result of @sec:dtl-detectability.
Completing Hypothesis 1's diagonal claim requires applying `estimate_shift_2d`
(@sec:meth-phaseplane) to the same baseline/monitor pairs and recovering the
joint $(#Dz, #Dx)$ estimate, which should satisfy the known $#Dx = 2 #Dz$
ratio of @tab:dtl-scenarios if the fit is working correctly --- this is the
single most important piece of evidence still missing from this chapter.]

Hypothesis 1.5, exploring a local, trace-based alternative to the global
phase-plane fit, is developed separately in @sec:hyp1-h15.
