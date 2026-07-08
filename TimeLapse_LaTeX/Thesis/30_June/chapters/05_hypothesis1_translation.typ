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
the baseline survey and displaced laterally by each of eight scenarios in the
monitor survey, from $2 lambda$ down to $1 \/ 32 lambda$ (@fig:tl-setup),
using the same gprMax domain and grid as @sec:meth-resolution.

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
    img("TL_006_Effect_of_Tapering_and_t0_Shift__2λ_dataset_single_trace.png"),
    img("TL_007_B-scan_effect_of_tapering_and_t0_shift__2λ_dataset.png"),
  ),
  caption: [Effect of tapering and the $t_0$ shift on the $2 lambda$ lateral
    displacement dataset: (a) a single trace; (b) the full B-scan.],
) <fig:tl-taper>

@fig:tl-kirchhoff, @fig:tl-gazdag, and @fig:tl-backprop show, respectively,
the Kirchhoff-, Gazdag-, and back-propagation-migrated images for all
displacement scenarios, and the time-lapse difference (migrated monitor minus
migrated baseline) that isolates the moved scatterer.

#figure(
  subfigs(cols: 1,
    img("TL_009_Kirchhoff_Migration_zoomed____f_c15_GHz____aperture40.png"),
    img("TL_011_Kirchhoff_Migration__TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Kirchhoff migration of the lateral time-lapse study ($f_c = 1.5 "GHz"$,
    aperture $= 40$), zoomed around the scatterers: (a) migrated image for all
    scenarios; (b) migrated-monitor-minus-migrated-baseline difference.],
) <fig:tl-kirchhoff>

#figure(
  subfigs(cols: 1,
    img("TL_013_Gazdag_Phase-Shift_Migration_zoomed____f_c15_GHz.png"),
    img("TL_015_Gazdag_Migration__TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Gazdag phase-shift migration of the lateral time-lapse study
    ($f_c = 1.5 "GHz"$), zoomed around the scatterers: (a) migrated image;
    (b) time-lapse difference.],
) <fig:tl-gazdag>

#figure(
  subfigs(cols: 1,
    img("TL_017_Back-Propagation_E_zoomed____focus_at_1906_ns.png"),
    img("TL_019_Back-Propagation_Ez_zoomed____focus_at_1906_ns.png"),
    img("TL_021_Back-Propagation__TimeLapse_Differences_Ez_zoomed.png"),
  ),
  caption: [Time-reversal back-propagation migration of the lateral time-lapse
    study, focused at $t = 19.06 "ns"$ and zoomed around the scatterers:
    (a) field magnitude $||bold(E)||$; (b) $E_z$; (c) $E_z$ time-lapse
    difference.],
) <fig:tl-backprop>

=== Amplitude-Based Detectability Summary <sec:tl-detectability>

@fig:tl-summary overlays the signed time-lapse-difference amplitude from all
three algorithms and plots the normalised lateral point-spread function of
that difference at the true scatterer depth, as a function of displacement.

#figure(
  subfigs(cols: 2,
    img("TL_022_TimeLapse_Migration_Comparison__Signed_Amplitude____f_c15_GH.png"),
    img("TL_023_Normalised_Lateral_PSF__TimeLapse_Difference_at_True_Scatter.png"),
  ),
  caption: [(a) Signed time-lapse-difference amplitude for all three migration
    algorithms at $f_c = 1.5 "GHz"$; (b) the normalised lateral PSF of the
    difference image at the true scatterer depth, swept across lateral
    displacements from $2 lambda$ to $1 \/ 32 lambda$.],
) <fig:tl-summary>

#draftnote[state the smallest displacement at which the difference image in
@fig:tl-summary (b) still shows a clear, unambiguous peak above the numerical
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
and is displaced _downward_ from its baseline depth across seven scenarios,
from $1 lambda$ down to $1 \/ 32 lambda$ (@fig:vtl-setup), using the same
domain and grid as @sec:meth-resolution and @sec:hyp1-lateral.

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
    img("VTL_005_Effect_of_Tapering_and_t0_Shift__1λ_dataset_single_trace.png"),
    img("VTL_006_B-scan_effect_of_tapering_and_t0_shift__1λ_dataset.png"),
  ),
  caption: [Effect of tapering and the $t_0$ shift on the $1 lambda$ vertical
    displacement dataset: (a) a single trace; (b) the full B-scan.],
) <fig:vtl-taper>

@fig:vtl-kirchhoff, @fig:vtl-gazdag, and @fig:vtl-backprop repeat the same
three-algorithm migration-and-differencing analysis for the vertical
displacement scenarios.

#figure(
  subfigs(cols: 1,
    img("VTL_008_Kirchhoff_Migration_zoomed____f_c15_GHz____aperture40.png"),
    img("VTL_010_Kirchhoff_Migration__TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Kirchhoff migration of the vertical time-lapse study
    ($f_c = 1.5 "GHz"$, aperture $= 40$), zoomed around the scatterer:
    (a) migrated image for all scenarios; (b) time-lapse difference.],
) <fig:vtl-kirchhoff>

#figure(
  subfigs(cols: 1,
    img("VTL_012_Gazdag_Phase-Shift_Migration_zoomed____f_cf_c_GHz_GHz.png"),
    img("VTL_014_Gazdag_Migration__TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Gazdag phase-shift migration of the vertical time-lapse study
    ($f_c = 1.5 "GHz"$), zoomed around the scatterer: (a) migrated image;
    (b) time-lapse difference.],
) <fig:vtl-gazdag>

#figure(
  subfigs(cols: 1,
    img("VTL_016_Back-Propagation_E_zoomed____focus_at_1906_ns.png"),
    img("VTL_018_Back-Propagation_Ez_zoomed____focus_at_1906_ns.png"),
    img("VTL_020_Back-Propagation__TimeLapse_Differences_Ez_zoomed.png"),
  ),
  caption: [Time-reversal back-propagation migration of the vertical time-lapse
    study, focused at $t = 19.06 "ns"$ and zoomed around the scatterer:
    (a) field magnitude $||bold(E)||$; (b) $E_z$; (c) $E_z$ time-lapse
    difference.],
) <fig:vtl-backprop>

=== Vertical Detectability Summary <sec:vtl-detectability>

@fig:vtl-summary overlays the signed time-lapse-difference amplitude from all
three algorithms and plots the normalised _vertical_ point-spread function of
that difference at $x = 2.0 "m"$, as a function of vertical displacement.

#figure(
  subfigs(cols: 2,
    img("VTL_021_Vertical_TimeLapse_Migration_Comparison__Signed_Amplitude.png"),
    img("VTL_022_Normalised_Vertical_PSF__TimeLapse_Difference_at_x__20_m.png"),
  ),
  caption: [(a) Signed time-lapse-difference amplitude for all three migration
    algorithms; (b) the normalised vertical PSF of the difference image at
    $x = 2.0 "m"$, swept across vertical displacements from $1 lambda$ to
    $1 \/ 32 lambda$.],
) <fig:vtl-summary>

#draftnote[compare the smallest reliably-detected vertical displacement in
@fig:vtl-summary (b) against the lateral result of @fig:tl-summary (b)
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

@fig:dtl-kirchhoff, @fig:dtl-gazdag, and @fig:dtl-backprop repeat the same
three-algorithm migration-and-differencing analysis for the diagonal
displacement scenarios.

#figure(
  subfigs(cols: 1,
    img("DTL_010_Kirchhoff_Migration_zoomed____f_c15_GHz____aperture40.png"),
    img("DTL_012_Kirchhoff_Migration__TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Kirchhoff migration of the diagonal time-lapse study
    ($f_c = 1.5 "GHz"$, aperture $= 40$), zoomed around the scatterer:
    (a) migrated image for all scenarios; (b) time-lapse difference.],
) <fig:dtl-kirchhoff>

#figure(
  subfigs(cols: 1,
    img("DTL_014_Gazdag_Phase-Shift_Migration_zoomed____f_c15_GHz.png"),
    img("DTL_016_Gazdag_Migration__TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Gazdag phase-shift migration of the diagonal time-lapse study
    ($f_c = 1.5 "GHz"$), zoomed around the scatterer: (a) migrated image;
    (b) time-lapse difference.],
) <fig:dtl-gazdag>

#figure(
  subfigs(cols: 1,
    img("DTL_018_Back-Propagation_E_zoomed____focus_at_1906_ns.png"),
    img("DTL_020_Back-Propagation_Ez_zoomed____focus_at_1906_ns.png"),
    img("DTL_022_Back-Propagation__TimeLapse_Differences_Ez_zoomed.png"),
  ),
  caption: [Time-reversal back-propagation migration of the diagonal time-lapse
    study, focused at $t = 19.06 "ns"$ and zoomed around the scatterer:
    (a) field magnitude $||bold(E)||$; (b) $E_z$; (c) $E_z$ time-lapse
    difference.],
) <fig:dtl-backprop>

=== Diagonal Detectability Summary <sec:dtl-detectability>

Because every scenario moves the scatterer along the same $2:1$ diagonal line
through the baseline, the detectability analysis samples each
time-lapse-difference image _along that line_ rather than along a single
Cartesian axis: @fig:dtl-summary (a) overlays the signed, diagonally-sampled
amplitude from all three algorithms, and @fig:dtl-summary (b) plots the
normalised diagonal point-spread function (signed amplitude and Hilbert
envelope) as a function of signed distance along the motion direction.

#figure(
  subfigs(cols: 2,
    img("DTL_023_Diagonal_TimeLapse_Migration_Comparison__Signed_Amplitude.png"),
    img("DTL_024_Normalised_Diagonal_PSF__TimeLapse_Difference_Along_Motion_D.png"),
  ),
  caption: [(a) Signed time-lapse-difference amplitude for all three migration
    algorithms, sampled along the diagonal motion direction; (b) the normalised
    diagonal PSF of the difference image along that same direction, swept across
    the five diagonal scenarios of @tab:dtl-scenarios.],
) <fig:dtl-summary>

#draftnote[state the smallest diagonal displacement at which @fig:dtl-summary
(b) still shows a clear, unambiguous peak above background, and compare it
with the lateral and vertical floors of @sec:tl-detectability and
@sec:vtl-detectability --- since the diagonal step combines a lateral and a
vertical component of different magnitude ($#Dx = 2 #Dz$), state whether
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
  subfigs(cols: 1,
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
    img("TLP_010_VerticalTimeLapse__Phase-plane_shift_estimation__Kirchhoff.png"),
    img("TLP_011_VerticalTimeLapse__Phase-plane_shift_estimation__Gazdag____B.png"),
    img("TLP_012_VerticalTimeLapse__Phase-plane_shift_estimation__Back-prop.png"),
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

== Hypothesis 1.5 (Optional): Local Phase-Gradient Methods <sec:hyp1-h15>

The phase-plane fit of @sec:hyp1-phaseplane is a _global_ operation: it acts
on the whole 2D spectrum at once. This section explores a _local_ alternative
--- the spatial instantaneous phase, extracted trace-by-trace with a
Riesz/Hilbert transform, and its time-frequency relatives --- to test
Hypothesis 1.5: whether the local phase gradients $partial phi \/ partial x$,
$partial phi \/ partial y$ offer an equivalent, simpler route to the same
displacement estimate.

#draftnote[this section is exploratory: the underlying
Riesz-transform/monogenic-signal machinery is implemented in
`CWT_playground.ipynb`, but it has not yet been validated quantitatively
against the global WLS fit of @sec:hyp1-phaseplane on the same datasets ---
treat the results below as a qualitative complement to, not a replacement for,
Hypothesis 1's primary evidence.]

=== Instantaneous Phase Imaging --- Lateral Displacement <sec:tlp-instphase-horizontal>

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

=== Instantaneous Phase Imaging --- Vertical Displacement <sec:tlp-instphase-vertical>

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

=== Spectral-Line Phase Analysis <sec:tlp-spectral-line>

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

=== Cross-Phase Spectrograms <sec:tlp-spectrogram>

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

=== Localised Fourier Shift via Short-Time Fourier Transform <sec:tlp-stft>

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
