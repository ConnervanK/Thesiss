#import "../template.typ": *

= Hypothesis 2: Fluid-Front Tracking and Material-Change Decoupling <ch:hyp2>

@ch:hyp1 used an idealised target: a point-like PEC cylinder, chosen because
its response is easy to interpret and ground-truth. This chapter applies the
same migration and phase-plane machinery to a more field-realistic scenario
--- a spatially _distributed_ fluid front advancing laterally through a
sub-wavelength fracture, following the thin-layer reflectivity model of
@sec:th-material-change. Unlike a point scatterer, a fluid front is not
expected to move as a rigid body, so this chapter is also the first practical
test of the geometry/material-change decoupling derived in
@sec:th-material-change.

#para-head[Hypothesis 2.] Multi-dimensional phase-plane regression allows
tracking the subwavelength translation of a fluid front and inferring the
associated material change in a subwavelength thin fracture.

== Setup: Raw B-Scans and Signal Conditioning

The fluid-front model reuses the domain, grid, and Ricker source of
@ch:methodology; the fracture is represented as a thin, sub-wavelength layer
whose air-filled and water-filled portions are separated by a front that
advances laterally between the baseline and each of seven monitor scenarios,
from $2 lambda$ down to $1 \/ 32 lambda$. @fig:ff-bscans shows the raw and
background-subtracted B-scans, and @fig:ff-taper the effect of the standard
tapering and $t_0$-shift conditioning (@sec:meth-conditioning) on the
$2 lambda$ scenario.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("FF_001_GPR_B-Scans__Background_Baseline_and_TimeLapsed_Models.png"),
    img("FF_002_GPR_B-Scans__Background_Subtracted.png"),
  ),
  caption: [Raw B-scans for the fluid-front study: (a) background and
    time-lapsed models; (b) background-subtracted.],
) <fig:ff-bscans>

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("FF_003_Effect_of_Tapering_and_t0_Shift__2λ_dataset_single_trace.png"),
    img("FF_004_B-scan_effect_of_tapering_and_t0_shift__2λ_dataset.png"),
  ),
  caption: [Effect of tapering and the $t_0$ shift on the $2 lambda$
    fluid-front dataset: (a) a single trace; (b) the full B-scan.],
) <fig:ff-taper>

== Kirchhoff Migration of the Fluid Front

@fig:ff-kirchhoff shows the Kirchhoff-migrated image for all eight fluid-front
scenarios and the resulting time-lapse difference.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("FF_005_Kirchhoff_Migration__All_8_Datasets____f_c15_GHz____aperture.png"),
    img("FF_006_Kirchhoff_Migration_zoomed____f_c15_GHz____aperture40.png"),
    img("FF_007_Kirchhoff_Migration__TimeLapse_Differences_migrated__migrate.png"),
    img("FF_008_Kirchhoff_Migration__TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Kirchhoff migration of the fluid-front study ($f_c = 15 "GHz"$,
    aperture $= 40$): (a, b) migrated image for all scenarios and zoomed;
    (c, d) time-lapse difference and zoomed.],
) <fig:ff-kirchhoff>

== Gazdag Migration of the Fluid Front

@fig:ff-gazdag repeats the analysis with Gazdag phase-shift migration.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("FF_009_Gazdag_Phase-Shift_Migration__All_7_Datasets____f_c15_GHz.png"),
    img("FF_010_Gazdag_Phase-Shift_Migration_zoomed____f_c15_GHz.png"),
    img("FF_011_Gazdag_Migration__TimeLapse_Differences_migrated__migrated_b.png"),
    img("FF_012_Gazdag_Migration__TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Gazdag phase-shift migration of the fluid-front study
    ($f_c = 15 "GHz"$): (a, b) migrated image and zoomed; (c, d) time-lapse
    difference and zoomed.],
) <fig:ff-gazdag>

== Back-Propagation Migration of the Fluid Front

@fig:ff-backprop shows the back-propagation result for both the
field-magnitude and $E_z$ images, focused at $t = 1906 "ns"$, and the $E_z$
time-lapse difference.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("FF_013_Back-Propagation_E__All_7_Datasets____focus_at_1906_ns.png"),
    img("FF_014_Back-Propagation_E_zoomed____focus_at_1906_ns.png"),
    img("FF_015_Back-Propagation_Ez__All_7_Datasets____focus_at_1906_ns.png"),
    img("FF_016_Back-Propagation_Ez_zoomed____focus_at_1906_ns.png"),
    img("FF_017_Back-Propagation__TimeLapse_Differences_Ez_Ez__Ez_baseline.png"),
    img("FF_018_Back-Propagation__TimeLapse_Differences_Ez_zoomed.png"),
  ),
  caption: [Time-reversal back-propagation migration of the fluid-front study,
    focused at $t = 1906 "ns"$: field-magnitude (a, b) and $E_z$ (c, d)
    images, and the $E_z$ time-lapse difference (e, f).],
) <fig:ff-backprop>

== Amplitude-Based Detectability Summary <sec:ff-detectability>

@fig:ff-summary overlays the signed time-lapse-difference amplitude from all
three algorithms and plots the normalised lateral PSF of that difference at
the true front location, exactly as in @sec:tl-detectability for the
point-scatterer case.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("FF_019_TimeLapse_Migration_Comparison__Signed_Amplitude____f_c15_GH.png"),
    img("FF_020_Normalised_Lateral_PSF__TimeLapse_Difference_at_True_Scatter.png"),
  ),
  caption: [(a) Signed time-lapse-difference amplitude for all three migration
    algorithms; (b) normalised lateral PSF of the difference image at the true
    front location.],
) <fig:ff-summary>

== Phase-Plane Fit Applied to the Fluid Front <sec:ff-phaseplane>

The fluid front is not a rigid point scatterer, so the phase-plane fit here
tests something qualitatively different from @ch:hyp1: rather than recovering
a single $(#Dz, #Dx)$ displacement, the goal is to confirm that the fit
correctly attributes the front's advance to the _material-change_ channel
$c = #Dtheta$ of @sec:th-material-change rather than to a spurious geometric
shift, since the fracture walls themselves do not move. @fig:ff-phaseplane
applies the estimator, cropped to $plus.minus 25 lambda$ around the front
position located from the smallest-shift scenario's difference-envelope peak,
to all three migration methods.

#figure(
  grid(columns: (1fr, 1fr, 1fr), gutter: 0.8em,
    img("TLP_013_FluidFlow__Kirchhoff____base_vs_mon_real____crop_25λ__baseli.png"),
    img("TLP_014_FluidFlow__Gazdag____base_vs_mon_real____crop_25λ__baseline.png"),
    img("TLP_015_FluidFlow__Back-prop____base_vs_mon_real____crop_25λ__baseli.png"),
  ),
  caption: [2D phase-plane fit applied to the fluid-front baseline/monitor
    pairs (ROI crop $plus.minus 25 lambda$), for all three migration methods:
    (a) Kirchhoff; (b) Gazdag; (c) back-propagation.],
) <fig:ff-phaseplane>

#draftnote[state, for each migration method in @fig:ff-phaseplane, whether
the fitted slopes $(#Dz, #Dx)$ remain near zero while the intercept $c$
tracks the front advance, as predicted by @eq:intercept-material --- this is
the key result of the fluid-flow application and should be stated explicitly
rather than left implicit in the figure. Recall also the centroid correction
noted in the source notebook: because the newly-flooded region spans
$[x_0, x_0 + Delta x_"true"]$, its centroid sits at $x_0 + Delta x_"true" \/ 2$, so
a raw centroid-shift readout must be doubled, $Delta x_"true" = 2 Delta x_"raw"$, to
recover the true front displacement --- confirm whether this correction is
already baked into the values plotted in @fig:ff-phaseplane or must be
applied separately. Noisy fluid-front data, where it exists, is presented in
@ch:hyp3 rather than here.]
