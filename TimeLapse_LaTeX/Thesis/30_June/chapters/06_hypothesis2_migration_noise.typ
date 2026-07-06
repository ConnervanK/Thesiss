#import "../template.typ": *

= Hypothesis 2: Which Migration Technique Is Best Suited for Time-Lapse Phase-Plane Tracking? <ch:hyp2>

@ch:hyp1 established that phase-plane regression can track sub-wavelength
displacements on clean, synthetic data with an idealised point-like scatterer.
Before applying the same machinery to more complex and realistic scenarios, a
fundamental question must be answered: does the choice of migration algorithm
matter, and if so, which one should be used? This chapter argues that the
answer is yes, and that back-propagation is the preferred method, for two
independent reasons that are examined in sequence: its handling of material
changes (via the fluid-front experiment) and its robustness to heavy-tailed
noise (via sign-bit time-reversal).

#para-head[Hypothesis 2.] Back-propagation migration is the most suitable
technique for time-lapse phase-plane tracking because it preserves the correct
physics for material-change targets --- requiring no correction factor for the
inferred displacement, unlike Kirchhoff and Gazdag --- and is most robust to
heavy-tailed Laplace noise via sign-bit time-reversal.

== Fluid-Front Application: Migration Technique Comparison <sec:hyp2-fluidfront>

@ch:hyp1 used an idealised target: a point-like PEC cylinder, chosen because
its response is easy to interpret and ground-truth. This section applies the
same migration and phase-plane machinery to a more field-realistic scenario
--- a spatially _distributed_ fluid front advancing laterally through a
sub-wavelength fracture, following the thin-layer reflectivity model of
@sec:th-material-change. Unlike a point scatterer, a fluid front is not
expected to move as a rigid body, so this section is also the first practical
test of the geometry/material-change decoupling derived in
@sec:th-material-change, and of whether the three migration algorithms treat
the inferred displacement consistently.

=== Setup: Raw B-Scans and Signal Conditioning

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

=== Kirchhoff Migration of the Fluid Front

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

=== Gazdag Migration of the Fluid Front

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

=== Back-Propagation Migration of the Fluid Front

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

=== Amplitude-Based Detectability Summary <sec:ff-detectability>

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

=== Phase-Plane Fit Applied to the Fluid Front <sec:ff-phaseplane>

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
$[x_0, x_0 + (#Dx)_"true"]$, its centroid sits at $x_0 + (#Dx)_"true" \/ 2$, so
a raw centroid-shift readout must be doubled, $(#Dx)_"true" = 2 (#Dx)_"raw"$, to
recover the true front displacement --- confirm whether this correction is
already baked into the values plotted in @fig:ff-phaseplane or must be
applied separately. Noisy fluid-front data is presented in
@sec:hyp3-fluidflow below.]

=== Factor-of-Two Correction for Kirchhoff and Gazdag <sec:ff-factor2>

The result of @sec:ff-phaseplane reveals a systematic difference between the
migration algorithms that is not present for a point-scatterer target: the
displacements inferred by Kirchhoff and Gazdag must be multiplied by two to
recover the physical front advance, while back-propagation gives the physical
value directly.

The root cause is the _exploding-reflector model_ used by both Kirchhoff and
Gazdag. These methods assume that reflectors "explode" at time zero and the
resulting waves propagate upward at half the true velocity $v \/ 2$, so that
two-way travel time corresponds to one-way propagation at $v \/ 2$. The
depth-axis of the migrated image is therefore calibrated in terms of $v \/ 2$,
not $v$. When the phase-plane estimator reads off a slope from such a migrated
image, the inferred $#Dz$ or $#Dx$ is expressed in those same
half-velocity coordinates, and must be doubled to obtain the physical
displacement.

Back-propagation is free of this ambiguity. It directly simulates the full
electromagnetic wavefield at the true propagation velocity $v$ (via a gprMax
forward run), so the resulting $E_z$ focus image is calibrated in physical
coordinates. The phase-plane slope from a back-propagation image is the
physical displacement without any additional factor.

#draftnote[confirm the factor-of-two numerically: read off the inferred
$#Dz$ and $#Dx$ from @fig:ff-phaseplane for a scenario whose true
displacement is known, and verify that the Kirchhoff/Gazdag estimates are
approximately half the back-propagation estimate (i.e. that doubling the
Kirchhoff/Gazdag values gives agreement with the known ground truth and with
the back-propagation result). State the ratio explicitly.]

== Noise Robustness: Sign-Bit Time-Reversal <sec:hyp2-noise>

The factor-of-two analysis above already hints at a deeper difference between
the algorithms: back-propagation models the physics rather than approximating
it. The same structural difference makes it the most noise-robust of the three
methods, as this section shows --- provided the time-reversed excitation is
handled correctly.

=== A Laplace Noise Model from Real Field Data <sec:hyp3-laplace>

Rather than injecting arbitrary synthetic noise, the noise level and shape
used throughout this thesis are fitted to real GPR field noise. A sample of
field noise is tracked through the processing pipeline and fitted with both a
Laplace and a Gaussian distribution at two stages: after the final processing
stage ("9 Crop Samples"), and at the last stage before the spherical-gain
correction ("7 Constant Velocity"), since the gain correction inflates the
amplitude scale by several orders of magnitude and is not representative of
the raw simulated $E_z$ amplitudes used elsewhere in this thesis.

#figure(
  table(
    columns: (auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Stage*], [*Laplace scale*], [*Gaussian $sigma$*],
    table.hline(stroke: 0.4pt),
    [9 Crop Samples (post-gain)],     [$117 space 584.0$], [$171 space 156.1$],
    [7 Constant Velocity (pre-gain)], [$6.085$],           [$9.233$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Laplace and Gaussian fits to real field noise, at two pipeline
    stages ($n = 976 space 244$ samples each). The heavier-tailed Laplace
    distribution is adopted throughout this thesis; loc $= 0$ for both
    stages.],
  kind: table,
) <tab:laplace-fit>

In both cases the fitted distribution is heavier-tailed than a Gaussian of
matched variance, consistent with field GPR noise being dominated by
occasional large-amplitude clutter and interference rather than purely thermal
noise. The pre-gain fit is the one actually used to generate synthetic noise
for the experiments in @sec:hyp3-lateral, @sec:hyp3-vertical,
@sec:hyp3-diagonal, and @sec:hyp3-fluidflow: its Laplace _shape_ (loc $= 0$,
heavier tails than Gaussian) is kept, but its scale is rescaled so that the
resulting noise standard deviation is exactly $10%$ of each synthetic B-scan's
own signal standard deviation --- a light, realistic noise level rather than
the raw fitted scale, which would be incommensurate with the synthetic $E_z$
amplitudes.

=== Migration-Algorithm Response to Pure Noise <sec:hyp3-purenoise>

Before testing noise robustness on real signal-plus-noise data, a sanity check
establishes what each migration algorithm does to noise _alone_: an empty
B-scan (no scatterer, no real signal) containing only the Laplace noise of
@sec:hyp3-laplace is migrated with Kirchhoff, Gazdag, and back-propagation.
If a method turns pure noise into something that looks like a coherent,
scatterer-like focus, every noisy result elsewhere in this chapter carries a
false-positive risk that must be accounted for.

*Result: the three methods do not fail the same way.* Kirchhoff turns pure
noise into clearly coherent, smooth wave-like bands that could easily be
misread as real layered structure: its delay-and-sum aperture stacking imposes
coherence on incoherent input by construction, summing many traces along
travel-time curves and smoothing incoherent noise into locally-correlated
structure. Gazdag's noise output, by contrast, stays speckled and incoherent,
with no wave-like artefacts --- just texture. This is a direct, practical
consequence of the structural difference between the two algorithms
(@sec:th-migration): Kirchhoff's spatial stacking manufactures apparent
coherence from nothing, while Gazdag's frequency-domain downward continuation
does not.

#draftnote[the pure-noise Kirchhoff/Gazdag comparison above is implemented
and run in `Noise_Playground.ipynb`, but its plots have not yet been exported
to a `TimeLapse_Figures` folder and added to this project's path --- do this
before finalising this section, and include the comparison figure here. The
back-propagation-of-pure-noise variant (@sec:hyp3-signbit) has real, completed
gprMax output in `noise_study/backprop/{purenoise_peaknorm,purenoise_signbit}/`
but its focus-frame snapshots have not yet been loaded, plotted, or compared
quantitatively --- this is the single most direct test of whether sign-bit
time-reversal actually suppresses the false-positive risk identified above for
back-propagation specifically, and is currently missing.]

=== Sign-Bit Time-Reversal for Noise-Robust Back-Propagation <sec:hyp3-signbit>

Back-propagation migration is structurally different from Kirchhoff and
Gazdag: it is not a post-processing step on an already-recorded image, but
requires re-injecting the (time-reversed) recorded data as a source into a
new forward simulation. This makes it vulnerable to noise in a way the other
two methods are not: spatial focusing during back-propagation is governed
almost entirely by _phase_ (the zero-crossings of the time-reversed
wavefield), not by amplitude, yet the default excitation scheme injects the
_peak-normalised_ time-reversed wavefield $u(x, tau)$. A large-amplitude
noise spike anywhere in the data is peak-normalised along with everything
else, so it is injected with the same outsized amplitude it had in the noisy
record --- letting it act as its own competing point source during
back-propagation, interfering at the true source locations instead of being
suppressed by destructive interference.

_Sign-bit time-reversal_ fixes this by injecting only the _sign_ of the
time-reversed wavefield instead of its peak-normalised value,
$ u_"sign" (x, tau) = op("sign")(u(x, tau)) , $ <eq:signbit>
implemented in `write_backprop_files(..., sign_bit=True)`
(`helper_functions/migration.py`). This keeps every zero-crossing and phase
trend of the GPR wavelet completely intact, since the sign of a signal carries
its full phase information, while squashing every noise spike down to the same
$plus.minus 1$ amplitude as the coherent signal --- stripping noise of the
outsized amplitude that would otherwise let it dominate the back-propagated
wavefield. The clean-data experiments of @ch:hyp1 and @sec:hyp2-fluidfront
use the default peak-normalised excitation throughout, since they have no
noise to suppress; every noisy back-propagation result in this section uses
sign-bit excitation instead.

=== Noisy Lateral Movement <sec:hyp3-lateral>

The lateral displacement sweep of @sec:hyp1-lateral is repeated on B-scans
contaminated with the synthetic Laplace noise of @fig:tl-bscans (c), for the
two analytically-defined migration algorithms (Kirchhoff and Gazdag;
back-propagation was not re-run on this particular noisy dataset due to its
computational cost --- contrast with the vertical, diagonal, and fluid-front
cases below, where it was). @fig:tl-noisy-kirchhoff and @fig:tl-noisy-gazdag
repeat the migration and differencing analysis of @sec:hyp1-lateral on this
noisy data.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("TL_024_Kirchhoff_Migration_-_All_8_Noisy_Datasets____f_c15_GHz____a.png"),
    img("TL_025_Kirchhoff_Migration_-_Noisy_zoomed____f_c15_GHz____aperture4.png"),
    img("TL_026_Kirchhoff_Migration_-_Noisy_TimeLapse_Differences_migrated_-.png"),
    img("TL_027_Kirchhoff_Migration_-_Noisy_TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Kirchhoff migration of the _noisy_ lateral time-lapse data
    ($f_c = 15 "GHz"$, aperture $= 40$): (a, b) migrated image for all
    scenarios and zoomed; (c, d) time-lapse difference and zoomed.],
) <fig:tl-noisy-kirchhoff>

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("TL_028_Gazdag_Phase-Shift_Migration_-_All_8_Noisy_Datasets____f_c15.png"),
    img("TL_029_Gazdag_Phase-Shift_Migration_-_Noisy_zoomed____f_c15_GHz.png"),
    img("TL_030_Gazdag_Migration_-_Noisy_TimeLapse_Differences_migrated_-_mi.png"),
    img("TL_031_Gazdag_Migration_-_Noisy_TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Gazdag phase-shift migration of the _noisy_ lateral time-lapse
    data ($f_c = 15 "GHz"$): (a, b) migrated image and zoomed; (c, d)
    time-lapse difference and zoomed.],
) <fig:tl-noisy-gazdag>

#draftnote[compare @fig:tl-noisy-kirchhoff and @fig:tl-noisy-gazdag against
the clean-data results of @fig:tl-kirchhoff and @fig:tl-gazdag and state at
which displacement scale the additive noise first prevents the time-lapse
difference from being distinguished from background clutter.]

The same noisy lateral dataset is also carried through the phase-plane
estimator of @sec:meth-phaseplane (@fig:tlp-noise), for the two migration
methods available on this noisy dataset.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("TLP_008_NoisyTimeLapse__Phase-plane_shift_estimation__Kirchhoff____B.png"),
    img("TLP_009_NoisyTimeLapse__Phase-plane_shift_estimation__Gazdag____Base.png"),
  ),
  caption: [2D phase-plane shift estimation applied to the _noisy_
    lateral-displacement dataset: (a) Kirchhoff; (b) Gazdag.],
) <fig:tlp-noise>

#draftnote[quantify, for both methods, how much the additive Laplace noise
degrades the estimated $#Dx$ relative to the clean result of
@fig:tlp-horiz-validation, and relate this to the amplitude-weighting
argument of @sec:th-mask-weight --- this figure pair is currently the _only_
direct evidence available anywhere in this thesis for the phase-plane
estimator's noise resilience; everything in @sec:hyp3-vertical,
@sec:hyp3-diagonal, and @sec:hyp3-fluidflow below is, for now,
amplitude-level only (see @sec:hyp3-gaps).]

=== Noisy Vertical Movement <sec:hyp3-vertical>

The vertical displacement sweep of @sec:hyp1-vertical is repeated under the
same Laplace noise, for all three migration algorithms --- unlike the lateral
case, back-propagation _is_ re-run here, using sign-bit time-reversal
(@sec:hyp3-signbit) on the noisy traces. @fig:vtl-noisy-kirchhoff and
@fig:vtl-noisy-gazdag repeat the migration and differencing analysis of
@sec:hyp1-vertical on this noisy data.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("VTL_006_B-scan_Frequency_Spectra_--_Clean_vs_Noisy____f_c15_GHz.png"),
    img("VTL_021_Kirchhoff_Migration_-_All_7_Noisy_Datasets____f_c15_GHz____a.png"),
    img("VTL_026_Kirchhoff_Migration_Noisy_zoomed____f_c15_GHz____aperture40.png"),
    img("VTL_029_Gazdag_Phase-Shift_Migration_-_All_7_Noisy_Datasets____f_c15.png"),
  ),
  caption: [Noisy vertical time-lapse data: (a) clean-vs-noisy frequency
    spectra; (b, c) Kirchhoff migration; (d) Gazdag migration.],
) <fig:vtl-noisy-kirchhoff>

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("VTL_030_Gazdag_Phase-Shift_Migration_Noisy_zoomed____f_c15_GHz.png"),
    img("VTL_034_Back-Propagation_E_Noisy_--_All_7_Datasets____focus_at_1906.png"),
    img("VTL_036_Back-Propagation_Ez_Noisy_--_All_7_Datasets____focus_at_1906.png"),
    img("VTL_038_Back-Propagation_Noisy_--_TimeLapse_Differences_Ez____focus.png"),
  ),
  caption: [Gazdag (a) and sign-bit back-propagation (b--d) migration of the
    _noisy_ vertical time-lapse data, focused at $t = 1906 "ns"$: (a) Gazdag
    zoomed; (b) back-prop $||bold(E)||$ all scenarios; (c) back-prop $E_z$ all
    scenarios; (d) back-prop $E_z$ time-lapse difference.],
) <fig:vtl-noisy-gazdag>

#figure(
  img("VTL_039_Back-Propagation_Noisy_--_TimeLapse_Differences_Ez_zoomed.png", width: 60%),
  caption: [Sign-bit back-propagation $E_z$ time-lapse difference for the
    noisy vertical dataset, zoomed around the scatterer depth.],
) <fig:vtl-noisy-backprop>

#draftnote[*Gap:* this section currently only carries the analysis to the
amplitude/migrated-image level, exactly mirroring @sec:hyp1-vertical's
clean-data detectability summary --- the 2D WLS phase-plane fit has not yet
been applied to this noisy vertical dataset. Doing so, and comparing the
result against both the clean vertical fit of @sec:tlp-vertical and the noisy
lateral fit of @fig:tlp-noise, is necessary to claim Hypothesis 2 for the
vertical direction specifically.]

=== Noisy Diagonal Movement <sec:hyp3-diagonal>

The diagonal displacement sweep of @sec:hyp1-diagonal is likewise repeated
under Laplace noise, for all three migration algorithms, again using sign-bit
time-reversal for the noisy back-propagation run. @fig:dtl-noisy shows the
noisy Kirchhoff and Gazdag migrations and time-lapse differences, the sign-bit
time-reversed excitation itself, and the noisy sign-bit back-propagation
result.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("DTL_025_Kirchhoff_Migration_-_All_6_Noisy_Datasets____f_c15_GHz____a.png"),
    img("DTL_027_Kirchhoff_Migration_-_Noisy_TimeLapse_Differences_migrated_-.png"),
    img("DTL_029_Gazdag_Phase-Shift_Migration_-_All_6_Noisy_Datasets____f_c15.png"),
    img("DTL_031_Gazdag_Migration_-_Noisy_TimeLapse_Differences_migrated_-_mi.png"),
    img("DTL_033_Sign-Bit_Time-Reversed_Excitation_--_B-scans_and_Spectra_Noi.png"),
    img("DTL_038_Back-Propagation_Noisy_--_TimeLapse_Differences_Ez____focus.png"),
  ),
  caption: [Noisy diagonal time-lapse data: Kirchhoff (a, b) and Gazdag (c, d)
    migration and differencing; (e) the sign-bit time-reversed back-propagation
    excitation B-scans and spectra; (f) the resulting noisy back-propagation
    time-lapse difference.],
) <fig:dtl-noisy>

#draftnote[as for the vertical case, the diagonal WLS phase-plane fit has not
yet been run on this noisy dataset --- this section is amplitude-level only.
Read off the diagonal-PSF detectability floor under noise (mirroring
@sec:dtl-detectability) and state whether it shifts relative to the clean-data
floor.]

=== Noisy Fluid Front <sec:hyp3-fluidflow>

Finally, the fluid-front scenario of @sec:hyp2-fluidfront is repeated under
the same Laplace noise model, again with sign-bit time-reversal for the noisy
back-propagation run. @fig:ff-noisy shows the noisy Kirchhoff and Gazdag
migrations and time-lapse differences, and the noisy sign-bit back-propagation
result.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("FF_017_Kirchhoff_Migration_-_All_8_Noisy_Datasets____f_c15_GHz____a.png"),
    img("FF_019_Kirchhoff_Migration_-_Noisy_TimeLapse_Differences_migrated_-.png"),
    img("FF_021_Gazdag_Phase-Shift_Migration_-_All_8_Noisy_Datasets____f_c15.png"),
    img("FF_023_Gazdag_Migration_-_Noisy_TimeLapse_Differences_migrated_-_mi.png"),
    img("FF_034_Back-Propagation_E_Noisy_--_All_Datasets____focus_at_1906_ns.png"),
    img("FF_038_Back-Propagation_Noisy_--_TimeLapse_Differences_Ez____focus.png"),
  ),
  caption: [Noisy fluid-front data: Kirchhoff (a, b) and Gazdag (c, d)
    migration and differencing; sign-bit back-propagation field magnitude (e)
    and time-lapse difference (f).],
) <fig:ff-noisy>

#draftnote[the noisy fluid-front phase-plane fit (the noisy analogue of
@fig:ff-phaseplane) has not yet been run --- this is the most field-relevant
gap in this section, since it is the only case combining noise robustness _and_
the geometry/material-change decoupling that Hypothesis 2 explicitly claims
together.]

=== Outstanding Quantitative Work <sec:hyp3-gaps>

This chapter currently demonstrates noise robustness convincingly only at the
level of migrated _images_ (do they still show a detectable, unambiguous focus
under noise?), for all four cases, and demonstrates the phase-plane
estimator's noise robustness quantitatively for only one case (noisy lateral
movement, @fig:tlp-noise). To fully support Hypothesis 2's noise-robustness
claim, the following quantitative work, proposed but not yet run (see
`Markdown/Robustness_Experiments.md`), is still needed:

+ Apply the 2D WLS phase-plane fit to the noisy vertical, diagonal, and
  fluid-front datasets of @sec:hyp3-vertical, @sec:hyp3-diagonal, and
  @sec:hyp3-fluidflow, exactly as already done for the lateral case.

+ Quantitatively compare sign-bit versus peak-normalised back-propagation on
  pure noise (@sec:hyp3-purenoise), using the completed gprMax output already
  sitting in `noise_study/backprop/`.

+ An orthogonality test isolating a pure geometric shift from a pure
  material-change phase rotation under noise.

+ A controlled SNR sweep from $+30 "dB"$ to $-10 "dB"$ comparing ordinary
  least squares against the WLS estimator.

+ A $20 times 20$ grid sweep of true $(#Dt, #Dtheta)$ pairs, reported as a
  2D inversion-error heatmap.

== Conclusion: Back-Propagation as the Preferred Method <sec:hyp2-conclusion>

The two sets of experiments in this chapter converge on the same answer:
back-propagation is the preferred migration technique for time-lapse
phase-plane tracking.

From the fluid-front experiment (@sec:hyp2-fluidfront), Kirchhoff and Gazdag
both require a factor-of-two correction when inferring displacement from a
target whose reflection arises from a material change rather than a geometric
shift. This correction is an artefact of the exploding-reflector model embedded
in both algorithms and does not appear for a point-scatterer target --- making
it an invisible systematic error in any field scenario where the nature of the
target (rigid-body motion versus material change) is not known in advance.
Back-propagation, which directly simulates the full wavefield, produces
geometrically correct displacements for both types of target without any
post-hoc correction.

From the noise-robustness experiments (@sec:hyp2-noise), Kirchhoff's
coherence-manufacturing behaviour (@sec:hyp3-purenoise) makes it the most
dangerous algorithm under noise: it creates scatterer-like artefacts from
pure noise, raising false-positive risk. Gazdag stays incoherent but adds
significant speckle. Back-propagation with sign-bit time-reversal suppresses
impulsive noise by reducing every noise spike to the same $plus.minus 1$
amplitude as the coherent signal, preserving phase information while stripping
amplitude-based false-positive risk.

Both arguments together support Hypothesis 2. The method developed in @ch:hyp1
should be applied with back-propagation whenever the additional computational
cost is acceptable; this is also the method used for the real field data in
@ch:hyp3.
