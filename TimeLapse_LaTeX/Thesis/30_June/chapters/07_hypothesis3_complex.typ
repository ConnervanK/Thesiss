#import "../template.typ": *

= Hypothesis 3: Generalisation to Complex Scenes and Real Field Data <ch:hyp3>

@ch:hyp1 and @ch:hyp2 established that the time-lapse phase-plane approach
works on idealised synthetic data with a single known target and a controlled
displacement. This chapter asks whether that generalises: first to synthetic
scenes with multiple independently-moving scatterers under noise, and second
to real borehole GPR field data where neither the target geometry nor the
displacement is controlled.

#para-head[Hypothesis 3.] The time-lapse phase-plane approach generalises
beyond idealised single-scatterer synthetic models to complex scenes with
multiple independently-moving scatterers and to real borehole GPR field data.

== Complex Synthetic Models <sec:hyp3-complex>

#draftnote[Placeholder --- synthetic multi-scatterer data with noise are still
being generated and will be inserted here once available. The planned
experiment contains multiple PEC cylinders displaced simultaneously in
different lateral and vertical directions, tested across Kirchhoff, Gazdag,
and back-propagation (with sign-bit time-reversal for the noisy variants),
with the same Laplace noise model of @sec:hyp3-laplace. The key question is
whether the phase-plane fit --- which assumes a single dominant displacement
within its ROI --- can be applied locally enough to separate multiple
independent events, and whether back-propagation's advantages from @ch:hyp2
carry over to this more complex geometry.]

== Real Field Data: Borehole GPR <sec:hyp3-fielddata>

The most direct test of Hypothesis 3 is whether the phase-plane pipeline,
developed entirely on synthetic gprMax data, produces physically interpretable
results when applied to real borehole GPR measurements. This section presents
the application to a controlled fluid-injection field experiment.

=== Dataset and Operational Stages <sec:hyp3-fd-setup>

The dataset consists of 38 zero-offset borehole GPR profiles acquired on
6 June 2016, with the following parameters: propagation velocity
$v = 0.10 "m/ns"$, centre frequency $f_0 = 0.10 "GHz"$, trace spacing
$d_L = 0.05 "m"$, and maximum depth $85 "m"$. The experiment was conducted in
four operational stages, each with a distinct expected fluid behaviour:

/ Push (profiles 1--4): fluid is actively injected; a large downward displacement and lateral spread are expected.
/ Chase (profiles 5--9): continued injection at a lower rate; smaller incremental advance expected.
/ Wait (profiles 10--20): injection paused; near-zero net displacement expected.
/ Pull (profiles 21--38): fluid is extracted; partial reversal of the push displacement expected.

The raw profiles are pre-processed using an 11-step pipeline: bandpass filter
($0.02$--$0.20 "GHz"$), DC and direct-wave removal, trace alignment (×5
upsampling), SVD rank-1 direct-wave suppression, time-lapse differencing
against the reference profile, linear spherical-gain correction, Tukey
tapering, f-k_z dip filter, excitation low-pass filter, and edge zeroing.

=== Migration of Field B-Scans <sec:hyp3-fd-migration>

Kirchhoff and Gazdag migrations are available for all 37 processed profiles
(profiles 2--38 differenced against profile 1 as the fixed baseline);
back-propagation is available for 14 profiles (profiles 1--5, 7--10, 13, 16,
20, 21, 38), the remainder still pending gprMax forward runs.

@fig:fd-profile-grid compares five representative profiles (1, 3, 8, 20 and
38, spanning the four operational stages) across all four representations:
the processed B-scan, the Kirchhoff-BP migration, the Gazdag migration, and
the back-propagation $E_z$ focus frame. All three migration techniques agree
on a single dominant reflector at approximately $70$--$80 "m"$ depth and
$4$--$8 "m"$ radial distance, which sharpens progressively from profile 1 to
profile 8 and then remains essentially stationary through profiles 20 and 38
-- a first visual indication of the push/chase/wait/pull kinematics discussed
below. The back-propagation column reproduces the same reflector but retains more
residual energy near the borehole (small radial distance) -- a structural
injection-halo artefact of single-sided time-reversal from a borehole source
array, not a pre-processing shortcoming.

#figure(
  cimg("FD_profile_migration_grid.png"),
  caption: [Processed B-scans and their migrated / back-propagated
    counterparts for profiles 1, 3, 8, 20 and 38. Columns: processed B-scan,
    Kirchhoff-BP migration, Gazdag migration, back-propagation $E_z$ focus
    frame (offset 28 snapshots from the nominal focus time).],
) <fig:fd-profile-grid>

=== Phase-Plane Displacement Tracking <sec:hyp3-fd-phaseplane>

The cross-spectrum phase-plane estimator of @sec:meth-phaseplane is applied
to the Gazdag-migrated images, using a fixed ROI of depth $70$--$79 "m"$ and
radial distance $4.5$--$7.5 "m"$. @fig:fd-stage-push to @fig:fd-stage-pull show,
for each of the four operational stages, the raw ingredients that feed this
estimator -- the time-lapse difference image and its monogenic envelope,
with the ROI overlaid -- computed independently for all three migration
techniques (Kirchhoff-BP, Gazdag, back-propagation). The envelope highlights
a single coherent patch inside the ROI in every stage and every technique,
confirming that the ROI is well-placed and that the "single dominant
displacement" assumption behind the phase-plane fit (@sec:meth-phaseplane)
holds throughout the experiment; the patch is largest and best-defined during
Push and shrinks progressively through Chase, Wait and Pull, mirroring the
weakening signal expected as the fluid front decelerates.

#figure(
  cimg("FD_stage_diff_envelope_pushing.png"),
  caption: [Push stage (profiles 1→4): time-lapse difference and monogenic
    envelope, ROI overlaid, for Kirchhoff-BP, Gazdag and back-propagation.],
) <fig:fd-stage-push>

#figure(
  cimg("FD_stage_diff_envelope_chasing.png"),
  caption: [Chase stage (profiles 5→9): time-lapse difference and monogenic
    envelope, ROI overlaid, for Kirchhoff-BP, Gazdag and back-propagation.],
) <fig:fd-stage-chase>

#figure(
  cimg("FD_stage_diff_envelope_waiting.png"),
  caption: [Wait stage (profiles 10→20): time-lapse difference and monogenic
    envelope, ROI overlaid, for Kirchhoff-BP, Gazdag and back-propagation.],
) <fig:fd-stage-wait>

#figure(
  cimg("FD_stage_diff_envelope_pulling.png"),
  caption: [Pull stage (profiles 21→38): time-lapse difference and monogenic
    envelope, ROI overlaid, for Kirchhoff-BP, Gazdag and back-propagation.],
) <fig:fd-stage-pull>

Three tracking strategies are then compared:

/ Strategy 1 (consecutive): each profile pair $(n, n+1)$ is fit independently;
  the cumulative trajectory is the running sum of the incremental estimates.
/ Strategy 2 (fixed global baseline): each profile $n$ is compared directly
  against profile 1; the result is a direct (non-cumulative) displacement
  estimate.
/ Strategy 3 (stage-anchored, recommended): one pair per stage is chosen to
  span the full stage (e.g. profiles 1→4 for Push, 5→9 for Chase), giving the
  highest SNR estimate of the total stage displacement before stitching into a
  global trajectory.

#draftnote[insert the three strategy summary plots here (3-panel line plots of
cumulative/direct $#Dx$, $#Dz$, and $phi_0$ vs profile number from cells 26,
28, 30, 32 of `FieldData_Playground.ipynb`). Export them to a
`TimeLapse_Figures/FieldData/` directory and reference via `img()` if not
already done.]

The Strategy 3 (stage-anchored) results, which give the highest-SNR estimate
of the total displacement per stage, are summarised in @tab:fielddata-stages.

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Stage*], [*Profiles*], [*$#Dz$ (depth)*], [*$#Dx$ (radial)*],
    table.hline(stroke: 0.4pt),
    [Push],  [1→4],   [$+1.41 "m"$], [$-0.26 "m"$],
    [Chase], [5→9],   [$+0.18 "m"$], [$-0.04 "m"$],
    [Wait],  [10→20], [$-0.11 "m"$], [$+0.03 "m"$],
    [Pull],  [21→38], [$-0.26 "m"$], [$+0.05 "m"$],
    table.hline(stroke: 0.4pt),
    [*Net (prof 1→38)*], [], [$+0.94 "m"$], [$-0.18 "m"$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Stage-anchored (Strategy 3) phase-plane displacement estimates
    from the borehole GPR field dataset. Positive $#Dz$ is downward (along
    the borehole); negative $#Dx$ is away from the borehole. Net values are
    the reconstructed global trajectory after stitching the four stage
    estimates; the Pull stage did not fully reverse the Push.],
  kind: table,
) <tab:fielddata-stages>

The dominant signal is a downward displacement of approximately $1.41 "m"$
during the Push stage, accompanied by a lateral spread of $0.26 "m"$ away
from the borehole, consistent with fluid being injected downward and outward.
The Chase stage adds a smaller increment in the same direction, the Wait stage
is near-zero (as expected for paused injection), and the Pull stage partially
reverses the Push but leaves a net residual of $#Dz approx +0.94 "m"$,
$#Dx approx -0.18 "m"$ at profile 38 relative to profile 1.

@fig:fd-disp-kirchhoff to @fig:fd-disp-backprop cross-check this result across
migration techniques: the same WLS cross-spectrum phase-plane fit
(@sec:meth-phaseplane), applied to the *same* four direct stage-boundary
pairs (1→4, 5→9, 10→20, 21→38) and the *same* ROI, is shown for Kirchhoff-BP,
Gazdag and back-propagation side by side, with the full five-panel
diagnostic (cross-spectrum phase, cross-spectrum energy with the WLS
amplitude-threshold contour, the fitted plane, and the 1-D $#kz$ and $#kx$
slices with their fits) for every stage. Unlike @tab:fielddata-stages --
which uses Strategy 3's intra-stage chaining with a stage-local reference --
these figures fit each stage-boundary pair directly, so the two sets of
numbers are not directly comparable; they are a consistency check on
*direction and relative magnitude*, not a replacement for the table.
Kirchhoff-BP and Gazdag agree closely on both the sign and the relative size
of the estimate in every stage (largest during Push, near-zero during Wait),
which is expected since both operate on the same underlying B-scan data.
Back-propagation's estimates are noticeably noisier -- visible in the more
scattered 1-D slices -- and disagree with Kirchhoff/Gazdag on the sign of
$#Dz$; this is addressed as an open question in @sec:hyp3-fd-interpretation.

#figure(
  cimg("FD_displacement_diagnostics_kirchhoff_bp.png"),
  caption: [Kirchhoff-BP: WLS cross-spectrum phase-plane diagnostics for the
    four stage-boundary pairs (rows) -- cross-spectrum phase, cross-spectrum
    energy with WLS threshold contour, fitted plane, 1-D $#kz$ slice, 1-D
    $#kx$ slice (columns).],
) <fig:fd-disp-kirchhoff>

#figure(
  cimg("FD_displacement_diagnostics_gazdag.png"),
  caption: [Gazdag: WLS cross-spectrum phase-plane diagnostics for the four
    stage-boundary pairs (rows), same panel layout as @fig:fd-disp-kirchhoff.],
) <fig:fd-disp-gazdag>

#figure(
  cimg("FD_displacement_diagnostics_backprop.png"),
  caption: [Back-propagation: WLS cross-spectrum phase-plane diagnostics for
    the four stage-boundary pairs (rows), same panel layout as
    @fig:fd-disp-kirchhoff.],
) <fig:fd-disp-backprop>

=== Interpretation <sec:hyp3-fd-interpretation>

#draftnote[Placeholder]

// #draftnote[fill in the physical interpretation once the borehole geometry,
// injection depth, and fluid-injection parameters are confirmed from the field
// survey metadata. Key questions to address: (1) do the inferred $#Dz$ and $#Dx$
// values agree with the known injection depth and the expected lateral spread for
// the given fracture geometry? (2) @fig:fd-disp-backprop shows back-propagation
// disagreeing with Kirchhoff/Gazdag (@fig:fd-disp-kirchhoff, @fig:fd-disp-gazdag)
// on the sign of $#Dz$, despite all three sharing the same ROI and stage pairs
// -- is this a genuine sign-convention difference between the two coordinate
// systems (the back-propagation depth axis runs in the opposite direction to
// the Kirchhoff/Gazdag depth array; @fig:fd-stage-push to @fig:fd-stage-pull show
// the same reflector location in both, so the ROI itself is not the issue), a
// consequence of back-propagation's lower SNR (14 profiles vs. 37, and visibly
// noisier 1-D slices), or evidence that the phase-plane fit is less reliable on
// this technique in the field, consistent with @ch:hyp2? (3) does the field
// result agree with the synthetic fluid-front experiment? State these
// explicitly once the field context is available.]
