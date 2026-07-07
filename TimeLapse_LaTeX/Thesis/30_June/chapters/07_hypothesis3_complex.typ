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
back-propagation is available for 11 profiles (profiles 1--5, 7--10, 13, 16),
the remainder still pending gprMax forward runs.

#draftnote[insert the Kirchhoff and Gazdag migrated images for representative
profiles here (e.g. one per stage: profs 3, 7, 15, 30), using the saved
`.npy` files in `TimeLapse_Notebooks/fielddata/output/migrated/` and the
corresponding difference PNGs in `.../difference_from_ref/`. Add back-prop
Ez focus-frame images from `.../backprop_snapshots/` for the available 11
profiles. Check whether the images are already in a `TimeLapse_Figures/`
subdirectory accessible via `img()`; if not, export them there first.]

=== Phase-Plane Displacement Tracking <sec:hyp3-fd-phaseplane>

The cross-spectrum phase-plane estimator of @sec:meth-phaseplane is applied
to the Gazdag-migrated images, using a fixed ROI of depth $70$--$79 "m"$ and
radial distance $4.5$--$7.5 "m"$. Three tracking strategies are compared:

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

=== Interpretation <sec:hyp3-fd-interpretation>

#draftnote[fill in the physical interpretation once the borehole geometry,
injection depth, and fluid-injection parameters are confirmed from the field
survey metadata. Key questions to address: (1) do the inferred $#Dz$ and $#Dx$
values agree with the known injection depth and the expected lateral spread for
the given fracture geometry? (2) does the back-prop phase-plane result
(available for 11 profiles) agree with the Gazdag result, as expected from the
noise-robustness comparison of @ch:hyp2, and which estimate is more reliable
given the field noise level? (3) does the synthetic fluid-front
experiment? State these explicitly once the field context is available.]
