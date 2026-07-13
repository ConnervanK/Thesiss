#import "../template.typ": *

= Hypothesis 3: Generalisation to Complex Scenes and Real Field Data <ch:hyp3>

@ch:hyp1 and @ch:hyp2 established that the time-lapse phase-plane approach
works on idealised synthetic data with a single known target and a controlled
displacement. This chapter asks whether that generalises to real borehole GPR
field data, where neither the target geometry nor the displacement is
controlled, and --- because a real reflector never comes with a
pre-labelled region of interest --- what happens to the estimate as the
method used to *choose* that region is varied.

#para-head[Hypothesis 3.] The time-lapse phase-plane approach generalises
beyond idealised single-scatterer synthetic models to real borehole GPR field
data.

Generalisation to complex synthetic scenes with multiple independently-moving
scatterers under noise (several PEC cylinders displaced simultaneously in
different lateral and vertical directions, tested across Kirchhoff, Gazdag
and back-propagation with the Laplace noise model of @sec:hyp3-laplace) is
left as future work and is not covered further in this chapter. <sec:hyp3-complex>

== Field Data Explanation <sec:hyp3-fielddata>

The dataset consists of 38 zero-offset borehole GPR profiles acquired on
6 June 2016, with the following parameters: propagation velocity
$v = 0.10 "m/ns"$, centre frequency $f_0 = 0.10 "GHz"$, trace spacing
$d_L = 0.05 "m"$, and maximum depth $85 "m"$. The experiment was conducted in
four operational stages, each with a distinct expected fluid behaviour:

/ Push (profiles 1--4): fluid is actively injected; a large downward displacement and lateral spread are expected.
/ Chase (profiles 5--9): continued injection at a lower rate; smaller incremental advance expected.
/ Wait (profiles 10--20): injection paused; near-zero net displacement expected.
/ Pull (profiles 21--38): fluid is extracted; partial reversal of the push displacement expected.

Each profile is a raw depth/two-way-travel-time B-scan along the borehole,
with no migration, differencing or gain correction applied yet --- that is
the subject of @sec:hyp3-fd-processing. The raw B-scans for the profiles used
throughout this chapter are shown, alongside their processed and migrated
counterparts, in @fig:fd-profile-grid.

== Processing <sec:hyp3-fd-processing>

The raw profiles are pre-processed using an 11-step pipeline: bandpass filter
($0.02$--$0.20 "GHz"$), DC and direct-wave removal, trace alignment (×5
upsampling), SVD rank-1 direct-wave suppression, time-lapse differencing
against the reference profile, linear spherical-gain correction, Tukey
tapering, f-k_z dip filter, excitation low-pass filter, and edge zeroing.
Kirchhoff and Gazdag migrations are available for all 37 processed profiles
(profiles 2--38 differenced against profile 1 as the fixed baseline);
back-propagation is available for 14 profiles (profiles 1--5, 7--10, 13, 16,
20, 21, 38), the remainder still pending gprMax forward runs.

Back-propagation additionally requires suppressing an *injection halo*: every
gprMax source in the time-reversal simulation is injected at the borehole
wall ($y approx 0 "m"$), and the destructive interference among sources that
should cancel the field away from the true reflector is never fully complete
near $y = 0$, leaving residual energy at small radial distance regardless of
how the input is pre-processed. This is a structural property of
single-sided time-reversal from a borehole source array, not a
pre-processing shortcoming, and is suppressed for display by zeroing the
first $2.0 "m"$ of the radial axis. One candidate fix was tested and
reverted: normalising each receiver trace by its RMS value amplifies
low-SNR traces far from the fluid front, and these noisy traces then
back-propagate *coherently* toward the borehole axis, making the halo worse
rather than better.

@fig:fd-profile-grid compares five representative profiles (1, 3, 8, 20 and
38, spanning the four operational stages) across all four representations:
the processed B-scan, the Kirchhoff-BP migration, the Gazdag migration, and
the back-propagation $E_z$ focus frame. All three migration techniques agree
on a single dominant reflector at approximately $70$--$80 "m"$ depth and
$4$--$8 "m"$ radial distance, which sharpens progressively from profile 1 to
profile 8 and then remains essentially stationary through profiles 20 and 38
-- a first visual indication of the push/chase/wait/pull kinematics discussed
below. The back-propagation column reproduces the same reflector but retains
more residual energy near the borehole -- the injection halo described above.

#figure(
  cimg("FD_profile_migration_grid.png"),
  caption: [Processed B-scans and their migrated / back-propagated
    counterparts for profiles 1, 3, 8, 20 and 38. Columns: processed B-scan,
    Kirchhoff-BP migration, Gazdag migration, back-propagation $E_z$ focus
    frame (offset 28 snapshots from the nominal focus time).],
) <fig:fd-profile-grid>

== Region of Influence Workflow <sec:hyp3-fd-roi>

The phase-plane estimator of @sec:meth-phaseplane fits a plane to the
cross-spectrum phase of two migrated images, but only inside a region of
influence (ROI): a subset of $(k_z, k_x)$ cells (or, equivalently, a spatial
window before the FFT) chosen to contain the reflector's coherent energy and
exclude noise. Every result in this chapter depends on that choice, so three
increasingly automatic ways of making it are compared here, all applied to
the same Gazdag-migrated images and the same four representative pairs
(1→3, 3→8, 8→20, 20→38 -- approximating the Push, Chase, Wait and Pull
stages respectively).

=== Rectangular Window <sec:hyp3-fd-roi-rect>

The baseline approach is a hand-picked axis-aligned box in depth and radial
distance, tuned per pair against the monogenic envelope of the time-lapse
difference image (@sec:meth-phaseplane), with a fixed cross-spectrum band and
amplitude threshold selecting the $(k_z, k_x)$ cells inside it. Using a fixed
ROI of depth $70$--$79 "m"$ and radial distance $4.5$--$7.5 "m"$,
@fig:fd-stage-push to @fig:fd-stage-pull show, for each stage, the time-lapse
difference image and its monogenic envelope with the ROI overlaid, computed
independently for all three migration techniques (Kirchhoff-BP, Gazdag,
back-propagation). The envelope highlights a single coherent patch inside
the ROI in every stage and every technique, confirming that the ROI is
well-placed and that the "single dominant displacement" assumption behind
the phase-plane fit holds throughout the experiment; the patch is largest
and best-defined during Push and shrinks progressively through Chase, Wait
and Pull, mirroring the weakening signal expected as the fluid front
decelerates.

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

The rectangular window is simple and, as @sec:hyp3-fd-phaseplane shows,
produces physically sensible results, but it is hand-tuned per pair and says
nothing about *where else* in the domain the same fit would or would not be
reliable. The next two subsections replace the fixed box with, respectively,
a systematic scan and a manually-drawn selection.

=== Sliding Window <sec:hyp3-fd-roi-sliding>

Rather than picking one ROI by hand, a fixed-size window can be slid across
the entire difference B-scan, with the phase-plane fit repeated independently
in every window position to build up a 2-D map of $#Dz$ and $#Dx$. The window
size is set as a multiple of the dominant wavelength $lambda = v \/ f_0 =
1.0 "m"$, with independent multiples for the depth and radial extent: the
sub-horizontal reflector geometry means coherent structure in depth extends
over many wavelengths ($k_z approx 0$), while the two energy lobes that carry
the displacement information sit at $k_x approx plus.minus 8 "rad/m"$ (a
radial wavelength of only $approx 0.8 "m"$), so a window narrower than about
one radial wavelength cannot resolve a stable phase gradient in $#Dx$ and
simply fits noise. The window advances by a fixed number of pixels each
step -- large enough to keep the scan fast, with an option to shrink it once
the window size itself has been validated.

A second gate is needed beyond the ROI-window's own amplitude threshold:
because that threshold is relative to *each window's own* spectral peak, a
window sitting entirely in noise can still return a "confident"-looking but
meaningless fit. An additional absolute-energy gate -- requiring a window's
mean $|Delta"amplitude"|$ to reach a fraction of the pair's global
98^"th"-percentile amplitude -- discards these before they reach the map.

@fig:fd-sliding-qc shows the resulting window grid for the Push-stage pair
(window $3.5 "m"$ deep × $0.96 "m"$ radial, $4424$ overlapping positions),
overlaid on the difference B-scan with one window highlighted at full size
for scale. @fig:fd-sliding-maps shows the resulting $#Dz$ / $#Dx$ maps for
all four pairs. In every pair, the surviving (non-blanked) region collapses
onto a single spatially-coherent patch at approximately $68$--$83 "m"$ depth
and $4$--$8 "m"$ radial distance -- the same location as the hand-picked ROI
of @sec:hyp3-fd-roi-rect, recovered here with no prior knowledge of where the
reflector was. The Push and Chase maps (1→3, 3→8) show the largest-magnitude,
best-defined patches, matching the strongest expected displacement; the Wait
map (8→20) is both smaller in magnitude (colour range roughly a third of
Push/Chase) and less coherent, consistent with near-zero net movement during
the paused-injection stage. This is an independent cross-check of both the
reflector's *location* and the qualitative *stage-to-stage trend* established
by the rectangular-window results, obtained by a method that never had the
ROI told to it.

#figure(
  cimg("FD_sliding_window_qc.png"),
  caption: [Sliding-window scan grid for the Push-stage pair (1→3), Gazdag
    migration: every scanned window overlaid (thin yellow) on the time-lapse
    difference B-scan, with one window highlighted (green) at full size.],
) <fig:fd-sliding-qc>

#figure(
  subfigs(cols: 1,
    cimg("FD_sliding_window_maps_1_to_3.png"),
    cimg("FD_sliding_window_maps_3_to_8.png"),
    cimg("FD_sliding_window_maps_8_to_20.png"),
    cimg("FD_sliding_window_maps_20_to_38.png"),
  ),
  caption: [Sliding-window $#Dz$ / $#Dx$ maps, Gazdag migration, for the four
    representative pairs: (a) 1→3 (Push), (b) 3→8 (Chase), (c) 8→20 (Wait),
    (d) 20→38 (Pull). Windows failing either the WLS mask-count or the
    absolute-energy gate are left blank.],
) <fig:fd-sliding-maps>

=== Picking Pixels <sec:hyp3-fd-roi-picking>

Both preceding approaches still rely on an automatic amplitude/band gate to
decide which cross-spectrum cells enter the fit inside whatever spatial
region is under consideration. The last variant removes that gate as well:
using the #link("https://napari.org")[napari] interactive image viewer, the
$(k_z, k_x)$ cells (or, in a second variant, the depth/radial pixels) that
should enter the fit are selected by hand, directly on top of the amplitude
image or the difference B-scan.

Picking the cross-spectrum directly lets the fit be restricted to exactly the
two energy lobes visible in the display, without a relative amplitude
threshold or a fixed $(k_z, k_x)$ band. For the Chase-stage pair (3→8),
@fig:fd-napari-kspace compares a tight manual pick around the lobe peaks
($194$ pixels) against the automatic band-plus-amplitude gate ($1300$
pixels): $#Dz = +1.147 "m"$, $#Dx = -0.129 "m"$ (manual) versus
$#Dz = +1.046 "m"$, $#Dx = -0.113 "m"$ (automatic) -- close agreement despite
using $15$ times fewer pixels, which is a useful sanity check that the
automatic gate is not being pulled off the true peak by the peripheral
cells @sec:meth-phaseplane's amplitude-cubed weighting is designed to
suppress. A naive wide, unweighted pick that includes only one of the two
lobes was tried during development and returned a badly *sign-reversed*
estimate; picking both lobes, or picking a narrow high-amplitude region
around the peak of one, both recover the correct sign and magnitude, so the
failure is a leakage artefact of an imprecise selection rather than a flaw in
the mapping itself.

#figure(
  cimg("FD_napari_kspace_manual_vs_auto.png"),
  caption: [Manually-picked versus automatic $(k_z, k_x)$ selection for pair
    3→8 (Gazdag), shown on the cross-spectrum phase. Green outlines mark the
    pixels entering each fit.],
) <fig:fd-napari-kspace>

Picking directly on the difference B-scan instead lets the ROI trace the
*actual shape* of the reflection rather than an axis-aligned box. Because a
hand-painted mask has hard edges that would otherwise leak spectral energy
the same way the naive k-space pick above did, the painted mask is
Gaussian-softened before being used as the spatial window for the
cross-spectrum, in place of the rectangular window's edge-tapered box.
@fig:fd-napari-bscan compares this against the rectangular ROI for the same
pair (3→8): the hand-traced selection gives $#Dz = +1.273 "m"$,
$#Dx = -0.152 "m"$, somewhat larger than the rectangular ROI's
$#Dz = +1.046 "m"$, $#Dx = -0.113 "m"$ (identical, as it should be, to the
automatic k-space result above -- both apply the same automatic gate to the
same rectangular window). Tracing the reflection tightly excludes
lower-amplitude pixels the box includes near its edges, and here that shifts
the estimate upward by roughly $20$--$35%$ -- a visible, if modest,
sensitivity of the final number to exactly which pixels are allowed to
contribute.

#figure(
  cimg("FD_napari_bscan_manual_vs_rect.png"),
  caption: [Manually-painted versus rectangular ROI for pair 3→8 (Gazdag),
    shown on the difference B-scan. Green outline: painted selection;
    yellow box: rectangular ROI.],
) <fig:fd-napari-bscan>

Taken together, the three ROI strategies agree on where the signal is and on
the sign and rough magnitude of the displacement, while disagreeing at the
$10$--$30%$ level on the exact number -- a useful indication of the
estimate's sensitivity to ROI choice that a single hand-picked box alone
would not have revealed.

== Phase-Plane Fit Workflow <sec:hyp3-fd-phaseplane>

With the rectangular ROI of @sec:hyp3-fd-roi-rect fixed, three strategies for
chaining the pairwise phase-plane fit into a displacement trajectory across
all 38 profiles are compared:

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
$#Dx approx -0.18 "m"$ at profile 38 relative to profile 1. These stage-level
magnitudes are consistent with the Push/Chase/Wait/Pull ordering already seen
qualitatively in the sliding-window maps of @fig:fd-sliding-maps and the
single-pair estimates of @sec:hyp3-fd-roi-picking, which used a different
(shorter) set of representative pairs and, for Chase, a somewhat larger
$#Dz$ estimate -- the two are not directly comparable (different pairs,
different ROI-selection method) but agree on sign and order of magnitude.

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

== Interpretation <sec:hyp3-fd-interpretation>

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
// result agree with the synthetic fluid-front experiment? (4) the ROI-selection
// comparison of @sec:hyp3-fd-roi shows a 10-30% spread in the displacement
// estimate depending on whether the region is chosen by hand, swept
// automatically, or painted directly onto the data -- does this spread bound a
// meaningful "ROI-choice uncertainty" that should be reported alongside the
// Strategy-3 numbers in @tab:fielddata-stages? State these explicitly once the
// field context is available.]
