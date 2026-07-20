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

The raw profiles share a common five-step conditioning chain -- bandpass
filter ($0.02$--$0.20 "GHz"$), DC and direct-wave removal, trace alignment
(×5 upsampling), SVD rank-1 direct-wave suppression, and time-lapse
differencing against the reference profile -- after which Kirchhoff/Gazdag
migration and back-propagation diverge. Kirchhoff and Gazdag apply a linear
spherical-gain correction (compensating $1\/r$ geometric spreading) directly
to the differenced B-scan before migrating. Back-propagation instead applies
a physics-based 3D-to-2D Green's-function correction -- a $sqrt(r)$
amplitude gain together with a $1\/sqrt(omega) dot e^(i phi)$ frequency-domain
phase filter -- because the survey records a 3D field but gprMax
back-propagates it through a 2D grid, and a 3D point source and a 2D line
source have different Green's functions; without this correction the
resulting amplitude and phase mismatch propagates directly into the
back-propagated image. The corrected B-scan is then spatially and temporally
Tukey-tapered, passed through an $f$-$k_z$ dip filter to remove evanescent
energy, time-reversed and peak-normalised (range $[-1,1]$), low-pass
filtered to gprMax's numerical-dispersion limit, and edge-zeroed before
injection. @fig:fd-pipeline summarises the full pipeline end-to-end,
including the downstream region-of-influence selection (@sec:hyp3-fd-roi)
and phase-plane fitting (@sec:hyp3-fd-phaseplane) stages described later in
this chapter.

#let _pipe-box(title, content, fill: rgb("#EFEFEF")) = block(
  fill: fill,
  stroke: 0.6pt + black,
  radius: 3pt,
  inset: (x: 0.8em, y: 0.6em),
  width: 100%,
  [
    #text(weight: "bold", size: 9.5pt)[#title]
    #v(0.3em)
    #text(size: 8.5pt)[#content]
  ]
)
#let _pipe-arrow = align(center, text(size: 15pt)[↓])

#figure(
  block(width: 100%)[
    #align(center)[#block(width: 65%)[
      #_pipe-box("Raw B-scans", [38 zero-offset borehole profiles, plus a
        pre-injection reference profile (prof. 0), acquired 6 June 2016])
    ]]
    #_pipe-arrow
    #align(center)[#block(width: 92%)[
      #_pipe-box("Shared Pre-Processing (every profile)", [
        1) Band-pass filter ($0.02$–$0.20 "GHz"$)\
        2) DC / direct-wave removal (`remove_mean`)\
        3) Trace alignment (`align_traces`, ×5 upsample)\
        4) SVD rank-0–1 direct-wave suppression (`remove_svd`)\
        5) Difference against the reference profile
      ])
    ]]
    #_pipe-arrow
    #grid(
      columns: (1fr, 1fr),
      gutter: 1.2em,
      _pipe-box("Kirchhoff / Gazdag Migration", [
        Linear ($1\/r$) gain correction\
        Kirchhoff migration\
        Gazdag migration
      ], fill: rgb("#E8F0FE")),
      _pipe-box("Back-Propagation (gprMax)", [
        $sqrt(r)$ gain $+$ $1\/sqrt(omega) dot e^(i phi)$ phase correction (3D-to-2D)\
        Spatial + temporal Tukey tapers\
        $f$-$k_z$ dip filter (evanescent-energy removal)\
        Time-reversal, peak-normalise, dispersion-limit low-pass, edge-zero\
        Explicit borehole geometry ($d_x = 0.02 "m"$, scaled permittivity)\
        gprMax FDTD run → $E_z$ focus-frame extraction (fixed snapshot offset)\
        Post-imaging: dip (fan) filter + amplitude gate + near-borehole radial taper
      ], fill: rgb("#FCE8E6")),
    )
    #_pipe-arrow
    #align(center)[#block(width: 92%)[
      #_pipe-box("Time-Lapse Differencing", [Consecutive (profile $n$ minus
        profile $n-1$) / against a fixed reference profile / by operational
        phase (Push / Chase / Wait / Pull)])
    ]]
    #_pipe-arrow
    #align(center)[#block(width: 92%)[
      #_pipe-box("Region-of-Influence Selection", [Rectangular window /
        sliding window / manual (napari) picking, on either the
        cross-spectrum or the difference B-scan])
    ]]
    #_pipe-arrow
    #align(center)[#block(width: 65%)[
      #_pipe-box("WLS Cross-Spectrum Phase-Plane Fit",
        [Displacement estimate: #Dz (depth), #Dx (radial)],
        fill: rgb("#E6F4EA"))
    ]]
  ],
  caption: [Field-data processing pipeline: shared pre-processing, the
    Kirchhoff/Gazdag and back-propagation branches, time-lapse differencing,
    region-of-influence selection and the WLS phase-plane displacement fit.],
) <fig:fd-pipeline>

#draftnote[the sign of the $1\/sqrt(omega)$ filter's $pi\/4$ phase term, and
whether a water-level regularisation should be added near DC, were both left
open pending a dedicated validation experiment (single synthetic reflector,
compared against Kirchhoff/Gazdag as ground truth). State the resolved
choice here once confirmed.]

*Borehole geometry.* Back-propagation now additionally models the borehole
itself, rather than treating the medium around the source array as a single
homogeneous material: a real single-hole survey has both antennas inside a
fluid-filled channel whose permittivity ($approx 81$ for water) differs
sharply from the surrounding rock/ice ($approx 9$), a discontinuity the
homogeneous model omitted entirely. The gprMax domain now includes an
explicit $10 "cm"$-wide, water-filled rectangle spanning the full modelled
depth range, with the source/receiver on its lateral centreline, $1 "m"$ of
background medium to its left (radial direction), and at least $12 "m"$ of
imaging clearance to its right. Both the background and the borehole water
use the same $times 4$ permittivity scaling ($v -> #vmig = v \/ 2$) that
implements the exploding-reflector convention shared with Kirchhoff and
Gazdag, so that the one-way/two-way travel-time equivalence the scaling
relies on holds consistently across every material in the model, not only
the background.

Modelling the borehole this way initially made results worse, not better,
than the homogeneous-domain baseline. At the grid spacing used elsewhere in
this study ($d_x = 0.05 "m"$), the $10 "cm"$ channel is resolved by only two
grid cells -- too coarse to represent a sharp, high-contrast material
boundary without introducing numerical stair-casing artefacts on top of the
genuine physics. A diagnostic on profile 1, isolating grid resolution from
the permittivity-scaling choice, confirmed that refining to $d_x = 0.02 "m"$
(five cells across the channel) removes most of the excess clutter and
recovers a reflector consistent in position with Kirchhoff and Gazdag; using
the borehole's true, unscaled permittivity instead of the scaled value did
not help and introduced a small position shift, so the scaled value is
retained. An alternative amplitude normalisation (rescaling each
time-reversed trace to $[0,1]$ rather than $[-1,1]$ around zero) was also
tested and rejected: it removes the physical zero baseline from a bipolar
field, and the resulting non-zero DC component injected by every current
source overwhelmed the simulation with spurious energy -- far worse than the
coarse-grid clutter it was meant to fix.

The refined ($d_x = 0.02 "m"$, explicit borehole) pipeline described above
has since been extended from profile 1 to the full five-profile set used
throughout the remainder of this chapter (1, 3, 8, 20, 38); the two
cross-profile registration issues that surfaced once time-lapse differences
were taken across this extended set, and how they were fixed, are described
in @sec:hyp3-fd-crossprofile below. The 14-profile back-propagation results
referenced below (@fig:fd-profile-grid onward, @fig:fd-disp-backprop) still
use the earlier homogeneous-domain model at $d_x = 0.05 "m"$ with no explicit
borehole material and predate the fixes of @sec:hyp3-fd-crossprofile; they
have not been re-run under the refined geometry.

#draftnote[source positions no longer sit at the domain edge ($y approx 0$)
under the refined geometry but on the borehole's centreline, so the
injection-halo masking distance discussed below (currently derived for the
homogeneous-domain, edge-sourced model) should be re-derived for the refined
geometry if the 14-profile homogeneous-domain results are ever re-run under
it.]

Kirchhoff and Gazdag migrations are available for all 37 processed profiles
(profiles 2--38 differenced against profile 1 as the fixed baseline);
back-propagation (homogeneous-domain model) is available for 14 profiles
(profiles 1--5, 7--10, 13, 16, 20, 21, 38), the remainder still pending
gprMax forward runs.

Back-propagation additionally requires suppressing an *injection halo*: every
gprMax source in the time-reversal simulation is injected at the borehole
wall ($y approx 0 "m"$ in this homogeneous-domain model), and the
destructive interference among sources that should cancel the field away
from the true reflector is never fully complete near $y = 0$, leaving
residual energy at small radial distance regardless of how the input is
pre-processed. This is a structural property of single-sided time-reversal
from a borehole source array, not a pre-processing shortcoming, and is
suppressed for display by zeroing the first $2.0 "m"$ of the radial axis.
One candidate fix was tested and reverted: normalising each receiver trace
by its RMS value amplifies low-SNR traces far from the fluid front, and
these noisy traces then back-propagate *coherently* toward the borehole
axis, making the halo worse rather than better.

*Post-imaging clutter removal.* The steps above all act before or during
injection into gprMax; a further round of processing was explored on the
completed back-propagation focus image itself (profile 1 only), following
two published time-reversal/migration post-processing schemes. Gaussian
smoothing followed by decomposition into low-rank and sparse components via
Robust Principal Component Analysis (RPCA), after Li and Yan (2021), was
tested and rejected: RPCA separates an image by amplitude sparsity, which
suits a compact, point-like target, but this reflector's response is
spatially extended along a consistent dip, so it behaves as low-rank rather
than sparse and RPCA could not isolate it from the low-rank background even
after sweeping the sparsity weight. An $f$-$k$ dip (fan) filter was adopted
instead, separating the reflector from clutter by orientation rather than
amplitude: the reflector's dip is estimated directly from the image (a
per-column peak-amplitude pick over a sub-window containing the primary
reflection, fitted with a straight line) rather than assumed, and a
cosine-tapered wedge of $plus.minus 20 "deg"$ around that orientation in
the $(k_z, k_x)$ domain is kept, with the taper avoiding the Gibbs ringing a
hard-edged wedge produces. Widening the fan enough to preserve the reflector
introduced a new artefact: incoherent noise sharing the target's orientation
was reconstructed as coherent diagonal streaks throughout the image, since
the filter is global and translation-invariant with no notion of proximity
to the true reflector. An amplitude gate built from the raw image's own
local energy envelope removes this cleanly; a windowed (spatially localised)
version of the fan filter was also tested and rejected, since restricting
the FFT to a local window does not stop locally-oriented noise from being
reconstructed within that window -- the artefact is a property of what
survives the dip criterion, not of the transform's spatial support. A final
cosine-ramped taper suppresses residual clutter between $1.0$ and
$3.0 "m"$ radial distance (widened from an initial $2.5 "m"$ after visual
inspection showed clutter persisting slightly past that point), beyond the
hard injection-halo mask already applied above.

This recipe is now applied identically across the five profiles used
throughout this chapter (1, 3, 8, 20, 38): the dip fit, fan filter, gate and
taper are all recomputed independently per profile, since each profile's
reflector sits at a different position. Profiles 3, 8, 20 and 38 additionally
required new $d_x = 0.02 "m"$ gprMax back-propagation runs (profile 1's
already existed from the grid-resolution diagnostic above) with the same
scaled-water-permittivity, peak-normalised settings.

This post-imaging recipe (20-degree tapered fan filter, amplitude gate and
near-borehole radial taper) is applied purely to the output focus image; it
does not change the gprMax injection file, which still only uses the
peak-normalisation choice described above.

=== Cross-Profile Registration and Amplitude Normalisation <sec:hyp3-fd-crossprofile>

Applying the recipe above independently to all five profiles and taking
consecutive time-lapse differences (profile $n$ minus profile $n-1$) at
first produced differences that looked nothing like Kirchhoff/Gazdag's
compact, localised residuals (@fig:fd-stage-push to @fig:fd-stage-pull):
instead, each difference showed near-complete smearing across the full
extent of the reflector, as if the two profiles being compared disagreed
almost everywhere rather than only where the fluid front had actually moved.
Isolating the cause -- by differencing raw, unfiltered back-propagation
frames at matched snapshot times and comparing against the fan-filtered
result -- ruled out the dip filter itself and traced the problem to two
independent, structural bugs in how the back-propagation focus frames were
selected and normalised.

*Snapshot-timing misalignment.* The focus frame for each profile was
originally selected by an independent per-profile search for the gprMax
snapshot with the highest masked peak amplitude -- a reasonable choice for
viewing one profile in isolation, but one that let each profile lock onto a
different snapshot *time* (indices $25$--$28$ out of the saved snapshot
sequence, a spread of $approx 26 "ns"$). Differencing two focus frames taken
at different times makes a still-converging, not-yet-focused wavefront look
like it moved, even with zero real displacement -- exactly the near-full
extent smearing observed. The fix mirrors the fixed-time convention already
used for the homogeneous-domain model's own consecutive differencing
(`FOCUS_IDX_OFFSET`, @sec:hyp3-fd-processing): every profile now uses the
same nominal focus time plus the same fixed snapshot-index offset, calibrated
once against profile 1's previously-validated best-focus index.

*Per-trace amplitude normalisation.* The fixed-timing frames alone did not
resolve the smearing -- differencing raw, unfiltered frames at matched times
still showed the same near-full-reflector residual -- which pointed at
amplitude rather than timing. The gprMax injection preparation
(`write_borehole_backprop_files`, @sec:hyp3-fd-processing) normalised every
trace to its own individual peak amplitude before injection. Measured
directly across the dataset: profile 3's real, pre-normalisation RMS
amplitude is $approx 1.8$ times profile 1's -- consistent with the reflector
genuinely sharpening through the push stage, exactly the kind of change a
time-lapse study is trying to detect -- and per-trace rescaling erased that
signal entirely, along with within-profile amplitude structure (individual
trace peaks varied $80$--$160 times$ before normalisation, so weak,
mostly-noise traces were boosted to the same injected amplitude as the
strongest genuine reflections). The fix replaces the per-trace scalar with a
single amplitude shared across all five profiles -- the peak found across
all five profiles' excitation data -- so every trace still divides by a
constant, as gprMax's numerical stability requires, but the same constant for
every profile, preserving both the within-profile and cross-profile
amplitude structure.

With both fixes applied and all five profiles re-run through gprMax,
@fig:fd-borehole-final-grid shows the raw focus frame, the final (fan +
gate + taper) processed frame, and the consecutive time-lapse difference for
each of the five profiles: the differences are now compact and localised,
qualitatively matching the character of the Kirchhoff/Gazdag residuals shown
earlier in this chapter, a marked change from the near-full-reflector
smearing the two bugs above had produced. One caveat carried forward
transparently: the per-column automatic dip fit that seeds the fan filter's
orientation, reliable for profile 1 alone, did not generalise once the
timing and normalisation fixes changed the underlying frames -- it also
began locking onto near-source clutter for profiles 1 and 3, not only
8/20/38 as before the fixes. All five profiles shown in
@fig:fd-borehole-final-grid therefore use the same manually-specified dip
override ($m_0 = 2.2$) rather than five independent automatic fits; the
robustness of the automatic dip estimator itself is not addressed further
here.

#figure(
  cimg("FD_borehole_final_processed_all_profiles.png"),
  caption: [Corrected back-propagation pipeline (fixed snapshot timing +
    shared global amplitude normalisation), all five profiles: raw focus
    frame (top row), final fan+gate+taper-processed frame (middle row), and
    consecutive time-lapse difference (bottom row, profile $n$ minus profile
    $n-1$).],
) <fig:fd-borehole-final-grid>

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
increasingly automatic ways of making it are compared here, applied to the
same Gazdag-migrated images and, from @sec:hyp3-fd-roi-sliding onward, the
same four representative pairs (1→3, 3→8, 8→20, 20→38 -- approximating the
Push, Chase, Wait and Pull stages respectively).
#draftnote[@fig:fd-stage-push to @fig:fd-stage-pull below (generated by the
separate Thesis Figure Compilations pipeline, out of scope for this
session's pair/sign corrections) still use the older stage-boundary pairs
(1→4, 5→9, 10→20, 21→38); regenerate them with the representative pairs, or
adjust this paragraph, once that pipeline is revisited.]

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
threshold or a fixed $(k_z, k_x)$ band. The $#Dz$ figures quoted throughout
this subsection are the *raw* WLS fit output, not yet converted to the
"positive = downward" convention established later in @sec:hyp3-fd-phaseplane
(@tab:fielddata-stages) -- the internal, relative comparisons below (manual
vs. automatic vs. painted) are unaffected either way, since all three share
the same raw convention consistently. For the Chase-stage pair (3→8),
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
  span the full stage (e.g. profiles 1→3 for Push, 3→8 for Chase), giving the
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
    [Push],  [1→3],   [$-1.32 "m"$], [$-0.26 "m"$],
    [Chase], [3→8],   [$-0.76 "m"$], [$-0.06 "m"$],
    [Wait],  [8→20],  [$+0.30 "m"$], [$+0.03 "m"$],
    [Pull],  [20→38], [$+0.28 "m"$], [$+0.15 "m"$],
    table.hline(stroke: 0.4pt),
    [*Net (prof 1→38)*], [], [$-1.50 "m"$], [$-0.15 "m"$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Stage-anchored (Strategy 3b) phase-plane displacement estimates
    from the borehole GPR field dataset, using the same four representative
    consecutive pairs as the rest of this chapter (1→3, 3→8, 8→20, 20→38).
    Positive $#Dz$ is downward (along the borehole, increasing depth);
    positive $#Dx$ is away from the borehole. Net values are the
    reconstructed global trajectory after stitching the four stage
    estimates -- a plain cumulative sum here, since consecutive pairs need no
    separate bridging step.],
  kind: table,
) <tab:fielddata-stages>

#draftnote[two corrections were made to this table this session, independent
of each other: (1) the sign of both columns -- verified via a synthetic test
with a known, directly injected shift, run against the exact WLS fit function
and depth/radial axis conventions used to produce this table -- confirmed the
raw fit output equals $-#Dz$ under the "positive = downward" convention (an
artefact of the depth axis decreasing with row index) but equals $+#Dx$
directly under "positive = away" (no flip needed there); (2) the profile
pairs themselves -- the table previously used stage-*boundary* pairs (1→4,
5→9, 10→20, 21→38) with three additional bridging pairs to chain them, which
had drifted out of sync with the representative pairs (1→3, 3→8, 8→20,
20→38) used everywhere else in this chapter and were producing numbers that
no longer matched any other part of the analysis. Both fixes are applied to
the `strat3b_cd` notebook cell and verified by direct execution against the
cached Gazdag migrations.]

The dominant signal is an *upward* displacement of approximately $1.32 "m"$
during the Push stage, accompanied by a lateral shift of $0.26 "m"$ *toward*
the borehole -- the opposite direction, on both axes, from what a naive
"downward and outward" injection picture would predict. The Chase stage adds
a smaller increment in the same (upward, inward) direction, the Wait stage
partially reverses it (downward, outward, but smaller in magnitude than
either Push or Chase), and the Pull stage continues in that reversing
direction, leaving a net residual of $#Dz approx -1.50 "m"$ (upward),
$#Dx approx -0.15 "m"$ (toward the borehole) at profile 38 relative to
profile 1. @sec:hyp3-fd-interpretation returns to whether this reversed
direction is physically plausible for this experiment -- it is reported here
as a corrected, sign-verified number, not yet as a settled physical
interpretation. These stage-level magnitudes are consistent in *sign* with
the Push/Chase/Wait/Pull ordering already seen qualitatively in the
sliding-window maps of @fig:fd-sliding-maps, the single-pair estimates of
@sec:hyp3-fd-roi-picking, and the corrected back-propagation estimates of
@tab:fielddata-bp-corrected -- all four now use the *same* representative
pairs, so this is a genuine cross-check rather than a qualitative one.
Magnitudes still differ meaningfully between this table and
@sec:hyp3-fd-roi-picking's single-pair estimate for the same Chase pair
(3→8): $-0.76 "m"$ here against $-1.05$ to $-1.27 "m"$ there, a larger
relative spread than the $10$--$30%$ already noted for ROI-selection choice
alone -- because this table's fit (`_estimate_shift_2d`, a fixed
padding/band/threshold configuration shared with Strategies 1 and 2) was
never re-tuned to match the rectangular-window recipe's per-technique
defaults established in @sec:hyp3-fd-roi. Both are legitimate estimates of
the same underlying displacement under different fit configurations, not a
contradiction, but the spread is a reminder that the *fit configuration*, not
only the ROI, is a real source of estimate variability.

@fig:fd-disp-kirchhoff to @fig:fd-disp-backprop cross-check this result across
migration techniques: the same WLS cross-spectrum phase-plane fit
(@sec:meth-phaseplane), applied to *stage-boundary* pairs (1→4, 5→9, 10→20,
21→38 -- distinct from @tab:fielddata-stages's representative pairs above)
and a shared ROI, is shown for Kirchhoff-BP, Gazdag and back-propagation side
by side, with the full five-panel diagnostic (cross-spectrum phase,
cross-spectrum energy with the WLS amplitude-threshold contour, the fitted
plane, and the 1-D $#kz$ and $#kx$ slices with their fits) for every stage.
These figures predate this session's pair/sign corrections and were not
regenerated with them; since they use different profile pairs from
@tab:fielddata-stages entirely (not just a different fit configuration, as
with @sec:hyp3-fd-roi-picking above), the two sets of numbers are not
directly comparable -- they are a consistency check on *direction and
relative magnitude*, not a replacement for the table.
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

== Corrected Back-Propagation Displacement Re-Estimation <sec:hyp3-fd-bp-corrected>

With the timing and normalisation fixes of @sec:hyp3-fd-crossprofile applied,
the manually-painted-ROI phase-plane fit of @sec:hyp3-fd-roi-picking was
re-run directly on the five corrected, final-processed focus frames, for the
same four representative consecutive pairs used throughout
@sec:hyp3-fd-roi (1→3, 3→8, 8→20, 20→38). Painting is done directly on the
difference image rather than the cross-spectrum -- the same B-scan-painting
variant of @fig:fd-napari-bscan -- since a Gaussian-softened painted mask
avoids the hard-edge spectral leakage a rectangular window would introduce,
and the reflector's shape (and, in two of the four pairs, a visible
side-lobe not present in the rectangular ROI) is easier to trace by eye on
the spatial difference than on the cross-spectrum phase.

A third bug surfaced while validating this re-estimation, independent of the
two in @sec:hyp3-fd-crossprofile: the corrected back-propagation frames'
depth axis increases with row index (row 0 $approx 60 "m"$, shallow; the
last row $approx 85 "m"$, deep), the opposite convention to the
Kirchhoff/Gazdag depth array used everywhere else in this chapter, which
*decreases* with row index (row 0 $= 85 "m"$, deep). The phase-plane fit
measures displacement along the increasing-row-index direction, so feeding
it a plain positive row spacing -- as @sec:meth-phaseplane's fit does
uniformly -- would silently report $+#Dz$ (raw fit output) for a shift
toward *shallower* depth under the Kirchhoff/Gazdag row-index convention,
but $+#Dz$ for a shift toward *deeper* depth under the back-propagation
row-index convention: the same physical event, opposite-signed raw number.
The fix negates the row-spacing constant used to build the fit's $k_z$ axis
for the back-propagation frames only, so the *raw* fit output means the same
physical row-index direction under both conventions -- the values below are
already reported after that fix, and additionally converted to the
"positive $#Dz$ = downward" convention established for
@tab:fielddata-stages (@sec:hyp3-fd-phaseplane), i.e. negated once more
relative to the raw fit output, consistently with that table.

@tab:fielddata-bp-corrected reports the corrected, sign-converted estimates.
Direction now agrees with the corrected @tab:fielddata-stages throughout:
upward and toward the borehole during Push and Chase, downward and away
during Wait and Pull -- the same reversed-from-naive-expectation pattern
discussed above, reproduced independently by a technique that shares no
processing steps with Kirchhoff/Gazdag downstream of the raw B-scans (one
migrates a differenced B-scan; the other back-propagates a gprMax simulation
through an explicit borehole geometry). Magnitudes are consistent in order
with the established Gazdag results for the same representative pairs, once
those are read under the same corrected convention (@sec:hyp3-fd-roi-picking):
the Chase-stage pair (3→8) gives $#Dz = -1.37 "m"$ here against
$-1.05$ to $-1.27 "m"$ across the three Gazdag ROI-selection variants, and
the Wait-stage pair (8→20) gives the smallest magnitude of the four here
($#Dz = +0.54 "m"$), matching the near-zero-net expectation of the
paused-injection stage already established qualitatively in
@fig:fd-sliding-maps. Two independently-processed techniques agreeing on
both sign and order of magnitude in every stage is a stronger check than
either alone, but it also means the reversed direction is very unlikely to
be a processing artefact specific to one pipeline -- @sec:hyp3-fd-interpretation
returns to what it might mean instead.

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Pair*], [*Stage*], [*$#Dz$ (depth)*], [*$#Dx$ (radial)*],
    table.hline(stroke: 0.4pt),
    [1→3],   [Push],  [$-1.91 "m"$], [$-0.34 "m"$],
    [3→8],   [Chase], [$-1.37 "m"$], [$-0.21 "m"$],
    [8→20],  [Wait],  [$+0.54 "m"$], [$+0.08 "m"$],
    [20→38], [Pull],  [$+0.75 "m"$], [$+0.14 "m"$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Corrected back-propagation pipeline: manually-painted-ROI
    phase-plane displacement estimates for the four representative
    consecutive pairs, after the timing, normalisation and depth-axis-sign
    fixes of @sec:hyp3-fd-crossprofile and @sec:hyp3-fd-bp-corrected.
    Positive $#Dz$ is downward and positive $#Dx$ is away from the borehole,
    matching the corrected convention of @tab:fielddata-stages.],
  kind: table,
) <tab:fielddata-bp-corrected>

#figure(
  subfigs(cols: 2,
    cimg("FD_borehole_napari_wls_1_to_3.png"),
    cimg("FD_borehole_napari_wls_3_to_8.png"),
    cimg("FD_borehole_napari_wls_8_to_20.png"),
    cimg("FD_borehole_napari_wls_20_to_38.png"),
  ),
  caption: [Corrected back-propagation pipeline, painted-ROI phase-plane fit
    for the four representative pairs: (a) 1→3 (Push), (b) 3→8 (Chase), (c)
    8→20 (Wait), (d) 20→38 (Pull). Green outline marks the painted region on
    the time-lapse difference image. The $#Dz$/$#Dx$ annotations in each
    panel are the *raw* fit output (row-index-direction convention, fixed
    for the depth-axis bug described above but *not* yet converted to the
    "positive = downward" convention of @tab:fielddata-stages) -- compare
    @tab:fielddata-bp-corrected for the sign-converted values used in the
    text.],
) <fig:fd-borehole-napari-pairs>

== Interpretation <sec:hyp3-fd-interpretation>

*On the back-propagation/Kirchhoff-Gazdag sign disagreement.* The
14-profile, homogeneous-domain back-propagation result of
@fig:fd-disp-backprop disagreeing with Kirchhoff/Gazdag on the sign of
$#Dz$ (noted above, @sec:hyp3-fd-phaseplane) was previously an open
question: a genuine sign-convention difference between the two coordinate
systems, back-propagation's lower SNR, or a real technique-dependent
reliability issue consistent with @ch:hyp2. @sec:hyp3-fd-bp-corrected's
independent re-derivation resolves this in favour of the first explanation
for the *corrected* ($d_x = 0.02 "m"$, five-profile, fixed-timing,
globally-normalised) back-propagation pipeline: once the back-propagation
depth axis's opposite row-index direction is accounted for
(@sec:hyp3-fd-bp-corrected), the two techniques agree on sign in every
stage, despite sharing no processing steps downstream of the raw B-scans.
This does not by itself explain the disagreement seen in the *original*
14-profile homogeneous-domain figures (@fig:fd-disp-backprop), which predate
all three fixes of @sec:hyp3-fd-crossprofile and @sec:hyp3-fd-bp-corrected
and have not been re-run under them; the SNR and technique-reliability
explanations remain open for that specific figure.

#draftnote[fill in the remaining physical interpretation once the borehole
geometry, injection depth, and fluid-injection parameters are confirmed from
the field survey metadata. Key questions to address: (1) do the inferred
$#Dz$ and $#Dx$ values agree with the known injection depth and the expected
lateral spread for the given fracture geometry -- and, now that
@tab:fielddata-stages and @tab:fielddata-bp-corrected have been corrected to
a verified sign convention (@sec:hyp3-fd-phaseplane), both agree that Push
and Chase move *upward and toward the borehole* rather than the
"downward and outward" direction a naive injection picture would predict:
is there a physical reason (fracture geometry, the specific reflector being
tracked vs. the injection point, a coordinate-system offset in how "depth"
and "radial distance" map onto the actual borehole/fracture geometry) this
reversed direction should be expected, or does it warrant re-checking against
the field survey log before being reported as a finding? (2) does the field
result agree with the synthetic fluid-front experiment? (3) the ROI-selection
comparison of @sec:hyp3-fd-roi shows a 10-30% spread in the displacement
estimate depending on whether the region is chosen by hand, swept
automatically, or painted directly onto the data -- does this spread bound a
meaningful "ROI-choice uncertainty" that should be reported alongside the
Strategy-3 numbers in @tab:fielddata-stages? State these explicitly once the
field context is available.]
