#import "../template.typ": *

= Hypothesis 3: Generalisation to Real Field Data <ch:hyp3>

The time-lapse phase-plane approach recovers a sub-wavelength displacement by
fitting a plane to the cross-spectrum phase of a baseline and a monitor
migrated image (@ch:theory). @ch:hyp1 and @ch:hyp2 established that it works on
idealised synthetic data --- a single known target, a controlled displacement,
and, under @ch:hyp2's noise model, added field noise. Real surveys offer none
of those guarantees. This chapter asks whether the approach generalises to real
borehole GPR field data, where neither the target geometry nor the true
displacement is known, and --- because a real reflector never comes with a
pre-labelled region of interest --- what happens to the estimate as the method
used to *choose* that region is varied.

#para-head[Hypothesis 3.] The time-lapse phase-plane approach generalises
beyond idealised single-scatterer synthetic models to real borehole GPR field
data.

The complementary generalisation to complex synthetic scenes with multiple
independently-moving scatterers is scoped as future work rather than tested
here; this chapter is concerned solely with the step from controlled synthetic
targets to uncontrolled real data.

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

The refined geometry requires a finer grid than the rest of the study: at
$d_x = 0.05 "m"$ the $10 "cm"$ channel spans only two cells and stair-cases
the sharp water/rock boundary, whereas $d_x = 0.02 "m"$ (five cells) removes
most of the excess clutter and recovers a reflector consistent in position
with Kirchhoff and Gazdag. The scaled ($times 4$) permittivity is retained for
the borehole water as well -- using its true unscaled value did not help and
introduced a small position shift -- and rescaling each time-reversed trace to
$[0,1]$ rather than $[-1,1]$ was rejected, since it injects a spurious DC
component that overwhelms the bipolar field.

The refined ($d_x = 0.02 "m"$, explicit borehole) pipeline described above
has since been extended from profile 1 to the full five-profile set used
throughout the remainder of this chapter (1, 3, 8, 20, 38); the two
cross-profile registration issues that surfaced once time-lapse differences
were taken across this extended set, and how they were fixed, are described
in @sec:hyp3-fd-crossprofile below. The 14-profile back-propagation results
referenced below (@fig:fd-profile-grid onward, and the back-propagation
diagnostic in the Supplementary Material, §S4.4.2) still
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
wall, and the destructive interference that should cancel the field away from
the true reflector is never complete near $y = 0$, leaving residual energy at
small radial distance. This is a structural property of single-sided
time-reversal from a borehole source array, not a pre-processing shortcoming,
and is suppressed for display by zeroing the first $2.0 "m"$ of the radial
axis. (Per-trace RMS normalisation was tried as a fix and reverted: it
amplifies low-SNR far traces that then back-propagate coherently toward the
axis, worsening the halo.)

*Post-imaging clutter removal.* A final round of processing acts on the
completed back-propagation focus image. The adopted recipe is a per-profile
$f$-$k$ dip (fan) filter -- the reflector's dip is estimated directly from the
image and a cosine-tapered $plus.minus 20 "deg"$ wedge around that orientation
in the $(k_z, k_x)$ domain is kept -- followed by an amplitude gate built from
the image's own local energy envelope (which removes the coherent diagonal
streaks a wider fan otherwise reconstructs from equally-oriented noise) and a
cosine-ramped near-borehole taper suppressing residual clutter out to
$3.0 "m"$ radial distance. Two alternatives were rejected: Gaussian smoothing
plus Robust PCA decomposition (after Li and Yan, 2021), because this
reflector's response is spatially extended along a consistent dip and so
behaves as low-rank rather than sparse; and a spatially-windowed fan filter,
because restricting the FFT to a local window does not stop locally-oriented
noise from being reconstructed within it. The recipe is recomputed
independently per profile (each reflector sits at a different position) and
acts purely on the output focus image -- it does not change the gprMax
injection file.

=== Cross-Profile Registration and Amplitude Normalisation <sec:hyp3-fd-crossprofile>

Taking consecutive time-lapse differences across the five-profile set at first
produced differences that smeared across the full extent of the reflector
rather than the compact residuals Kirchhoff/Gazdag produce. Isolating the cause
(differencing raw, unfiltered frames at matched snapshot times) ruled out the
dip filter and traced it to two structural issues in how the back-propagation
focus frames were selected and normalised, both since fixed:

+ *Snapshot timing.* Selecting each profile's focus frame by an independent
  per-profile peak-amplitude search let different profiles lock onto different
  snapshot _times_ (a spread of $approx 26 "ns"$); differencing frames at
  different times makes a still-converging wavefront look displaced. Every
  profile now uses the same nominal focus time plus a fixed snapshot-index
  offset (`FOCUS_IDX_OFFSET`), calibrated once against profile 1.

+ *Amplitude normalisation.* Normalising every trace to its own peak before
  injection erased both the genuine cross-profile amplitude growth (profile 3's
  pre-normalisation RMS is $approx 1.8 times$ profile 1's -- exactly the change
  a time-lapse study targets) and within-profile structure. A single amplitude
  shared across all five profiles -- the peak over all five profiles'
  excitation data -- replaces the per-trace scalar, preserving both.

With both fixes applied and all five profiles re-run, @fig:fd-borehole-final-grid
shows compact, localised differences matching the Kirchhoff/Gazdag character, a
marked change from the earlier smearing. (One caveat: the automatic per-column
dip fit did not generalise across the corrected frames, so all five profiles
use a single manual dip override, $m_0 = 2.2$.)

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
increasingly automatic ways of making it are compared here, all applied to the
same Gazdag-migrated images. The rectangular-window examples in this section
are shown at the four operational-stage boundaries (profiles 1→4, 5→9, 10→20,
21→38); from @sec:hyp3-fd-roi-sliding onward the analysis uses the four
representative pairs (1→3, 3→8, 8→20, 20→38 -- approximating the Push, Chase,
Wait and Pull stages) that carry through the rest of the chapter.

=== Rectangular Window <sec:hyp3-fd-roi-rect>

The baseline approach is a hand-picked axis-aligned box in depth and radial
distance, tuned per pair against the monogenic envelope of the time-lapse
difference image (@sec:meth-phaseplane), with a fixed cross-spectrum band and
amplitude threshold selecting the $(k_z, k_x)$ cells inside it. Using a fixed
ROI of depth $70$--$79 "m"$ and radial distance $4.5$--$7.5 "m"$,
@fig:fd-stage-push shows, for the Push stage, the time-lapse difference image
and its monogenic envelope with the ROI overlaid, computed independently for
all three migration techniques (Kirchhoff-BP, Gazdag, back-propagation). The
envelope highlights a single coherent patch inside the ROI, confirming that
the ROI is well-placed and that the "single dominant displacement" assumption
behind the phase-plane fit holds. The same construction applied to Chase,
Wait, and Pull shows the patch shrinking progressively through those stages,
mirroring the weakening signal expected as the fluid front decelerates.

#figure(
  cimg("FD_stage_diff_envelope_pushing.png"),
  caption: [Push stage (profiles 1→4): time-lapse difference and monogenic
    envelope, ROI overlaid, for Kirchhoff-BP, Gazdag and back-propagation.],
) <fig:fd-stage-push>

#supp-note[The Chase, Wait, and Pull stage envelope figures are provided in
the Supplementary Material, §S4.3.1.]

The rectangular window is simple and, as @sec:hyp3-fd-phaseplane shows,
produces physically sensible results, but it is hand-tuned per pair and says
nothing about *where else* in the domain the same fit would or would not be
reliable. The next two subsections replace the fixed box with, respectively,
a systematic scan and a manually-drawn selection.

=== Sliding Window <sec:hyp3-fd-roi-sliding>

Rather than picking one ROI by hand, a fixed-size window can instead be slid
across the entire difference B-scan, with the phase-plane fit repeated
independently in every window position to build up a 2-D map of $#Dz$ and
$#Dx$. In practice this scan proved less useful than the rectangular window
above as a primary result: it is far more expensive to compute and, once its
absolute-energy gate is tuned to suppress noise-only windows returning
"confident"-looking but meaningless fits, adds little beyond confirming what
the rectangular ROI already shows. Its value is instead as an independent
cross-check that requires no prior knowledge of where the reflector is: in
every pair, the surviving (non-blanked) region collapses onto a single
spatially-coherent patch at approximately $68$--$83 "m"$ depth and
$4$--$8 "m"$ radial distance -- the same location as the hand-picked ROI of
@sec:hyp3-fd-roi-rect -- with the same Push/Chase-largest, Wait-weakest
stage-to-stage trend seen there.

#supp-note[The sliding-window scan grid and the resulting $#Dz$/$#Dx$ maps
for all four pairs are provided in the Supplementary Material, §S4.3.2.]

=== Picking Pixels <sec:hyp3-fd-roi-picking>

Both preceding approaches still rely on an automatic amplitude/band gate to
decide which cross-spectrum cells enter the fit inside whatever spatial
region is under consideration. The last variant removes that gate as well:
using the #link("https://napari.org")[napari] interactive image viewer, the
$(k_z, k_x)$ cells (or, in a second variant, the depth/radial pixels) that
should enter the fit are selected by hand, directly on top of the amplitude
image or the difference B-scan. @fig:fd-roi-sensitivity compares both
variants, for the Chase-stage pair (3→8): (a) manual versus automatic
picking in $(k_z, k_x)$ space; (b) manual painting versus the rectangular
ROI on the difference B-scan itself. Gazdag is the worked example here; the
same painted-ROI workflow applied to the corrected back-propagation frames is
presented in @sec:hyp3-fd-bp-corrected.

Picking the cross-spectrum directly lets the fit be restricted to exactly the
two energy lobes visible in the display, without a relative amplitude
threshold or a fixed $(k_z, k_x)$ band. The $#Dz$ figures quoted throughout
this subsection are the *raw* WLS fit output, not yet converted to the
"positive = downward" convention established later in @sec:hyp3-fd-phaseplane
(@tab:fielddata-stages) -- the internal, relative comparisons below (manual
vs. automatic vs. painted) are unaffected either way, since all three share
the same raw convention consistently. For the Chase-stage pair (3→8),
@fig:fd-roi-sensitivity (a) compares a tight manual pick around the lobe peaks
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

Picking directly on the difference B-scan instead lets the ROI trace the
*actual shape* of the reflection rather than an axis-aligned box. Because a
hand-painted mask has hard edges that would otherwise leak spectral energy
the same way the naive k-space pick above did, the painted mask is
Gaussian-softened before being used as the spatial window for the
cross-spectrum, in place of the rectangular window's edge-tapered box.
@fig:fd-roi-sensitivity (b) compares this against the rectangular ROI for the
same pair (3→8): the hand-traced selection gives $#Dz = +1.273 "m"$,
$#Dx = -0.152 "m"$, somewhat larger than the rectangular ROI's
$#Dz = +1.046 "m"$, $#Dx = -0.113 "m"$ (identical, as it should be, to the
automatic k-space result above -- both apply the same automatic gate to the
same rectangular window). Tracing the reflection tightly excludes
lower-amplitude pixels the box includes near its edges, and here that shifts
the estimate upward by roughly $20$--$35%$ -- a visible, if modest,
sensitivity of the final number to exactly which pixels are allowed to
contribute.

#figure(
  subfigs(cols: 1,
    cimg("FD_napari_kspace_manual_vs_auto.png"),
    cimg("FD_napari_bscan_manual_vs_rect.png"),
  ),
  caption: [ROI-selection sensitivity for pair 3→8 (Gazdag): (a)
    manually-picked versus automatic $(k_z, k_x)$ selection, shown on the
    cross-spectrum phase; (b) manually-painted versus rectangular ROI, shown
    on the difference B-scan. Green outlines mark the pixels/region entering
    each fit; the yellow box in (b) marks the rectangular ROI.],
) <fig:fd-roi-sensitivity>

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

@fig:fd-strategies shows all three strategies applied across the 38 profiles.
The stage-anchored Strategy 3, which gives the highest-SNR estimate of the
displacement per stage, is the one carried forward into @tab:fielddata-stages.

#figure(
  subfigs(cols: 1,
    img("FD_strat1_consecutive_summary.png"),
    img("FD_strat2_baseline_summary.png"),
    img("FD_strat3_stage_anchored_summary.png"),
  ),
  caption: [Displacement trajectories over all 38 profiles under the three
    chaining strategies (Gazdag migration): (a) consecutive increments summed
    cumulatively; (b) each profile against the fixed profile-1 baseline; (c)
    stage-anchored. Each panel shows $#Dx$, $#Dz$ and the phase offset $phi_0$
    versus profile number.],
) <fig:fd-strategies>

The Strategy 3 (stage-anchored) results, which give the highest-SNR estimate
of the total displacement per stage, are summarised alongside the
independently-processed, corrected back-propagation cross-check in
@tab:fielddata-stages (@sec:hyp3-fd-bp-corrected), once that cross-check's
own methodology has been introduced below.

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
sliding-window maps (Supplementary Material, §S4.3.2), the single-pair estimates of
@sec:hyp3-fd-roi-picking, and the corrected back-propagation columns of
@tab:fielddata-stages below -- all four now use the *same* representative
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

@fig:fd-disp-gazdag cross-checks this result across migration techniques:
the same WLS cross-spectrum phase-plane fit (@sec:meth-phaseplane), applied
to *stage-boundary* pairs (1→4, 5→9, 10→20, 21→38 -- distinct from
@tab:fielddata-stages's representative pairs) and a shared ROI, is shown for
Kirchhoff-BP, Gazdag and back-propagation side by side, with the full
five-panel diagnostic (cross-spectrum phase, cross-spectrum energy with the
WLS amplitude-threshold contour, the fitted plane, and the 1-D $#kz$ and
$#kx$ slices with their fits) for every stage. These figures predate this
session's pair/sign corrections and were not regenerated with them; since
they use different profile pairs from @tab:fielddata-stages entirely (not
just a different fit configuration, as with @sec:hyp3-fd-roi-picking above),
the two sets of numbers are not directly comparable -- they are a
consistency check on *direction and relative magnitude*, not a replacement
for the table.
Kirchhoff-BP and Gazdag agree closely on both the sign and the relative size
of the estimate in every stage (largest during Push, near-zero during Wait),
which is expected since both operate on the same underlying B-scan data.
Back-propagation's estimates are noticeably noisier -- visible in the more
scattered 1-D slices -- and disagree with Kirchhoff/Gazdag on the sign of
$#Dz$; this is addressed as an open question in @sec:hyp3-fd-interpretation.

#figure(
  cimg("FD_displacement_diagnostics_gazdag.png"),
  caption: [Gazdag: WLS cross-spectrum phase-plane diagnostics for the four
    stage-boundary pairs (rows) -- cross-spectrum phase, cross-spectrum
    energy with WLS threshold contour, fitted plane, 1-D $#kz$ slice, 1-D
    $#kx$ slice (columns).],
) <fig:fd-disp-gazdag>

#supp-note[The equivalent Kirchhoff-BP and back-propagation five-panel
diagnostics, same panel layout as @fig:fd-disp-gazdag, are provided in the
Supplementary Material, §S4.4.1 and §S4.4.2.]

=== Back-Propagation Cross-Check and Stage Displacements <sec:hyp3-fd-bp-corrected>

The painted-ROI workflow of @sec:hyp3-fd-roi-picking is applied to the five
corrected, final-processed back-propagation focus frames, for the same four
representative consecutive pairs (1→3, 3→8, 8→20, 20→38), painting on the
difference image so a Gaussian-softened mask can trace the reflector's shape
(including, in two pairs, a side-lobe the rectangular ROI misses).

One back-propagation-specific correction is required: its focus frames' depth
axis increases with row index, opposite to the Kirchhoff/Gazdag arrays used
elsewhere in this chapter. Since the fit measures displacement along increasing
row index, the $k_z$-axis row spacing is negated for the back-propagation
frames only; the raw fit output then denotes the same physical direction under
both conventions, and the values quoted are additionally converted to the
"positive $#Dz$ = downward" convention of @tab:fielddata-stages.

@tab:fielddata-stages's BP columns report the corrected, sign-converted
estimates. Direction now agrees with the Gazdag columns throughout:
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
paused-injection stage already established qualitatively in the
sliding-window maps (Supplementary Material, §S4.3.2). Two independently-processed techniques agreeing on
both sign and order of magnitude in every stage is a stronger check than
either alone, but it also means the reversed direction is very unlikely to
be a processing artefact specific to one pipeline -- @sec:hyp3-fd-interpretation
returns to what it might mean instead.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.7em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Stage*], [*Pair*], [*$#Dz$ (Gazdag)*], [*$#Dx$ (Gazdag)*],
    [*$#Dz$ (BP, corrected)*], [*$#Dx$ (BP, corrected)*],
    table.hline(stroke: 0.4pt),
    [Push],  [1→3],   [$-1.32 "m"$], [$-0.26 "m"$], [$-1.91 "m"$], [$-0.34 "m"$],
    [Chase], [3→8],   [$-0.76 "m"$], [$-0.06 "m"$], [$-1.37 "m"$], [$-0.21 "m"$],
    [Wait],  [8→20],  [$+0.30 "m"$], [$+0.03 "m"$], [$+0.54 "m"$], [$+0.08 "m"$],
    [Pull],  [20→38], [$+0.28 "m"$], [$+0.15 "m"$], [$+0.75 "m"$], [$+0.14 "m"$],
    table.hline(stroke: 0.4pt),
    [*Net (prof 1→38)*], [], [$-1.50 "m"$], [$-0.15 "m"$], [--], [--],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Stage-anchored phase-plane displacement estimates from the
    borehole GPR field dataset, using the same four representative
    consecutive pairs throughout this chapter (1→3, 3→8, 8→20, 20→38): the
    Strategy-3b (Gazdag-migrated) estimate of @sec:hyp3-fd-phaseplane
    alongside the independently-processed, corrected back-propagation
    pipeline of this section (manually-painted ROI, after the timing,
    normalisation, and depth-axis-sign fixes of @sec:hyp3-fd-crossprofile).
    Positive $#Dz$ is downward (along the borehole, increasing depth);
    positive $#Dx$ is away from the borehole. Net values are the
    reconstructed global trajectory after stitching the four Gazdag stage
    estimates -- a plain cumulative sum, since consecutive pairs need no
    separate bridging step; back-propagation has no corresponding net row
    since it was re-estimated only for the four representative pairs shown.],
  kind: table,
) <tab:fielddata-stages>

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
    @tab:fielddata-stages's BP columns for the sign-converted values used in
    the text.],
) <fig:fd-borehole-napari-pairs>

=== WLS vs RANSAC <sec:hyp3-fd-ransac>

The amplitude weighting of the WLS fit (@sec:th-wls) down-weights low-energy
bins but still lets every masked cell contribute, so a coherent band of
phase-wrapped or noise-dominated cells inside the mask can still bias the
plane. As a robustness check, the same fit is repeated with a RANSAC (random
sample consensus) estimator: it fits the plane to random cell subsets, keeps
the largest consensus set of inliers (phase residual $< 0.35 "rad"$, over
$2000$ iterations), and refits WLS on those inliers alone. The comparison is
run for all four representative pairs and all three migration techniques, in
both ROI picking domains of @sec:hyp3-fd-roi-picking.

When the $(k_z, k_x)$ cells are painted by hand (@fig:fd-ransac-kspace), RANSAC
flags no outliers -- every painted cell is already a consensus inlier -- and so
reproduces the WLS estimate exactly in all twelve cases: a tight hand-pick
around the coherent lobes needs no further robustification. When the cells are
instead populated by the automatic band-plus-amplitude gate over a hand-painted
_spatial_ ROI (@fig:fd-ransac-bscan), RANSAC rejects anywhere from a negligible
fraction up to about half of the gated cells (e.g. $53%$ for back-propagation,
Pull) and shifts the estimate accordingly -- usually by only a few percent, but
by up to $23%$ in the worst case (Kirchhoff-BP, Chase: $#Dz = -1.07 "m"$ under
WLS versus $-0.83 "m"$ under RANSAC). Direction and stage-to-stage ordering are
preserved in every case. The practical reading is that the fitting _method_
matters only once the cross-spectrum cell population is left to an automatic
gate; a curated k-space pick makes WLS and RANSAC interchangeable.

#figure(
  cimg("FD_ransac_vs_wls_displacement_summary.png"),
  caption: [WLS versus RANSAC displacement estimates for hand-painted
    $(k_z, k_x)$ picking, all four representative pairs and three migration
    techniques (hatched = WLS, solid = RANSAC). Every painted cell is a RANSAC
    inlier, so the two estimates coincide.],
) <fig:fd-ransac-kspace>

#figure(
  cimg("FD_ransac_vs_wls_bscan_displacement_summary.png"),
  caption: [WLS versus RANSAC displacement estimates for amplitude-domain
    (difference B-scan) picking, where an automatic gate selects the
    $(k_z, k_x)$ cells inside a hand-painted spatial ROI. RANSAC rejects up to
    roughly half the gated cells; the estimate shifts by up to $23%$
    (Kirchhoff-BP, Chase) but never changes sign or stage ordering.],
) <fig:fd-ransac-bscan>

#supp-note[The per-pick WLS/RANSAC inlier--outlier phase panels for every stage
and technique, in both picking domains, are provided in the Supplementary
Material, §S4.4.3.]

== Interpretation <sec:hyp3-fd-interpretation>

*On the back-propagation/Kirchhoff-Gazdag sign disagreement.* The
14-profile, homogeneous-domain back-propagation result (Supplementary
Material, §S4.4.2) disagreeing with Kirchhoff/Gazdag on the sign of
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
14-profile homogeneous-domain figures (Supplementary Material, §S4.4.2), which predate
all three fixes of @sec:hyp3-fd-crossprofile and @sec:hyp3-fd-bp-corrected
and have not been re-run under them; the SNR and technique-reliability
explanations remain open for that specific figure.

#draftnote[fill in the remaining physical interpretation once the borehole
geometry, injection depth, and fluid-injection parameters are confirmed from
the field survey metadata. Key questions to address: (1) do the inferred
$#Dz$ and $#Dx$ values agree with the known injection depth and the expected
lateral spread for the given fracture geometry -- and, now that
@tab:fielddata-stages's Gazdag and BP columns have both been corrected to
a verified sign convention (@sec:hyp3-fd-phaseplane), they agree that Push
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
