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
pre-labelled region of interest --- how that region is chosen.

#linebreak()

#para-head[Hypothesis 3.] The time-lapse phase-plane approach generalises
beyond idealised single-scatterer synthetic models to real borehole GPR field
data.

== Field Data Explanation <sec:hyp3-fielddata>

The dataset consists of 38 zero-offset borehole GPR profiles acquired on
6 June 2016, from the same push-pull, single-hole GPR tracer-test field
campaign in fractured rock described in @shakas2016, with the following
parameters: medium propagation velocity
$v = 0.10 "m/ns"$, centre frequency $f_0 = 0.10 "GHz"$, trace spacing
$d_L = 0.05 "m"$, and maximum depth $85 "m"$. The experiment was conducted in
four operational stages, each with a distinct expected fluid behaviour:

/ Push (profiles 1--4): fluid is actively injected at $78.7 "m"$ depth; a large upward displacement and lateral spread are expected.
/ Chase (profiles 5--9): continued injection at a lower rate; smaller incremental advance expected.
/ Wait (profiles 10--20): injection paused; near-zero net displacement expected.
/ Pull (profiles 21--38): fluid is extracted; partial reversal of the push displacement expected.

Each profile is a raw depth/two-way-travel-time B-scan along the borehole,
with no migration, differencing or gain correction applied yet --- that is
the subject of @sec:hyp3-fd-processing.

Back-propagation migration (@sec:hyp3-fd-processing) additionally requires an
explicit model of the borehole itself: a real single-hole survey has both
antennas inside a fluid-filled channel whose relative permittivity ($approx 81$ for
water) differs sharply from the surrounding rock/ice ($approx 9$).
@fig:fd-borehole-schematic shows this domain: a $10 "cm"$-wide, water-filled
borehole spanning the full modelled depth range, with the source/receiver on
its lateral centreline, background medium to its left (radial direction), and
imaging buffer to its right, at the fine grid spacing ($d_x = 0.02 "m"$)
needed to resolve the borehole without causing numerical dispersion.

#figure(
  cimg("FD_borehole_domain_schematic.png", width: 110%),
  caption: [Borehole domain used in the back-propagation migration: the
    explicit water-filled channel (permittivity discontinuity against the
    background medium), source/receiver centreline, and imaging aperture. Grid
    spacing and permittivity scaling are detailed in @sec:hyp3-fd-processing.],
) <fig:fd-borehole-schematic>

== Processing and Migration <sec:hyp3-fd-processing>

Every profile passes through a shared five-step conditioning chain --
bandpass filter ($0.02$--$0.20 "GHz"$), DC and direct-wave removal, trace
alignment (×5 upsampling), SVD rank-1 direct-wave suppression, and time-lapse
differencing against the reference profile -- after which Kirchhoff, Gazdag
and back-propagation migration diverge, and then reconverge: *every*
migrated image receives the same post-imaging
clean-up before it is used downstream. An amplitude gate
and near-borehole taper (described in Post-imaging processing) are generic operations on a migrated image. Without post-imaging processing, back-propagation, and Kirchhoff and Gazdag images would carry the same class of migration-artefact clutter.

#linebreak()

*Migration.* Kirchhoff and Gazdag apply a linear spherical-gain correction directly to the differenced B-scan before migrating; by compensating for $1\/r$ geometric spreading, this step simultaneously acts as the requisite 3D-to-2D conversion. Kirchhoff migration itself --- "Kirchhoff-BP" throughout this chapter --- is built with a delta (spike) wavelet rather than the Ricker wavelet used elsewhere in this thesis: because the Pylops created Kirchhoff adjoint operator correlates the image against the wavelet a second time, a delta input avoids a redundant autocorrelation step, giving a sharper point-spread function that matches the one already produced by Gazdag's $t=0$ imaging condition, so the two techniques' native resolution stays directly comparable. Back-propagation instead applies a rigorous physics-based 3D-to-2D Green’s-function correction—an amplitude gain of $sqrt(r)$ together
with a $1\/sqrt(omega) dot e^(i phi)$ frequency-domain phase filter. This ensures dimensional consistency, because while the survey records a 3D field, gprMax back-propagates it through a 2D grid, and the Green's function governing the dipole response of a 3D medium is fundamentally different from that of a 2D medium. The corrected B-scan is then spatially and temporally Tukey-window-tapered, passed through a $f$-$k_z$ dip filter to remove evanescent energy, time-reversed and peak-normalised, low-pass filtered to gprMax’s numerical-dispersion limit, and edge-zeroed before injection into the explicit borehole geometry of @fig:fd-borehole-schematic ($d_x = 0.02 "m"$ and
$times 4$ permittivity scaling is needed to agree with the exploding-reflector
convention $#vmig = v\/2$. This is applied consistently with all three migration techniques for both the background and the borehole water). Back-propagation
additionally requires suppressing *injection artefacts*: every gprMax source in
the time-reversal simulation is injected along the borehole wall rather than
at a single point, so the injected wavefields only destructively interfere
and cancel out well away from the true reflector; close to the borehole this
cancellation is incomplete, leaving residual energy that the post-imaging
taper below is tuned to remove.

#linebreak()

*Post-imaging processing.* A final round of processing
now acts uniformly on every completed migrated image: an amplitude gate built
from the image's own local energy envelope, and, for back-propagation
specifically, a cosine-ramped near-borehole taper suppressing the injection
halo out to $3.0 "m"$ radial distance. The recipe is recomputed
independently per profile and per technique, and acts purely on the migrated
output, not on the underlying B-scan or gprMax injection file.

=== Cross-Profile Registration and Amplitude Normalisation <sec:hyp3-fd-crossprofile>

Back-propagation's focus frames must additionally be registered across
profiles before differencing: a fixed snapshot-index, calibrated once
against profile 1, was chosen by hand so that the reflection's depth and
radial distance matched its position in the Kirchhoff and Gazdag images as
closely as possible. A single peak amplitude, found across all five profiles'
excitation data combined, is used to normalise every profile's traces, in
place of each profile normalising to its own peak. With both fixes
applied, time-lapse differences (e.g. profile 8 minus 3) across the five-profile set (1, 3, 8, 20, 38) are
compact and localised, matching the Kirchhoff/Gazdag character, rather than
smeared across the reflector as occurred before the fix.

Kirchhoff, Gazdag and back-propagation migrations are all available for the
five representative profiles used throughout this chapter (1, 3, 8, 20, 38),
back-propagation under the explicit-borehole, $d_x = 0.02 "m"$ pipeline
described above.

#linebreak()

@fig:fd-profile-grid compares these five profiles across all four
representations: the processed B-scan, and the post-imaging-processed
Kirchhoff-BP, Gazdag, and back-propagation migrations. All three migration
techniques agree on a single dominant reflector at approximately
$70$--$80 "m"$ depth and $4$--$8 "m"$ radial distance, which sharpens
progressively from profile 1 to profile 8 and then remains essentially
stationary through profiles 20 and 38 -- a first visual indication of the
push/chase/wait/pull kinematics discussed in @sec:hyp3-fd-phaseplane. The
back-propagation column retains more residual energy near the borehole than
Kirchhoff or Gazdag even after the post-imaging taper, consistent with the
injection artefacts being a structural property of single-sided time-reversal
rather than a processing shortcoming.

#figure(
  cimg("FD_profile_migration_grid.png", width: 110%),
  caption: [Processed B-scans and their migrated counterparts for profiles 1,
    3, 8, 20 and 38, all post-imaging-processed: processed B-scan,
    Kirchhoff-BP migration, Gazdag migration, back-propagation $E_z$ focus
    frame.],
) <fig:fd-profile-grid>

== Timelapse Differencing and Region of Influence Selection <sec:hyp3-fd-roi>

Every result in this chapter is computed on one of four representative
consecutive profile pairs -- 1→3, 3→8, 8→20, 20→38 -- chosen to approximate
the Push, Chase, Wait and Pull stages respectively while spanning the full
38-profile survey with only four differenced images. For a given pair, the
time-lapse difference is simply the later migrated (and post-imaging
processed) image minus the earlier one.

The phase-plane estimator of @sec:meth-phaseplane fits a plane to the
cross-spectrum phase of the two images, but only inside a region of
influence (ROI): a subset of $(k_z, k_x)$ cells, or equivalently a spatial
window before the FFT, chosen to contain the reflector's coherent energy and
exclude noise. ROI selection in this chapter is done exclusively by hand,
using the #link("https://napari.org")[napari] interactive image viewer,
directly on top of the data -- not by an automatic rectangular or
sliding-window search. Two variants are used together: the $(k_z, k_x)$
cells can be painted directly on the cross-spectrum phase or amplitude
display, or the spatial region can be painted on the difference B-scan
itself, letting the ROI trace the reflector's actual shape rather than having a rectangular bounding box around the main reflection.

#linebreak()

Because the two migrated images being compared are real-valued, their cross-spectrum is exactly Hermitian, $#XS (-#kz, -#kx) = #XS (#kz, #kx)^*$, so every genuine $(#kz, #kx)$ point has a phase-negated mirror duplicate at $(-#kz, -#kx)$. While fitting a plane through both halves simultaneously does not bias an ideal synthetic dataset, field data inherently contains physical non-linearities—such as wavelet dispersion or deviations from a perfectly straight borehole trajectory. These real-world effects can introduce subtle asymmetries during phase unwrapping, causing a full-spectrum plane fit to become biased. To avoid this, the fit must be restricted to one half of the spectrum. For this borehole/VRP geometry, the antenna radiates along the radial
direction $x$ as it is lowered down the borehole, so $x$ --- not $z$ --- is
the carrier/wavelet axis, and it is this mirror duplicate that sits at
$-#kx$. The $#kz$ axis is simply the spatial direction the antenna traverses
along the borehole, and a real reflector's diffraction response has genuine,
independent structure on both its $+#kz$ (antenna receding) and $-#kz$
(antenna approaching) flanks --- restricting $#kz$ instead would discard half
of every real target's structure while leaving its Hermitian mirror
untouched. Both pixel picking variants are therefore restricted to
$#kx >= 0$, with $#kz$ left unrestricted: the wavenumber-domain display only
shows, and only allows pixel picking on, $#kx >= 0$, and the amplitude-domain
automatic gate that resolves the $(#kz, #kx)$ population behind a
spatially-painted B-scan region (@fig:fd-wls-workflow) applies the same restriction.

#linebreak()

Picking the cross-spectrum directly restricts the fit to exactly the
coherent energy lobes visible in the display, with no relative amplitude
threshold or fixed $(k_z, k_x)$ band needed (@fig:fd-wls-workflow). Picking on the difference
B-scan instead lets the ROI follow the reflection's shape. Because a
hand-painted mask has hard edges that would otherwise leak spectral energy,
the painted mask is Gaussian-softened before being used as the spatial
window for the cross-spectrum.

== WLS Cross-Spectrum Phase Plane Fit <sec:hyp3-fd-phaseplane>

For each representative pair, a plane is fit by weighted least squares (WLS)
to the cross-spectrum phase inside the pixel-painted ROI of
@sec:hyp3-fd-roi, following the pipeline of @sec:meth-phaseplane: crop, taper,
FFT, form the cross-spectrum, weight every retained cell by its amplitude,
and solve for $(#Dz, #Dx, c)$. The four pairs are chained into a single
global trajectory across all profiles by a stage-anchored construction, where each stage represents the Push, Chase, Wait, and Pull phase respectively. The global displacement
at each stage boundary is simply the cumulative sum of the four independent
per-pair (e.g. between profile 1 and 3) estimates.

@fig:fd-wls-workflow shows this workflow for the Push-stage pair (profiles 1→3,
Gazdag): the two migrated images being compared, and the two pixel picks --
one on the amplitude-domain difference B-scan, one on the wavenumber-domain
cross-spectrum -- that together define the ROI for the WLS fit. The resulting $(#Dz, #Dx)$ estimate is then chained with the other three representative pairs to form a single global displacement trajectory across all 38 profiles.

#figure(
  cimg("FD_wls_workflow_1_to_3.png"),
  caption: [WLS cross-spectrum phase-plane workflow, Push-stage pair (1→3),
    Gazdag migration. Top row: the two migrated images being compared
    (profiles 1 and 3). Bottom row: the pixel-painted ROI on the
    amplitude-domain difference B-scan (left) and on the wavenumber-domain
    cross-spectrum (right).],
) <fig:fd-wls-workflow>

=== Back-Propagation Cross-Check and Stage Displacements <sec:hyp3-fd-bp-corrected>

The same pixel-picked workflow is applied to the back-propagation
focus frames for the same four representative pairs. Throughout this chapter, positive $#Dz$ corresponds to
downward (increasing depth) movement, and negative $#Dz$ to upward
(decreasing depth) movement, whereas a positive $#Dx$ corresponds to movement away from the borehole, and negative $#Dx$ to movement toward it.

#linebreak()

@tab:fielddata-stages summarises the resulting stage-by-stage displacement
estimates for Gazdag and back-propagation, with the corrected $#kx > 0$ Hermitian-mirror restriction
(@sec:hyp3-fd-roi) in place. The dominant signal is an *upward* displacement
of approximately $1.66 "m"$ during the Push stage, accompanied by a smaller
lateral shift of $0.17 "m"$ *toward* the borehole. The Chase stage adds a
comparable increment in the same (upward) direction, and the Wait stage
substantially reverses it ($+0.68 "m"$, i.e. downward); both methods continue a
smaller reversal during Pull ($+0.14 "m"$ for Gazdag, $+0.24 "m"$ for
back-propagation), leaving a net residual of $#Dz approx -2.44 "m"$
(upward), $#Dx approx +0.22 "m"$ (away from
the borehole) at profile 38 relative to profile 1 for Gazdag. Back-propagation
agrees with Gazdag on direction in every stage -- upward and
toward the borehole during Push, upward and away during Chase, downward
during Wait and Pull, despite
the results coming from two different migration techniques. @sec:hyp3-fd-interpretation
returns to what this displacement pattern means physically.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.7em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Stage*], [*Pair*], [*$#Dz$ (Gazdag)*], [*$#Dx$ (Gazdag)*],
    [*$#Dz$ (BP, corrected)*], [*$#Dx$ (BP, corrected)*],
    table.hline(stroke: 0.4pt),
    [Push],  [1→3],   [$-1.66 "m"$], [$-0.17 "m"$], [$-2.06 "m"$], [$-0.25 "m"$],
    [Chase], [3→8],   [$-1.60 "m"$], [$+0.07 "m"$], [$-1.65 "m"$], [$+0.06 "m"$],
    [Wait],  [8→20],  [$+0.68 "m"$], [$+0.01 "m"$], [$+0.28 "m"$], [$+0.03 "m"$],
    [Pull],  [20→38], [$+0.14 "m"$], [$+0.31 "m"$], [$+0.24 "m"$], [$-0.14 "m"$],
    table.hline(stroke: 0.4pt),
    [*Net (prof 1→38)*], [], [$-2.44 "m"$], [$+0.22 "m"$], [$-3.19 "m"$], [$-0.30 "m"$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Stage-anchored phase-plane displacement estimates from the
    borehole GPR field dataset, for the four representative pairs (1→3,
    3→8, 8→20, 20→38): the Gazdag-migrated estimate alongside the
    independently-processed, corrected back-propagation cross-check
    (pixel-painted ROI). Positive $#Dz$ is downward (along the borehole,
    increasing depth); positive $#Dx$ is away from the borehole. Net values
    are the cumulative sum of the four Gazdag and back-propagation stage estimates.],
  kind: table,
) <tab:fielddata-stages>

== Interpretation (WLS vs. RANSAC) <sec:hyp3-fd-interpretation>

Back-propagation and Gazdag agree on the sign of the displacement in every
stage (@tab:fielddata-stages), despite sharing no processing steps
downstream of the raw B-scans. This cross-method agreement is evidence that the recovered $#Dz$
pattern -- upward during Push and Chase, reversing downward during Wait and
continuing through Pull -- is
a genuine feature of this dataset rather than a processing or
sign-convention artefact specific to one migration technique. The dominant
vertical displacement is consistent with the injection point itself: fluid is
injected at $78.7 "m"$ depth, within the same $70$--$80 "m"$ band as the
tracked reflector, so an upward-moving front while fluid is actively pumped
in (Push, Chase) and a downward settling once pumping stops (Wait, Pull) is
the physically expected pattern. The lateral component is less
straightforward, moving *toward* the borehole during Push before reversing
*away* from it in the later stages. This could be explained by the orientation of the fracture network, being not perfectly parallel to the borehole.

#linebreak()

*WLS vs. RANSAC.* The amplitude weighting of the WLS fit (@sec:th-wls)
down-weights low-energy bins but still lets every masked cell contribute, so
a coherent band of phase-wrapped or noise-dominated cells inside a pixel-painted
mask could in principle still bias the plane. As a robustness check, each
fit is repeated with RANSAC #cite(<zhu2025>, form: "prose"): it fits the plane to random cell subsets, keeps
the largest consensus set of inliers (phase residual $< 0.35 "rad"$, over
$2000$ iterations), and refits WLS on those inliers alone. The comparison is
run for all four representative pairs and all three migration techniques, in
both the wavenumber-domain and amplitude-domain pixel-painted picking variants of
@sec:hyp3-fd-roi. @fig:fd-ransac-phase-slices shows a representative example
(Gazdag, Chase-stage pair): the cross-spectrum phase and amplitude with the
fit population highlighted, and the $k_x$-/$k_z$-direction 1-D fits, for WLS
and RANSAC side by side.

#linebreak()

@fig:fd-ransac-summary summarises the comparison across all four
representative pairs and three migration techniques. When the $(k_z, k_x)$
cells are painted by hand directly on the
cross-spectrum (@fig:fd-ransac-summary (a)), RANSAC flags few or no outliers and closely reproduces the
WLS estimate in every case. When the
cells are instead populated by gating a hand-painted amplitude-domain
(difference B-scan) region (@fig:fd-ransac-summary (b)), RANSAC rejects anywhere from a negligible
fraction up to about half of the gated cells and shifts the estimate
accordingly. Direction is preserved in every case, and stage-to-stage ordering in most.
The takeaway is the same as for a tight k-space pick: which fitting *method*
is used matters far more when the cross-spectrum cell population comes from
an amplitude-threshold gate over a painted B-scan region than when it comes
from a direct, hand-picked selection in the wavenumber domain.

#figure(
  cimg("FD_ransac_vs_wls_phase_amp_fit_gazdag_chasing_3_to_8.png", width: 120%),
  caption: [WLS (top row) versus RANSAC (bottom row) phase-plane fit,
    Chase-stage pair (3→8), Gazdag migration. Columns: cross-spectrum phase,
    cross-spectrum amplitude, $x$-direction ($#kx$) fit, $z$-direction
    ($#kz$) fit. ROI within the wavenumber band and amplitude threshold is
    highlighted on the phase and amplitude panels.],
) <fig:fd-ransac-phase-slices>

#figure(
  subfigs(cols: 1,
    cimg("FD_ransac_vs_wls_displacement_summary.png", width: 120%),
    cimg("FD_ransac_vs_wls_bscan_displacement_summary.png", width: 120%),
  ),
  caption: [WLS versus RANSAC displacement estimates, all four representative
    pairs and three migration techniques (hatched = WLS, solid = RANSAC): (a)
    wavenumber-domain (hand-painted $(k_z,k_x)$) picking, where nearly every
    painted cell is a RANSAC inlier so the two estimates closely coincide; (b)
    amplitude-domain (difference B-scan) picking, where an automatic gate
    populates the $(k_z,k_x)$ cells inside the painted spatial ROI and RANSAC
    rejects a larger fraction of them.],
) <fig:fd-ransac-summary>

#linebreak()

#supp-note[A representative selection of the per-pick inlier/outlier phase
panels underlying the two aggregate summaries above -- spanning both picking
domains and the full range of WLS/RANSAC agreement seen, from near-identical
estimates to the largest disagreements -- is provided in the Supplementary
Material, §S4.1.]
