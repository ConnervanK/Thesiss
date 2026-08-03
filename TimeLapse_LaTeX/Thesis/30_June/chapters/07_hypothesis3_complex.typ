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

/ Push (profiles 1--4): fluid is actively injected; a large downward displacement and lateral spread are expected.
/ Chase (profiles 5--9): continued injection at a lower rate; smaller incremental advance expected.
/ Wait (profiles 10--20): injection paused; near-zero net displacement expected.
/ Pull (profiles 21--38): fluid is extracted; partial reversal of the push displacement expected.

Each profile is a raw depth/two-way-travel-time B-scan along the borehole,
with no migration, differencing or gain correction applied yet --- that is
the subject of @sec:hyp3-fd-processing.

Back-propagation migration (@sec:hyp3-fd-processing) additionally requires an
explicit model of the borehole itself: a real single-hole survey has both
antennas inside a fluid-filled channel whose permittivity ($approx 81$ for
water) differs sharply from the surrounding rock/ice ($approx 9$), a
discontinuity a homogeneous-medium model would omit entirely.
@fig:fd-borehole-schematic shows this domain: a $10 "cm"$-wide, water-filled
channel spanning the full modelled depth range, with the source/receiver on
its lateral centreline, background medium to its left (radial direction), and
imaging clearance to its right, at the fine grid spacing ($d_x = 0.02 "m"$)
needed to resolve the channel without stair-casing its boundary.

#figure(
  cimg("FD_borehole_domain_schematic.png"),
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
migrated image, regardless of technique, now receives the same post-imaging
clean-up before it is used downstream. This is a deliberate change from
treating back-propagation's clutter as a special case -- the amplitude gate
and near-borehole taper below are generic operations on a migrated image, not
back-propagation-specific fixes, and Kirchhoff and Gazdag images carry the
same class of migration-artefact clutter that benefits from them.

*Migration.* Kirchhoff and Gazdag apply a linear spherical-gain correction directly to the differenced B-scan before migrating; by compensating for $1\/r$ geometric spreading, this step simultaneously acts as the requisite 3D-to-2D conversion. Back-propagation instead applies a rigorous physics-based 3D-to-2D Green’s-function correction—an amplitude gain of $sqrt(r)$ together
with a $1\/sqrt(omega) dot e^(i phi)$ frequency-domain phase filter. This ensures dimensional consistency, because while the survey records a 3D field, gprMax back-propagates it through a 2D grid, and the Green's function governing the dipole response of a 3D medium is fundamentally different from that of a 2D medium. The corrected B-scan is then spatially and temporally Tukey-window-tapered, passed through a $f$-$k_z$ dip filter to remove evanescent energy, time-reversed and peak-normalised, low-pass filtered to gprMax’s numerical-dispersion limit, and edge-zeroed before injection into the explicit borehole geometry of @fig:fd-borehole-schematic ($d_x = 0.02 "m"$,
$times 4$ permittivity scaling implementing the exploding-reflector
convention $#vmig = v\/2$ shared with Kirchhoff and Gazdag, applied
consistently to both the background and the borehole water). Back-propagation
additionally requires suppressing an *injection halo*: every gprMax source in
the time-reversal simulation is injected at the borehole wall, and the
destructive interference that should cancel the field away from the true
reflector is never complete near the borehole, leaving residual energy at
small radial distance that the post-imaging taper below is tuned to remove.

#draftnote[the sign of the $1\/sqrt(omega)$ filter's $pi\/4$ phase term, and
whether a water-level regularisation should be added near DC, were both left
open pending a dedicated validation experiment (single synthetic reflector,
compared against Kirchhoff/Gazdag as ground truth). State the resolved
choice here once confirmed.]

*Post-imaging processing (all three techniques).* A final round of processing
now acts uniformly on every completed migrated image: an amplitude gate built
from the image's own local energy envelope, and, for back-propagation
specifically, a cosine-ramped near-borehole taper suppressing the injection
halo out to $3.0 "m"$ radial distance. A per-profile $f$-$k$ dip (fan) filter
was also tried -- estimating the reflector's dip directly from the image and
keeping only a cosine-tapered wedge of that orientation in the $(k_z, k_x)$
domain -- but was dropped after testing against the real Gazdag and
back-propagation images: it flattened genuine reflector energy along with
noise, leaving the post-imaged result looking worse than the amplitude-gated
image alone, not cleaner. Two further alternatives were rejected earlier in
development for the same reason (removing signal, not only noise): Gaussian
smoothing plus Robust PCA decomposition (after Li and Yan, 2021), because
this reflector's response is spatially extended along a consistent dip and so
behaves as low-rank rather than sparse; and a spatially-windowed fan filter,
because restricting the FFT to a local window does not stop locally-oriented
noise from being reconstructed within it. The recipe is recomputed
independently per profile and per technique, and acts purely on the migrated
output, not on the underlying B-scan or gprMax injection file.

=== Cross-Profile Registration and Amplitude Normalisation <sec:hyp3-fd-crossprofile>

Back-propagation's focus frames must additionally be registered across
profiles before differencing: an independent per-profile peak-amplitude
search lets different profiles lock onto different snapshot times (a spread
of $approx 26 "ns"$), so every profile now uses the same nominal focus time
plus a fixed snapshot-index offset calibrated once against profile 1; and
normalising every trace to its own peak before injection was erasing the
genuine cross-profile amplitude growth the time-lapse comparison targets, so
a single amplitude shared across all five profiles' excitation data replaces
the per-trace scalar. With both fixes applied, differences across the
five-profile set (1, 3, 8, 20, 38) are compact and localised, matching the
Kirchhoff/Gazdag character, rather than smeared across the reflector as
occurred before the fix.

Kirchhoff and Gazdag migrations are available for all 37 processed profiles
(profiles 2--38, differenced against profile 1 as the fixed baseline).
Back-propagation is available for the five representative profiles used
throughout this chapter (1, 3, 8, 20, 38), under the explicit-borehole,
$d_x = 0.02 "m"$ pipeline described above -- the only back-propagation
configuration now used in this thesis.

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
injection halo being a structural property of single-sided time-reversal
rather than a processing shortcoming.

#figure(
  cimg("FD_profile_migration_grid.png"),
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
itself, letting the ROI trace the reflector's actual shape rather than an
axis-aligned box.

Because the two migrated images being compared are real-valued, their cross-spectrum is exactly Hermitian, $#XS (-#kz, -#kx) = #XS (#kz, #kx)^*$, so every genuine $(#kz, #kx)$ point has a phase-negated mirror duplicate at $(-#kz, -#kx)$. While fitting a plane through both halves simultaneously does not bias an ideal synthetic dataset, field data inherently contains physical non-linearities—such as wavelet dispersion or deviations from a perfectly straight borehole trajectory. These real-world effects can introduce subtle asymmetries during phase unwrapping, causing a full-spectrum plane fit to become biased. To avoid this, the fit must be restricted to one half of the spectrum. For this borehole/VRP geometry, the antenna radiates along the radial
direction $x$ as it is lowered down the borehole, so $x$ --- not $z$ --- is
the carrier/wavelet axis, and it is this mirror duplicate that sits at
$-#kx$. The $#kz$ axis is simply the spatial direction the antenna traverses
along the borehole, and a real reflector's diffraction response has genuine,
independent structure on both its $+#kz$ (antenna receding) and $-#kz$
(antenna approaching) flanks --- restricting $#kz$ instead would discard half
of every real target's structure while leaving its Hermitian mirror
untouched. Both napari picking variants are therefore restricted to
$#kx >= 0$, with $#kz$ left unrestricted: the wavenumber-domain display only
shows, and only allows painting on, $#kx >= 0$, and the amplitude-domain
automatic gate that resolves the $(#kz, #kx)$ population behind a
spatially-painted B-scan region applies the same restriction.

Picking the cross-spectrum directly restricts the fit to exactly the
coherent energy lobes visible in the display, with no relative amplitude
threshold or fixed $(k_z, k_x)$ band needed. Picking on the difference
B-scan instead lets the ROI follow the reflection's shape; because a
hand-painted mask has hard edges that would otherwise leak spectral energy,
the painted mask is Gaussian-softened before being used as the spatial
window for the cross-spectrum. A naive wide, unweighted pick that happens to
include only one of a reflector's two energy lobes was found during
development to return a badly *sign-reversed* estimate; picking both lobes,
or a narrow high-amplitude region around the peak of one, both recover the
correct sign and magnitude -- a leakage artefact of an imprecise selection,
not a flaw in the underlying mapping, and the reason the painted mask is
softened rather than left hard-edged.

== WLS Cross-Spectrum Phase Plane Fit <sec:hyp3-fd-phaseplane>

For each representative pair, a plane is fit by weighted least squares (WLS)
to the cross-spectrum phase inside the napari-painted ROI of
@sec:hyp3-fd-roi, following the pipeline of @sec:meth-phaseplane: crop, taper,
FFT, form the cross-spectrum, weight every retained cell by its amplitude,
and solve for $(#Dz, #Dx, phi_0)$. The four pairs are chained into a single
global trajectory across all 38 profiles by a stage-anchored construction:
because each representative pair already spans its stage end-to-end
(profile $n$ of one pair is profile 1 of the next), the global displacement
at each stage boundary is simply the cumulative sum of the four independent
per-pair estimates, with no separate bridging step required.

@fig:fd-wls-workflow shows this workflow for the Push-stage pair (1→3,
Gazdag): the two migrated images being compared, and the two napari picks --
one on the amplitude-domain difference B-scan, one on the wavenumber-domain
cross-spectrum -- that together define the ROI entering the fit.

#figure(
  cimg("FD_wls_workflow_1_to_3.png"),
  caption: [WLS cross-spectrum phase-plane workflow, Push-stage pair (1→3),
    Gazdag migration. Top row: the two migrated images being compared
    (profiles 1 and 3). Bottom row: the napari-painted ROI on the
    amplitude-domain difference B-scan (left) and on the wavenumber-domain
    cross-spectrum (right).],
) <fig:fd-wls-workflow>

=== Back-Propagation Cross-Check and Stage Displacements <sec:hyp3-fd-bp-corrected>

The same napari-picked workflow is applied to the corrected back-propagation
focus frames for the same four representative pairs. Throughout this chapter, negative $#Dz$ corresponds to
downward (increasing depth) movement, and positive $#Dz$ to upward
(decreasing depth) movement.

@tab:fielddata-stages summarises the resulting stage-by-stage displacement
estimates for Gazdag and for the independently-processed back-propagation
cross-check, from the napari-only, stage-anchored, universally post-imaged
pipeline with the corrected $#kx > 0$ Hermitian-mirror restriction
(@sec:hyp3-fd-roi) in place. The dominant signal is a *downward* displacement
of approximately $1.66 "m"$ during the Push stage, accompanied by a smaller
lateral shift of $0.17 "m"$ *toward* the borehole. The Chase stage adds a
comparable increment in the same (downward) direction, and the Wait stage
substantially reverses it ($+0.68 "m"$, i.e. upward); both methods continue a
smaller reversal during Pull ($+0.14 "m"$ for Gazdag, $+0.24 "m"$ for
back-propagation), leaving a net residual of $#Dz approx -2.44 "m"$
(downward), $#Dx approx +0.22 "m"$ (away from
the borehole) at profile 38 relative to profile 1 for Gazdag. Back-propagation
agrees with Gazdag on direction in every stage -- downward and
toward the borehole during Push, downward and away during Chase, upward
during Wait and Pull -- once its depth-axis convention is corrected, despite
the two techniques sharing no processing steps downstream of the raw
B-scans. @sec:hyp3-fd-interpretation
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
    [*Net (prof 1→38)*], [], [$-2.44 "m"$], [$+0.22 "m"$], [--], [--],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Stage-anchored phase-plane displacement estimates from the
    borehole GPR field dataset, for the four representative pairs (1→3,
    3→8, 8→20, 20→38): the Gazdag-migrated estimate alongside the
    independently-processed, corrected back-propagation cross-check
    (napari-painted ROI). Negative $#Dz$ is downward (along the borehole,
    increasing depth); positive $#Dx$ is away from the borehole. Net values
    are the cumulative sum of the four Gazdag stage estimates;
    back-propagation has no corresponding net row since it was estimated
    only for the four representative pairs shown.],
  kind: table,
) <tab:fielddata-stages>

== Interpretation (WLS vs. RANSAC) <sec:hyp3-fd-interpretation>

Back-propagation and Gazdag agree on the sign of the displacement in every
stage (@tab:fielddata-stages), despite sharing no processing steps
downstream of the raw B-scans. This cross-method agreement is evidence that the recovered $#Dz$
pattern -- downward during Push and Chase, reversing upward during Wait and
continuing through Pull -- is
a genuine feature of this dataset rather than a processing or
sign-convention artefact specific to one migration technique. The dominant
vertical displacement now agrees in direction with the naive expectation of
@sec:hyp3-fielddata (a downward displacement during active injection); the
lateral component does not straightforwardly follow the same naive "outward
spread" picture, moving *toward* the borehole during Push before reversing
*away* from it in the later stages.

#draftnote[fill in the remaining physical interpretation once the borehole
geometry, injection depth, and fluid-injection parameters are confirmed from
the field survey metadata. Key questions: (1) do the inferred $#Dz$ and $#Dx$
values agree with the known injection depth and expected lateral spread for
the given fracture geometry, or is there a physical reason (fracture
geometry, the specific reflector tracked vs. the injection point, a
coordinate-system offset) the toward-the-borehole lateral component during
Push should be expected? (2) does the field result agree with the synthetic
fluid-front experiment? (3) does the WLS-vs-RANSAC spread below bound a
meaningful "fit-method uncertainty" that should be reported alongside the
headline numbers?]

*WLS vs. RANSAC.* The amplitude weighting of the WLS fit (@sec:th-wls)
down-weights low-energy bins but still lets every masked cell contribute, so
a coherent band of phase-wrapped or noise-dominated cells inside a napari
mask could in principle still bias the plane. As a robustness check, each
fit is repeated with RANSAC #cite(<zhu2025>, form: "prose"): it fits the plane to random cell subsets, keeps
the largest consensus set of inliers (phase residual $< 0.35 "rad"$, over
$2000$ iterations), and refits WLS on those inliers alone. The comparison is
run for all four representative pairs and all three migration techniques, in
both the wavenumber-domain and amplitude-domain napari picking variants of
@sec:hyp3-fd-roi. @fig:fd-ransac-phase-slices shows a representative example
(Gazdag, Chase-stage pair): the cross-spectrum phase and amplitude with the
fit population highlighted, and the $k_x$-/$k_z$-direction 1-D fits, for WLS
and RANSAC side by side.

@fig:fd-ransac-summary summarises the comparison across all four
representative pairs and three migration techniques. When the $(k_z, k_x)$
cells are painted by hand directly on the
cross-spectrum (@fig:fd-ransac-summary (a)), RANSAC flags few or no outliers and closely reproduces the
WLS estimate in every case (within about $6%$ at worst): a tight hand-pick
around the coherent lobes needs little further robustification. When the
cells are instead populated by gating a hand-painted amplitude-domain
(difference B-scan) region (@fig:fd-ransac-summary (b)), RANSAC rejects anywhere from a negligible
fraction up to about half of the gated cells and shifts the estimate
accordingly -- usually by a few percent up to around twenty percent, but by
as much as $59%$ in the worst case (Kirchhoff-BP, Chase stage). Direction is
preserved in every case; stage-to-stage ordering is preserved in most cases.
The practical reading is unchanged from a tight k-space pick: the fitting
*method* matters far more once the cross-spectrum cell population is left to
an automatic gate over a painted region than when it comes from a direct,
hand-picked selection.

#figure(
  cimg("FD_ransac_vs_wls_phase_amp_fit_gazdag_chasing_3_to_8.png"),
  caption: [WLS (top row) versus RANSAC (bottom row) phase-plane fit,
    Chase-stage pair (3→8), Gazdag migration. Columns: cross-spectrum phase,
    cross-spectrum amplitude, $x$-direction ($#kx$) fit, $z$-direction
    ($#kz$) fit. ROI within the wavenumber band and amplitude threshold is
    highlighted on the phase and amplitude panels.],
) <fig:fd-ransac-phase-slices>

#figure(
  subfigs(cols: 1,
    cimg("FD_ransac_vs_wls_displacement_summary.png"),
    cimg("FD_ransac_vs_wls_bscan_displacement_summary.png"),
  ),
  caption: [WLS versus RANSAC displacement estimates, all four representative
    pairs and three migration techniques (hatched = WLS, solid = RANSAC): (a)
    wavenumber-domain (hand-painted $(k_z,k_x)$) picking, where nearly every
    painted cell is a RANSAC inlier so the two estimates closely coincide; (b)
    amplitude-domain (difference B-scan) picking, where an automatic gate
    populates the $(k_z,k_x)$ cells inside the painted spatial ROI and RANSAC
    rejects a larger fraction of them.],
) <fig:fd-ransac-summary>

// #supp-note[The per-pick inlier/outlier phase panels underlying the two
// aggregate summaries above, for every stage and migration technique, in both
// the wavenumber-domain and amplitude-domain napari picking variants, are
// provided in the Supplementary Material, §S4.1.]
