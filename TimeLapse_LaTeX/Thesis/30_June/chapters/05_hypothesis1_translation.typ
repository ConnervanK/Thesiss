#import "../template.typ": *

= Hypothesis 1: Can We Infer Subwavelength Movement from Phase Changes? <ch:hyp1>

@sec:meth-resolution established that amplitude-based detection of a
_stationary_ point scatterer fails below roughly half a wavelength of
separation. This chapter scales that finding to the time-lapse setting: if a
target moves by a sub-wavelength amount between a baseline and a monitor
survey, can that motion be inferred? Four synthetic gprMax experiments probe
this question, each sweeping a target displacement from $2 lambda$ down to
$1 \/ 32 lambda$ ($lambda approx 112.6 "mm"$ at the $1.5 "GHz"$ centre
frequency and pure-ice velocity used throughout). Two are developed in full:
*lateral* translation of a point scatterer (@sec:hyp1-lateral), the worked
example through which the diagnostics are introduced, and a graded
*fluid-flow* front (@sec:hyp1-fluidflow), the more realistic target relevant
to the field data of @ch:hyp3. Two further directions --- *vertical*
(@sec:hyp1-vertical) and *diagonal* (@sec:hyp1-diagonal) translation of the
same point scatterer --- confirm the finding and are reported alongside the
cross-movement comparison in @sec:hyp1-summary.

For every experiment, the conclusion is consistent: simple amplitude differencing fails below the resolution floor found in @sec:meth-resolution. While windowed cross-correlation of migrated amplitudes can offer an intermediate solution by indirectly exploiting phase, directly examining the _phase change_ in the two-dimensional Fourier domain of the migrated images (@ch:theory) recovers the displacement accurately down to the smallest scale tested. @sec:hyp1-workflow defines the shared processing and analysis pipeline used identically across all four experiments.

#linebreak()

#para-head[Hypothesis 1.] Can multi-dimensional phase-plane regression infer
lateral, vertical, and diagonal subwavelength displacements from time-lapse
migrated GPR images, at scales where amplitude differencing has already
failed?

// An optional, local trace-based alternative to the global phase-plane fit used
// in this chapter (Hypothesis 1.5: tracking the local phase gradients
// $partial phi \/ partial x$ and $partial phi \/ partial y$ directly) is
// explored separately in the Supplementary Material, §S5, and is not required
// for the results below.

== Workflow Defined <sec:hyp1-workflow>

Every experiment in this chapter is analysed by the same five-step pipeline
(@fig:h1-workflow), applied identically to Lateral, Vertical, Diagonal, and
FluidFlow --- only the gprMax domain; the moving scatterer or fluid front; and the sweep
direction of the moving scatterer differ between them.

#figure(
  img("H1_040_Hypothesis_1_--_Shared_Analysis_Pipeline.png"),
  caption: [The shared five-step pipeline applied identically to every
    experiment in this chapter, parameterised only by the gprMax domain, the
    moving target, and the sweep direction.],
) <fig:h1-workflow>

+ *Model set up* (@sec:hyp1-lat-setup, @sec:hyp1-ff-setup). \
  The gprMax domain and grid are plotted with the
  baseline target position, followed by a second view marking every
  scenario's target position across the full displacement sweep.
+ *Raw and processed B-scans* (@sec:hyp1-lat-bscans, §S2.4). \
  The background-subtracted B-scan is shown for every scenario, using the
  tapering and $t_0$-shift conditioning of @sec:meth-conditioning.
+ *Migration results* (@sec:hyp1-lat-bscans, @sec:hyp1-fluidflow). \
  Every scenario is migrated with Kirchhoff,
  Gazdag, and back-propagation, and the signed time-lapse-difference
  amplitude (monitor-minus-baseline) is compared across all three methods in
  one zoomed grid.
+ *Amplitude test* (@sec:hyp1-lat-amplitude, @sec:hyp1-fluidflow). \
  A Rayleigh-criterion argument is built from the
  baseline and monitor migrated images, not from a time-lapse
  difference image. The migrated baseline and all monitor images are each single-lobed
  point-spread functions. We measure the FWHM value of the baseline PSF. The peak of the monitor PSF is then located, and the distance between the two peaks is measured. This distance is then divided by the baseline FWHM to give a separation/FWHM ratio. A ratio above $1$ means the two peaks are resolvable, while a ratio below $1$ means the two peaks cannot be distinguished from amplitude alone.
+ *Phase test* (@sec:hyp1-lat-phase, @sec:hyp1-fluidflow). \
  The 2D weighted-least-squares
  phase-plane fit is applied to the baseline/monitor cross-spectrum, using a
  window cropped around the target apex, and recovers the sub-wavelength
  displacement $(#Dz, #Dx)$ directly from the cross-spectrum's phase plane.
  The error metric is simply this estimate minus the true (grid-rounded)
  displacement, reported in millimetres and as a percentage of the true
  displacement, for every scenario and migration method.

=== The Phase-Plane Fit (WLS) <sec:hyp1-phaseplane>

The 2D weighted-least-squares (WLS) phase-plane fit derived in
@sec:th-fourier-shift and @sec:th-wls, and
implemented as described in @sec:meth-phaseplane, is applied to the
baseline/monitor cross-spectrum. The target apex is first localised by
taking the time-lapse difference (Monitor minus Baseline) of the migrated
images and computing its envelope (Hilbert magnitude) to find the
high-amplitude region where the two images differ most, then refining to
the nearby peak of the Baseline envelope itself, which anchors the crop
window precisely on the true, undisplaced target position; both baseline and monitor images are
then cropped to a $plus.minus 2.5 lambda$ window around that
apex, and the WLS fit recovers the sub-wavelength displacement
$(#Dz, #Dx)$ directly from the cropped cross-spectrum's phase plane, using
the central wavenumber $k_(z,c) = 2 pi \/ lambda$. 
A diagnostic image
(containing the cross-spectrum phase,
cross-spectrum energy, and the 1D plane-fit panels) illustrates this process for one representative scenario, like
in the example of @ch:theory (@fig:phaseplane-schematic). 

#linebreak()

The full set of per-scenario,
per-method diagnostics is instead provided in the Supplementary Material,
with the numeric results collected in this chapter's error tables (e.g.
@tab:h1-lat-phase). Because the fitted plane $phi = k_z #Dz + k_x #Dx +
phi_0$ lives in the 3-D space $(k_z, k_x, phi)$, it cannot be checked by
eye from a single 2-D plot; the two plane-fit panels instead give two 1-D
cross-sections through it. In the $k_x$ panel, the already-fitted $k_z
#Dz$ contribution is subtracted from every measured phase value ($phi -
k_z #Dz$) and the residual is plotted against $k_x$ alone, so a correct
fit collapses onto a straight line of slope $#Dx$; the $k_z$ panel does the
same with the roles of $k_z$ and $k_x$ reversed. A good fit therefore shows
the scattered per-bin phase measurements clustering tightly along that
fitted line in both panels, rather than scattering around it. A worked example of the
full six-panel diagnostic --- Lateral study, Gazdag migration, Baseline
versus each of the seven displacement scenarios --- is given in the
Supplementary Material, §S2.1.2.

#linebreak()

Each point is
coloured by the cross-spectrum magnitude $|#XS|$, i.e. the weight that bin
actually received in the WLS fit, making visible which measurements drove the
result versus which were downweighted as noise.

#linebreak()

One further subtlety affects the smallest scenarios tested: because the
FDTD grid cell is $1 "mm"$, the *nominal* fraction-of-$lambda$ displacement
requested of gprMax (e.g. $1\/32 lambda approx 3.52 "mm"$) is rounded to the
nearest grid cell before the simulation is run. The *true* displacement used
throughout this chapter's error tables is this grid-rounded value.

== Lateral Movement <sec:hyp1-lateral>

A single PEC cylinder ($r = 28 "mm"$, baseline depth $0.676 "m"$) is held
fixed in the baseline survey and displaced laterally by each of seven
scenarios (@fig:h1-lat-setup), from $2 lambda$ down to $1 \/ 32 lambda$,
in a $4.0 times 1.0 "m"$ pure-ice domain ($eps_r = 3.15$).

#pagebreak(weak: true)

=== Model Set Up <sec:hyp1-lat-setup>

#figure(
  subfigs(cols: 1,
    img("H1_001_Lateral_Movement_--_Model_Set_Up.png"),
    img("H1_002_Lateral_--_target_position_all_8_scenarios.png", width: 70%),
    table(
      columns: (auto, auto, auto),
      stroke: none,
      inset: (x: 0.8em, y: 0.3em),
      table.hline(stroke: 0.7pt),
      [*Scenario*], [*$#Dx$*], [*$#Dz$*],
      table.hline(stroke: 0.4pt),
      [Baseline], [$0$], [$0$],
      [$2 lambda$], [$22.5 "cm"$], [$0$],
      [$1 lambda$], [$11.3 "cm"$], [$0$],
      [$1\/2 lambda$], [$5.6 "cm"$], [$0$],
      [$1\/4 lambda$], [$2.8 "cm"$], [$0$],
      [$1\/8 lambda$], [$1.4 "cm"$], [$0$],
      [$1\/16 lambda$], [$0.7 "cm"$], [$0$],
      [$1\/32 lambda$], [$0.4 "cm"$], [$0$],
      table.hline(stroke: 0.7pt),
    ),
  ),
  caption: [Forward-model setup for the lateral time-lapse study: (a) the
    gprMax domain and grid, with the baseline scatterer position and (where
    present) the fixed second scatterer used for @sec:hyp1-lat-amplitude; (b)
    the target position for all eight scenarios (Baseline plus seven
    displacements), colour-coded from $2 lambda$ to $1 \/ 32 lambda$; (c) the
    eight displacement scenarios: the scatterer moves laterally ($#Dx$) while
    depth is fixed ($#Dz = 0$), swept from $2 lambda$ down to $1\/32 lambda$.],
) <fig:h1-lat-setup>

=== Migration Results <sec:hyp1-lat-bscans>

The raw and background-subtracted B-scans follow the standard conditioning
pipeline of @sec:meth-conditioning, already illustrated in @ch:methodology
(@fig:res-bscans); @fig:tl-summary-amp shows the resulting time-lapse
(migrated monitor image minus migrated baseline image) comparison directly.

#page(flipped: false)[
#figure(
  img("H1_004_Lateral_--_TimeLapse_Migration_Comparison_Clean_--_Signed_Am.png", width: 95%),
  caption: [Signed time-lapse-difference amplitude (monitor-minus-baseline)
    for all three migration algorithms, clean data, all eight scenarios.],
) <fig:tl-summary-amp>
]

Individual per-method migrated images and time-lapse differences, zoomed
around the scatterer, are given in the Supplementary Material, §S7.

=== Amplitude Test <sec:hyp1-lat-amplitude>

The amplitude test asks whether the Baseline and Monitor point-spread
functions remain distinguishable as two separate peaks after migration,
using the classical Rayleigh criterion: two peaks are resolved when their
separation is at least the Baseline peak's Full Width at Half Maximum
(FWHM), i.e. a separation-to-FWHM ratio of $1$ or more; below that, the two
peaks blur into one and amplitude alone can no longer tell Baseline and
Monitor apart. 

#linebreak()

Concretely, a 1-D profile through the *raw* migrated image
is extracted along the movement axis at the true target depth for both
Baseline and each Monitor scenario, windowed to $plus.minus 2.5 lambda$
around the baseline position so that unrelated features (e.g. noise for later scenarios) are ignored. The separation is
measured as the distance between the two profiles' peak positions, while
the FWHM is measured once, from the Baseline profile alone; the reported
ratio is simply that separation divided by the Baseline FWHM. @tab:h1-lat-amp
reports this ratio for every scenario and migration method.

#supp-note[The zoomed Baseline-versus-Monitor PSF comparison for this
scenario set is provided in the Supplementary Material, §S2.1.1.]

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$2 lambda$],     [7.306], [4.403], [3.656],
    [$1 lambda$],     [3.494], [2.106], [1.908],
    [$1\/2 lambda$],  [1.906], [1.149], [0.954],
    [$1\/4 lambda$],  [0.953], [0.574], [0.477],
    [$1\/8 lambda$],  [0.318], [0.383], [0.318],
    [$1\/16 lambda$], [0.318], [0.191], [0.159],
    [$1\/32 lambda$], [0.0],   [0.191], [0.159],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Lateral Rayleigh-criterion ratio (Baseline--Monitor peak
    separation / Baseline FWHM). A ratio well below $1$ means the two peaks
    cannot be distinguished from amplitude alone.],
  kind: table,
) <tab:h1-lat-amp>

Every method's ratio drops below $1$ at or before $1\/4 lambda$ (Kirchhoff
already at $1\/2 lambda$, ratio $0.954$) and falls below $0.4$ by
$1\/8 lambda$: from $1\/4 lambda$ downward, amplitude differencing alone
cannot resolve the lateral displacement for any of the three migration
algorithms. Back-propagation's ratio at $1\/32 lambda$ is exactly $0$ rather
than a small positive value like the other two methods at that scale: the
true displacement there ($0.4 "cm"$) falls below the spatial sampling of the
peak-position estimate, so the Baseline and Monitor peaks land on the
identical grid sample -- a finite-grid-size artefact of the FDTD
discretisation, not a qualitatively different failure mode.

=== Phase Test <sec:hyp1-lat-phase>

#supp-note[The Kirchhoff-, Gazdag-, and back-propagation-migrated
phase-plane shift-estimation diagnostics for this scenario set are provided
in the Supplementary Material, §S2.1.2.]

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$2 lambda$],     [$-221.66$], [$-201.58$], [$-212.45$],
    [$1 lambda$],     [$-79.14$],  [$-39.65$],  [$-88.31$],
    [$1\/2 lambda$],  [$-92.44$],  [$-0.00$],   [$+0.04$],
    [$1\/4 lambda$],  [$-0.35$],   [$-0.00$],   [$+0.04$],
    [$1\/8 lambda$],  [$-0.17$],   [$-0.00$],   [$+0.02$],
    [$1\/16 lambda$], [$-0.10$],   [$-0.00$],   [$+0.01$],
    [$1\/32 lambda$], [$-0.06$],   [$-0.00$],   [$+0.01$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Lateral phase-plane WLS displacement error, $#Dx$
    (estimated $-$ true), millimetres.],
  kind: table,
) <tab:h1-lat-phase>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$2 lambda$],     [$-98.5%$],  [$-89.6%$], [$-94.4%$],
    [$1 lambda$],     [$-70.0%$],  [$-35.1%$], [$-78.2%$],
    [$1\/2 lambda$],  [$-165.1%$], [$-0.0%$],  [$+0.1%$],
    [$1\/4 lambda$],  [$-1.2%$],   [$-0.0%$],  [$+0.1%$],
    [$1\/8 lambda$],  [$-1.2%$],   [$-0.0%$],  [$+0.1%$],
    [$1\/16 lambda$], [$-1.5%$],   [$-0.0%$],  [$+0.1%$],
    [$1\/32 lambda$], [$-1.6%$],   [$-0.0%$],  [$+0.1%$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Lateral phase-plane WLS displacement error, $#Dx$, as a
    percentage of the true displacement.],
  kind: table,
) <tab:h1-lat-phase-pct>

From $1\/4 lambda$ down to $1\/32 lambda$ --- exactly the regime
(@tab:h1-lat-amp) in which amplitude differencing has already collapsed to an
unresolvable single lobe for every method --- every method recovers $#Dx$ to
within $0.35 "mm"$: the central empirical claim of Hypothesis 1. The two
largest scales behave as expected outside the estimator's unambiguous range:
at $2 lambda$ and $1 lambda$ the cross-spectrum has wrapped and every method
is off by tens to hundreds of millimetres. In between, Gazdag and Kirchhoff
are already accurate at $1\/2 lambda$ (within $0.04 "mm"$), while
back-propagation's fit degrades at that single scale ($-92.44 "mm"$) before
also becoming sub-millimetre-accurate from $1\/4 lambda$ downward.

== Fluid Flow <sec:hyp1-fluidflow>

Where @sec:hyp1-lateral moved a rigid point scatterer, this second experiment
repeats the comparison for the more realistic target that motivates the field
study of @ch:hyp3 --- a _graded wetting zone_ rather than a discrete PEC
cylinder. The domain, grid,
and centre frequency match @sec:hyp1-lateral exactly ($4.0 times 1.0 "m"$,
$f_c = 1.5 "GHz"$, $lambda = 112.6 "mm"$), but the moving target is a
$7$-step graded permittivity transition (box width $11.3 "mm" approx
lambda \/ 10$, total transition $78.8 "mm" approx 0.7 lambda$) between
water ($eps_r = 80$) and air ($eps_r = 1$) at depth $0.676 "m"$, stepping
down in $eps_r$ by $10$ per box ($70, 60, dots, 10$) before the final,
slightly smaller step to air, standing
in for a fluid front advancing through a horizontal fracture. As in
@sec:hyp1-lateral, the front is swept laterally across the same seven
scenarios, from $2 lambda$ down to $1 \/ 32 lambda$ relative to its baseline
position (@fig:h1-ff-setup).

#pagebreak(weak: true)

=== Model Set Up <sec:hyp1-ff-setup>

#figure(
  subfigs(cols: 1,
    img("H1_028_FluidFlow_--_Model_Set_Up_baseline_graded_zone_true_scale.png"),
    img("H1_029_FluidFlow_--_graded_wetting_zone_all_8_scenarios_box_width11.png", width: 90%),
    table(
      columns: (auto, auto),
      stroke: none,
      inset: (x: 0.8em, y: 0.3em),
      table.hline(stroke: 0.7pt),
      [*Scenario*], [*Front displacement $#Dx$*],
      table.hline(stroke: 0.4pt),
      [Baseline], [$0$],
      [$2 lambda$], [$22.5 "cm"$],
      [$1 lambda$], [$11.3 "cm"$],
      [$1\/2 lambda$], [$5.6 "cm"$],
      [$1\/4 lambda$], [$2.8 "cm"$],
      [$1\/8 lambda$], [$1.4 "cm"$],
      [$1\/16 lambda$], [$0.7 "cm"$],
      [$1\/32 lambda$], [$0.4 "cm"$],
      table.hline(stroke: 0.7pt),
    ),
  ),
  caption: [Forward-model setup for the fluid-flow time-lapse study: (a) the
    gprMax domain with the baseline graded wetting zone at true scale; (b)
    the graded zone's position for all eight scenarios (Baseline plus seven
    displacements); (c) the eight displacement scenarios: the graded wetting
    front's centroid moves laterally, swept from $2 lambda$ down to
    $1\/32 lambda$ --- identical sweep to @sec:hyp1-lateral.],
) <fig:h1-ff-setup>

#linebreak()

The graded front's broader intrinsic point-spread function (visible in the
zoomed fluid-flow migration results and PSF's: Supplementary Material, §S2.4.2, §S2.4.1) makes amplitude
differencing far harder than for a discrete point scatterer, even at large
scales: Gazdag's FWHM ratio is already below $1$ at $2 lambda$, and Kirchhoff ---
the only method still resolvable ($1.237$) at $1\/2 lambda$ --- fails by
$1\/4 lambda$; every method has collapsed to $ratio <= 0.09$ by
$1\/16 lambda$. The phase-plane fit is unaffected by this earlier collapse:
Gazdag and Kirchhoff are already accurate (within $1.2 "mm"$) from
$1 lambda$ downward, and all three methods stay within $1 "mm"$ displacement
estimate error from
$1\/2 lambda$ down to $1\/32 lambda$ --- the entire regime in which
amplitude differencing has already failed. The graded, spatially-extended
front is therefore recovered by the same phase-plane approach used for the
point scatterers above, despite its inherently broader PSF.

#linebreak()

#supp-note[Background-subtracted B-scans, the migration-comparison figure,
the Rayleigh-criterion ratio table and PSF comparison, and the phase-plane
WLS front-displacement-error table and shift-estimation diagnostics for
this scenario set are all provided in the Supplementary Material, §S2.4.]

== Results across Movement Types <sec:hyp1-summary>

The lateral worked example (@sec:hyp1-lateral) and the fluid-flow front
(@sec:hyp1-fluidflow) are joined by two confirmatory directions --- vertical
and diagonal translation of the same point scatterer --- before all four are
compared side by side.

=== Vertical <sec:hyp1-vertical>

Displacing the same point scatterer in depth rather than laterally confirms the
phase-plane finding, with one direction-specific nuance predicted by
@sec:th-duality. Because a migrated image does not separate _vertical_ spatial
frequencies the way it separates lateral ones, vertical amplitude differencing
stays resolvable a full octave lower --- every method's Rayleigh ratio holds
above $1$ down to $1\/4 lambda$, only failing between $1\/4 lambda$ and
$1\/8 lambda$, against $1\/2 lambda$ for the lateral case. Because amplitude
differencing already works down to a smaller scale here, the regime where
phase estimation offers a genuine advantage over it is correspondingly
narrower; even so, within that regime the fit remains accurate to
$<=0.03 "mm"$ from $1\/4 lambda$ downward (Supplementary Material, §S2.2).

=== Diagonal <sec:hyp1-diagonal>

Moving the scatterer along a combined $2$:$1$ lateral-and-vertical path
($#Dx = 2 #Dz$) confirms that the 2D fit separates two simultaneous components
correctly: from Scenario 3 downward the recovered $#Dz$/$#Dx$ ratio matches the
known $#Dx = 2 #Dz$ geometry to within the reported error, with Gazdag and
Kirchhoff accurate to $0.02 "mm"$ on both axes. The combined displacement
collapses amplitude at roughly the lateral scale, not the looser vertical one
(Supplementary Material, §S2.3).

=== Cross-Movement Comparison

@tab:h1-lat-phase and @tab:h1-lat-phase-pct give the full per-scenario
phase-plane error for the lateral case; the equivalent per-scenario tables for Vertical, Diagonal, and
FluidFlow are provided in the Supplementary Material (§S2.2.3, §S2.3.3,
§S2.4.3). @fig:h1-mae-summary condenses the mean absolute error (MAE) for
all four movement types and three migration methods into one figure,
restricted to the sub-half-wavelength regime ($1\/4 lambda$ down to
$1\/32 lambda$, plus any nominal $1\/2 lambda$ scenario that the $1 "mm"$
FDTD grid rounds to just under $0.5 lambda$, @sec:hyp1-workflow) --- the
regime in which every amplitude test above has already collapsed; exact
values are given in @tab:h1-mae below.

#figure(
  img("H1_039_Hypothesis_1_--_MAE_Summary_Across_Movement_Types.png", width: 85%),
  caption: [Mean absolute phase-plane displacement error by movement type
    and migration method, clean data, sub-half-wavelength regime only.
    Background shading groups rows by movement type.],
) <fig:h1-mae-summary>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Movement*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [Lateral],   [23.645], [0.000], [0.027],
    [Vertical],  [11.950], [16.079], [16.799],
    [Diagonal],  [0.683],  [0.003],  [0.016],
    [FluidFlow], [0.967],  [0.172],  [0.335],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Mean absolute phase-plane displacement error [mm] underlying
    @fig:h1-mae-summary.],
  kind: table,
) <tab:h1-mae>

#page(flipped: true)[
#figure(
  img("H1_037_Hypothesis_1_--_Detectability_Map_Amplitude_vs_Phase.png"),
  caption: [Detectability map, clean data: top row, Rayleigh-criterion
    amplitude ratio (threshold $1$); bottom row, absolute phase-plane WLS
    displacement error as a percentage of the true displacement (log scale,
    threshold $5%$); for all three migration methods, across all four
    movement types. The full per-scenario numeric tables underlying this
    figure are given in @tab:h1-lat-amp/@tab:h1-lat-phase/@tab:h1-lat-phase-pct
    (Lateral, in the main text) and the Supplementary Material, §S2.2--§S2.4 (Vertical,
    Diagonal, FluidFlow).],
) <fig:h1-detectability>
]

In every column of @fig:h1-detectability, the bottom-row phase-error curves
drop below their $5%$ threshold at or before the top-row amplitude-ratio
curves cross below $1$: phase overtakes amplitude precisely where amplitude
differencing gives out, for every movement type tested. For every
point-scatterer geometry (Lateral, Vertical, Diagonal), at least one method
recovers the true displacement to a few hundredths of a millimetre MAE
(Gazdag on Lateral and Diagonal; Kirchhoff on Lateral) --- several orders of
magnitude below the wavelength scale. The larger MAE values visible in
@fig:h1-mae-summary are not evidence of systematic sub-wavelength failure:
each is driven almost entirely by a single poorly-resolved scenario per
(movement, method) pair, close to the amplitude/phase crossover itself
(back-propagation's Lateral fit at $1\/2 lambda$; the shared Vertical
outlier at $1\/2 lambda$ for Gazdag and Kirchhoff), with every scale below
that already accurate to a few tenths of a millimetre or better. FluidFlow's
graded, spatially-extended front is the hardest case tested (MAE around
$1 "mm"$) but is still two to three orders of magnitude below the
wavelength scale, confirming that the method's advantage holds not only for
a moving point scatterer but also for a moving fluid front, and is not
simply an artefact of the idealised point-scatterer geometry used elsewhere
in this chapter.

#linebreak()

Across all four experiments, from $1\/4 lambda$ down to $1\/32 lambda$, the
phase-plane fit remains accurate to a few tenths of a millimetre or better
for the large majority of (movement, method) pairs, exactly where amplitude
differencing has already failed --- supporting Hypothesis 1's central claim
that sub-wavelength displacement is recoverable from phase, not amplitude,
information. The one systematic weak point is back-propagation's fit at
$1\/2 lambda$ (Lateral) and at $1 lambda$/$1\/2 lambda$ (Vertical, shared
with the two analytic methods); whether this is a genuine
phase wrapping effect or noise-free numerical sensitivity of the WLS
fit is revisited in @ch:hyp2 once Laplace noise is introduced.
