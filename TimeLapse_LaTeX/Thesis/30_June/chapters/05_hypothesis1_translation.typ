#import "../template.typ": *

= Hypothesis 1: Can We Infer Subwavelength Movement from Phase Changes? <ch:hyp1>

@sec:meth-resolution established that amplitude-based detection of a
_stationary_ point scatterer fails below roughly half a wavelength of
separation. This chapter scales that finding to the time-lapse setting: if a
target moves by a sub-wavelength amount between a baseline and a monitor
survey, can that motion be inferred? Four synthetic gprMax experiments probe
this question from different angles, each sweeping a target displacement from
$2 lambda$ down to $1 \/ 32 lambda$ ($lambda approx 112.6 "mm"$ at the
$1.5 "GHz"$ centre frequency and pure-ice velocity used throughout):

+ *@sec:hyp1-lateral* --- a point scatterer shifted along the survey line.
+ *@sec:hyp1-vertical* --- the same point scatterer shifted deeper into the ice.
+ *@sec:hyp1-diagonal* --- the point scatterer shifted along a combined
  lateral-and-vertical ($2$:$1$) path.
+ *@sec:hyp1-fluidflow* --- a graded, spatially-extended wetting-front
  reflector rather than a point target, testing whether the phase-based
  method generalises beyond an idealised point scatterer.

For every experiment, the answer is the same: amplitude differencing fails
below the same resolution floor found in @sec:meth-resolution, but examining
the _phase change_ in the two-dimensional Fourier domain of the migrated
images (@ch:theory) recovers the displacement accurately down to the smallest
scale tested. @sec:hyp1-workflow defines the shared processing and analysis
pipeline used identically across all four experiments; @sec:hyp1-summary then
collects the resulting displacement errors into one master comparison.

#para-head[Hypothesis 1.] Can multi-dimensional phase-plane regression infer
lateral, vertical, and diagonal subwavelength displacements from time-lapse
migrated GPR images, at scales where amplitude differencing has already
failed?

An optional, local trace-based alternative to the global phase-plane fit used
in this chapter (Hypothesis 1.5: tracking the local phase gradients
$partial phi \/ partial x$ and $partial phi \/ partial y$ directly) is
explored separately in @sec:hyp1-h15 and is not required for the results
below.

== Workflow Defined <sec:hyp1-workflow>

Every experiment in this chapter is analysed by the same five-step pipeline,
applied identically to Lateral, Vertical, Diagonal, and FluidFlow --- only the
gprMax domain, the moving target, and the sweep direction differ between
them. No gprMax forward model or migration is re-run to produce this
chapter's figures and tables: every B-scan and migrated image is loaded from
the `.npz` caches already written by the corresponding
`*_TimeLapse_Playground.ipynb` / `FluidFlow_Playground.ipynb` notebooks.

+ *Model set up* (§X.1). The gprMax domain and grid are plotted with the
  baseline target position, followed by a second view marking every
  scenario's target position across the full displacement sweep.
+ *Raw and processed B-scans* (§X.2). The background-subtracted B-scan is
  shown for every scenario, using the tapering and $t_0$-shift conditioning
  of @sec:meth-conditioning.
+ *Migration results* (§X.3). Every scenario is migrated with Kirchhoff,
  Gazdag, and back-propagation, and the signed time-lapse-difference
  amplitude (monitor-minus-baseline) is compared across all three methods in
  one zoomed grid.
+ *Amplitude test* (§X.4). A Rayleigh-criterion argument is built from the
  *raw* (non-difference) migrated images, not the bipolar time-lapse
  difference: the Baseline and each Monitor scenario are each single-lobed
  point-spread functions, so the Baseline-to-Monitor peak-to-peak separation,
  measured against the Baseline PSF's own full width at half maximum (FWHM),
  is the standard two-point resolution argument. A separation/FWHM ratio well
  below $1$ means the two peaks cannot be distinguished from amplitude alone.
+ *Phase test* (§X.5, @sec:hyp1-phaseplane below). The 2D weighted-least-squares
  phase-plane fit is applied to the baseline/monitor cross-spectrum, using a
  window cropped around the target apex, and recovers the sub-wavelength
  displacement $(#Dz, #Dx)$ directly from the cross-spectrum's phase plane.

=== The Phase-Plane Fit (WLS) <sec:hyp1-phaseplane>

The 2D weighted-least-squares (WLS) phase-plane fit derived in
@sec:th-fourier-shift, @sec:th-wls, and @sec:th-material-change, and
implemented as described in @sec:meth-phaseplane, is applied to the
baseline/monitor cross-spectrum. The target apex is first localised from the
peak of the Baseline envelope nearest to where the two images differ most,
both images are then cropped to a $plus.minus 2.5 lambda$ window around that
apex, and the WLS fit recovers the sub-wavelength displacement
$(#Dz, #Dx)$ directly from the cropped cross-spectrum's phase plane, using
the central wavenumber $k_(z,c) = 2 pi \/ lambda$. A six-panel diagnostic
(cropped difference image with the apex marked, cross-spectrum phase,
cross-spectrum energy, a numeric True/Estimated/Error summary, and two
plane-fit panels) is shown for every scenario and every migration method, in
every §X.5 section below. The two plane-fit panels isolate the fitted plane
$phi = k_z #Dz + k_x #Dx + phi_0$ along each wavenumber axis separately, by
subtracting the *other* axis's fitted contribution from the measured
cross-spectrum phase ($phi - k_z #Dz$ plotted against $k_x$, and
$phi - k_x #Dx$ plotted against $k_z$): a good fit collapses the scattered,
per-bin phase measurements onto the fitted line in both panels. Each point is
coloured by the cross-spectrum magnitude $|X S|$, i.e. the weight that bin
actually received in the WLS fit, making visible which measurements drove the
result versus which were downweighted as noise.

One further subtlety affects the smallest scenarios tested: because the
FDTD grid cell is $1 "mm"$, the *nominal* fraction-of-$lambda$ displacement
requested of gprMax (e.g. $1\/32 lambda approx 3.52 "mm"$) is rounded to the
nearest grid cell before the simulation is run. The *true* displacement used
throughout this chapter's error tables is this grid-rounded value, not the
raw continuous fraction --- at the smallest scale tested this is a
$approx 14%$ correction, and it is why a handful of nominal-$1\/2 lambda$
scenarios below fall just inside a "$< 1\/2 lambda$" regime rather than
exactly on its boundary.

== Lateral Movement <sec:hyp1-lateral>

A single PEC cylinder ($r = 28 "mm"$, baseline depth $0.676 "m"$) is held
fixed in the baseline survey and displaced laterally by each of seven
scenarios (@fig:h1-lat-setup), from $2 lambda$ down to $1 \/ 32 lambda$,
in a $4.0 times 1.0 "m"$ pure-ice domain ($eps_r = 3.15$).

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
      [$2 lambda$], [$2 lambda$], [$0$],
      [$1 lambda$], [$1 lambda$], [$0$],
      [$1\/2 lambda$], [$1\/2 lambda$], [$0$],
      [$1\/4 lambda$], [$1\/4 lambda$], [$0$],
      [$1\/8 lambda$], [$1\/8 lambda$], [$0$],
      [$1\/16 lambda$], [$1\/16 lambda$], [$0$],
      [$1\/32 lambda$], [$1\/32 lambda$], [$0$],
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

=== Raw and Processed B-scans, and Migration Results <sec:hyp1-lat-bscans>

#figure(
  subfigs(cols: 1,
    img("H1_003_Lateral_Movement_--_Background-Subtracted_B-Scans.png"),
    img("H1_004_Lateral_--_TimeLapse_Migration_Comparison_Clean_--_Signed_Am.png", width: 95%),
  ),
  caption: [(a) Background-subtracted B-scans for the lateral time-lapse
    study, all eight scenarios; (b) signed time-lapse-difference amplitude
    (monitor-minus-baseline) for all three migration algorithms, clean
    data.],
) <fig:tl-summary-amp>

Individual per-method migrated images and time-lapse differences, zoomed
around the scatterer, are given in @app:hyp1-methods.

=== Amplitude Test <sec:hyp1-lat-amplitude>

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
algorithms.

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
    [$2 lambda$],     [$-221.66 (-98.5%)$],  [$-201.58 (-89.6%)$], [$-212.45 (-94.4%)$],
    [$1 lambda$],     [$-79.14 (-70.0%)$],   [$-39.65 (-35.1%)$],  [$-88.31 (-78.2%)$],
    [$1\/2 lambda$],  [$-92.44 (-165.1%)$],  [$-0.00 (-0.0%)$],    [$+0.04 (+0.1%)$],
    [$1\/4 lambda$],  [$-0.35 (-1.2%)$],     [$-0.00 (-0.0%)$],    [$+0.04 (+0.1%)$],
    [$1\/8 lambda$],  [$-0.17 (-1.2%)$],     [$-0.00 (-0.0%)$],    [$+0.02 (+0.1%)$],
    [$1\/16 lambda$], [$-0.10 (-1.5%)$],     [$-0.00 (-0.0%)$],    [$+0.01 (+0.1%)$],
    [$1\/32 lambda$], [$-0.06 (-1.6%)$],     [$-0.00 (-0.0%)$],    [$+0.01 (+0.1%)$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Lateral phase-plane WLS displacement error, $#Dx$
    (estimated $-$ true), in millimetres, with the equivalent percentage of
    the true displacement in parentheses.],
  kind: table,
) <tab:h1-lat-phase>

At $2 lambda$ and $1 lambda$ the cross-spectrum has wrapped and every method
is off by tens to hundreds of millimetres, as expected outside the
estimator's unambiguous range. From $1\/4 lambda$ down to $1\/32 lambda$,
however, every method recovers $#Dx$ to within $0.35 "mm"$ --- Gazdag and
Kirchhoff are already accurate at $1\/2 lambda$ (within $0.04 "mm"$), while
back-propagation's fit degrades at that single scale ($-92.44 "mm"$) before
also becoming sub-millimetre-accurate from $1\/4 lambda$ downward. This is
exactly the regime (@tab:h1-lat-amp) in which amplitude differencing has
already collapsed to an unresolvable single lobe for every method --- the
central empirical claim of Hypothesis 1.

== Vertical Movement <sec:hyp1-vertical>

This section mirrors @sec:hyp1-lateral, replacing lateral displacement with
_vertical_ (depth) displacement. The two directions are not expected to
behave identically: @sec:th-duality showed that a migrated image separates
lateral spatial frequencies like a prism but does not separate vertical
ones, so the results below are an important point of comparison for
@sec:hyp1-lateral.

A single PEC cylinder ($r = 28 "mm"$) sits at lateral position $x = 2.0 "m"$
and is displaced _downward_ from its baseline depth ($0.676 "m"$) across six
scenarios, from $1 lambda$ down to $1 \/ 32 lambda$, using the same domain
and grid as @sec:hyp1-lateral.

#supp-note[The forward-model setup (domain, scenario table, target depths),
background-subtracted B-scans, migration-comparison figure, Rayleigh-criterion
ratio table and PSF comparison, and phase-plane WLS displacement-error table
and shift-estimation diagnostics for this scenario set are all provided in
the Supplementary Material, §S2.2. Individual per-method migrated images and
time-lapse differences, zoomed around the scatterer, are given in
@app:hyp1-methods.]

Unlike the lateral case, amplitude differencing here stays resolvable a full
octave lower: every method's Rayleigh ratio remains above $1$ down to
$1\/4 lambda$ (minimum $1.590$), only dropping below $1$ at $1\/8 lambda$,
before collapsing ($<0.33$ for all methods) by $1\/32 lambda$ (compare
@tab:h1-lat-amp). The phase-plane fit only overtakes amplitude within this
narrower failure band: it is already accurate ($<=0.03 "mm"$ error) from
$1\/4 lambda$ downward, but at $1 lambda$ and $1\/2 lambda$ --- where
amplitude differencing remains comfortably resolvable --- the phase fit is
still wrapped and inaccurate. Vertical movement therefore shows the phase
estimator's advantage is concentrated specifically where amplitude has
already failed, not a blanket improvement across the whole sweep.

== Diagonal Movement <sec:hyp1-diagonal>

The third translation direction combines the previous two: the scatterer
moves simultaneously laterally and vertically, along a fixed $2$:$1$
diagonal ($#Dx = 2 thin #Dz$ in every scenario), so that all scenarios lie on
the same line through the baseline position. A single PEC cylinder at
$x = 2.0 "m"$, baseline depth $0.676 "m"$, is displaced diagonally across
five scenarios, using the same domain and grid as
@sec:hyp1-lateral and @sec:hyp1-vertical.

#supp-note[The forward-model setup (domain, scenario table, target
positions), background-subtracted B-scans, migration-comparison figure,
Rayleigh-criterion ratio table and PSF comparison, and phase-plane WLS
displacement-error tables ($#Dz$ and $#Dx$) and shift-estimation
diagnostics for this scenario set are all provided in the Supplementary
Material, §S2.3. Individual per-method migrated images and time-lapse
differences, zoomed around the scatterer, are given in @app:hyp1-methods.]

The diagonal amplitude ratio crosses below $1$ at Scenario 3
($1\/4 lambda_x$, $1\/8 lambda_z$) for Gazdag and Kirchhoff, and remains just
above $1$ for back-propagation at that scenario ($1.266$) before also
collapsing by Scenario 4 --- roughly the same combined-displacement failure
scale as the lateral case (@tab:h1-lat-amp), not the looser vertical floor
of @sec:hyp1-vertical. The phase-plane fit tracks this failure point
closely: Scenarios 1 and 2 remain wrapped for every method, as in the
lateral and vertical cases, but from Scenario 3 downward Gazdag and
Kirchhoff recover both $#Dz$ and $#Dx$ to within $0.02 "mm"$ (back-propagation
follows from Scenario 4, with a single $-1.51 "mm"$ outlier in $#Dz$ at
Scenario 3). The recovered $#Dz$/$#Dx$ ratio at Scenarios 3--5 matches the
known $#Dx = 2 #Dz$ scenario geometry to well within the reported error,
confirming the 2D WLS fit correctly separates the two simultaneous
displacement components.

== Fluid Flow <sec:hyp1-fluidflow>

The three translation studies above all move a rigid point scatterer; this
final experiment instead repeats the clean-data comparison for a target that
is directly relevant to the real fluid-injection field data of @ch:hyp3 ---
a _graded wetting zone_ rather than a discrete PEC cylinder. The domain, grid,
and centre frequency match @sec:hyp1-lateral exactly ($4.0 times 1.0 "m"$,
$f_c = 1.5 "GHz"$, $lambda = 112.6 "mm"$), but the moving target is a
$7$-step graded permittivity transition (box width $11.3 "mm"$, total
transition $78.8 "mm"$) between water and ice at depth $0.676 "m"$, standing
in for a fluid front advancing through a horizontal fracture. As in
@sec:hyp1-lateral, the front is swept laterally across the same seven
scenarios, from $2 lambda$ down to $1 \/ 32 lambda$ relative to its baseline
position (@fig:h1-ff-setup).

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
      [$2 lambda$], [$2 lambda$],
      [$1 lambda$], [$1 lambda$],
      [$1\/2 lambda$], [$1\/2 lambda$],
      [$1\/4 lambda$], [$1\/4 lambda$],
      [$1\/8 lambda$], [$1\/8 lambda$],
      [$1\/16 lambda$], [$1\/16 lambda$],
      [$1\/32 lambda$], [$1\/32 lambda$],
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

#supp-note[Background-subtracted B-scans, the migration-comparison figure,
the Rayleigh-criterion ratio table and PSF comparison, and the phase-plane
WLS front-displacement-error table and shift-estimation diagnostics for
this scenario set are all provided in the Supplementary Material, §S2.4.]

The graded front's broader intrinsic point-spread function makes amplitude
differencing far harder than for a discrete point scatterer, even at large
scales: Gazdag's ratio is already below $1$ at $2 lambda$, and Kirchhoff ---
the only method still resolvable ($1.237$) at $1\/2 lambda$ --- fails by
$1\/4 lambda$; every method has collapsed to $ratio <= 0.09$ by
$1\/16 lambda$. The phase-plane fit is unaffected by this earlier collapse:
Gazdag and Kirchhoff are already accurate (within $1.2 "mm"$) from
$1 lambda$ downward, and all three methods stay within $1 "mm"$ from
$1\/2 lambda$ down to $1\/32 lambda$ --- the entire regime in which
amplitude differencing has already failed. The graded, spatially-extended
front is therefore recovered by the same phase-plane approach used for the
point scatterers above, despite its inherently broader PSF.

== Summary of the Results <sec:hyp1-summary>

@tab:h1-lat-phase gives the full per-scenario phase-plane error for the
lateral case, in millimetres with the equivalent percentage of the true
displacement alongside in parentheses; the equivalent per-scenario tables
for Vertical, Diagonal, and FluidFlow are provided in the Supplementary
Material (§S2.2.3, §S2.3.3, §S2.4.3). @tab:h1-mae condenses the millimetre
errors for all four movement types into one mean absolute error (MAE)
per movement type and migration method, restricted to the sub-half-wavelength
regime ($1\/4 lambda$ down to $1\/32 lambda$, plus any nominal $1\/2 lambda$
scenario that the $1 "mm"$ FDTD grid rounds to just under $0.5 lambda$,
@sec:hyp1-workflow) --- the regime in which every amplitude test above has
already collapsed.

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
  caption: [Mean absolute phase-plane displacement error [mm], clean data,
    sub-half-wavelength regime only.],
  kind: table,
) <tab:h1-mae>

@fig:h1-detectability puts every movement type's Rayleigh-criterion ratio
and phase-plane error on one shared displacement axis, making @tab:h1-mae's
numbers directly comparable across Lateral, Vertical, Diagonal, and
FluidFlow at a glance: in every column, the bottom-row phase-error curves
drop below their $5%$ threshold at or before the top-row amplitude-ratio
curves cross below $1$, visually confirming that phase overtakes amplitude
precisely where amplitude differencing gives out, for every movement type
tested.

#figure(
  img("H1_037_Hypothesis_1_--_Detectability_Map_Amplitude_vs_Phase.png"),
  caption: [Detectability map, clean data: top row, Rayleigh-criterion
    amplitude ratio (threshold $1$); bottom row, absolute phase-plane WLS
    displacement error as a percentage of the true displacement (log scale,
    threshold $5%$); for all three migration methods, across all four
    movement types. The full per-scenario numeric tables underlying this
    figure are given in @tab:h1-lat-amp/@tab:h1-lat-phase (Lateral, in the
    main text) and the Supplementary Material, §S2.2--§S2.4 (Vertical,
    Diagonal, FluidFlow).],
) <fig:h1-detectability>

Two patterns stand out. First, for every point-scatterer geometry
(Lateral, Vertical, Diagonal), at least one migration method recovers the
true displacement to a few hundredths of a millimetre mean absolute error
(Gazdag on Lateral and Diagonal; Kirchhoff on Lateral) --- several orders of
magnitude below the wavelength scale, at displacements where
@sec:hyp1-lat-amplitude, @sec:hyp1-vertical, and
@sec:hyp1-diagonal showed amplitude differencing has already failed.
Second, the larger MAE values in @tab:h1-mae are driven almost entirely by a
*single* remaining scenario per (movement, method) pair, not by a systematic
sub-wavelength failure: back-propagation's $23.6 "mm"$ Lateral MAE is $99%$
attributable to its one poorly-resolved $1\/2 lambda$ scenario
(@tab:h1-lat-phase); Gazdag's and Kirchhoff's $16$--$17 "mm"$ Vertical MAE is
likewise dominated by their shared $1\/2 lambda$ scenario (Supplementary
Material, §S2.2.3), while $1\/4 lambda$ and below are already accurate to
$0.03 "mm"$ for all three methods. Diagonal, which combines both axes, is
the easiest case for every method once past its two largest scenarios.
FluidFlow's graded, spatially-extended front is recovered to within roughly
$1 "mm"$ mean absolute error by every method --- harder than the sharpest
point-scatterer results, but still two to three orders of magnitude below
the wavelength scale, and, per @sec:hyp1-fluidflow, at displacement scales
where amplitude differencing has already collapsed for all three algorithms.

Across all four experiments, from $1\/4 lambda$ down to $1\/32 lambda$, the
phase-plane fit remains accurate to a few tenths of a millimetre or better
for the large majority of (movement, method) pairs, exactly where amplitude
differencing has already failed --- supporting Hypothesis 1's central claim
that sub-wavelength displacement is recoverable from phase, not amplitude,
information. The one systematic weak point is back-propagation's fit at
$1\/2 lambda$ (Lateral) and at $1 lambda$/$1\/2 lambda$ (Vertical, shared
with the two analytic methods): whether this reflects a genuine
displacement-scale effect specific to back-propagation's excitation scheme,
or noise-free numerical sensitivity of the WLS fit at that particular scale,
is revisited in @ch:hyp2 once Laplace noise is introduced.
