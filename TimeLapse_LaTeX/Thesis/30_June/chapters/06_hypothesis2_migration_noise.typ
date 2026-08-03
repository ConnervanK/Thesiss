#import "../template.typ": *

= Hypothesis 2: Which Migration Technique Is Best Suited for Noise-Robust Phase-Plane Tracking? <ch:hyp2>

The time-lapse phase-plane method of this thesis recovers a sub-wavelength
displacement by fitting a plane to the phase of the cross-spectrum of a
baseline and a monitor migrated image (@ch:theory), the plane's two slopes
giving the vertical and lateral shift. @ch:hyp1 showed this works on clean
synthetic data across four experiments --- lateral, vertical, and diagonal
translation of a point scatterer, and a graded fluid front --- recovering
displacements well below the amplitude-based resolution floor. Real GPR data,
though, are never clean. This chapter therefore repeats all four experiments
under a realistic Laplace noise model fitted to real GPR field data, to answer
a question the clean study could not: does the choice of migration algorithm
matter once the data are noisy? Migration choice turns out to matter
substantially, but no single method is unconditionally "most noise-robust":
Kirchhoff achieves the lowest aggregate displacement error of the three,
back-propagation with sign-bit time-reversal is the only method that does not
manufacture spurious coherent structure from noise alone, and Gazdag --- while
substantially improved once the target-localisation problem of
@sec:hyp3-groundtruth-apex is corrected --- remains the least accurate
overall.

#linebreak()

#para-head[Hypothesis 2.] Which migration algorithm is best suited to
noise-robust time-lapse phase-plane tracking under a realistic, heavy-tailed
Laplace noise model, and does the answer depend on what "noise-robust" is
taken to mean --- lowest aggregate displacement error, or avoiding
false-positive structure manufactured from noise alone?

== Noise Creation <sec:hyp3-laplace>

Rather than injecting arbitrary synthetic noise, the noise level and shape
used throughout this thesis are fitted to real GPR field noise. A sample of
field noise is tracked through eleven stages of the processing pipeline and
fitted with both a Laplace and a Gaussian distribution at the last stage
before the spherical-gain correction, since the gain correction inflates the
amplitude scale by several orders of magnitude and is not representative of
the raw simulated $E_z$ amplitudes used elsewhere in this thesis.

#supp-note[The noise amplitude distribution at every one of the eleven
tracked processing stages, each with a Gaussian reference overlay, is
provided in the Supplementary Material, §S3.5.]

#figure(
  subfigs(cols: 1,
    img("H2_031_Hypothesis_2_--_Noise_window_signal_zone_reference.png", width: 90%),
    img("H2_002_Noise_distribution_at_7_Constant_Velocity_pre-gain_with_Lapl.png", width: 70%),
    table(
      columns: (auto, auto),
      stroke: none,
      inset: (x: 0.8em, y: 0.3em),
      table.hline(stroke: 0.7pt),
      [*Laplace scale*], [*Gaussian $sigma$*],
      table.hline(stroke: 0.4pt),
      [$6.085$], [$9.233$],
      table.hline(stroke: 0.7pt),
    ),
  ),
  caption: [(a) B-scan showing the signal zone the noise sample is drawn
    from; (b) noise distribution at the pipeline stage actually sampled to
    generate every noisy dataset used from @sec:hyp3-purenoise onward (the
    last stage before the spherical-gain correction), with both Laplace and
    Gaussian fits overlaid; (c) the fitted Laplace scale and Gaussian
    $sigma$ at that stage; loc $= 0$ ($n = 976 space 244$ samples). The
    heavier-tailed Laplace distribution is adopted throughout this thesis.],
) <fig:h2-noise-laplace-fit>

The fitted distribution is heavier-tailed than a Gaussian of
matched variance, consistent with field GPR noise being dominated by
occasional large-amplitude clutter and interference rather than purely
thermal noise. The pre-gain fit is the one actually used to generate
synthetic noise for every experiment in @sec:hyp3-lateral,
@sec:hyp3-vertical, @sec:hyp3-diagonal, and @sec:hyp3-fluidflow --- its
Laplace _shape_ (loc $= 0$, heavier tails than Gaussian) is kept, but its
scale is
rescaled so that the resulting noise standard deviation is exactly $10%$ of
each synthetic B-scan's own signal standard deviation --- a light, realistic
noise level rather than the raw fitted scale, which would be incommensurate
with the synthetic $E_z$ amplitudes.

== Migrating Noise <sec:hyp3-purenoise>

Before testing noise robustness on real signal-plus-noise data, a sanity
check establishes what each migration algorithm does to noise _alone_: an
empty B-scan (no scatterer, no real signal) containing only the Laplace
noise of @sec:hyp3-laplace is migrated with Kirchhoff, Gazdag, and
back-propagation. If a method turns pure noise into something that looks
like a coherent, scatterer-like focus, every noisy result elsewhere in this
chapter carries a false-positive risk that must be accounted for.

#figure(
  img("H2_003_Migrating_Pure_Noise_no_scatterers_no_signal_--_Kirchhoff_Ga.png"),
  caption: [Migrating pure noise (no scatterers, no signal): the input
    pure-noise B-scan, its Kirchhoff migration, its Gazdag migration, and the
    focus-time snapshot of its back-propagated wavefield (default,
    peak-normalised excitation).],
) <fig:h2-purenoise-kg>

*Result: the three methods do not fail the same way.* Kirchhoff turns pure
noise into clearly coherent, smooth wave-like bands that could easily be
misread as real layered structure: its delay-and-sum aperture stacking
imposes coherence on incoherent input by construction, summing many traces
along travel-time curves and smoothing incoherent noise into
locally-correlated structure. Gazdag's noise output, by contrast, stays
speckled and incoherent, with no wave-like artefacts --- just texture. This
is a direct, practical consequence of the structural difference between the
two algorithms (@sec:th-migration): Kirchhoff's spatial stacking manufactures
apparent coherence from nothing, while Gazdag's frequency-domain downward
continuation does not.

With no scatterer present, there is also no true location for either
back-propagation excitation scheme to focus on, so the back-propagated
wavefield (rightmost panel) stays diffuse speckle throughout the domain
rather than collapsing into an obvious spurious bright spot --- unlike
Kirchhoff's coherent bands, there is no clearly visible artefact here to
point at. What does stand out is the amplitude scale:
peak-normalised excitation reaches only $approx 116 space 600 "V/m"$, since
only each trace's single largest sample is normalised to $plus.minus 1$ and
every other noise sample stays small. @sec:hyp3-signbit revisits this with
sign-bit excitation, where every sample --- not just the peak --- is forced
to $plus.minus 1$.

The image-domain comparison above shows Kirchhoff manufacturing spatially
coherent bands from noise, but says nothing about whether that coherence is
tied to a specific frequency. A complementary spectral view --- each
method's pure-noise output averaged into one frequency spectrum, with every
domain's natural sample axis (time, depth, or back-propagation snapshot
depth) converted to an equivalent frequency via $f = v_"ice" k$ so all four
curves share one physically comparable axis --- checks this directly in @fig:h2-purenoise-spectral.

#figure(
  img("H2_004_Migrating_Pure_Noise_--_Spectral_Content_by_Method__axes_res.png", width: 80%),
  caption: [Spectral content of migrated pure noise by method: the input
    noise floor, Kirchhoff migration, Gazdag migration, and peak-normalised
    back-propagation, each averaged into one frequency spectrum (axes
    rescaled to an equivalent frequency via $f = v_"ice" k$).],
) <fig:h2-purenoise-spectral>

The input noise floor, Gazdag, and back-propagation all stay close to flat
across the full $0$--$5 "GHz"$ range shown, with no pronounced peak anywhere
near the Ricker centre frequency ($f_c = 1.5 "GHz"$) --- consistent with
Gazdag's and back-propagation's incoherent, speckled image-domain output
above. Kirchhoff is the outlier: its spectrum is suppressed below
$approx 1 "GHz"$, rises steeply through $1$--$2.5 "GHz"$, peaks around
$2.5$--$3 "GHz"$ (*above*, not at, $f_c$), and then rolls off toward
$5 "GHz"$ --- a pronounced band-pass shape none of the other three methods
share. This does not contradict the coherent-bands finding above; rather, it
localises the mechanism: Kirchhoff is not manufacturing energy at the
GPR wavelet's own frequency, so its false-coherence risk is a
*spatial*-stacking effect (the delay-and-sum aperture's geometry, not
its frequency response) rather than a narrowband spectral one. This is a
purely spectral view, though: it discards the phase/spatial-alignment
information that produces the low frequency speckle visible in
@fig:h2-purenoise-kg, so it cannot by itself explain *why* Kirchhoff's
aperture stacking prefers this particular frequency band.

== Making Back-Propagation Noise-Robust <sec:hyp3-signbit>

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

While percentile normalization is often used to clip such extreme outliers, it still preserves relative amplitude variations. Because focusing is driven by phase, this thesis bypasses amplitude-capping entirely in favor of _sign-bit time-reversal_:

By injecting only the _sign_ of the
time-reversed wavefield instead of its peak-normalised value,
$ u_"sign" (x, tau) = op("sign")(u(x, tau)) , $ <eq:signbit> every zero-crossing and phase trend is kept completely intact, since the sign of a signal
carries its full phase information, while squashing every noise spike down
to the same $plus.minus 1$ amplitude as the coherent signal --- stripping noise of the outsized amplitude that would otherwise dominate the back-propagated wavefield. The clean-data experiments of @ch:hyp1 use the
default peak-normalised excitation throughout, since they have no noise to
suppress; every noisy back-propagation result in this chapter uses sign-bit
excitation instead.

#linebreak()

In @fig:h2-signbit-excitation (b) sign-bit excitation reaches $approx 15 space 100 "V/m"$, about $25 times$
larger than peak-normalised's, since forcing every sample (not just each
trace's single peak) to $plus.minus 1$ injects far more total energy into
the medium. Unlike the idealised, zero-signal test of @sec:hyp3-purenoise,
both excitation schemes here have a real target --- the Baseline scatterer
--- to focus on, so @fig:h2-signbit-excitation (b) is the direct,
signal-bearing analogue of that pure-noise comparison: compare how tightly
each panel's energy collapses onto the true scatterer position rather than
staying diffuse artefact, as it did for pure noise. 
#draftnote[Describe what
the regenerated focus-frame comparison actually shows once the peak-norm
back-propagation gprMax run has completed --- see the run instructions
printed by the corresponding Hypothesis_2.ipynb cell.] 
The actual,
quantitative evidence that sign-bit back-propagation is noise-robust comes
from @sec:hyp3-summary's master MAE table (@tab:h2-mae) on the real noisy studies below,
where back-propagation is the second most accurate method overall.

#figure(
  subfigs(cols: 1,
    img("H2_005_Sign-Bit_Time-Reversed_Excitation_--_B-scans_and_Spectra_Lat.png", width: 100%),
    img("H2_006_Back-Propagation_of_a_Signal-Bearing_B-scan_Lateral_Baseline.png", width: 100%),
  ),
  caption: [Sign-bit time-reversal, Lateral, Noisy: (a) original noisy
    B-scans, sign-bit B-scans, and their respective frequency spectra, for
    every displacement scenario; (b) back-propagation of a real,
    signal-bearing B-scan (Baseline scatterer, noisy) --- default
    peak-normalised excitation versus sign-bit excitation, focus-time
    snapshots side by side.],
) <fig:h2-signbit-excitation>



== Extra Processing Steps in the Phase Domain to Mitigate Noise <sec:hyp3-phase-denoise>

Every phase-plane fit so far weights each masked cross-spectrum bin by its
magnitude $|X S|$ (WLS, @sec:hyp1-phaseplane), so that noisy, low-power bins
contribute less than high-power, high-SNR ones. This section tests that
weighting choice directly, by re-running the identical fit with every masked
bin given *equal* weight instead (OLS) and comparing the two on the hardest
realistic case in this chapter: back-propagation-migrated *diagonal*
movement, the smallest sub-wavelength shift tested ($1\/16 lambda_x$,
$1\/32 lambda_y$) against Baseline --- the scenario with the weakest signal
relative to the noise floor.

#figure(
  img("H2_007_Chapter_54_--_OLS_vs_WLS_Phase-Plane_Fitting_on_Noisy_Data.png"),
  caption: [OLS (unweighted, top row) versus WLS (weighted, bottom row)
    phase-plane fitting on the same noisy cross-spectrum, diagonal movement:
    cross-spectrum phase and cross-spectrum power, followed by two 1D
    cross-sections through the fitted plane along $k_x$ and $k_z$, with the
    same masked bins scattered and coloured by their (raw) weight $|X S|$ in
    both rows.],
) <fig:h2-ols-vs-wls>

Both fits see exactly the same masked cross-spectrum bins (the scatter
points in the two cross-section panels are identical between rows); only
the weight each bin is given during the least-squares solve differs. In
this diagonal-movement example the true shift along $x$ is $+7.000 "mm"$.
OLS recovers $+6.867 "mm"$ ($-0.133 "mm"$ error) and WLS recovers $+6.945
"mm"$ ($-0.055 "mm"$ error) --- WLS's error is under half of OLS's here. The
cross-section panels make the
*mechanism* visible: the bright
(high-weight) points visibly cluster closer to the fitted line than the dim
(low-weight, noise-dominated) points in both rows, confirming that $|X
S|$-weighting does discount incoherent bins as intended.
@sec:hyp3-summary's master MAE table (computed with WLS throughout this
chapter, consistent with @ch:hyp1) remains the relevant comparison for
overall noise robustness across methods, not this single scenario's
OLS/WLS pair.

== Locating the Target Under Noise <sec:hyp3-groundtruth-apex>

Every phase-plane fit in this chapter needs a crop window centred on the
target before the WLS fit of @sec:hyp1-phaseplane can run
(@sec:hyp1-phaseplane). On clean data (@ch:hyp1) that centre is found by
locating the peak of the Baseline envelope nearest to where Baseline and
Monitor differ most --- a search that works because the clean signal is, by
construction, the dominant feature in the image. Under Laplace noise that
assumption breaks down: a noisy-envelope search, cropped only along $x$ (keeping the full
depth range), locks onto a noise-driven false peak instead
of the true target --- especially for Gazdag, whose incoherent speckle
(@sec:hyp3-purenoise) both dominates the envelope search itself and, once a
wrong window is cropped, floods the WLS fit with off-target energy.

Because every experiment in this chapter is synthetic FDTD data, the true
target position is always known exactly, in both $x$ and $z$ --- a shortcut
that would not be available on real field data, but is legitimate here. Every
result below therefore crops directly around that known position (a
$plus.minus 2.5 lambda$ window in *both* $x$ and $z$, not $x$ alone), reused
identically across every scenario and all three migration methods, rather
than re-localising it from the noisy image. Restricting the crop in depth as
well as laterally matters specifically because Gazdag's noise speckle is not
confined to the target's depth: it fills the migrated image at every depth,
so a crop that is tight in $x$ but left open in $z$ still lets speckle from
other depths dominate the fit's weighted bins. This localisation fix is
applied uniformly to Kirchhoff, Gazdag, and back-propagation alike, so it
cannot by itself explain any remaining *difference* between methods below
--- but it substantially changes the absolute accuracy each one achieves, most visibly for Gazdag (@sec:hyp3-summary).

== Noisy Lateral Movement <sec:hyp3-lateral>

The lateral displacement sweep of @sec:hyp1-lateral (@fig:h1-lat-setup)
is repeated on B-scans contaminated with the Laplace noise of
@sec:hyp3-laplace, for all three migration algorithms, using sign-bit
time-reversal (@sec:hyp3-signbit) for back-propagation.

=== Migration Results <sec:hyp3-lat-migration>

#figure(
  img("H2_008_Lateral_--_TimeLapse_Migration_Comparison_Noisy_--_Signed_Am.png", width: 95%),
  caption: [Signed time-lapse-difference amplitude (monitor-minus-baseline)
    for all three migration algorithms, lateral time-lapse study, _noisy_
    data.],
) <fig:h2-lat-summary-amp>

=== Amplitude Test <sec:hyp3-lat-amplitude>

#supp-note[The zoomed Baseline-versus-Monitor PSF comparison for this
_noisy_ scenario set is provided in the Supplementary Material, §S3.1.1.]

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$2 lambda$],     [7.839], [4.152], [3.656],
    [$1 lambda$],     [3.749], [2.076], [1.907],
    [$1\/2 lambda$],  [2.045], [0.944], [0.954],
    [$1\/4 lambda$],  [1.023], [0.566], [0.477],
    [$1\/8 lambda$],  [0.341], [0.377], [0.318],
    [$1\/16 lambda$], [0.341], [0.189], [0.159],
    [$1\/32 lambda$], [0.341], [0.0],   [0.159],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Lateral Rayleigh-criterion ratio, _noisy_ data (compare
    @tab:h1-lat-amp).],
  kind: table,
) <tab:h2-lat-amp>

Noise barely changes the amplitude-detectability floor: the ratio still
crosses below $1$ between $1\/2 lambda$ and $1\/4 lambda$ for every method,
matching the clean-data result almost exactly (@tab:h1-lat-amp) --- the
Rayleigh criterion is a property of the raw PSF width, which $10%$ noise
perturbs only slightly.

=== Phase Test <sec:hyp3-lat-phase>

#supp-note[The Kirchhoff-, Gazdag-, and back-propagation-migrated (sign-bit)
phase-plane shift-estimation diagnostics for this _noisy_ scenario set are
provided in the Supplementary Material, §S3.1.2.]

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$2 lambda$],     [$-214.51$], [$-216.16$], [$-208.96$],
    [$1 lambda$],     [$-102.61$], [$-73.20$],  [$-83.94$],
    [$1\/2 lambda$],  [$-67.94$],  [$-51.30$],  [$+0.23$],
    [$1\/4 lambda$],  [$-0.05$],   [$+6.24$],   [$+0.19$],
    [$1\/8 lambda$],  [$+0.20$],   [$-32.87$],  [$+0.24$],
    [$1\/16 lambda$], [$+0.03$],   [$-10.58$],  [$+0.05$],
    [$1\/32 lambda$], [$+0.18$],   [$-8.59$],   [$+0.28$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Lateral phase-plane WLS displacement error, $#Dx$
    (estimated $-$ true), millimetres, _noisy_ data (compare
    @tab:h1-lat-phase).],
  kind: table,
) <tab:h2-lat-phase>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$2 lambda$],     [$-95.3%$],  [$-96.1%$],  [$-92.9%$],
    [$1 lambda$],     [$-90.8%$],  [$-64.8%$],  [$-74.3%$],
    [$1\/2 lambda$],  [$-121.3%$], [$-91.6%$],  [$+0.4%$],
    [$1\/4 lambda$],  [$-0.2%$],   [$+22.3%$],  [$+0.7%$],
    [$1\/8 lambda$],  [$+1.4%$],   [$-234.8%$], [$+1.7%$],
    [$1\/16 lambda$], [$+0.4%$],   [$-151.1%$], [$+0.7%$],
    [$1\/32 lambda$], [$+4.5%$],   [$-214.8%$], [$+7.1%$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Lateral phase-plane WLS displacement error, $#Dx$, as a
    percentage of the true displacement, _noisy_ data (compare
    @tab:h1-lat-phase).],
  kind: table,
) <tab:h2-lat-phase-pct>

Noise degrades every method relative to the clean-data result, but with the
crop window now centred on the *known* scatterer position in both $x$ and
$z$ (@sec:hyp3-groundtruth-apex) rather than an apex hunted for in the noisy
envelope, Kirchhoff and back-propagation both stay within a third of a
millimetre from $1\/4 lambda$ downward (Kirchhoff $<=0.28 "mm"$,
back-propagation $<=0.20 "mm"$) --- noise-driven scatter around zero rather
than a systematic bias. Gazdag remains qualitatively worse throughout the
sub-half-wavelength regime ($+6.2$ to $-32.9 "mm"$ from $1\/4 lambda$ to
$1\/32 lambda$), though see @sec:hyp3-summary for how its aggregate accuracy
compares once every movement type is combined.

== Noisy Fluid Flow <sec:hyp3-fluidflow>

The clean-data fluid-flow experiment of @sec:hyp1-fluidflow
(@fig:h1-ff-setup) is repeated under the same Laplace noise, for a
target directly relevant to the real fluid-injection field data of
@ch:hyp3 --- a graded wetting zone rather than a discrete point scatterer,
swept laterally across the same seven scenarios.

#supp-note[The migration-comparison figure, the Rayleigh-criterion ratio
table and zoomed PSF comparison (compare the clean-data version in
Supplementary Material §S2.4), and the per-scenario phase-plane WLS
front-displacement-error table and shift-estimation diagnostics for this
_noisy_ scenario set are all provided in the Supplementary Material, §S3.4.]

The amplitude ratio behaves similarly to the clean case. Gazdag improves
the most dramatically of any result in this chapter: where
the noisy-envelope-search version of this table left it flat at $-100%$
error from $1\/2 lambda$ down to $1\/32 lambda$ (never resolving any front
displacement at all), it now tracks the front to within $1.50 "mm"$ from
$1\/4 lambda$ downward --- closely matching Kirchhoff. Kirchhoff remains the
single most accurate method overall, within $0.63 "mm"$ from
$1\/2 lambda$ downward. Back-propagation is, for the first time in this
chapter, the *least* accurate of the three from $1\/2 lambda$ down to
$1\/32 lambda$ ($0.52$--$1.27 "mm"$) --- still sub-millimetre-to-low-single-digit
accurate in absolute terms, but the only movement type where both Kirchhoff
and Gazdag now outperform it.

== Results across Movement Types <sec:hyp3-summary>

The noisy lateral worked example (@sec:hyp3-lateral) and the fluid-flow front
(@sec:hyp3-fluidflow) are joined by the vertical and diagonal directions ---
which, under noise, are where the three migration methods separate most clearly
--- before all four are compared side by side.

=== Vertical <sec:hyp3-vertical>

Vertical is the one movement type where the method ranking inverts. Repeating
the depth sweep under noise, Kirchhoff is the most consistent
scenario-by-scenario below $1\/4 lambda$ (within $0.08 "mm"$) with
back-propagation close behind (within $0.82 "mm"$), yet in aggregate mean
absolute error back-propagation is the single most accurate method here
($13.4 "mm"$ against Kirchhoff's $17.5 "mm"$) --- the only movement type for
which it beats Kirchhoff. Gazdag is markedly worse ($-2.93 "mm"$ to a
$+6.95 "mm"$ outlier at $1\/16 lambda$), the cross-axis leakage it assigns to
purely-vertical motion being its defining weakness under noise (Supplementary
Material, §S3.2).

=== Diagonal <sec:hyp3-diagonal>

The combined $2$:$1$ path under noise reproduces the clean-data ranking:
Kirchhoff is the most accurate from Scenario 3 onward (within $0.13 "mm"$ on
both axes), back-propagation close behind (within $1.05 "mm"$), and Gazdag
again carries a persistent tens-of-millimetre error through Scenarios 3--4 (up
to $-41.00 "mm"$ in $#Dx$), only approaching the others at Scenario 5
(Supplementary Material, §S3.3).

=== Cross-Movement Comparison

@tab:h2-lat-phase and @tab:h2-lat-phase-pct give the full per-scenario noisy
phase-plane error for the lateral case; the equivalent per-scenario tables
for Vertical, Diagonal, and FluidFlow are provided in the Supplementary
Material (§S3.2.3, §S3.3.3, §S3.4.3). @fig:h2-detectability is this
chapter's noisy counterpart to @fig:h1-detectability: in every column the
bottom-row phase-error curves still drop below their $5%$ threshold at or
before the top-row amplitude-ratio curves cross below $1$, though noise
pushes both crossing points later, and Gazdag's phase-error curve sits
noticeably higher than the other two methods' throughout.

#page(flipped: true)[
#figure(
  img("H2_028_Hypothesis_2_--_Detectability_Map_Noisy.png"),
  caption: [Detectability map, _noisy_ data (compare @fig:h1-detectability):
    top row, Rayleigh-criterion amplitude ratio (threshold $1$); bottom row,
    absolute phase-plane WLS displacement error as a percentage of the true
    displacement (log scale, threshold $5%$); for all three migration
    methods, across all four movement types. The full per-scenario numeric
    tables underlying this figure are given in
    @tab:h2-lat-amp/@tab:h2-lat-phase (Lateral, in the main text) and the
    Supplementary Material, §S3.2--§S3.4 (Vertical, Diagonal, FluidFlow).],
) <fig:h2-detectability>
]

@fig:h2-dumbbell puts every (movement, method) pair's clean-data and
noisy-data MAE, as a percentage of the true displacement, on one shared log
axis, so @tab:h1-mae and @tab:h2-mae --- and the whole of @ch:hyp1 versus
this chapter --- can be compared directly in one picture. Background
shading groups the twelve rows by movement type.

#page(flipped: true)[
#figure(
  img("H2_029_Hypothesis_2_--_Clean_vs_Noisy_MAE_Dumbbell_Plot.png"),
  caption: [Clean-versus-noisy mean absolute phase-plane displacement error,
    as a percentage of the true displacement, for every movement type and
    migration method: circle marker = clean data (@tab:h1-mae, @ch:hyp1),
    square marker = noisy data (@tab:h2-mae, this chapter), connecting line
    shows the resulting noise-driven degradation. Background colour bands
    group rows by movement type (Lateral, Vertical, Diagonal, FluidFlow).
    Log-scaled horizontal axis.],
) <fig:h2-dumbbell>
]

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Movement*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [Lateral],   [14.381], [22.254], [0.236],
    [Vertical],  [13.410], [34.476], [17.485],
    [Diagonal],  [0.685],  [23.308], [0.065],
    [FluidFlow], [15.758], [0.980],  [0.234],
    table.hline(stroke: 0.4pt),
    [*Mean*],    [*11.059*], [*20.254*], [*4.505*],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Mean absolute phase-plane displacement error [mm] underlying
    @fig:h2-dumbbell's noisy-data series (compare @tab:h1-mae), sub-half-wavelength
    regime only.],
  kind: table,
) <tab:h2-mae>

With the localisation fix of @sec:hyp3-groundtruth-apex applied uniformly to
all three methods, Kirchhoff is now the most accurate method overall
($4.5 "mm"$ mean MAE) --- roughly $2.5 times$ below sign-bit
back-propagation's ($11.1 "mm"$) and $4.5 times$ below Gazdag's
($20.3 "mm"$) --- and the most accurate method for three of the four
movement types (Lateral, Diagonal, FluidFlow); back-propagation remains
most accurate only for Vertical ($13.4 "mm"$ against Kirchhoff's
$17.5 "mm"$).

Locating the target correctly was necessary but not sufficient for Gazdag:
the fix removed most of its excess error (mean MAE $28.9 -> 20.3 "mm"$) but
left a real, movement-specific weakness behind. FluidFlow improves most
dramatically ($25.6 -> 0.98 "mm"$, @sec:hyp3-fluidflow), consistent with a
genuine localisation failure being the dominant error source there;
Vertical, conversely, gets *worse* ($26.0 -> 34.5 "mm"$), driven by spurious
lateral ($#Dx$) error the WLS fit assigns even though Vertical's true
$#Dx = 0$ by construction (Supplementary Material, §S3.2.3) --- a
cross-axis-leakage artefact that a better crop window does not fix, and
whose exact numerical cause remains unresolved (left as outstanding future
work).

Answering Hypothesis 2's question directly: no single method is
unconditionally "most noise-robust" here --- the answer depends on what
"robust" is taken to mean. By raw aggregate accuracy, Kirchhoff wins clearly
--- but @sec:hyp3-purenoise showed Kirchhoff is also the only method that
turns pure noise into coherent, wave-like bands that could be misread as
real structure, a false-positive risk this chapter's MAE metric cannot see
because every scenario it is computed on contains a genuine target.
Back-propagation with sign-bit time-reversal trades some of that raw
accuracy (a factor of $2$--$3$ worse than Kirchhoff on three of four
movement types, though still the best method for Vertical) for staying
diffuse, incoherent speckle on pure noise (@sec:hyp3-signbit) rather than
manufacturing false structure --- the more conservative choice where false
positives, not raw displacement accuracy, are the primary concern. Gazdag,
despite its substantial improvement once correctly localised, remains the
weakest choice on both counts and is not recommended for noisy time-lapse
phase-plane tracking by either criterion. Underlying all three verdicts, the
qualitative conclusion of @ch:hyp1 survives the introduction of noise:
comparing @tab:h2-mae against the clean-data @tab:h1-mae, noise degrades
every method's accuracy by roughly one to two orders of magnitude in the
sub-half-wavelength regime, yet at least one method still recovers every
movement type to a few tenths of a millimetre or better, at displacement
scales where amplitude differencing has already collapsed.
