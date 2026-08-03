#import "../template.typ": *

= Discussion <ch:discussion>

This chapter returns to the Research Question of @ch:introduction --- can
time-lapse GPR accurately track subwavelength movement by analysing phase
changes in migrated images? --- and integrates the results of Hypotheses
1--3 to answer it. Hypothesis 1 is supported without qualification: the
phase-plane fit recovers displacement, on synthetic data, well below the
amplitude-based resolution floor established in @sec:meth-resolution.
Hypothesis 2's original framing --- that a single migration technique is
unconditionally "most noise-robust" --- turns out to be the wrong question;
@ch:hyp2 shows the answer instead depends on what "robust" is taken to mean.
Hypothesis 3 is supported procedurally rather than conclusively: the pipeline
generalises to real, noisy borehole field data and produces stable,
cross-technique-corroborated displacement estimates, but whether those
estimates are physically correct is still open, and generalisation to
multiple simultaneously-moving scatterers was not tested.

== Hypothesis 1: Phase Changes Reveal Sub-Resolution Movement

Across all four @ch:hyp1 experiments --- Lateral, Vertical, Diagonal, and
FluidFlow --- the same pattern holds: amplitude differencing of the migrated
time-lapse images fails once the true displacement drops below roughly
$1\/4$--$1\/2 lambda$ (@tab:h1-lat-amp and its Vertical/Diagonal/FluidFlow
counterparts), while the 2D weighted-least-squares phase-plane fit of
@sec:hyp1-phaseplane remains accurate from that same failure point down to
the smallest scale tested, $1\/32 lambda$. @tab:h1-mae and
@fig:h1-detectability make this directly visible: in every column, the
phase-error curve drops below its $5%$ threshold at or before the
amplitude-ratio curve crosses below $1$, for all three migration methods and
all four movement types. Where the fit is accurate it is accurate to a few
hundredths of a millimetre in the best case (Gazdag and Kirchhoff, Lateral
and Diagonal) --- several orders of magnitude below the wavelength scale.
The larger mean-absolute-error values in @tab:h1-mae are not evidence of a
systematic sub-wavelength failure; they are driven almost entirely by a
single poorly-resolved scenario per (movement, method) pair, close to the
amplitude/phase crossover itself (back-propagation's Lateral fit at
$1\/2 lambda$; the shared Vertical outlier at $1\/2 lambda$ for Gazdag and
Kirchhoff), with every scale below that already accurate to a few tenths of a
millimetre or better. The FluidFlow experiment extends this conclusion beyond
an idealised point scatterer: despite a graded wetting-front target whose
broader intrinsic point-spread function makes amplitude differencing fail
even earlier than for a discrete PEC cylinder (Gazdag's ratio is already
below $1$ at $2 lambda$), the phase-plane fit recovers the front's
displacement to within about a millimetre from $1\/2 lambda$ downward,
confirming that the method's advantage is not an artefact of the idealised
point-scatterer geometry used in the other three experiments.

=== Lateral/Vertical Asymmetry

@sec:th-duality predicted, from first principles, that a migrated image
separates lateral spatial frequency like a prism (@eq:kx-inst) but leaves
vertical spatial frequency at a fixed carrier wavenumber regardless of depth
(@eq:kz-inst) --- so a lateral displacement should collapse the amplitude
image's resolvability sooner than an equivalent vertical displacement. The
amplitude tests bear this out directly: the Lateral Rayleigh ratio
(@tab:h1-lat-amp) already crosses below $1$ for Kirchhoff at $1\/2 lambda$,
whereas the Vertical ratio (@sec:hyp1-vertical) stays above $1$ a full octave
lower, only failing between $1\/4 lambda$ and $1\/8 lambda$. The Diagonal
experiment (@sec:hyp1-diagonal), which combines both components on a fixed
$2$:$1$ path, collapses at roughly the same scale as the pure Lateral case
(Scenario 3, $1\/4 lambda_x$/$1\/8 lambda_z$) rather than at the looser
Vertical floor --- i.e. the lateral component of a combined displacement,
not the vertical one, sets when amplitude differencing gives out. For a
field deployment, this means the direction of an expected displacement
matters for how soon phase-based inference becomes necessary: a process that
is known to be purely vertical (e.g. compaction or settling) remains
amplitude-detectable over a wider range of sub-wavelength scales than one
with any lateral component.

Importantly, this asymmetry is specific to the amplitude-detectability
_floor_, not to the phase estimator's accuracy once that floor is crossed:
@sec:hyp1-vertical notes explicitly that the phase fit is already accurate
($<=0.03 "mm"$) from $1\/4 lambda$ downward in both the lateral and vertical
cases, but at $1 lambda$ and $1\/2 lambda$ --- where amplitude differencing
is still comfortably resolvable in the vertical case --- the phase fit
itself remains wrapped and unreliable regardless of direction. The phase
estimator's advantage is therefore concentrated specifically where amplitude
has already failed, and the lateral/vertical asymmetry predicted by
@sec:th-duality shows up in _where that failure point sits_, not in how well
phase performs once it is reached.

=== Hypothesis 1.5: Local Phase-Gradient Methods

The local, trace-based alternative explored in the Supplementary Material,
§S5 --- tracking the instantaneous phase gradient directly, rather than
fitting a global 2D plane to the cross-spectrum --- can only be reported as
exploratory. Where it has actually been computed, the results are
qualitatively consistent with the global fit and with @sec:th-duality's
predicted asymmetry: the lateral instantaneous-phase cross-section produces
a slope proportional to displacement, while the vertical cross-section
produces a plateau rather than a slope, exactly mirroring @eq:dphi-lateral
and @eq:dphi-vertical. However, three of the four diagnostic views described
alongside it in the Supplementary Material --- the spectral-line fit, the
cross-phase spectrogram, and polar vector rotation --- currently have no
corresponding figures at all, because the notebook sections that generate
them are disabled, and no direct quantitative comparison between any local
method and the global WLS fit has been run on the same dataset. Hypothesis
1.5 is therefore neither confirmed nor refuted by the evidence collected so
far: the one local method with results available (instantaneous-phase
imaging) is consistent with the global fit wherever both exist, but the
claim that local phase gradients offer an "equivalent, simpler" route to the
same displacement estimate remains untested.

== Hypothesis 2: Accuracy and False-Positive Avoidance Favour Different Methods

@ch:hyp2's own framing of Hypothesis 2 already anticipates the result: which
migration algorithm is "best suited" to noise-robust tracking depends on
whether robustness is measured as lowest aggregate displacement error, or as
avoiding manufactured false-positive structure. The two criteria pick
different winners. On raw accuracy, @tab:h2-mae is unambiguous: Kirchhoff's
mean absolute error across all four movement types ($4.5 "mm"$) is
roughly $2.5 times$ below sign-bit back-propagation's ($11.1 "mm"$) and
$4.5 times$ below Gazdag's ($20.3 "mm"$), and Kirchhoff is individually the
most accurate method for three of the four movement types (Lateral,
Diagonal, FluidFlow); back-propagation is more accurate only for Vertical.
On false-positive risk, the pure-noise sanity check of @sec:hyp3-purenoise
reverses the ranking: Kirchhoff's delay-and-sum aperture stacking turns pure
Laplace noise into smooth, wave-like coherent bands that could plausibly be
misread as real layered structure, while Gazdag's frequency-domain
continuation and back-propagation's phase-governed focusing both leave pure
noise as incoherent speckle with no comparable artefact. Sign-bit
time-reversal (@sec:hyp3-signbit) is what makes back-propagation usable
under noise at all: by injecting only the sign of the time-reversed
wavefield rather than its peak-normalised amplitude, it keeps every
zero-crossing of the true signal intact while clamping noise spikes to the
same $plus.minus 1$ amplitude as genuine reflections, removing the outsized
amplitude that would otherwise let a single spike compete with the real
source during back-propagation.

Gazdag is the clear loser on both criteria, though not for a fixed reason
across movement types. Once the target-localisation problem of
@sec:hyp3-groundtruth-apex is corrected (cropping the WLS fit around the
_known_ target position in both $x$ and $z$, rather than an apex hunted for
in the noisy envelope), Gazdag's mean MAE falls from $28.9 "mm"$ to
$20.3 "mm"$, and its FluidFlow error in particular improves by more than an
order of magnitude ($25.6 -> 0.98 "mm"$) --- consistent with localisation,
not a fundamental weakness in the phase-shift operator itself, having been
the dominant error source there. Vertical movement is the exception: Gazdag's
error there gets _worse_ after the same localisation fix ($26.0 ->
34.5 "mm"$), traced to spurious lateral ($#Dx$) error the WLS fit assigns
even though the true $#Dx = 0$ by construction --- a cross-axis leakage that
a better crop window does not resolve and that remains an open,
Gazdag-specific weakness (@sec:hyp3-summary).

=== Implications for Migration Choice in Practice

The recommendation is therefore conditional, not absolute. Where raw
displacement accuracy is the primary objective and the survey volume is
large enough that back-propagation's added cost (a full gprMax forward
simulation per profile pair, versus a comparatively cheap post-processing
step for Kirchhoff or Gazdag) is a real constraint, Kirchhoff's phase-plane
fit is both the cheapest and the most accurate of the three methods tested,
for three of the four movement types --- provided its false-positive risk is
managed by some independent check (e.g. cross-referencing a suspicious
detection against Gazdag or back-propagation, or against the pure-noise
signature of @fig:h2-purenoise-kg) rather than trusted blindly. Where
avoiding a manufactured false detection matters more than raw accuracy ---
for instance, a first-pass anomaly screen in a monitoring context where a
false positive is costly to chase down --- sign-bit back-propagation is the
more conservative choice, and it is the single most accurate method
available specifically for vertical-only displacement monitoring. Gazdag is
not recommended under noise by either criterion: its aggregate error remains
the worst of the three even after the localisation fix, and its unresolved
cross-axis leakage on purely vertical motion is a specific, uncorrected
failure mode. This recommendation is based on the single noise level actually
tested here ($10%$ of each B-scan's own signal standard deviation,
@sec:hyp3-laplace); a systematic sweep across multiple noise levels, which
would locate any crossover between Kirchhoff's accuracy advantage and
back-propagation's false-positive-avoidance advantage, has not been run and
is identified as future work below.

== Hypothesis 3: Generalisation to Field Data

@sec:hyp3-fielddata applied the full pipeline to 38 real, zero-offset
borehole GPR profiles from a controlled fluid-injection experiment, where ---
unlike every @ch:hyp1/@ch:hyp2 experiment --- neither the target geometry nor
the true displacement is known in advance. The region of influence is
therefore chosen entirely by hand with the napari image viewer
(@sec:hyp3-fd-roi) --- painting on the wavenumber-domain cross-spectrum or
the spatial difference B-scan --- rather than by an automatic rectangular or
sliding-window search. A RANSAC check confirms this is largely trustworthy:
a direct cross-spectrum pick reproduces the WLS estimate to within about
$6%$, while an amplitude-gated pick from a painted B-scan region shifts it
by up to $59%$ in the worst case, though direction is always preserved ---
a real source of estimate uncertainty absent from the synthetic studies,
where the crop is instead centred on known ground truth.

With the ROI fixed, the stage-anchored estimates of @tab:fielddata-stages
(@sec:hyp3-fd-phaseplane) tell a consistent story, but only partly the naive
one: the dominant signal is a *downward* displacement of $1.66 "m"$
(Gazdag) / $2.08 "m"$ (back-propagation) during Push, with a smaller lateral
shift *toward* the borehole ($-0.17$/$-0.25 "m"$) --- matching the naive
vertical expectation for active injection, but not the lateral one. Chase
adds a comparable downward increment ($+1.62$/$+1.65 "m"$) while the lateral
shift reverses to move away from the borehole ($+0.07$/$+0.06 "m"$); Wait
substantially reverses the accumulated downward movement ($-0.67$/
$-0.28 "m"$, i.e. upward); and Pull is where the two techniques part ways,
Gazdag settling at the noise floor ($+0.00 "m"$) while back-propagation
continues a smaller upward reversal ($-0.24 "m"$), leaving a net residual of
$#Dz approx +2.61 "m"$ (downward) and $#Dx approx -0.10 "m"$ (toward the
borehole) at profile 38. Back-propagation, sharing no processing steps with
the Gazdag branch, agrees with it on direction and roughly on magnitude
during Push, Chase and Wait, diverging only in Pull, where Gazdag sits at
the noise floor rather than truly disagreeing in sign --- evidence that
neither the vertical reversal nor the toward-the-borehole lateral shift is a
processing artefact.

Whether that lateral shift is _physically_ expected, however, is a question
this thesis cannot answer: @sec:hyp3-fd-interpretation leaves it open
pending independent ground truth (injection depth, volume, fracture
geometry) not available during processing. Hypothesis 3 is therefore best
read as partially, not fully, supported: the pipeline generalises to real,
noisy, geometry-unknown field data procedurally, producing stable,
cross-technique-corroborated estimates, but whether those estimates are
physically correct remains unvalidated. The complementary test proposed
alongside it --- a complex synthetic scene with multiple
independently-moving scatterers under noise --- was not run and is left to
future work.

== Limitations

+ Hypothesis 1.5's local phase-gradient methods (Supplementary Material,
  §S5) remain exploratory: three of the four diagnostic views have no corresponding
  results because the notebook cells that generate them are currently
  disabled, and no quantitative comparison against the global WLS fit has
  been run on shared data.

+ The clean-versus-noisy comparison of @ch:hyp1/@ch:hyp2 relies on cropping
  the phase-plane fit around the _known_ target position under noise
  (@sec:hyp3-groundtruth-apex) --- legitimate for synthetic data with exact
  ground truth, but unavailable on real field data, where hand-painted ROI
  selection (@sec:hyp3-fd-roi) was needed instead: a RANSAC check shows this
  introduces an estimate spread of a few percent up to $59%$ in the worst
  case depending on how the ROI is picked, with no equivalent uncertainty in
  the synthetic study.

+ The migration-technique comparison of @ch:hyp2 is based on a single fixed
  noise level ($10%$ of each B-scan's own signal standard deviation,
  @sec:hyp3-laplace); no sweep across multiple noise levels has been run, so
  it is not known whether, or at what noise level, Kirchhoff's raw-accuracy
  advantage and back-propagation's false-positive-avoidance advantage trade
  places.

+ Field-data back-propagation coverage is partial: Kirchhoff and Gazdag are
  available for all 37 processed profiles, but the explicit-borehole-geometry
  pipeline ($d_x = 0.02 "m"$, @sec:hyp3-fd-crossprofile) --- the only
  back-propagation configuration now used in this thesis --- has been
  validated on only the five representative profiles (1, 3, 8, 20, 38) used
  throughout @ch:hyp3, not the full dataset.

+ The complex synthetic scene with multiple, independently-moving scatterers
  under noise was not run, so simultaneous or interacting displacements ---
  plausibly closer to a real field scenario than any single-scatterer
  experiment in this thesis --- remain untested even synthetically.

+ The field-data displacement direction reported in @tab:fielddata-stages is
  internally consistent and corroborated across two independently-processed
  migration techniques, but its physical interpretation is pending
  confirmation against independent field-survey metadata not available
  during this work (@sec:hyp3-fd-interpretation).

+ Gazdag's noise sensitivity is only partially diagnosed: correcting the
  target-localisation crop window resolved most of its excess error, but a
  residual cross-axis leakage on purely vertical motion remains unexplained,
  and an earlier hypothesis attributing a related streaking artefact to the
  phase-shift depth-stepping operator itself has not been re-tested since the
  localisation fix (@sec:hyp3-summary).

+ The velocity-recalibration problem flagged in @sec:meth-phaseplane --- a
  genuine change in soil or ice moisture between baseline and monitor
  surveys is, without recalibrating against a known static reflector,
  indistinguishable from a spurious vertical shift $#Dz$ --- is not tested by
  any experiment in this thesis, synthetic or field.

+ The automatic dip-filter orientation estimator used in the back-propagation
  field pipeline's post-imaging clutter removal did not generalise once the
  timing and normalisation fixes of @sec:hyp3-fd-crossprofile were applied,
  and was replaced by a single manually-specified value shared across all
  five profiles; the estimator's robustness on a larger profile set is not
  otherwise established.

== Outlook

+ Extend the corrected, explicit-borehole-geometry back-propagation pipeline
  from the five representative profiles used in this thesis to the full
  38-profile field dataset, and confirm that it continues to agree with the
  Gazdag trajectory over the complete Push/Chase/Wait/Pull sequence rather
  than only the four representative pairs checked here.

+ Run a complex synthetic experiment to test whether
  the phase-plane fit, or a spatially-localised variant of it, can separate
  multiple simultaneously-moving scatterers under noise --- the natural next
  step between the idealised single-scatterer studies of @ch:hyp1/@ch:hyp2
  and the uncontrolled real field data of @sec:hyp3-fielddata.

+ Obtain independent field metadata (injection depth and volume, fracture
  and borehole geometry) and cross-check it against the field-data
  displacement estimates of @tab:fielddata-stages, to resolve whether the
  observed toward-the-borehole lateral shift during Push is physically
  expected for this experiment or points to a coordinate-convention issue not
  yet identified.

+ Run a systematic sweep across multiple noise levels, rather than the
  single $10%$ level tested in @ch:hyp2, to locate any crossover between
  Kirchhoff's raw-accuracy advantage and back-propagation's
  false-positive-avoidance advantage, and to test whether Gazdag's
  cross-axis leakage on vertical motion is the same numerical artefact
  suspected earlier in the phase-shift operator or a distinct effect.

+ Complete the quantitative validation of Hypothesis 1.5 by re-enabling the
  disabled processing steps that generate the missing spectral-line,
  cross-phase-spectrogram, and localised-STFT figures, and run a direct,
  same-dataset comparison between the local phase-gradient methods and the
  global WLS fit.

+ Quantify the velocity-recalibration procedure needed when soil or ice
  moisture genuinely changes between a baseline and monitor survey, so that
  a true velocity change can be distinguished from a spurious vertical shift
  $#Dz$ in field deployment.

+ Propagate the ROI/fitting-method sensitivity found in the field-data study
  (up to $59%$ in the worst case, @sec:hyp3-fd-roi) into a reported
  uncertainty band alongside the headline displacement numbers of
  @tab:fielddata-stages, rather than a single point estimate per stage.
