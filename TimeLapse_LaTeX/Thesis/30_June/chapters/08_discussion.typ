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
cross-migration-technique consistent displacement estimates.

== Discussion Regarding Hypothesis 1: Phase Changes Reveal Sub-Resolution Movement

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
even earlier than for a discrete PEC cylinder (Gazdag's FWHM-ratio is already
below $1$ at $2 lambda$), the phase-plane fit recovers the front's
displacement to within about a millimetre from $1\/2 lambda$ downward,
confirming that the phase-plane fit's sub-wavelength accuracy is not a
special case of the discrete point-scatterer targets used in the Lateral,
Vertical and Diagonal experiments, but extends to a spatially-extended,
graded target, like a fluid front, as well.

#pagebreak()

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

// === Hypothesis 1.5: Local Phase-Gradient Methods

// The local, trace-based alternative explored in the Supplementary Material,
// §S5 --- tracking the instantaneous phase gradient directly, rather than
// fitting a global 2D plane to the cross-spectrum --- can only be reported as
// exploratory. Where it has actually been computed, the results are
// qualitatively consistent with the global fit and with @sec:th-duality's
// predicted asymmetry: the lateral instantaneous-phase cross-section produces
// a slope proportional to displacement, while the vertical cross-section
// produces a plateau rather than a slope, exactly mirroring @eq:dphi-lateral
// and @eq:dphi-vertical. However, three of the four diagnostic views described
// alongside it in the Supplementary Material --- the spectral-line fit, the
// cross-phase spectrogram, and polar vector rotation --- currently have no
// corresponding figures at all, because the notebook sections that generate
// them are disabled, and no direct quantitative comparison between any local
// method and the global WLS fit has been run on the same dataset. Hypothesis
// 1.5 is therefore neither confirmed nor refuted by the evidence collected so
// far: the one local method with results available (instantaneous-phase
// imaging) is consistent with the global fit wherever both exist, but the
// claim that local phase gradients offer an "equivalent, simpler" route to the
// same displacement estimate remains untested.

== Discussion Regarding Hypothesis 2: Accuracy and False-Positive Avoidance Favour Different Methods

@ch:hyp2's own framing of Hypothesis 2 already anticipates the result: which
migration algorithm is "best suited" to noise-robust tracking depends on
whether robustness is measured as lowest displacement error, or as
avoiding migration artefacts that mimic a genuine, high-amplitude reflection
where none exists. The two criteria pick different winners. 

#linebreak()

On raw accuracy, @tab:h2-mae is unambiguous: Kirchhoff's
mean absolute error across all four movement types ($4.5 "mm"$) is
roughly $2.5 times$ below sign-bit back-propagation's ($11.1 "mm"$) and
$4.5 times$ below Gazdag's ($20.3 "mm"$), and Kirchhoff is individually the
most accurate method for three of the four movement types (Lateral,
Diagonal, FluidFlow); back-propagation is more accurate only for Vertical.
On false-positive risk, the pure-noise sanity check of @sec:hyp3-purenoise
reverses the ranking: Kirchhoff's delay-and-sum aperture stacking turns pure
Laplace noise into exactly this kind of artefact --- clustered noise spikes
of high amplitude that could plausibly be misread as a genuine reflection ---
while Gazdag's frequency-domain
continuation and back-propagation's sign-bit phase-governed focusing both leave pure
noise as incoherent speckle with no comparable artefact. Sign-bit
time-reversal (@sec:hyp3-signbit) is what makes back-propagation usable
under noise at all: by injecting only the sign of the time-reversed
wavefield rather than its peak-normalised amplitude, it keeps every
zero-crossing of the true signal intact while clamping noise spikes to the
same $plus.minus 1$ amplitude as genuine reflections, removing the outsized
amplitude that would otherwise let a single spike compete with the real
source during back-propagation.

#linebreak()

Gazdag is the clear loser on both criteria. @fig:h2-dumbbell visualises this
directly: plotting every (movement, method) pair's clean- and noisy-data MAE
on the same log axis, Gazdag shows by far the largest clean-to-noisy
degradation of the three methods --- most dramatically for Lateral and
Diagonal, where its clean-data error is essentially zero ($0.000 "mm"$ and
$0.003 "mm"$ respectively, @tab:h1-mae) but its noisy-data error exceeds
$20 "mm"$ on both (@tab:h2-mae), the longest connecting lines in the figure.
Vertical is Gazdag's single worst case in absolute terms in either regime
($16.1 -> 34.5 "mm"$), traced to spurious lateral ($#Dx$) error the WLS fit
assigns even though the true $#Dx = 0$ by construction --- a cross-axis
leakage that remains an open, Gazdag-specific weakness (@sec:hyp3-summary).

=== Implications for Migration Choice in Practice

The recommendation is therefore conditional, not absolute. Where raw
displacement accuracy is the primary objective, and back-propagation's added
cost (a full gprMax forward simulation per profile pair, versus a
comparatively cheap post-processing step for Kirchhoff or Gazdag) is a real
constraint, Kirchhoff's phase-plane fit is the cheapest and most accurate
method for three of the four movement types. Its false-positive risk ---
noise clusters that can look like genuine reflections --- should be checked
independently, for example against Gazdag, back-propagation, or the
pure-noise signature of @fig:h2-purenoise-kg, rather than trusted on its own.
Where avoiding a false detection matters more than raw accuracy, sign-bit
back-propagation is the safer choice; it is also the most accurate method
for vertical-only displacement. Gazdag is
not recommended under noise by either criterion: its MAE remains
the worst of the three, and its unresolved
cross-axis leakage on purely vertical motion is a specific, uncorrected
failure mode. This recommendation is based on the single noise level actually
tested here ($10%$ of each B-scan's own signal standard deviation,
@sec:hyp3-laplace); a systematic sweep across multiple noise levels, which
would locate any crossover between Kirchhoff's accuracy advantage and
back-propagation's false-positive-avoidance advantage, has not been run and
is identified as future work below.

== Discussion Regarding Hypothesis 3: Generalisation to Field Data

@sec:hyp3-fielddata applied the full pipeline to real, zero-offset
borehole GPR profiles from a controlled fluid-injection experiment, where ---
unlike every @ch:hyp1/@ch:hyp2 experiment --- neither the target geometry nor
the true displacement is known in advance. The region of influence is
therefore chosen entirely by hand with the napari image viewer
(@sec:hyp3-fd-roi) --- painting on the wavenumber-domain cross-spectrum or
the spatial difference B-scan --- rather than by an automatic rectangular or
sliding-window search. A RANSAC check confirms this is largely trustworthy:
refitting each estimate with RANSAC, which discards outlier cells before
re-running WLS on the inliers alone, changes a direct cross-spectrum pick's
displacement estimate by at most about $6%$, while it changes an
amplitude-gated pick's (from a painted B-scan region) by up to $59%$ in the
worst case --- direction is preserved either way, but the size of the shift
shows how sensitive each picking strategy is to noisy or wrapped cells
slipping into the WLS mask. This is a real source of estimate uncertainty
absent from the synthetic studies, where the crop is instead centred on
known ground truth.

#linebreak()

With the ROI fixed, the stage-anchored estimates of @tab:fielddata-stages
(@sec:hyp3-fd-phaseplane) tell a consistent story that matches the injection
geometry once it is accounted for: fluid is injected at $78.7 "m"$ depth,
inside the same $70$--$80 "m"$ band as the tracked reflector, so an
*upward*-moving front while fluid is actively pumped in (Push, Chase) and a
downward settling once pumping stops (Wait, Pull) is the physically expected
pattern. The dominant signal is an upward displacement of $1.66 "m"$
(Gazdag) / $2.06 "m"$ (back-propagation) during Push, with a smaller lateral
shift *toward* the borehole ($0.17$/$0.25 "m"$). Chase adds a comparable
upward increment ($1.60$/$1.65 "m"$) while the lateral shift reverses to
move away from the borehole ($0.07$/$0.06 "m"$); Wait substantially reverses
the accumulated upward movement ($+0.68$/$+0.28 "m"$, i.e. downward); and
Pull continues that smaller downward reversal ($+0.14 "m"$ for Gazdag,
$+0.24 "m"$ for back-propagation), leaving a net residual at profile 38 of
$#Dz approx -2.44 "m"$ (upward) and $#Dx approx +0.22 "m"$ (away from the
borehole) for Gazdag, and $#Dz approx -3.19 "m"$, $#Dx approx -0.30 "m"$ for
back-propagation. Back-propagation, sharing no processing steps with the
Gazdag branch, agrees with it on the sign of the displacement in every stage
--- evidence that neither the vertical reversal nor the toward-the-borehole
lateral shift during Push is a processing artefact, but a real feature of the field data. The two methods disagree on the lateral direction of the final residual, but the magnitude is small and the sign is sensitive to the painted ROI, so this is not a strong contradiction.

#linebreak()

The lateral shift is the one part of this picture without a direct
geometric explanation from the injection depth alone; @sec:hyp3-fd-interpretation
suggests it may instead reflect the fracture network's orientation, which
need not run perfectly parallel to the borehole, but this remains
unconfirmed against independent ground truth. Hypothesis 3 is therefore
best read as partially, not fully, supported: the pipeline generalises to
real, noisy, geometry-unknown field data procedurally, producing stable,
cross-migration-technique confirmed displacement estimates whose dominant vertical
component now has a physical explanation, but whose lateral component's
physical origin remains unconfirmed.

== Limitations

*Unexplored extensions*:
Several extensions of the core method remain only partially explored.
Hypothesis 1.5's local phase-gradient methods (Supplementary Material, §S5)
remain exploratory, and the migration-technique comparison of @ch:hyp2 is
based on a single fixed noise level ($10%$ of each B-scan's own signal
standard deviation, @sec:hyp3-laplace); no sweep across multiple noise
levels has been run, so it is not known whether, or at what noise level,
Kirchhoff's raw-accuracy advantage and back-propagation's
false-positive-avoidance advantage trade places. The velocity-recalibration
problem flagged in @sec:meth-phaseplane --- a genuine change in soil or ice
moisture between baseline and monitor surveys is, without recalibrating
against a known static reflector, indistinguishable from a spurious
vertical shift $#Dz$ --- is likewise untested by any experiment in this
thesis, synthetic or field.

#linebreak()

*Region-of-influence limitations*:
A second group of limitations centres on how the region of influence is
chosen. The clean-versus-noisy comparison of @ch:hyp1/@ch:hyp2 relies on
cropping the phase-plane fit around the _known_ target position under noise
(@sec:hyp3-groundtruth-apex) --- legitimate for synthetic data with exact
ground truth, but unavailable on real field data. The complex synthetic
scene with multiple, independently-moving scatterers under noise was not
run either. Simultaneous displacements --- plausibly
closer to a real field scenario than any single-scatterer experiment --- remain untested even synthetically. In the wavenumber
domain, each scatterer imprints its own linear phase ramp on the
cross-spectrum (@eq:phase-plane), and where two or more scatterers move in
different directions, their overlapping spectral energy superimposes these
differently-tilted ramps rather than adding cleanly. A global WLS fit
would return some noise-sensitive, energy-weighted blend of the true
displacement vectors unless the region of influence is first restricted. This should be done via the envelope of the difference B-scan, picking one coherent event at a time.
This underlines how much a good region of influence matters in practice.
Gazdag's noise sensitivity is only partially diagnosed, with a residual
cross-axis leakage on purely vertical motion left unexplained, and more
generally the envelope of the difference B-scan becomes harder to read
cleanly under noise. An
automatic rule that simply picks the
highest-amplitude region is not trustworthy. @sec:hyp3-purenoise shows
Kirchhoff can turn pure noise into a coherent, high-amplitude cluster that
such a rule would pick as a genuine target. Better, more
robust automatic ROI selection under noise remains an open problem for
future work.

#linebreak()

*Field-data trajectory granularity*:
Finally, the field-data trajectory (@ch:hyp3) is estimated stage-by-stage
rather than step-by-step: each stage's displacement comes from a single WLS
fit between the two profiles at its endpoints (e.g. profile 8 to profile 20
for Wait, @sec:hyp3-fd-phaseplane), not from accumulating a separate WLS
estimate between every consecutive profile pair (1→2, 2→3, $dots$, 37→38)
into one global trajectory. The latter would resolve the fluid front's
motion at the full resolution of the survey rather than in four coarse
steps, but was not attempted here and is left to future work.


== Outlook

*Direct follow-up experiments*:
Several of these open questions have a direct next experiment. Completing
the quantitative validation of Hypothesis 1.5 means exploring the
spectral-line, cross-phase-spectrogram, and localised-STFT figures, and
running a direct, same-dataset comparison against the global WLS fit ---
not only to check that the local methods agree with it, but to test
whether they add any genuinely useable information the global fit does not
already provide (e.g. finer spatial localisation of the displacement, or
robustness in regions too small or irregular for a reliable plane fit). A
systematic sweep across multiple noise levels, rather than the single
$10%$ level tested in @ch:hyp2, would locate any crossover between
Kirchhoff's raw-accuracy advantage and back-propagation's advantage at
avoiding false positives --- back-propagation's time-reversal focusing
concentrates energy at the true scatterer location rather than clustering
noise into high-amplitude regions the way Kirchhoff's delay-and-sum
stacking does. The velocity-recalibration procedure needed when soil or ice moisture
genuinely changes between a baseline and monitor survey should also be
quantified, so that a true velocity change can be distinguished from a
spurious vertical shift $#Dz$ in field deployment.

#linebreak()

*Testing beyond a single scatterer*:
A second priority is testing the method beyond a single, isolated
scatterer. Running a complex synthetic experiment to test whether the
phase-plane fit, or a spatially-localised variant of it, can separate
multiple simultaneously-moving scatterers under noise is the natural next
step between the idealised single-scatterer studies of @ch:hyp1/@ch:hyp2
and the uncontrolled real field data of @sec:hyp3-fielddata. The
ROI/fitting-method sensitivity found in the field-data study should also be
turned into a properly reported uncertainty band: picking a different
region of influence --- whether in the wavenumber domain or on the
difference B-scan --- measurably shifts the recovered $#Dz$/$#Dx$ estimate
(@sec:hyp3-fd-interpretation), so each stage's headline number in
@tab:fielddata-stages should be reported as a sensitivity range around it,
not a single point estimate.

#linebreak()

*Field-data trajectory extension*:
Finally, the corrected, explicit-borehole-geometry back-propagation
pipeline should be extended from the five representative profiles used in
this thesis to the full 38-profile field dataset, confirming that it
continues to agree with the Gazdag trajectory over the complete
Push/Chase/Wait/Pull sequence rather than only the four representative
pairs checked here. Doing so would also let the trajectory be accumulated
stepwise from consecutive profile pairs (1→2, 2→3, $dots$, 37→38), rather than as a single
stage-anchored estimate per stage (push, chase, wait, pull), so the two approaches can be directly
compared.
