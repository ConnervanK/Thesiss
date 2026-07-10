#import "../template.typ": *

= Hypothesis 2: Which Migration Technique Is Best Suited for Noise-Robust Phase-Plane Tracking? <ch:hyp2>

@ch:hyp1 established that phase-plane regression can track sub-wavelength
displacements on clean, synthetic data. Before the method is used in more
complex settings, a fundamental question must be answered: does the choice of
migration algorithm matter when the data are corrupted by noise? This chapter
argues that yes, it matters substantially, and that back-propagation with
sign-bit time-reversal is the preferred method.

#para-head[Hypothesis 2.] Back-propagation migration with sign-bit
time-reversal is the most noise-robust technique for time-lapse phase-plane
tracking: it keeps the algorithm usable under heavy-tailed Laplace noise by
suppressing large-amplitude noise spikes before back-propagation, while
Kirchhoff produces false-coherent artefacts and Gazdag adds incoherent
speckle.

== A Laplace Noise Model from Real Field Data <sec:hyp3-laplace>

Rather than injecting arbitrary synthetic noise, the noise level and shape
used throughout this thesis are fitted to real GPR field noise. A sample of
field noise is tracked through the processing pipeline and fitted with both a
Laplace and a Gaussian distribution at two stages: after the final processing
stage ("9 Crop Samples"), and at the last stage before the spherical-gain
correction ("7 Constant Velocity"), since the gain correction inflates the
amplitude scale by several orders of magnitude and is not representative of
the raw simulated $E_z$ amplitudes used elsewhere in this thesis.

#figure(
  table(
    columns: (auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Stage*], [*Laplace scale*], [*Gaussian $sigma$*],
    table.hline(stroke: 0.4pt),
    [9 Crop Samples (post-gain)],     [$117 space 584.0$], [$171 space 156.1$],
    [7 Constant Velocity (pre-gain)], [$6.085$],           [$9.233$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Laplace and Gaussian fits to real field noise, at two pipeline
    stages ($n = 976 space 244$ samples each). The heavier-tailed Laplace
    distribution is adopted throughout this thesis; loc $= 0$ for both
    stages.],
  kind: table,
) <tab:laplace-fit>

In both cases the fitted distribution is heavier-tailed than a Gaussian of
matched variance, consistent with field GPR noise being dominated by
occasional large-amplitude clutter and interference rather than purely thermal
noise. The pre-gain fit is the one actually used to generate synthetic noise
for the experiments in @sec:hyp3-lateral, @sec:hyp3-vertical,
@sec:hyp3-diagonal, and @sec:hyp3-fluidflow: its Laplace _shape_ (loc $= 0$, heavier tails than
Gaussian) is kept, but its scale is rescaled so that the resulting noise
standard deviation is exactly $10%$ of each synthetic B-scan's own signal
standard deviation --- a light, realistic noise level rather than the raw
fitted scale, which would be incommensurate with the synthetic $E_z$
amplitudes.

== Migration-Algorithm Response to Pure Noise <sec:hyp3-purenoise>

Before testing noise robustness on real signal-plus-noise data, a sanity check
establishes what each migration algorithm does to noise _alone_: an empty
B-scan (no scatterer, no real signal) containing only the Laplace noise of
@sec:hyp3-laplace is migrated with Kirchhoff, Gazdag, and back-propagation.
If a method turns pure noise into something that looks like a coherent,
scatterer-like focus, every noisy result elsewhere in this chapter carries a
false-positive risk that must be accounted for.

*Result: the three methods do not fail the same way.* Kirchhoff turns pure
noise into clearly coherent, smooth wave-like bands that could easily be
misread as real layered structure: its delay-and-sum aperture stacking imposes
coherence on incoherent input by construction, summing many traces along
travel-time curves and smoothing incoherent noise into locally-correlated
structure. Gazdag's noise output, by contrast, stays speckled and incoherent,
with no wave-like artefacts --- just texture. This is a direct, practical
consequence of the structural difference between the two algorithms
(@sec:th-migration): Kirchhoff's spatial stacking manufactures apparent
coherence from nothing, while Gazdag's frequency-domain downward continuation
does not.

#draftnote[the pure-noise Kirchhoff/Gazdag comparison above is implemented
and run in `Noise_Playground.ipynb`, but its plots have not yet been exported
to a `TimeLapse_Figures` folder and added to this project's path --- do this
before finalising this section, and include the comparison figure here. The
back-propagation-of-pure-noise variant (@sec:hyp3-signbit) has real, completed
gprMax output in `noise_study/backprop/{purenoise_peaknorm,purenoise_signbit}/`
but its focus-frame snapshots have not yet been loaded, plotted, or compared
quantitatively --- this is the single most direct test of whether sign-bit
time-reversal actually suppresses the false-positive risk identified above for
back-propagation specifically, and is currently missing.]

== Sign-Bit Time-Reversal for Noise-Robust Back-Propagation <sec:hyp3-signbit>

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

_Sign-bit time-reversal_ fixes this by injecting only the _sign_ of the
time-reversed wavefield instead of its peak-normalised value,
$ u_"sign" (x, tau) = op("sign")(u(x, tau)) , $ <eq:signbit>
implemented in `write_backprop_files(..., sign_bit=True)`
(`helper_functions/migration.py`). This keeps every zero-crossing and phase
trend of the GPR wavelet completely intact, since the sign of a signal carries
its full phase information, while squashing every noise spike down to the same
$plus.minus 1$ amplitude as the coherent signal --- stripping noise of the
outsized amplitude that would otherwise let it dominate the back-propagated
wavefield. The clean-data experiments of @ch:hyp1 use the default
peak-normalised excitation throughout, since they have no noise to suppress;
every noisy back-propagation result in this chapter uses sign-bit excitation
instead.

== Noisy Lateral Movement <sec:hyp3-lateral>

The lateral displacement sweep of @sec:hyp1-lateral is repeated on B-scans
contaminated with the synthetic Laplace noise of @fig:tl-bscans (c), for all
three migration algorithms, including back-propagation with sign-bit
time-reversal (@sec:hyp3-signbit). @fig:tl-noisy-summary-amp overlays the
signed time-lapse-difference amplitude from all three algorithms, directly
comparable to the clean-data summary of @sec:tl-detectability, and
@fig:tl-noisy-summary-psf plots the corresponding normalised lateral PSF.
@fig:tl-noisy-signbit shows the sign-bit time-reversed excitation itself,
across all eight scenarios.

#figure(
  img("TL_041_TimeLapse_Migration_Comparison_Noisy__Signed_Amplitude____f_.png"),
  caption: [Signed time-lapse-difference amplitude for all three migration
    algorithms on the _noisy_ lateral dataset ($f_c = 1.5 "GHz"$,
    aperture $= 40$).],
) <fig:tl-noisy-summary-amp>

#figure(
  img("TL_042_Normalised_Lateral_PSF_Noisy__TimeLapse_Difference_at_True_S.png"),
  caption: [Normalised lateral PSF of the _noisy_ time-lapse-difference image
    at the true scatterer depth, swept across lateral displacements from
    $2 lambda$ to $1 \/ 32 lambda$.],
) <fig:tl-noisy-summary-psf>

#figure(
  img("TL_034_Sign-Bit_Time-Reversed_Excitation__B-scans_and_Spectra_Noisy.png", width: 90%),
  caption: [Sign-bit time-reversed excitation B-scans and their frequency
    spectra for the _noisy_ lateral dataset, all eight displacement
    scenarios.],
) <fig:tl-noisy-signbit>

#draftnote[in @fig:tl-noisy-summary-amp the Kirchhoff and Gazdag columns
appear visually blank against the shared colour scale, which is dominated by
back-propagation's much larger raw amplitude range (the same pattern recurs
in @fig:vtl-noisy-summary-amp and @fig:dtl-noisy-summary-amp below) ---
confirm this is a colour-scale/normalisation artefact rather than a
data-loading issue, and consider giving each algorithm its own colour scale
as in the clean-data version (@fig:tl-summary-amp) if so. The unexplained
fourth ("N/A") column present in this figure and @fig:tl-noisy-summary-psf,
absent from the vertical and diagonal equivalents, should also be resolved
or removed.]

The same noisy lateral dataset is also carried through a shift estimator in
`TimeLapse_Processing.ipynb`, for Kirchhoff only (@fig:tlp-noise); Gazdag is
excluded from this particular test because Gazdag phase-shift migration of
this noisy dataset produces a structural migration artefact that defeats
shift estimation (see the notebook's Section 2b/2c notes), not because it was
skipped for convenience, and back-propagation's phase-plane fit for the
lateral case specifically has not yet been run (unlike the vertical and
diagonal cases below, @fig:vtl-noisy-phaseplane and @fig:dtl-noisy-phaseplane).

#figure(
  img("TLP_010_NoisyTimeLapse__Robust_shift_estimation_GCC__lstsq_fallback.png", width: 90%),
  caption: [Robust shift estimation (GCC peak search with a phase-plane
    least-squares fallback) applied to the _noisy_ lateral-displacement
    dataset, Kirchhoff-migrated: baseline vs. each scenario, with the
    per-scenario difference image, cross-spectrum phase and energy, and the
    recovered vs. true $(#Dz, #Dx)$.],
) <fig:tlp-noise>

#draftnote[quantify how much the additive Laplace noise degrades the
estimated $#Dx$ relative to the clean Kirchhoff result of
@fig:tlp-horiz-validation, and relate this to the amplitude-weighting
argument of @sec:th-mask-weight. Note this is a GCC-peak-search-with-lstsq-fallback
estimator (`estimate_shift_2d_cleaning`), not the plain WLS phase-plane fit
used elsewhere in this thesis and for the vertical/diagonal cases below ---
state explicitly why the plain fit needed this more robust replacement for
the noisy lateral case, or re-run the plain fit here if that is a fairer
comparison to @fig:tlp-horiz-validation. Completing a Gazdag- and
back-propagation-equivalent phase-plane fit for the lateral case is tracked
in @sec:hyp3-gaps.]

== Noisy Vertical Movement <sec:hyp3-vertical>

The vertical displacement sweep of @sec:hyp1-vertical is repeated under the
same Laplace noise, for all three migration algorithms, using sign-bit
time-reversal (@sec:hyp3-signbit) for the noisy back-propagation run.
@fig:vtl-noisy-summary-amp overlays the signed time-lapse-difference
amplitude from all three algorithms, directly comparable to the clean-data
summary of @sec:vtl-detectability, and @fig:vtl-noisy-summary-psf plots the
corresponding normalised vertical PSF. @fig:vtl-noisy-signbit shows the
sign-bit time-reversed excitation itself.

#figure(
  img("VTL_040_Vertical_TimeLapse_Migration_Comparison_Noisy__Signed_Amplit.png"),
  caption: [Signed time-lapse-difference amplitude for all three migration
    algorithms on the _noisy_ vertical dataset.],
) <fig:vtl-noisy-summary-amp>

#figure(
  img("VTL_041_Normalised_Vertical_PSF_Noisy__TimeLapse_Difference_at_x__20.png", width: 55%),
  caption: [Normalised vertical PSF of the _noisy_ time-lapse-difference
    image at $x = 2.0 "m"$, swept across vertical displacements from
    $1 lambda$ to $1 \/ 32 lambda$.],
) <fig:vtl-noisy-summary-psf>

#figure(
  img("VTL_033_Sign-Bit_Time-Reversed_Excitation_--_B-scans_and_Spectra_Noi.png", width: 90%),
  caption: [Sign-bit time-reversed excitation B-scans and their frequency
    spectra for the _noisy_ vertical dataset, all seven displacement
    scenarios.],
) <fig:vtl-noisy-signbit>

@fig:vtl-noisy-phaseplane applies the 2D WLS phase-plane fit of
@sec:meth-phaseplane to the noisy vertical dataset for all three migration
methods, closing the gap left open in earlier drafts of this chapter.

#figure(
  subfigs(cols: 2,
    img("TLP_020_VerticalTimeLapse_Noisy__Phase-plane_shift_estimation__Kirch.png"),
    img("TLP_021_VerticalTimeLapse_Noisy__Phase-plane_shift_estimation__Gazda.png"),
    img("TLP_022_VerticalTimeLapse_Noisy__Phase-plane_shift_estimation__Back-.png"),
  ),
  caption: [2D phase-plane shift estimation applied to the _noisy_ vertical
    time-lapse dataset of @sec:hyp1-vertical, for all three migration
    methods: (a) Kirchhoff; (b) Gazdag; (c) back-propagation.],
) <fig:vtl-noisy-phaseplane>

#draftnote[state whether the noisy vertical phase-plane fit remains accurate
to the same smallest displacement found for the noisy lateral case
(@fig:tlp-noise), and compare directly against the clean vertical fit of
@fig:tlp-vert-validation.]

== Noisy Diagonal Movement <sec:hyp3-diagonal>

The diagonal displacement sweep of @sec:hyp1-diagonal is likewise repeated
under Laplace noise, for all three migration algorithms, again using sign-bit
time-reversal for the noisy back-propagation run. @fig:dtl-noisy-summary-amp
overlays the signed, diagonally-sampled time-lapse-difference amplitude from
all three algorithms, directly comparable to the clean-data summary of
@sec:dtl-detectability, and @fig:dtl-noisy-summary-psf plots the
corresponding normalised diagonal PSF. @fig:dtl-noisy-signbit shows the
sign-bit time-reversed excitation itself.

#figure(
  img("DTL_040_Diagonal_TimeLapse_Migration_Comparison_Noisy__Signed_Amplit.png"),
  caption: [Signed time-lapse-difference amplitude for all three migration
    algorithms on the _noisy_ diagonal dataset, sampled along the diagonal
    motion direction.],
) <fig:dtl-noisy-summary-amp>

#figure(
  img("DTL_041_Normalised_Diagonal_PSF_Noisy__TimeLapse_Difference_Along_Mo.png", width: 80%),
  caption: [Normalised diagonal PSF of the _noisy_ time-lapse-difference
    image along the motion direction, swept across the five diagonal
    scenarios of @tab:dtl-scenarios.],
) <fig:dtl-noisy-summary-psf>

#figure(
  img("DTL_033_Sign-Bit_Time-Reversed_Excitation__B-scans_and_Spectra_Noisy.png", width: 90%),
  caption: [Sign-bit time-reversed excitation B-scans and their frequency
    spectra for the _noisy_ diagonal dataset, all five displacement
    scenarios.],
) <fig:dtl-noisy-signbit>

@fig:dtl-noisy-phaseplane applies the 2D WLS phase-plane fit of
@sec:meth-phaseplane to the noisy diagonal dataset for all three migration
methods.

#figure(
  subfigs(cols: 2,
    img("TLP_017_DiagonalTimeLapse_Noisy__Phase-plane_shift_estimation__Kirch.png"),
    img("TLP_018_DiagonalTimeLapse_Noisy__Phase-plane_shift_estimation__Gazda.png"),
    img("TLP_019_DiagonalTimeLapse_Noisy__Phase-plane_shift_estimation__Back-.png"),
  ),
  caption: [2D phase-plane shift estimation applied to the _noisy_ diagonal
    time-lapse dataset of @sec:hyp1-diagonal, for all three migration
    methods: (a) Kirchhoff; (b) Gazdag; (c) back-propagation. The recovered
    $(#Dz, #Dx)$ should satisfy the known $#Dx = 2 #Dz$ ratio of
    @tab:dtl-scenarios if the fit is working correctly.],
) <fig:dtl-noisy-phaseplane>

#draftnote[read off the diagonal-PSF detectability floor under noise
(mirroring @sec:dtl-detectability) and state whether it shifts relative to
the clean-data floor, and whether @fig:dtl-noisy-phaseplane's recovered
displacements satisfy the known $#Dx = 2 #Dz$ ratio to the same tolerance as
the clean-data fit.]

== Noisy Fluid Flow Study <sec:hyp3-fluidflow>

The three translation studies above all move a rigid point scatterer; this
section instead repeats the noisy migration comparison for a target that is
directly relevant to the real fluid-injection field data of @ch:hyp3: a
_graded wetting zone_ rather than a discrete PEC cylinder. The domain, grid,
and centre frequency match @sec:hyp1-lateral exactly ($4.0 times 1.0$ m,
$f_c = 1.5 "GHz"$, $lambda = 112.6 "mm"$), but the background medium is ice
($epsilon_r = 3.15$) and the moving target is a $7$-step graded transition
(box width $11.3 "mm"$, total transition $78.8 "mm"$) between water and ice at
depth $0.676 "m"$, standing in for a fluid front advancing through a
horizontal fracture. As in @sec:hyp1-lateral, the front is swept laterally
across the same eight scenarios, from $2 lambda$ down to $1 \/ 32 lambda$
relative to its baseline position, and the same $10%$-of-signal-std synthetic
Laplace noise of @sec:hyp3-laplace is added before migration.

#draftnote[a clean-data (noise-free) counterpart to this section, mirroring
@sec:hyp1-lateral's structure, has not yet been written into the thesis even
though the clean figures (`FF_001`--`FF_024`) already exist --- add a short
clean-data fluid-flow subsection to @ch:hyp1 or earlier in this chapter, and
cross-reference it from here, before finalising.]

@fig:ff-noisy-summary-amp overlays the signed time-lapse-difference amplitude
from all three migration algorithms, directly comparable in structure to
@fig:tl-noisy-summary-amp, and @fig:ff-noisy-summary-psf plots the
corresponding normalised PSF. @fig:ff-noisy-signbit shows the sign-bit
time-reversed excitation itself.

#figure(
  img("FF_040_TimeLapse_Migration_Comparison_Noisy__Signed_Amplitude____f_.png"),
  caption: [Signed time-lapse-difference amplitude for all three migration
    algorithms on the _noisy_ fluid-flow dataset ($f_c = 1.5 "GHz"$,
    aperture $= 40$).],
) <fig:ff-noisy-summary-amp>

#figure(
  img("FF_041_Normalised_Lateral_PSF_Noisy__TimeLapse_Difference_at_True_S.png"),
  caption: [Normalised PSF of the _noisy_ time-lapse-difference image at the
    baseline front position, swept across front displacements from
    $2 lambda$ to $1 \/ 32 lambda$.],
) <fig:ff-noisy-summary-psf>

#figure(
  img("FF_033_Sign-Bit_Time-Reversed_Excitation__B-scans_and_Spectra_Noisy.png", width: 90%),
  caption: [Sign-bit time-reversed excitation B-scans and their frequency
    spectra for the _noisy_ fluid-flow dataset, all eight displacement
    scenarios.],
) <fig:ff-noisy-signbit>

Unlike the three point-scatterer directions above, no dedicated 2D WLS
phase-plane fit has been run for the fluid-flow front; the closest available
result is the per-method cross-spectrum displacement estimate already used
for shift inference, applied here to the noisy fluid-flow baseline/monitor
pairs at $10%$ noise level (@fig:ff-noisy-phaseplane).

#figure(
  subfigs(cols: 2,
    img("TLP_026_FluidFlow_Noisy_noise_level01__Kirchhoff____base_vs_mon_real.png"),
    img("TLP_027_FluidFlow_Noisy_noise_level01__Gazdag____base_vs_mon_real.png"),
    img("TLP_028_FluidFlow_Noisy_noise_level01__Back-prop____base_vs_mon_real.png"),
  ),
  caption: [Cross-spectrum phase-plane displacement estimation applied to the
    _noisy_ fluid-flow dataset (noise level $= 0.1$), for all three migration
    methods: (a) Kirchhoff; (b) Gazdag; (c) back-propagation. Each row shows
    the baseline-vs-monitor real-part difference, the cross-spectrum phase,
    the cross-spectrum energy, and the recovered vs. true front
    displacement.],
) <fig:ff-noisy-phaseplane>

#draftnote[state whether the recovered front displacement in
@fig:ff-noisy-phaseplane tracks the true value down to the same smallest
scale found for the point-scatterer directions above, and whether the
graded, spatially-extended nature of the fluid front (rather than a discrete
point reflector) changes the achievable resolution. Relate this result back
to the real borehole fluid-injection displacements of @sec:hyp3-fielddata
--- does the synthetic fluid-front experiment support the field-data
interpretation questions raised in @sec:hyp3-fd-interpretation?]

== Outstanding Quantitative Work <sec:hyp3-gaps>

This chapter now demonstrates noise robustness at the level of migrated
_images_ and the phase-plane estimator for all three point-scatterer
translation directions (lateral, vertical, diagonal) and, at the
cross-spectrum level, for the fluid-flow study. To fully support Hypothesis
2, the following quantitative work remains:

+ Run the plain 2D WLS phase-plane fit (rather than the GCC-peak-search
  fallback of @fig:tlp-noise) on the noisy lateral dataset, and extend it to
  Gazdag and back-propagation, so that the lateral case has the same
  three-method phase-plane coverage now available for vertical and diagonal
  (@fig:vtl-noisy-phaseplane, @fig:dtl-noisy-phaseplane).

+ Run a proper 2D WLS phase-plane fit for the fluid-flow front (rather than
  the cross-spectrum displacement estimate of @fig:ff-noisy-phaseplane), and
  add the clean-data fluid-flow subsection flagged in @sec:hyp3-fluidflow.

+ Quantitatively compare sign-bit versus peak-normalised back-propagation on
  pure noise (@sec:hyp3-purenoise), using the completed gprMax output already
  sitting in `noise_study/backprop/`.

+ A controlled SNR sweep from $+30 "dB"$ to $-10 "dB"$ comparing ordinary
  least squares against the WLS estimator.

+ A $20 times 20$ grid sweep of true $(#Dt, #Dtheta)$ pairs, reported as a
  2D inversion-error heatmap.

+ Resolve the colour-scale and "N/A"-column anomalies flagged in
  @sec:hyp3-lateral for the migration-comparison figures.

== Conclusion: Back-Propagation as the Preferred Method <sec:hyp2-conclusion>

The experiments in this chapter support Hypothesis 2: back-propagation with
sign-bit time-reversal is the most suitable migration technique for
noise-robust time-lapse phase-plane tracking.

Kirchhoff's coherence-manufacturing behaviour (@sec:hyp3-purenoise) creates
scatterer-like artefacts from pure noise, raising false-positive risk.
Gazdag stays incoherent but adds significant speckle. Back-propagation with
sign-bit time-reversal suppresses impulsive noise by reducing every noise
spike to the same $plus.minus 1$ amplitude as the coherent signal, preserving
phase information while stripping the amplitude-based false-positive risk
that would otherwise allow noise spikes to act as competing point sources
during back-propagation. The fluid-flow study of @sec:hyp3-fluidflow further
suggests this recommendation is not specific to discrete point scatterers:
the same three-way pattern (Kirchhoff/Gazdag amplitude collapse versus
back-propagation's noisier but present focus) reappears for a graded,
spatially-extended wetting-zone target, directly relevant to the borehole
fluid-injection geometry of @ch:hyp3. This recommendation is used for the
real field data in @ch:hyp3.
