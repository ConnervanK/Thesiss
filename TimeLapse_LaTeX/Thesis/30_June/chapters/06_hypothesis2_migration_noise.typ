#import "../template.typ": *

= Hypothesis 2: Which Migration Technique Is Best Suited for Noise-Robust Phase-Plane Tracking? <ch:hyp2>

@ch:hyp1 established that phase-plane regression can track sub-wavelength
displacements on clean, synthetic data. Before the method is used in more
complex settings, a fundamental question must be answered: does the choice of
migration algorithm matter when the data are corrupted by noise? This
chapter repeats all four @ch:hyp1 experiments (Lateral, Vertical, Diagonal,
FluidFlow) under a realistic Laplace noise model fitted to real GPR field
data, and argues that migration choice matters substantially: back-propagation
with sign-bit time-reversal is the most noise-robust technique overall.

#para-head[Hypothesis 2.] Back-propagation migration with sign-bit
time-reversal is the most noise-robust technique for time-lapse phase-plane
tracking: it keeps the algorithm usable under heavy-tailed Laplace noise by
suppressing large-amplitude noise spikes before back-propagation, while
Kirchhoff produces false-coherent artefacts and Gazdag adds incoherent
speckle.

== Noise Creation <sec:hyp3-laplace>

Rather than injecting arbitrary synthetic noise, the noise level and shape
used throughout this thesis are fitted to real GPR field noise. A sample of
field noise is tracked through eleven stages of the processing pipeline and
fitted with both a Laplace and a Gaussian distribution at two stages: after
the final processing stage ("9 Crop Samples"), and at the last stage before
the spherical-gain correction ("7 Constant Velocity"), since the gain
correction inflates the amplitude scale by several orders of magnitude and
is not representative of the raw simulated $E_z$ amplitudes used elsewhere in
this thesis.

#figure(
  img("H2_001_Noise_amplitude_distribution_by_processing_stage.png"),
  caption: [Noise amplitude distribution at every one of the eleven tracked
    processing stages, each with a Gaussian reference overlay.],
) <fig:h2-noise-stages>

#figure(
  img("H2_002_Noise_distribution_at_9_Crop_Samples_post-gain_with_Laplace.png"),
  caption: [Noise distribution at the final processing stage ("9 Crop
    Samples", post-gain), with both Laplace and Gaussian fits overlaid.],
) <fig:h2-noise-laplace-fit>

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
  img("H2_003_Migrating_Pure_Noise_no_scatterers_no_signal_--_Kirchhoff_vs.png"),
  caption: [Migrating pure noise (no scatterers, no signal): the input
    pure-noise B-scan, its Kirchhoff migration, and its Gazdag migration.],
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

#figure(
  img("H2_004_Back-propagation_of_pure_noise_--_peak-normalised_excitation.png"),
  caption: [Back-propagation of pure noise, default peak-normalised
    excitation: focus-time snapshot of the back-propagated wavefield.],
) <fig:h2-purenoise-backprop>

With no scatterer present, there is no true location for either excitation
scheme to focus on, so the back-propagated wavefield stays diffuse speckle
throughout the domain rather than collapsing into an obvious spurious bright
spot --- unlike Kirchhoff's coherent bands above, there is no clearly visible
artefact here to point at. What does stand out is the amplitude scale:
peak-normalised excitation reaches only $approx 116 space 600 "V/m"$, since
only each trace's single largest sample is normalised to $plus.minus 1$ and
every other noise sample stays small. @sec:hyp3-signbit revisits this with
sign-bit excitation, where every sample --- not just the peak --- is forced
to $plus.minus 1$.

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

_Sign-bit time-reversal_ fixes this by injecting only the _sign_ of the
time-reversed wavefield instead of its peak-normalised value,
$ u_"sign" (x, tau) = op("sign")(u(x, tau)) , $ <eq:signbit>
implemented in `write_backprop_files(..., sign_bit=True)`
(`helper_functions/migration.py`). This keeps every zero-crossing and phase
trend of the GPR wavelet completely intact, since the sign of a signal
carries its full phase information, while squashing every noise spike down
to the same $plus.minus 1$ amplitude as the coherent signal --- stripping
noise of the outsized amplitude that would otherwise let it dominate the
back-propagated wavefield. The clean-data experiments of @ch:hyp1 use the
default peak-normalised excitation throughout, since they have no noise to
suppress; every noisy back-propagation result in this chapter uses sign-bit
excitation instead.

#figure(
  img("H2_005_Sign-Bit_Time-Reversed_Excitation_--_B-scans_and_Spectra_Lat.png", width: 90%),
  caption: [Sign-bit time-reversed excitation for the noisy lateral dataset:
    B-scans (top) and their frequency spectra (bottom), every displacement
    scenario.],
) <fig:h2-signbit-excitation>

#figure(
  img("H2_006_Back-Propagation_of_Pure_Noise_--_Peak-Normalised_vs_Sign-Bi.png"),
  caption: [Back-propagation of pure noise, peak-normalised versus sign-bit
    excitation: focus-time snapshots side by side.],
) <fig:h2-signbit-purenoise>

Sign-bit excitation reaches $approx 713 space 100 "V/m"$, about $6 times$
larger than peak-normalised's, since forcing every sample (not just each
trace's single peak) to $plus.minus 1$ injects far more total energy into
the medium --- but the wavefield itself still stays diffuse speckle in both
panels of @fig:h2-signbit-purenoise, exactly as in @sec:hyp3-purenoise ---
pure noise has no true target to focus on, so this idealised zero-signal test
cannot visually demonstrate whether sign-bit suppresses spurious focusing
the way it can on real, signal-bearing data. The actual, quantitative
evidence that sign-bit back-propagation is noise-robust comes from
@sec:hyp3-summary's master MAE table on the real noisy studies below, where
back-propagation is one of the two most accurate methods overall.

== Extra Processing Steps in the Phase Domain to Remove Noise <sec:hyp3-phase-denoise>

#draftnote[Reserved for future work: extra phase-domain noise-suppression
processing applied to the migrated (post-Kirchhoff/Gazdag/back-propagation)
images, on top of the sign-bit time-reversal already applied at the
back-propagation excitation stage (@sec:hyp3-signbit). Nothing is
implemented here yet; a natural candidate, given the Gazdag streaking
artefact identified in @sec:hyp3-summary, is a targeted filter in the
vertical-wavenumber domain, though @sec:hyp3-summary shows a first attempt at
this (a raised-cosine low-$#kz$ taper) did not resolve that artefact.]

== Noisy Lateral Movement <sec:hyp3-lateral>

The lateral displacement sweep of @sec:hyp1-lateral (@tab:h1-lat-scenarios)
is repeated on B-scans contaminated with the Laplace noise of
@sec:hyp3-laplace, for all three migration algorithms, using sign-bit
time-reversal (@sec:hyp3-signbit) for back-propagation.

=== Migration Results <sec:hyp3-lat-migration>

#figure(
  img("H2_007_Lateral_--_TimeLapse_Migration_Comparison_Noisy_--_Signed_Am.png", width: 85%),
  caption: [Signed time-lapse-difference amplitude (monitor-minus-baseline)
    for all three migration algorithms, lateral time-lapse study, _noisy_
    data.],
) <fig:h2-lat-summary-amp>

=== Amplitude Test <sec:hyp3-lat-amplitude>

#figure(
  img("H2_008_Lateral_Noisy_--_Amplitude_PSF_zoom_Baseline_vs_Monitor_RAW.png", width: 75%),
  caption: [Zoomed Baseline-versus-Monitor point-spread-function comparison
    for the raw _noisy_ migrated lateral images, all scenarios and migration
    methods.],
) <fig:h2-lat-amp-zoom>

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

#figure(
  img("H2_009_Lateral_Noisy_--_Phase-plane_shift_estimation_--_Kirchhoff.png", width: 75%),
  caption: [Kirchhoff-migrated phase-plane shift estimation, _noisy_ lateral
    study: Baseline versus each of the seven displacement scenarios.],
) <fig:h2-lat-phase-kirchhoff>

#figure(
  img("H2_010_Lateral_Noisy_--_Phase-plane_shift_estimation_--_Gazdag____B.png", width: 75%),
  caption: [Gazdag-migrated phase-plane shift estimation, _noisy_ lateral
    study: Baseline versus each of the seven displacement scenarios.],
) <fig:h2-lat-phase-gazdag>

#figure(
  img("H2_011_Lateral_Noisy_--_Phase-plane_shift_estimation_--_Back-prop.png", width: 75%),
  caption: [Back-propagation-migrated (sign-bit) phase-plane shift
    estimation, _noisy_ lateral study: Baseline versus each of the seven
    displacement scenarios.],
) <fig:h2-lat-phase-backprop>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$2 lambda$],     [$-224.00$], [$-216.21$], [$-204.23$],
    [$1 lambda$],     [$-96.92$],  [$-113.50$], [$-81.57$],
    [$1\/2 lambda$],  [$-58.06$],  [$-73.69$],  [$+5.71$],
    [$1\/4 lambda$],  [$-0.05$],   [$-22.17$],  [$-0.29$],
    [$1\/8 lambda$],  [$+0.24$],   [$-29.01$],  [$-14.60$],
    [$1\/16 lambda$], [$+0.02$],   [$-11.86$],  [$+5.61$],
    [$1\/32 lambda$], [$+0.14$],   [$-1.79$],   [$+7.31$],
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
    [$2 lambda$],     [$-99.6%$],  [$-96.1%$],  [$-90.8%$],
    [$1 lambda$],     [$-85.8%$],  [$-100.4%$], [$-72.2%$],
    [$1\/2 lambda$],  [$-103.7%$], [$-131.6%$], [$+10.2%$],
    [$1\/4 lambda$],  [$-0.2%$],   [$-79.2%$],  [$-1.0%$],
    [$1\/8 lambda$],  [$+1.7%$],   [$-207.2%$], [$-104.3%$],
    [$1\/16 lambda$], [$+0.2%$],   [$-169.4%$], [$+80.2%$],
    [$1\/32 lambda$], [$+3.4%$],   [$-44.8%$],  [$+182.6%$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Lateral phase-plane WLS displacement error, $#Dx$, as a
    percentage of the true displacement, _noisy_ data (compare
    @tab:h1-lat-phase-pct).],
  kind: table,
) <tab:h2-lat-phase-pct>

Noise degrades every method relative to the clean-data result. Kirchhoff and
back-propagation stay within a few millimetres from $1\/4 lambda$ downward
(back-propagation in particular recovers to within $0.24 "mm"$), but
Kirchhoff still shows a $14.6 "mm"$ excursion at $1\/8 lambda$ and
$5$--$7 "mm"$ residual error at the two smallest scales --- noise-driven
scatter around zero rather than a systematic bias. Gazdag is qualitatively
worse throughout the sub-half-wavelength regime ($-1.8$ to $-29 "mm"$ from
$1\/4 lambda$ to $1\/32 lambda$), foreshadowing the structural artefact
identified in @sec:hyp3-summary.

== Noisy Vertical Movement <sec:hyp3-vertical>

The vertical displacement sweep of @sec:hyp1-vertical (@tab:h1-vert-scenarios)
is repeated under the same Laplace noise, for all three migration
algorithms, using sign-bit time-reversal for the noisy back-propagation run.

=== Migration Results <sec:hyp3-vert-migration>

#figure(
  img("H2_012_Vertical_--_TimeLapse_Migration_Comparison_Noisy_--_Signed_A.png", width: 85%),
  caption: [Signed time-lapse-difference amplitude (monitor-minus-baseline)
    for all three migration algorithms, vertical time-lapse study, _noisy_
    data.],
) <fig:h2-vert-summary-amp>

=== Amplitude Test <sec:hyp3-vert-amplitude>

#figure(
  img("H2_013_Vertical_Noisy_--_Amplitude_PSF_zoom_Baseline_vs_Monitor_RAW.png", width: 75%),
  caption: [Zoomed Baseline-versus-Monitor point-spread-function comparison
    for the raw _noisy_ migrated vertical images, all scenarios and
    migration methods.],
) <fig:h2-vert-amp-zoom>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$1 lambda$],     [6.508], [5.695], [6.528],
    [$1\/2 lambda$],  [3.977], [3.480], [3.406],
    [$1\/4 lambda$],  [2.169], [1.582], [1.703],
    [$1\/8 lambda$],  [1.085], [0.633], [0.851],
    [$1\/16 lambda$], [1.085], [0.316], [0.568],
    [$1\/32 lambda$], [0.362], [0.0],   [0.284],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Vertical Rayleigh-criterion ratio, _noisy_ data (compare
    @tab:h1-vert-amp).],
  kind: table,
) <tab:h2-vert-amp>

As in the clean case, vertical amplitude differencing stays resolvable
comfortably below $1\/4 lambda$ for every method, only crossing below $1$
between $1\/4 lambda$ and $1\/8 lambda$ --- noise has almost no effect on
this floor either.

=== Phase Test <sec:hyp3-vert-phase>

#figure(
  img("H2_014_Vertical_Noisy_--_Phase-plane_shift_estimation_--_Kirchhoff.png", width: 75%),
  caption: [Kirchhoff-migrated phase-plane shift estimation, _noisy_
    vertical study: Baseline versus each of the six displacement scenarios.],
) <fig:h2-vert-phase-kirchhoff>

#figure(
  img("H2_015_Vertical_Noisy_--_Phase-plane_shift_estimation_--_Gazdag.png", width: 75%),
  caption: [Gazdag-migrated phase-plane shift estimation, _noisy_ vertical
    study: Baseline versus each of the six displacement scenarios.],
) <fig:h2-vert-phase-gazdag>

#figure(
  img("H2_016_Vertical_Noisy_--_Phase-plane_shift_estimation_--_Back-prop.png", width: 75%),
  caption: [Back-propagation-migrated (sign-bit) phase-plane shift
    estimation, _noisy_ vertical study: Baseline versus each of the six
    displacement scenarios.],
) <fig:h2-vert-phase-backprop>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$1 lambda$],     [$-82.09$], [$-101.78$], [$-93.93$],
    [$1\/2 lambda$],  [$-58.33$], [$-50.71$],  [$-77.10$],
    [$1\/4 lambda$],  [$-0.80$],  [$-25.32$],  [$+0.67$],
    [$1\/8 lambda$],  [$-0.44$],  [$-24.74$],  [$-2.23$],
    [$1\/16 lambda$], [$-0.49$],  [$-2.34$],   [$+0.84$],
    [$1\/32 lambda$], [$-0.38$],  [$+0.22$],   [$-0.85$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Vertical phase-plane WLS displacement error, $#Dz$
    (estimated $-$ true), millimetres, _noisy_ data (compare
    @tab:h1-vert-phase).],
  kind: table,
) <tab:h2-vert-phase>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$1 lambda$],     [$-72.6%$],  [$-90.1%$],  [$-83.1%$],
    [$1\/2 lambda$],  [$-104.2%$], [$-90.6%$],  [$-137.7%$],
    [$1\/4 lambda$],  [$-2.8%$],   [$-90.4%$],  [$+2.4%$],
    [$1\/8 lambda$],  [$-3.1%$],   [$-176.7%$], [$-15.9%$],
    [$1\/16 lambda$], [$-7.0%$],   [$-33.5%$],  [$+12.0%$],
    [$1\/32 lambda$], [$-9.5%$],   [$+5.5%$],   [$-21.1%$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Vertical phase-plane WLS displacement error, $#Dz$, as a
    percentage of the true displacement, _noisy_ data (compare
    @tab:h1-vert-phase-pct).],
  kind: table,
) <tab:h2-vert-phase-pct>

Back-propagation is the most consistent method here, staying within
$0.8 "mm"$ from $1\/4 lambda$ downward. Kirchhoff is close behind but with
more scatter (up to $2.2 "mm"$). Gazdag is markedly worse, with a
$24$--$25 "mm"$ error persisting at $1\/4 lambda$ and $1\/8 lambda$ before
improving at the two smallest scales --- the same qualitative pattern seen
for Lateral.

== Noisy Diagonal Movement <sec:hyp3-diagonal>

The diagonal displacement sweep of @sec:hyp1-diagonal
(@tab:h1-diag-scenarios) is likewise repeated under Laplace noise, for all
three migration algorithms, again using sign-bit time-reversal for the noisy
back-propagation run.

=== Migration Results <sec:hyp3-diag-migration>

#figure(
  img("H2_017_Diagonal_--_TimeLapse_Migration_Comparison_Noisy_--_Signed_A.png", width: 85%),
  caption: [Signed time-lapse-difference amplitude (monitor-minus-baseline)
    for all three migration algorithms, diagonal time-lapse study, _noisy_
    data.],
) <fig:h2-diag-summary-amp>

=== Amplitude Test <sec:hyp3-diag-amplitude>

#figure(
  img("H2_018_Diagonal_Noisy_--_Amplitude_PSF_zoom_Baseline_vs_Monitor_RAW.png", width: 75%),
  caption: [Zoomed Baseline-versus-Monitor point-spread-function comparison
    for the raw _noisy_ migrated diagonal images, all scenarios and
    migration methods.],
) <fig:h2-diag-amp-zoom>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [1 ($1 lambda_x$, $1\/2 lambda_z$)],   [5.295], [4.157], [3.430],
    [2 ($1\/2 lambda_x$, $1\/4 lambda_z$)], [2.492], [2.019], [1.663],
    [3 ($1\/4 lambda_x$, $1\/8 lambda_z$)], [1.090], [1.069], [0.831],
    [4 ($1\/8 lambda_x$, $1\/16 lambda_z$)],[0.623], [0.594], [0.312],
    [5 ($1\/16 lambda_x$, $1\/32 lambda_z$)],[0.156],[0.238], [0.104],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Diagonal Rayleigh-criterion ratio, _noisy_ data (compare
    @tab:h1-diag-amp).],
  kind: table,
) <tab:h2-diag-amp>

=== Phase Test <sec:hyp3-diag-phase>

#figure(
  img("H2_019_Diagonal_Noisy_--_Phase-plane_shift_estimation_--_Kirchhoff.png", width: 75%),
  caption: [Kirchhoff-migrated phase-plane shift estimation, _noisy_
    diagonal study: Baseline versus each of the five displacement scenarios.],
) <fig:h2-diag-phase-kirchhoff>

#figure(
  img("H2_020_Diagonal_Noisy_--_Phase-plane_shift_estimation_--_Gazdag.png", width: 75%),
  caption: [Gazdag-migrated phase-plane shift estimation, _noisy_ diagonal
    study: Baseline versus each of the five displacement scenarios.],
) <fig:h2-diag-phase-gazdag>

#figure(
  img("H2_021_Diagonal_Noisy_--_Phase-plane_shift_estimation_--_Back-prop.png", width: 75%),
  caption: [Back-propagation-migrated (sign-bit) phase-plane shift
    estimation, _noisy_ diagonal study: Baseline versus each of the five
    displacement scenarios.],
) <fig:h2-diag-phase-backprop>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [1], [$-53.62$], [$-44.78$], [$-62.55$],
    [2], [$-37.96$], [$-22.71$], [$-15.05$],
    [3], [$-0.89$],  [$-11.32$], [$+0.02$],
    [4], [$-0.18$],  [$-17.74$], [$-1.75$],
    [5], [$-0.35$],  [$+0.66$],  [$+0.89$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Diagonal phase-plane WLS displacement error, $#Dz$
    (estimated $-$ true), millimetres, _noisy_ data (compare
    @tab:h1-diag-phase-dz).],
  kind: table,
) <tab:h2-diag-phase-dz>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [1], [$-114.43$], [$-108.25$], [$-135.49$],
    [2], [$-46.38$],  [$-54.99$],  [$-99.64$],
    [3], [$-0.91$],   [$-32.44$],  [$+0.34$],
    [4], [$+0.10$],   [$-26.34$],  [$+1.23$],
    [5], [$-0.07$],   [$-12.34$],  [$+1.74$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Diagonal phase-plane WLS displacement error, $#Dx$
    (estimated $-$ true), millimetres, _noisy_ data (compare
    @tab:h1-diag-phase-dx).],
  kind: table,
) <tab:h2-diag-phase-dx>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [1], [$-95.7%$],  [$-80.0%$],  [$-111.7%$],
    [2], [$-135.6%$], [$-81.1%$],  [$-53.7%$],
    [3], [$-6.3%$],   [$-80.9%$],  [$+0.2%$],
    [4], [$-2.5%$],   [$-253.4%$], [$-25.0%$],
    [5], [$-8.6%$],   [$+16.4%$],  [$+22.1%$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Diagonal phase-plane WLS displacement error, $#Dz$, as a
    percentage of the true displacement, _noisy_ data (compare
    @tab:h1-diag-phase-dz-pct).],
  kind: table,
) <tab:h2-diag-phase-dz-pct>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [1], [$-101.3%$], [$-95.8%$],  [$-119.9%$],
    [2], [$-82.8%$],  [$-98.2%$],  [$-177.9%$],
    [3], [$-3.3%$],   [$-115.9%$], [$+1.2%$],
    [4], [$+0.7%$],   [$-188.1%$], [$+8.8%$],
    [5], [$-1.0%$],   [$-176.3%$], [$+24.8%$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Diagonal phase-plane WLS displacement error, $#Dx$, as a
    percentage of the true displacement, _noisy_ data (compare
    @tab:h1-diag-phase-dx-pct).],
  kind: table,
) <tab:h2-diag-phase-dx-pct>

Back-propagation is again the most accurate method from Scenario 3 onward
(within $0.9 "mm"$ on both axes). Kirchhoff is comparable at Scenarios 3 and
5 but shows a larger excursion at Scenario 4 ($1.75$/$1.23 "mm"$). Gazdag
again carries a persistent tens-of-millimetre error through Scenarios 3 and
4, only improving at Scenario 5 --- consistent with the same structural
artefact seen in the lateral and vertical cases.

== Noisy Fluid Flow <sec:hyp3-fluidflow>

The clean-data fluid-flow experiment of @sec:hyp1-fluidflow
(@tab:h1-ff-scenarios) is repeated under the same Laplace noise, for a
target directly relevant to the real fluid-injection field data of
@ch:hyp3 --- a graded wetting zone rather than a discrete point scatterer,
swept laterally across the same seven scenarios.

=== Migration Results <sec:hyp3-ff-migration>

#figure(
  img("H2_022_FluidFlow_--_TimeLapse_Migration_Comparison_Noisy_--_Signed.png", width: 85%),
  caption: [Signed time-lapse-difference amplitude (monitor-minus-baseline)
    for all three migration algorithms, fluid-flow time-lapse study, _noisy_
    data.],
) <fig:h2-ff-summary-amp>

=== Amplitude Test <sec:hyp3-ff-amplitude>

#figure(
  img("H2_023_FluidFlow_Noisy_--_Amplitude_PSF_zoom_Baseline_vs_Monitor_RA.png", width: 75%),
  caption: [Zoomed Baseline-versus-Monitor point-spread-function comparison
    for the raw _noisy_ migrated fluid-flow images, all scenarios and
    migration methods.],
) <fig:h2-ff-amp-zoom>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$2 lambda$],     [4.853], [0.036], [5.450],
    [$1 lambda$],     [2.427], [0.214], [2.725],
    [$1\/2 lambda$],  [1.103], [0.605], [1.239],
    [$1\/4 lambda$],  [0.441], [0.107], [0.495],
    [$1\/8 lambda$],  [0.221], [0.071], [0.248],
    [$1\/16 lambda$], [0.0],   [0.0],   [0.0],
    [$1\/32 lambda$], [0.0],   [0.036], [0.0],
    table.hline(stroke: 0.7pt),
  ),
  caption: [FluidFlow Rayleigh-criterion ratio, _noisy_ data (compare
    @tab:h1-ff-amp).],
  kind: table,
) <tab:h2-ff-amp>

=== Phase Test <sec:hyp3-ff-phase>

#figure(
  img("H2_024_FluidFlow_Noisy_--_Phase-plane_shift_estimation_--_Kirchhoff.png", width: 75%),
  caption: [Kirchhoff-migrated phase-plane shift estimation, _noisy_
    fluid-flow study: Baseline versus each of the seven displacement
    scenarios.],
) <fig:h2-ff-phase-kirchhoff>

#figure(
  img("H2_025_FluidFlow_Noisy_--_Phase-plane_shift_estimation_--_Gazdag.png", width: 75%),
  caption: [Gazdag-migrated phase-plane shift estimation, _noisy_ fluid-flow
    study: Baseline versus each of the seven displacement scenarios.],
) <fig:h2-ff-phase-gazdag>

#figure(
  img("H2_026_FluidFlow_Noisy_--_Phase-plane_shift_estimation_--_Back-prop.png", width: 75%),
  caption: [Back-propagation-migrated (sign-bit) phase-plane shift
    estimation, _noisy_ fluid-flow study: Baseline versus each of the seven
    displacement scenarios.],
) <fig:h2-ff-phase-backprop>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$2 lambda$],     [$-218.10$], [$-225.00$], [$-225.00$],
    [$1 lambda$],     [$-104.24$], [$-113.00$], [$+1.11$],
    [$1\/2 lambda$],  [$-65.06$],  [$-56.00$],  [$+0.46$],
    [$1\/4 lambda$],  [$-0.52$],   [$-28.00$],  [$+0.35$],
    [$1\/8 lambda$],  [$-0.46$],  [$-14.00$],  [$+0.23$],
    [$1\/16 lambda$], [$-0.73$],   [$-7.00$],   [$-0.00$],
    [$1\/32 lambda$], [$-1.17$],   [$-4.00$],   [$-0.08$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [FluidFlow phase-plane WLS front-displacement error, $#Dx$
    (estimated $-$ true), millimetres, _noisy_ data (compare
    @tab:h1-ff-phase).],
  kind: table,
) <tab:h2-ff-phase>

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Scenario*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [$2 lambda$],     [$-96.9%$],  [$-100.0%$], [$-100.0%$],
    [$1 lambda$],     [$-92.3%$],  [$-100.0%$], [$+1.0%$],
    [$1\/2 lambda$],  [$-116.2%$], [$-100.0%$], [$+0.8%$],
    [$1\/4 lambda$],  [$-1.9%$],   [$-100.0%$], [$+1.2%$],
    [$1\/8 lambda$],  [$-3.3%$],   [$-100.0%$], [$+1.6%$],
    [$1\/16 lambda$], [$-10.4%$],  [$-100.0%$], [$-0.1%$],
    [$1\/32 lambda$], [$-29.3%$],  [$-100.0%$], [$-2.1%$],
    table.hline(stroke: 0.7pt),
  ),
  caption: [FluidFlow phase-plane WLS front-displacement error, $#Dx$, as a
    percentage of the true displacement, _noisy_ data (compare
    @tab:h1-ff-phase-pct).],
  kind: table,
) <tab:h2-ff-phase-pct>

Kirchhoff is by far the most accurate method for the noisy fluid-flow front:
already within $1.2 "mm"$ from $1 lambda$ downward, and within $0.5 "mm"$
from $1\/4 lambda$ downward. Back-propagation only becomes reliable from
$1\/4 lambda$ ($0.52 "mm"$), and Gazdag never recovers a small error --- its
error stays between $-4$ and $-28 "mm"$ across the entire sub-half-wavelength
regime, again matching the structural weakness identified in
@sec:hyp3-summary.

== Summary of the Results <sec:hyp3-summary>

@tab:h2-lat-phase, @tab:h2-vert-phase, @tab:h2-diag-phase-dz/@tab:h2-diag-phase-dx,
and @tab:h2-ff-phase give the full per-scenario noisy phase-plane error for
each experiment in millimetres, with the same errors expressed as a
percentage of the true displacement in @tab:h2-lat-phase-pct,
@tab:h2-vert-phase-pct, @tab:h2-diag-phase-dz-pct/@tab:h2-diag-phase-dx-pct,
and @tab:h2-ff-phase-pct. @tab:h2-mae condenses the millimetre errors into
one mean absolute error (MAE) per movement type and method, using the same
sub-half-wavelength regime as @tab:h1-mae, so that it is directly comparable
to the clean-data result of @ch:hyp1.

#figure(
  table(
    columns: (auto, auto, auto, auto),
    stroke: none,
    inset: (x: 0.8em, y: 0.3em),
    table.hline(stroke: 0.7pt),
    [*Movement*], [*Back-prop*], [*Gazdag*], [*Kirchhoff*],
    table.hline(stroke: 0.4pt),
    [Lateral],   [12.805], [37.900], [8.268],
    [Vertical],  [12.394], [25.979], [22.651],
    [Diagonal],  [0.609],  [26.157], [1.478],
    [FluidFlow], [15.971], [25.601], [0.254],
    table.hline(stroke: 0.4pt),
    [*Mean*],    [*10.445*], [*28.909*], [*8.162*],
    table.hline(stroke: 0.7pt),
  ),
  caption: [Mean absolute phase-plane displacement error [mm], _noisy_ data,
    sub-half-wavelength regime only (compare @tab:h1-mae).],
  kind: table,
) <tab:h2-mae>

The best-performing method differs by movement type: Kirchhoff is most
accurate for Lateral ($8.3 "mm"$) and FluidFlow ($0.3 "mm"$), while sign-bit
back-propagation is most accurate for Vertical ($12.4 "mm"$) and Diagonal
($0.6 "mm"$). Averaged across all four movement types, Kirchhoff
($8.2 "mm"$ mean MAE) and back-propagation ($10.4 "mm"$ mean MAE) are close
to each other and both far more noise-robust than Gazdag ($28.9 "mm"$ mean
MAE, worst for every single movement type by a wide margin).

Gazdag's poor showing is *not* an apex-finding problem: even with the same
robust envelope-based apex search used for every method, Gazdag's noisy
migrated images carry a structural streaking artefact (independently
documented via `TimeLapse_Processing.ipynb`, amplitude up to $approx 60
times$ the genuine scatterer signal) that dominates both the apex search and
the WLS phase-plane fit regardless of how the apex is located. The
underlying mechanism was tested directly: a raised-cosine taper suppressing
low vertical-wavenumber energy (up to half the central wavenumber) made no
measurable difference to the artefact's location or amplitude, and the raw
noisy input trace at the artefact location shows no anomalous amplitude
either --- ruling out a "noise concentrated near $#kz = 0$" explanation. The
true mechanism is most likely a numerical property of the phase-shift
depth-stepping operator itself, not the noise's spectral content; it remains
unidentified and is left as outstanding work (@sec:hyp3-phase-denoise).

Kirchhoff's aperture-stacking sums over many traces and partially averages
the noise down, which particularly helps recover Lateral's and FluidFlow's
smallest sub-wavelength shifts once the apex is correctly located.
Sign-bit back-propagation remains the most accurate for Vertical and
Diagonal, consistent with @sec:hyp3-signbit's motivation: clamping every
noise spike to $plus.minus 1$ before back-propagation keeps its
phase-governed focusing intact. Diagonal is the easiest case for both of
these methods (under $1.5 "mm"$ MAE), likely because its combined 2D
$(#Dz, #Dx)$ error norm partially cancels axis-wise noise scatter that would
otherwise show up as pure along-axis error in the Lateral or Vertical cases.

Comparing @tab:h2-mae directly against the clean-data @tab:h1-mae, noise
degrades every method's accuracy by roughly one to two orders of magnitude
in the sub-half-wavelength regime, yet the qualitative conclusion of
@ch:hyp1 survives: Kirchhoff and back-propagation both remain accurate to
single-digit millimetres or better for the majority of movement types, at
displacement scales where the corresponding amplitude tests
(@tab:h2-lat-amp, @tab:h2-vert-amp, @tab:h2-diag-amp, @tab:h2-ff-amp) show
amplitude differencing has already collapsed. Gazdag, in contrast, is
unsuitable for noisy time-lapse phase-plane tracking regardless of movement
type --- supporting Hypothesis 2's specific recommendation of
back-propagation with sign-bit time-reversal (matched or exceeded only by
Kirchhoff, and only for the two movement types, Lateral and FluidFlow, where
Kirchhoff's trace-stacking noise averaging happens to help most).
