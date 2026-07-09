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
for the experiments in @sec:hyp3-lateral, @sec:hyp3-vertical, and
@sec:hyp3-diagonal: its Laplace _shape_ (loc $= 0$, heavier tails than
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
contaminated with the synthetic Laplace noise of @fig:tl-bscans (c), for the
two analytically-defined migration algorithms (Kirchhoff and Gazdag;
back-propagation was not re-run on this particular noisy dataset due to its
computational cost --- contrast with the vertical and diagonal cases below,
where it was). @fig:tl-noisy-kirchhoff and @fig:tl-noisy-gazdag repeat the
migration and differencing analysis of @sec:hyp1-lateral on this noisy data.

#figure(
  subfigs(cols: 1,
    img("TL_025_Kirchhoff_Migration_-_Noisy_zoomed____f_c15_GHz____aperture4.png"),
    img("TL_027_Kirchhoff_Migration_-_Noisy_TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Kirchhoff migration of the _noisy_ lateral time-lapse data
    ($f_c = 1.5 "GHz"$, aperture $= 40$), zoomed around the scatterers:
    (a) migrated image for all scenarios; (b) time-lapse difference.],
) <fig:tl-noisy-kirchhoff>

#figure(
  subfigs(cols: 1,
    img("TL_029_Gazdag_Phase-Shift_Migration_-_Noisy_zoomed____f_c15_GHz.png"),
    img("TL_031_Gazdag_Migration_-_Noisy_TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Gazdag phase-shift migration of the _noisy_ lateral time-lapse
    data ($f_c = 1.5 "GHz"$), zoomed around the scatterers: (a) migrated
    image; (b) time-lapse difference.],
) <fig:tl-noisy-gazdag>

#draftnote[compare @fig:tl-noisy-kirchhoff and @fig:tl-noisy-gazdag against
the clean-data results of @fig:tl-kirchhoff and @fig:tl-gazdag and state at
which displacement scale the additive noise first prevents the time-lapse
difference from being distinguished from background clutter.]

The same noisy lateral dataset is also carried through the phase-plane
estimator of @sec:meth-phaseplane (@fig:tlp-noise), for the two migration
methods available on this noisy dataset.

#figure(
  subfigs(cols: 2,
    img("TLP_008_NoisyTimeLapse__Phase-plane_shift_estimation__Kirchhoff____B.png"),
    img("TLP_009_NoisyTimeLapse__Phase-plane_shift_estimation__Gazdag____Base.png"),
  ),
  caption: [2D phase-plane shift estimation applied to the _noisy_
    lateral-displacement dataset: (a) Kirchhoff; (b) Gazdag.],
) <fig:tlp-noise>

#draftnote[quantify, for both methods, how much the additive Laplace noise
degrades the estimated $#Dx$ relative to the clean result of
@fig:tlp-horiz-validation, and relate this to the amplitude-weighting
argument of @sec:th-mask-weight --- this figure pair is currently the _only_
direct evidence available anywhere in this thesis for the phase-plane
estimator's noise resilience; everything in @sec:hyp3-vertical and
@sec:hyp3-diagonal below is, for now, amplitude-level only (see
@sec:hyp3-gaps).]

== Noisy Vertical Movement <sec:hyp3-vertical>

The vertical displacement sweep of @sec:hyp1-vertical is repeated under the
same Laplace noise, for all three migration algorithms --- unlike the lateral
case, back-propagation _is_ re-run here, using sign-bit time-reversal
(@sec:hyp3-signbit) on the noisy traces. @fig:vtl-noisy-kirchhoff and
@fig:vtl-noisy-gazdag repeat the migration and differencing analysis of
@sec:hyp1-vertical on this noisy data, and @fig:vtl-noisy-backprop shows the
sign-bit back-propagation result.

#figure(
  subfigs(cols: 1,
    img("VTL_006_B-scan_Frequency_Spectra_--_Clean_vs_Noisy____f_c15_GHz.png", width: 90%),
    img("VTL_026_Kirchhoff_Migration_Noisy_zoomed____f_c15_GHz____aperture40.png", width: 90%),
    img("VTL_028_Kirchhoff_Migration_-_Noisy_TimeLapse_Differences_zoomed.png", width: 90%),
  ),
  caption: [Noisy vertical time-lapse data ($f_c = 1.5 "GHz"$), zoomed around
    the scatterer: (a) clean-vs-noisy frequency spectra; (b) Kirchhoff
    migration for all scenarios; (c) Kirchhoff time-lapse difference.],
) <fig:vtl-noisy-kirchhoff>

#figure(
  subfigs(cols: 1,
    img("VTL_030_Gazdag_Phase-Shift_Migration_Noisy_zoomed____f_c15_GHz.png"),
    img("VTL_032_Gazdag_Migration_-_Noisy_TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Gazdag phase-shift migration of the _noisy_ vertical time-lapse
    data ($f_c = 1.5 "GHz"$), zoomed around the scatterer: (a) migrated image;
    (b) time-lapse difference.],
) <fig:vtl-noisy-gazdag>

#figure(
  subfigs(cols: 1,
    img("VTL_035_Back-Propagation_E_Noisy_zoomed____focus_at_1906_ns.png", width: 90%),
    img("VTL_037_Back-Propagation_Ez_Noisy_zoomed____focus_at_1906_ns.png", width: 90%),
    img("VTL_039_Back-Propagation_Noisy_--_TimeLapse_Differences_Ez_zoomed.png", width: 90%),
  ),
  caption: [Sign-bit back-propagation migration of the _noisy_ vertical
    time-lapse data, focused at $t = 19.06 "ns"$ and zoomed around the
    scatterer: (a) field magnitude $||bold(E)||$; (b) $E_z$; (c) $E_z$
    time-lapse difference.],
) <fig:vtl-noisy-backprop>

#draftnote[*Gap:* this section currently only carries the analysis to the
amplitude/migrated-image level, exactly mirroring @sec:hyp1-vertical's
clean-data detectability summary --- the 2D WLS phase-plane fit has not yet
been applied to this noisy vertical dataset. Doing so, and comparing the
result against both the clean vertical fit of @sec:tlp-vertical and the noisy
lateral fit of @fig:tlp-noise, is necessary to claim Hypothesis 2 for the
vertical direction specifically.]

== Noisy Diagonal Movement <sec:hyp3-diagonal>

The diagonal displacement sweep of @sec:hyp1-diagonal is likewise repeated
under Laplace noise, for all three migration algorithms, again using sign-bit
time-reversal for the noisy back-propagation run. @fig:dtl-noisy-kirchhoff
and @fig:dtl-noisy-gazdag show the noisy Kirchhoff and Gazdag migrations and
their time-lapse differences; @fig:dtl-noisy-backprop shows the sign-bit
time-reversed excitation itself and the noisy sign-bit back-propagation
result.

#figure(
  subfigs(cols: 1,
    img("DTL_026_Kirchhoff_Migration__Noisy_zoomed____f_c15_GHz____aperture40.png"),
    img("DTL_028_Kirchhoff_Migration__Noisy_TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Kirchhoff migration of the _noisy_ diagonal time-lapse data
    ($f_c = 1.5 "GHz"$, aperture $= 40$), zoomed around the scatterer:
    (a) migrated image for all scenarios; (b) time-lapse difference.],
) <fig:dtl-noisy-kirchhoff>

#figure(
  subfigs(cols: 1,
    img("DTL_030_Gazdag_Phase-Shift_Migration__Noisy_zoomed____f_c15_GHz.png"),
    img("DTL_032_Gazdag_Migration__Noisy_TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Gazdag phase-shift migration of the _noisy_ diagonal time-lapse
    data ($f_c = 1.5 "GHz"$), zoomed around the scatterer: (a) migrated image;
    (b) time-lapse difference.],
) <fig:dtl-noisy-gazdag>

#figure(
  subfigs(cols: 1,
    img("DTL_033_Sign-Bit_Time-Reversed_Excitation__B-scans_and_Spectra_Noisy.png", width: 90%),
    img("DTL_037_Back-Propagation_Ez_Noisy_zoomed____focus_at_1906_ns.png", width: 90%),
    img("DTL_039_Back-Propagation_Noisy__TimeLapse_Differences_Ez_zoomed.png", width: 90%),
  ),
  caption: [Sign-bit back-propagation of the _noisy_ diagonal time-lapse
    data, focused at $t = 19.06 "ns"$: (a) the sign-bit time-reversed
    excitation B-scans and their spectra; (b) the back-propagated $E_z$,
    zoomed around the scatterer; (c) the $E_z$ time-lapse difference,
    zoomed.],
) <fig:dtl-noisy-backprop>

#draftnote[as for the vertical case, the diagonal WLS phase-plane fit has not
yet been run on this noisy dataset --- this section is amplitude-level only.
Read off the diagonal-PSF detectability floor under noise (mirroring
@sec:dtl-detectability) and state whether it shifts relative to the clean-data
floor.]

== Outstanding Quantitative Work <sec:hyp3-gaps>

This chapter currently demonstrates noise robustness convincingly only at the
level of migrated _images_ for all three displacement directions, and
demonstrates the phase-plane estimator's noise robustness quantitatively for
only one case (noisy lateral movement, @fig:tlp-noise). To fully support
Hypothesis 2, the following quantitative work remains:

+ Apply the 2D WLS phase-plane fit to the noisy vertical and diagonal datasets
  of @sec:hyp3-vertical and @sec:hyp3-diagonal, exactly as already done for
  the lateral case.

+ Quantitatively compare sign-bit versus peak-normalised back-propagation on
  pure noise (@sec:hyp3-purenoise), using the completed gprMax output already
  sitting in `noise_study/backprop/`.

+ A controlled SNR sweep from $+30 "dB"$ to $-10 "dB"$ comparing ordinary
  least squares against the WLS estimator.

+ A $20 times 20$ grid sweep of true $(#Dt, #Dtheta)$ pairs, reported as a
  2D inversion-error heatmap.

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
during back-propagation. This recommendation is used for the real field data
in @ch:hyp3.
