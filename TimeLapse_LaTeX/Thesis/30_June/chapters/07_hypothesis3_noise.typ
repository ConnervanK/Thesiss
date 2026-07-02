#import "../template.typ": *

= Hypothesis 3: Robustness to Heavy-Tailed Noise <ch:hyp3>

@ch:hyp1 and @ch:hyp2 established that the phase-plane estimator recovers
sub-wavelength displacement and material change on clean, synthetic data.
This chapter asks whether that still holds once the data are corrupted by
realistic, heavy-tailed noise, and --- critically --- shows that the answer
depends on _which_ migration technique is used.

#para-head[Hypothesis 3.] Multi-dimensional phase-plane regression can deal
with noisy data, given that the right migration technique is used.

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
@sec:hyp3-diagonal, and @sec:hyp3-fluidflow: its Laplace _shape_ (loc $= 0$,
heavier tails than Gaussian) is kept, but its scale is rescaled so that the
resulting noise standard deviation is exactly $10%$ of each synthetic B-scan's
own signal standard deviation --- a light, realistic noise level rather than
the raw fitted scale, which would be incommensurate with the synthetic $E_z$
amplitudes.

== Migration-Algorithm Response to Pure Noise <sec:hyp3-purenoise>

Before testing noise robustness on real signal-plus-noise data, a sanity check
establishes what each migration algorithm does to noise _alone_: an empty
B-scan (no scatterer, no real signal) containing only the Laplace noise of
@sec:hyp3-laplace is migrated with Kirchhoff, Gazdag, and back-propagation.
If a method turns pure noise into something that looks like a coherent,
scatterer-like focus, every noisy result elsewhere in this thesis carries a
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
wavefield. The clean-data experiments of @ch:hyp1 and @ch:hyp2 use the
default peak-normalised excitation throughout, since they have no noise to
suppress; every noisy back-propagation result in this chapter uses sign-bit
excitation instead.

== Noisy Lateral Movement <sec:hyp3-lateral>

The lateral displacement sweep of @sec:hyp1-lateral is repeated on B-scans
contaminated with the synthetic Laplace noise of @fig:tl-bscans (c), for the
two analytically-defined migration algorithms (Kirchhoff and Gazdag;
back-propagation was not re-run on this particular noisy dataset due to its
computational cost --- contrast with the vertical, diagonal, and fluid-front
cases below, where it was). @fig:tl-noisy-kirchhoff and @fig:tl-noisy-gazdag
repeat the migration and differencing analysis of @sec:hyp1-lateral on this
noisy data.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("TL_024_Kirchhoff_Migration_-_All_8_Noisy_Datasets____f_c15_GHz____a.png"),
    img("TL_025_Kirchhoff_Migration_-_Noisy_zoomed____f_c15_GHz____aperture4.png"),
    img("TL_026_Kirchhoff_Migration_-_Noisy_TimeLapse_Differences_migrated_-.png"),
    img("TL_027_Kirchhoff_Migration_-_Noisy_TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Kirchhoff migration of the _noisy_ lateral time-lapse data
    ($f_c = 15 "GHz"$, aperture $= 40$): (a, b) migrated image for all
    scenarios and zoomed; (c, d) time-lapse difference and zoomed.],
) <fig:tl-noisy-kirchhoff>

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("TL_028_Gazdag_Phase-Shift_Migration_-_All_8_Noisy_Datasets____f_c15.png"),
    img("TL_029_Gazdag_Phase-Shift_Migration_-_Noisy_zoomed____f_c15_GHz.png"),
    img("TL_030_Gazdag_Migration_-_Noisy_TimeLapse_Differences_migrated_-_mi.png"),
    img("TL_031_Gazdag_Migration_-_Noisy_TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Gazdag phase-shift migration of the _noisy_ lateral time-lapse
    data ($f_c = 15 "GHz"$): (a, b) migrated image and zoomed; (c, d)
    time-lapse difference and zoomed.],
) <fig:tl-noisy-gazdag>

#draftnote[compare @fig:tl-noisy-kirchhoff and @fig:tl-noisy-gazdag against
the clean-data results of @fig:tl-kirchhoff and @fig:tl-gazdag and state at
which displacement scale the additive noise first prevents the time-lapse
difference from being distinguished from background clutter.]

The same noisy lateral dataset is also carried through the phase-plane
estimator of @sec:meth-phaseplane (@fig:tlp-noise), for the two migration
methods available on this noisy dataset.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
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
estimator's noise resilience; everything in @sec:hyp3-vertical,
@sec:hyp3-diagonal, and @sec:hyp3-fluidflow below is, for now,
amplitude-level only (see @sec:hyp3-gaps).]

== Noisy Vertical Movement <sec:hyp3-vertical>

The vertical displacement sweep of @sec:hyp1-vertical is repeated under the
same Laplace noise, for all three migration algorithms --- unlike the lateral
case, back-propagation _is_ re-run here, using sign-bit time-reversal
(@sec:hyp3-signbit) on the noisy traces. @fig:vtl-noisy-kirchhoff and
@fig:vtl-noisy-gazdag repeat the migration and differencing analysis of
@sec:hyp1-vertical on this noisy data.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("VTL_006_B-scan_Frequency_Spectra_--_Clean_vs_Noisy____f_c15_GHz.png"),
    img("VTL_021_Kirchhoff_Migration_-_All_7_Noisy_Datasets____f_c15_GHz____a.png"),
    img("VTL_026_Kirchhoff_Migration_Noisy_zoomed____f_c15_GHz____aperture40.png"),
    img("VTL_029_Gazdag_Phase-Shift_Migration_-_All_7_Noisy_Datasets____f_c15.png"),
  ),
  caption: [Noisy vertical time-lapse data: (a) clean-vs-noisy frequency
    spectra; (b, c) Kirchhoff migration; (d) Gazdag migration.],
) <fig:vtl-noisy-kirchhoff>

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("VTL_030_Gazdag_Phase-Shift_Migration_Noisy_zoomed____f_c15_GHz.png"),
    img("VTL_034_Back-Propagation_E_Noisy_--_All_7_Datasets____focus_at_1906.png"),
    img("VTL_036_Back-Propagation_Ez_Noisy_--_All_7_Datasets____focus_at_1906.png"),
    img("VTL_038_Back-Propagation_Noisy_--_TimeLapse_Differences_Ez____focus.png"),
  ),
  caption: [Gazdag (a) and sign-bit back-propagation (b--d) migration of the
    _noisy_ vertical time-lapse data, focused at $t = 1906 "ns"$: (a) Gazdag
    zoomed; (b) back-prop $||bold(E)||$ all scenarios; (c) back-prop $E_z$ all
    scenarios; (d) back-prop $E_z$ time-lapse difference.],
) <fig:vtl-noisy-gazdag>

#figure(
  img("VTL_039_Back-Propagation_Noisy_--_TimeLapse_Differences_Ez_zoomed.png", width: 60%),
  caption: [Sign-bit back-propagation $E_z$ time-lapse difference for the
    noisy vertical dataset, zoomed around the scatterer depth.],
) <fig:vtl-noisy-backprop>

#draftnote[*Gap:* this section currently only carries the analysis to the
amplitude/migrated-image level, exactly mirroring @sec:hyp1-vertical's
clean-data detectability summary --- the 2D WLS phase-plane fit has not yet
been applied to this noisy vertical dataset. Doing so, and comparing the
result against both the clean vertical fit of @sec:tlp-vertical and the noisy
lateral fit of @fig:tlp-noise, is necessary to claim Hypothesis 3 for the
vertical direction specifically.]

== Noisy Diagonal Movement <sec:hyp3-diagonal>

The diagonal displacement sweep of @sec:hyp1-diagonal is likewise repeated
under Laplace noise, for all three migration algorithms, again using sign-bit
time-reversal for the noisy back-propagation run. @fig:dtl-noisy shows the
noisy Kirchhoff and Gazdag migrations and time-lapse differences, the sign-bit
time-reversed excitation itself, and the noisy sign-bit back-propagation
result.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("DTL_025_Kirchhoff_Migration_-_All_6_Noisy_Datasets____f_c15_GHz____a.png"),
    img("DTL_027_Kirchhoff_Migration_-_Noisy_TimeLapse_Differences_migrated_-.png"),
    img("DTL_029_Gazdag_Phase-Shift_Migration_-_All_6_Noisy_Datasets____f_c15.png"),
    img("DTL_031_Gazdag_Migration_-_Noisy_TimeLapse_Differences_migrated_-_mi.png"),
    img("DTL_033_Sign-Bit_Time-Reversed_Excitation_--_B-scans_and_Spectra_Noi.png"),
    img("DTL_038_Back-Propagation_Noisy_--_TimeLapse_Differences_Ez____focus.png"),
  ),
  caption: [Noisy diagonal time-lapse data: Kirchhoff (a, b) and Gazdag (c, d)
    migration and differencing; (e) the sign-bit time-reversed back-propagation
    excitation B-scans and spectra; (f) the resulting noisy back-propagation
    time-lapse difference.],
) <fig:dtl-noisy>

#draftnote[as for the vertical case, the diagonal WLS phase-plane fit has not
yet been run on this noisy dataset --- this section is amplitude-level only.
Read off the diagonal-PSF detectability floor under noise (mirroring
@sec:dtl-detectability) and state whether it shifts relative to the clean-data
floor.]

== Noisy Fluid Front <sec:hyp3-fluidflow>

Finally, the fluid-front scenario of @ch:hyp2 is repeated under the same
Laplace noise model, again with sign-bit time-reversal for the noisy
back-propagation run. @fig:ff-noisy shows the noisy Kirchhoff and Gazdag
migrations and time-lapse differences, and the noisy sign-bit back-propagation
result.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("FF_017_Kirchhoff_Migration_-_All_8_Noisy_Datasets____f_c15_GHz____a.png"),
    img("FF_019_Kirchhoff_Migration_-_Noisy_TimeLapse_Differences_migrated_-.png"),
    img("FF_021_Gazdag_Phase-Shift_Migration_-_All_8_Noisy_Datasets____f_c15.png"),
    img("FF_023_Gazdag_Migration_-_Noisy_TimeLapse_Differences_migrated_-_mi.png"),
    img("FF_034_Back-Propagation_E_Noisy_--_All_Datasets____focus_at_1906_ns.png"),
    img("FF_038_Back-Propagation_Noisy_--_TimeLapse_Differences_Ez____focus.png"),
  ),
  caption: [Noisy fluid-front data: Kirchhoff (a, b) and Gazdag (c, d)
    migration and differencing; sign-bit back-propagation field magnitude (e)
    and time-lapse difference (f).],
) <fig:ff-noisy>

#draftnote[the noisy fluid-front phase-plane fit (the noisy analogue of
@fig:ff-phaseplane) has not yet been run --- this is the most field-relevant
gap in this chapter, since it is the only case combining noise robustness _and_
the geometry/material-change decoupling that the unifying hypothesis explicitly
claims together.]

== Outstanding Quantitative Work <sec:hyp3-gaps>

This chapter currently demonstrates noise robustness convincingly only at the
level of migrated _images_ (do they still show a detectable, unambiguous focus
under noise?), for all four cases, and demonstrates the phase-plane
estimator's noise robustness quantitatively for only one case (noisy lateral
movement, @fig:tlp-noise). To fully support Hypothesis 3, the following
quantitative work, proposed but not yet run (see
`Markdown/Robustness_Experiments.md`), is still needed:

+ Apply the 2D WLS phase-plane fit to the noisy vertical, diagonal, and
  fluid-front datasets of @sec:hyp3-vertical, @sec:hyp3-diagonal, and
  @sec:hyp3-fluidflow, exactly as already done for the lateral case.

+ Quantitatively compare sign-bit versus peak-normalised back-propagation on
  pure noise (@sec:hyp3-purenoise), using the completed gprMax output already
  sitting in `noise_study/backprop/`.

+ An orthogonality test isolating a pure geometric shift from a pure
  material-change phase rotation under noise.

+ A controlled SNR sweep from $+30 "dB"$ to $-10 "dB"$ comparing ordinary
  least squares against the WLS estimator.

+ A $20 times 20$ grid sweep of true $(#Dt, #Dtheta)$ pairs, reported as a
  2D inversion-error heatmap.
