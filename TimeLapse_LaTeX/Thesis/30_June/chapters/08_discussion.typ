#import "../template.typ": *

= Discussion <ch:discussion>

#draftnote[this chapter is a placeholder skeleton built around the draft-note
prompts left in @ch:hyp1, @ch:hyp2, and @ch:hyp3. Resolve those prompts
first (they ask for specific numbers read off specific figures), then rewrite
this chapter as connected prose rather than a list.]

This chapter returns to the Research Question of @ch:introduction --- can
time-lapse GPR accurately track subwavelength movement by analysing phase
changes in migrated images? --- and integrates the results of Hypotheses 1--3
to answer it.

== Hypothesis 1: Phase Changes Reveal Sub-Resolution Movement

@sec:hyp1-lateral, @sec:hyp1-vertical, and @sec:hyp1-diagonal established an
amplitude-based detectability floor in all three translation directions, and
@sec:hyp1-phaseplane showed the phase-plane estimator remaining accurate at
lateral and vertical scales where that amplitude floor has already been
reached.

#draftnote[state the improvement factor explicitly: floor scale found in
@sec:tl-detectability, @sec:vtl-detectability, and @sec:dtl-detectability
versus the smallest scale at which @sec:tlp-horizontal-validation and
@sec:tlp-vertical still recover $(#Dz, #Dx)$ within tolerance. Note that the
diagonal phase-plane fit itself is still outstanding (@sec:tlp-diagonal), so
Hypothesis 1's diagonal claim is currently supported only at the amplitude
level --- flag this explicitly rather than overstating the result.]

=== Lateral/Vertical Asymmetry

@sec:th-duality predicted, from first principles, that a migrated image
separates lateral spatial frequency like a prism but not vertical frequency.

#draftnote[state whether this asymmetry was visible only in the
instantaneous-phase analysis (@fig:tlp-instphase-horiz-summary vs.
@fig:tlp-instphase-vert-summary) or also at the amplitude level
(@fig:tl-summary vs. @fig:vtl-summary), and what that implies for a field
deployment that cares more about one direction than the other. Discuss whether
the diagonal detectability floor of @sec:dtl-detectability sits closer to the
(tighter) lateral floor or the (looser) vertical one, as a further test of
the same asymmetry.]

=== Hypothesis 1.5: Local Phase-Gradient Methods

#draftnote[state whether the local, trace-based methods of @sec:hyp1-h15
(instantaneous phase, spectral-line fitting, cross-phase spectrograms,
localised STFT) gave results consistent with the global WLS fit wherever both
were computed, and whether Hypothesis 1.5 can be considered supported,
exploratory-but-promising, or unresolved given that no direct quantitative
comparison between the two has yet been run.]

== Hypothesis 2: Back-Propagation Is the Preferred Method

@sec:hyp2-fluidfront showed that the geometry/material-change decoupling of
@sec:th-material-change survives contact with a non-rigid, spatially
distributed target, and @sec:ff-factor2 identified a systematic factor-of-two
error that affects Kirchhoff and Gazdag but not back-propagation.
@sec:hyp2-noise then showed that back-propagation with sign-bit time-reversal
is the most noise-robust of the three algorithms tested.

#draftnote[state explicitly, for each migration method in @fig:ff-phaseplane,
whether the fitted slopes $(#Dz, #Dx)$ stayed near zero while the intercept
$c$ tracked the front advance as predicted by @eq:intercept-material. State
the factor-of-two ratio numerically (from @sec:ff-factor2). State whether the
noise robustness conclusion is supported quantitatively (currently only the
noisy lateral case, @fig:tlp-noise, has a phase-plane fit; the vertical,
diagonal, and fluid-front noisy fits are amplitude-level only --- flag this
as a limitation).]

=== Implications for Migration Choice in Practice

The recommendation to use back-propagation has a practical cost: it requires
a full gprMax forward simulation for each profile pair, which is substantially
more expensive than Kirchhoff or Gazdag. For field surveys with many profiles,
this cost must be weighed against the benefits identified in @ch:hyp2.

#draftnote[add a brief practical recommendation: under what conditions
(target type, noise level, available compute) should a practitioner use
Kirchhoff/Gazdag with the factor-of-two correction versus back-propagation
with sign-bit time-reversal? State whether the correction is reliable enough
to use Kirchhoff/Gazdag when back-propagation is too expensive.]

== Hypothesis 3: Generalisation to Field Data

@sec:hyp3-fielddata applied the full pipeline to real borehole GPR data and
obtained physically interpretable displacement estimates across four
operational stages of a fluid-injection experiment (@tab:fielddata-stages).
The Push stage produced a downward displacement of approximately $1.41 "m"$,
consistent with active injection; the Wait stage produced near-zero
displacement, as expected; and the Pull stage only partially reversed the Push,
leaving a net residual.

#draftnote[state whether the field-data estimates are consistent with any
independent ground-truth available from the field experiment (e.g. injection
volume, borehole depth, known fracture geometry), and whether the factor-of-two
correction of @sec:ff-factor2 was applied to the Gazdag estimates in
@tab:fielddata-stages. Note explicitly that the complex synthetic model
component of @sec:hyp3-complex is still pending, so Hypothesis 3 is only
partially supported at this stage.]

== Limitations

#draftnote[list concrete limitations: (1) synthetic gprMax data only for
Hypotheses 1--2, pending Hypothesis 3's complex synthetic results; (2)
back-propagation not run on the noisy lateral dataset due to computational
cost (@sec:hyp3-lateral); (3) the domain-size and STFT-scale-labelling
uncertainties flagged in @ch:methodology and @sec:tlp-stft; (4) the
fluid-front model uses a single idealised thin-layer geometry rather than a
swept range of fracture thicknesses or contrasts; (5) the diagonal experiment
(@sec:hyp1-diagonal) sweeps only a single fixed $2:1$ lateral-to-vertical
ratio rather than a range of diagonal angles; (6) the field-data back-prop
runs are incomplete (11/37 profiles), so the preferred method cannot yet be
fully applied to the real data.]

== Outlook

The field-data results of @sec:hyp3-fielddata demonstrate that the pipeline is
not limited to synthetic data, but several extensions are needed before field
deployment can be recommended without reservation:

#draftnote[expand with concrete next steps: (1) complete the outstanding
back-propagation runs for the remaining 26 field profiles and confirm that
the factor-of-two-corrected Gazdag estimates agree with the back-prop
estimates, validating the correction for field data; (2) complete the complex
synthetic experiment of @sec:hyp3-complex and test whether local application
of the phase-plane fit can separate multiple simultaneously-moving scatterers;
(3) extend the fluid-front model of @sec:hyp2-fluidfront to a swept range of
fracture thickness and fluid contrast, and to a front geometry that is not
perfectly straight; (4) quantify the velocity-recalibration procedure needed
when soil or ice moisture genuinely changes between baseline and monitor
surveys, since this is otherwise indistinguishable from a true vertical shift
$#Dz$ in any of the three migration algorithms.]
