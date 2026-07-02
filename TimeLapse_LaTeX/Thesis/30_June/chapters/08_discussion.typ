#import "../template.typ": *

= Discussion <ch:discussion>

#draftnote[this chapter is a placeholder skeleton built around the draft-note
prompts left in @ch:hyp1, @ch:hyp2, and @ch:hyp3. Resolve those prompts
first (they ask for specific numbers read off specific figures), then rewrite
this chapter as connected prose rather than a list.]

This chapter returns to the unifying hypothesis of @ch:introduction --- that
time-lapse GPR can decouple and accurately track subwavelength material
substitutions and mechanical fluid-front movements using noise-robust
multi-dimensional phase-plane regression and sign-bit time-reversal, even
under heavy-tailed Laplace noise --- and integrates the results of Hypotheses
1--3 to evaluate it directly.

== Hypothesis 1: Translation Tracking

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

== Hypothesis 2: Fluid-Front Tracking and Material Change

@sec:ff-phaseplane is the key test of whether the geometry/material-change
decoupling of @sec:th-material-change survives contact with a non-rigid,
spatially distributed target.

#draftnote[state explicitly, for each migration method, whether the fitted
slopes $(#Dz, #Dx)$ stayed near zero while the intercept $c$ tracked the
front advance as predicted --- this is the central claim of Hypothesis 2 and
should be stated as a number, not left implicit in a figure reference.]

== Hypothesis 3: Noise Robustness

@sec:hyp3-purenoise showed that Kirchhoff and Gazdag respond differently to
pure noise (false-coherent bands versus incoherent texture), and
@sec:hyp3-signbit introduced sign-bit time-reversal as the noise-robust
excitation scheme for back-propagation.

#draftnote[state whether the amplitude-level noisy results of
@sec:hyp3-lateral, @sec:hyp3-vertical, @sec:hyp3-diagonal, and
@sec:hyp3-fluidflow support the claim that "the right migration technique"
(per Hypothesis 3) makes noisy detection possible, and be explicit that ---
as flagged repeatedly in @ch:hyp3 --- the quantitative phase-plane-regression
evidence for this currently exists only for the lateral case
(@fig:tlp-noise); the vertical, diagonal, and fluid-front cases are
amplitude-level only until the outstanding work of @sec:hyp3-gaps is
completed. Hypothesis 3 should therefore be reported as _partially_ supported,
not fully supported, at the current state of the thesis.]

== Hypothesis 4: Future Work <sec:disc-hyp4>

Hypothesis 4 --- that the method generalises to complex synthetic scenes with
multiple, independently-moving scatterers and non-uniform fluid fronts, and to
real field data --- has not been tested in this thesis. Validating it would
require, at minimum: (1) a synthetic model containing several scatterers at
different depths moving in different directions and amounts simultaneously, to
test whether the phase-plane fit (which assumes a single dominant displacement
within its ROI) can be applied locally enough to separate multiple independent
events; (2) a fluid-front geometry that is not a straight line, to test
whether the front-tracking approach of @ch:hyp2 generalises beyond the
idealised straight front used there; and (3) a real zero-offset field survey,
ideally one with a fluid front advancing away from a borehole under at least
partially known conditions, to test the field pre-processing assumptions
(dewow, time-zero correction, trace re-binning) that the synthetic data in
this thesis sidesteps by construction. None of this exists yet in this
project; Hypothesis 4 is recorded here as the clearest direction for future
work, not as a result.

== Limitations

#draftnote[list concrete limitations: synthetic gprMax data only (no field
validation, pending Hypothesis 4); back-propagation migration not run on the
noisy lateral dataset due to computational cost (@sec:hyp3-lateral); the
domain-size and STFT-scale-labelling uncertainties flagged in @ch:methodology
and @sec:tlp-stft; the fluid-front model uses a single idealised thin-layer
geometry rather than a swept range of fracture thicknesses or contrasts; the
diagonal experiment (@sec:hyp1-diagonal) sweeps only a single fixed $2:1$
lateral-to-vertical ratio rather than a range of diagonal angles.]

== Outlook for Field Application

#draftnote[expand with concrete next steps, building on @sec:disc-hyp4's
Hypothesis 4 discussion: (1) validate the pipeline on a real zero-offset
field survey with a known, controlled sub-wavelength displacement (e.g. a
target on a calibrated micrometre stage) to test the field pre-processing
assumptions of @sec:meth-phaseplane that the synthetic data in this thesis
sidesteps by construction; (2) complete the outstanding quantitative
noise-robustness work of @sec:hyp3-gaps; (3) extend the fluid-front model of
@ch:hyp2 to a swept range of fracture thickness and fluid contrast, and to a
front geometry that is not perfectly straight; (4) quantify the
velocity-recalibration procedure needed when soil or ice moisture genuinely
changes between baseline and monitor surveys, since @sec:meth-phaseplane noted
that this is otherwise indistinguishable from a true vertical shift $#Dz$.]
