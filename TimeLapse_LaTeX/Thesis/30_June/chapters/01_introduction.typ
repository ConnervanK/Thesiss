#import "../template.typ": *

= Introduction <ch:introduction>

#draftnote[this chapter is a placeholder skeleton. Expand with a precise
statement of the research questions before submission.]

Ground-penetrating radar (GPR) images subsurface structure by emitting a
short electromagnetic pulse and recording its reflections from dielectric
contrasts in the ground. The achievable image resolution is fundamentally
limited by the wavelength of the probing pulse: two reflectors closer than
roughly half a wavelength apart, or a single reflector that moves by a small
fraction of a wavelength between two surveys, cannot be distinguished from
the migrated _amplitude_ image alone (@sec:meth-resolution). This
is a practical obstacle for time-lapse monitoring applications --- tracking
millimetre-scale ground movement, or the advance of a fluid front through a
sub-wavelength fracture --- where the displacement or material change of
interest is, by construction, far smaller than the wavelength of the radar
pulse used to image it.

This thesis develops and validates a _phase_-based alternative: rather
than reading displacement off the migrated amplitude image, a baseline and a
monitor survey are migrated and compared in the two-dimensional Fourier
domain, where the Fourier shift theorem turns any sub-wavelength translation
into a fully resolvable linear phase ramp (@ch:theory). A weighted
least-squares fit of this phase plane recovers the displacement with
sub-millimetre precision, and --- critically --- separates a purely
geometric shift from a phase rotation caused by a change in the dielectric
properties of the target itself, which is the signature of, for example, a
fracture filling with water. Where the back-propagation migration algorithm
is used, the method is further combined with _sign-bit time-reversal_,
a noise-robust excitation scheme that keeps the algorithm usable even when
the data are heavily corrupted by noise.

== Unifying Hypothesis

The series of experiments in this thesis tests a single unifying hypothesis,
which no individual experiment can test on its own:

#block(inset: (left: 2em, right: 2em), above: 1em, below: 1em)[
  _Time-lapse ground-penetrating radar can successfully decouple and
  accurately track subwavelength material substitutions and mechanical fluid
  front movements using (noise-robust) multi-dimensional phase-plane
  regression, and sign-bit time-reversal, even in environments heavily
  corrupted by heavy-tailed Laplace noise._
]

This is broken down into four specific, individually testable hypotheses,
each the subject of one experimental chapter:

/ Hypothesis 1 (@ch:hyp1): Multi-dimensional phase-plane
  regression allows tracking the subwavelength translation of a scatterer
  in the lateral, vertical, and diagonal directions. An optional
  Hypothesis 1.5 explores whether the local phase gradients
  $partial phi \/ partial x, thin partial phi \/ partial y$ give an
  equivalent, simpler alternative to the global plane fit.

/ Hypothesis 2 (@ch:hyp2): Multi-dimensional phase-plane
  regression allows tracking the subwavelength translation of a fluid
  front and inferring the associated material change in a subwavelength
  thin fracture.

/ Hypothesis 3 (@ch:hyp3): Multi-dimensional phase-plane
  regression can deal with noisy data, given that the right migration
  technique --- specifically, sign-bit time-reversal for back-propagation
  --- is used.

/ Hypothesis 4: Multi-dimensional phase-plane regression can deal with
  complex field data containing multiple scatterers moving in different
  directions, and non-uniform fluid fronts, inferring every event
  separately. This hypothesis has not yet been tested against real field
  data or complex synthetic scenes; it is stated here for completeness and
  revisited as future work in @ch:discussion.

== Outline

@ch:litreview reviews the literature this thesis builds on.
@ch:theory derives the migration and phase-plane theory used
throughout the thesis. @ch:methodology describes the shared simulation
and processing pipeline, including a validation of the amplitude resolution
floor of the migration algorithms used. @ch:hyp1 tests Hypothesis 1
(and the optional Hypothesis 1.5) on lateral, vertical, and diagonal
sub-wavelength translation of a point scatterer. @ch:hyp2 tests
Hypothesis 2 on a spatially distributed fluid front. @ch:hyp3 tests
Hypothesis 3 under heavy-tailed Laplace noise, across migration techniques.
@ch:discussion closes with a discussion of the results in light of the
unifying hypothesis above --- including Hypothesis 4, which remains
untested --- and the Summary restates the thesis's main conclusions.
