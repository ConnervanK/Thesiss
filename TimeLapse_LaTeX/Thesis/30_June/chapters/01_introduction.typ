#import "../template.typ": *

= Introduction <ch:introduction>

#draftnote[Placeholder]

// Ground-penetrating radar (GPR) images subsurface structure by emitting a
// short electromagnetic pulse and recording its reflections from dielectric
// contrasts in the ground. The achievable image resolution is fundamentally
// limited by the wavelength of the probing pulse: two reflectors closer than
// roughly half a wavelength apart, or a single reflector that moves by a small
// fraction of a wavelength between two surveys, cannot be distinguished from
// the migrated _amplitude_ image alone (@sec:meth-resolution). This
// is a practical obstacle for time-lapse monitoring applications --- tracking
// millimetre-scale ground movement, or the advance of a fluid front through a
// sub-wavelength fracture --- where the displacement or material change of
// interest is, by construction, far smaller than the wavelength of the radar
// pulse used to image it.

// This thesis develops and validates a _phase_-based alternative: rather
// than reading displacement off the migrated amplitude image, a baseline and a
// monitor survey are migrated and compared in the two-dimensional Fourier
// domain, where the Fourier shift theorem turns any sub-wavelength translation
// into a fully resolvable linear phase ramp (@ch:theory). A weighted
// least-squares fit of this phase plane recovers the displacement with
// sub-millimetre precision, and --- critically --- separates a purely
// geometric shift from a phase rotation caused by a change in the dielectric
// properties of the target itself, which is the signature of, for example, a
// fracture filling with water. Where the back-propagation migration algorithm
// is used, the method is further combined with _sign-bit time-reversal_,
// a noise-robust excitation scheme that keeps the algorithm usable even when
// the data are heavily corrupted by noise.

== Research Question

#para-head[Research Question.] Can time-lapse ground-penetrating radar
accurately track subwavelength movement --- achieving a form of
super-resolution --- by analysing _phase_ changes in migrated images rather
than their amplitude?

This thesis addresses this question through three specific, individually
testable hypotheses, each the subject of one experimental chapter:

/ Hypothesis 1 (@ch:hyp1): Phase changes in the two-dimensional
  Fourier domain of time-lapse migrated images allow inferring sub-wavelength
  displacements in the lateral, vertical, and diagonal directions, down to
  scales where amplitude differencing has already failed. An optional
  Hypothesis 1.5 explores whether local phase gradients
  $partial phi \/ partial x, thin partial phi \/ partial y$ give an
  equivalent, simpler alternative.

/ Hypothesis 2 (@ch:hyp2): Back-propagation migration with sign-bit
  time-reversal is the most noise-robust technique for time-lapse phase-plane
  tracking: it suppresses impulsive Laplace noise while Kirchhoff creates
  false-coherent artefacts and Gazdag adds incoherent speckle.

/ Hypothesis 3 (@ch:hyp3): The time-lapse phase-plane approach generalises
  beyond idealised single-scatterer synthetic models to complex scenes with
  multiple independently-moving scatterers and to real borehole GPR field data.

== Outline

@ch:litreview reviews the literature this thesis builds on.
@ch:theory covers the theoretical background (migration algorithms,
phase-plane shift estimation) and the shared simulation and processing
methodology, including a validation of the amplitude resolution floor.
@ch:hyp1 tests Hypothesis 1 (and the optional Hypothesis 1.5) on lateral,
vertical, and diagonal sub-wavelength translation of a point scatterer.
@ch:hyp2 tests Hypothesis 2: which migration technique is most robust to
heavy-tailed Laplace noise, concluding that back-propagation with sign-bit
time-reversal is the preferred method. @ch:hyp3 tests Hypothesis 3 by
applying the full pipeline to complex synthetic scenes and to real borehole
GPR field data. @ch:discussion integrates the results of all three
hypotheses into a final answer to the Research Question, and the Summary
restates the main conclusions.
