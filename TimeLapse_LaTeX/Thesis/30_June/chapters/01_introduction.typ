#import "../template.typ": *

= Introduction <ch:introduction>

Ground-penetrating radar (GPR) images the subsurface by emitting a short
electromagnetic pulse and recording its reflections from dielectric contrasts
in the ground. The resolution of that image is fundamentally limited by the
wavelength of the probing pulse: two reflectors closer than roughly half a
wavelength, or a single reflector that moves by a small fraction of a
wavelength between two surveys, cannot be told apart from the migrated
_amplitude_ image alone. This is a practical obstacle for time-lapse
monitoring --- tracking millimetre-scale ground movement, or the advance of a
fluid front through a thin fracture --- where the displacement or material
change of interest is, by construction, far smaller than the wavelength of the
radar pulse used to image it.

This thesis develops and validates a way around that limit that reads
displacement from the _phase_ of the migrated wavefield rather than its
amplitude. The idea rests on a short chain of statements, each supported in
turn by the review of @ch:litreview and the theory of @ch:theory:

// + *GPR resolution is wavelength-limited.* The migrated amplitude image cannot
//   separate two reflectors, or localise a moved one, below roughly half a
//   wavelength --- a floor set laterally by the Fresnel zone and vertically by
  // the pulse bandwidth (@ch:litreview).

+ *GPR resolution is generally wavelength-limited,* with the Fresnel zone and pulse bandwidth setting a lateral and vertical floor of roughly half a wavelength for separating two distinct reflectors. While a single migrated amplitude image cannot resolve features below this limit, sub-wavelength localisation of a _moved_ reflector is possible using windowed cross-correlation on two migrated amplitude images, which bypasses the amplitude resolution limit by indirectly leveraging phase shifts. (@ch:litreview).

+ *The changes this thesis targets are sub-wavelength in the regime that
  matters.* Depending on the GPR centre frequency and the scale of the target, dynamic changes are not always sub-wavelength. However, for the frequencies employed in this thesis, the targeted processes—such as the millimetre-scale aperture of a fluid-filled fracture or a slowly advancing infiltration front—are physically smaller than the dominant wavelength. Because GPR reflectivity is strongly governed by the permittivity contrast of water, these specific monitoring targets produce measurable signals that remain trapped below the classical resolution floor.

+ *Migration relocates energy but does not beat the resolution floor.*
  Focusing the raw B-scan repositions recorded energy toward its true
  location, but the result is still read off as an _amplitude_ image. Neither
  migration nor its amplitude-based refinements (deconvolution, least-squares
  and full-waveform inversion) step outside that amplitude imaging condition,
  so all remain bound by the same half-wavelength resolution floor identified
  above (@ch:litreview).

+ *Phase varies linearly with a sub-wavelength shift; amplitude does not.* The envelope of a 
  migrated image of a point scatterer is a rough amplitude lobe, so a shift
  much smaller than a wavelength leaves the envelope itself essentially unchanged ---
  the same resolution floor identified above. The phase of the wave
  oscillating beneath that envelope, however, advances continuously and linearly with position. So the same shift still produces a
  small but exactly proportional, and therefore measurable, phase change
  (@ch:theory).

+ *A time-lapse pair turns sub-wavelength change into a recoverable phase
  ramp.* Comparing a baseline and a monitor survey in the two-dimensional
  Fourier domain, the Fourier shift theorem turns any sub-wavelength
  translation between the two migrated images into an exactly linear phase ramp
  whose slope _is_ the displacement (@ch:theory) --- so a shift invisible in
  amplitude becomes, in principle, exactly recoverable from phase. This is the
  sense in which _super-resolution_ is used throughout this thesis: not
  sharpening the migrated image itself, but recovering a displacement smaller
  than its classical half-wavelength resolution limit from phase information
  the image already contains.

== The Unifying Hypothesis

Taken together, these statements make one hypothesis a sensible thing to test
--- the research question of this thesis:

#para-head[Research question.] Can time-lapse ground-penetrating radar
accurately track subwavelength movement --- achieving a form of
super-resolution monitoring --- by analysing _phase_ changes in migrated images rather
than their amplitude?

This unifying hypothesis cannot be settled by any single experiment: for it to
hold, phase-based inference must work across displacement direction and scale,
survive realistic field noise, and generalise from idealised synthetic models
to real data. It is therefore tested through three specific, individually
testable hypotheses, each the subject of one experimental chapter and each
asking one of those three questions in turn:

/ Hypothesis 1 --- _does it work?_ (@ch:hyp1): Phase changes in the
  two-dimensional Fourier domain of time-lapse migrated images allow inferring
  sub-wavelength displacements in the lateral, vertical, and diagonal
  directions, down to scales where amplitude differencing has already failed.
  An optional Hypothesis 1.5 explores whether local phase gradients
  $partial phi \/ partial x, thin partial phi \/ partial y$ give an equivalent,
  simpler alternative.

/ Hypothesis 2 --- _does it survive noise, and which migration is best?_ (@ch:hyp2): Back-propagation migration with sign-bit time-reversal is the
  most noise-robust technique for time-lapse phase-plane tracking: it
  suppresses impulsive Laplace noise whereas Kirchhoff creates false-coherent
  artefacts and Gazdag adds incoherent speckle.

/ Hypothesis 3 --- _does it generalise to the field?_ (@ch:hyp3): The
  time-lapse phase-plane approach generalises beyond idealised single-scatterer
  synthetic models to real borehole GPR field data.

== Outline

@ch:litreview reviews the literature this thesis builds on, developing the five
statements above into the arguments that motivate a phase-based approach.
@ch:theory derives the theoretical background --- the migration algorithms and
the phase-plane shift-estimation method. @ch:methodology describes the shared simulation and processing pipeline that applies it to data, and closes with a validation of the amplitude resolution floor (@sec:meth-resolution) against which the phase method is measured.
@ch:hyp1 tests Hypothesis 1 (and the optional Hypothesis 1.5) on lateral,
vertical, and diagonal sub-wavelength translation of a point scatterer.
@ch:hyp2 tests Hypothesis 2: which migration technique is most robust to
heavy-tailed Laplace noise. @ch:hyp3 tests Hypothesis 3 by applying the full
pipeline to real borehole GPR field data. @ch:discussion integrates the
results of all three hypotheses into a final answer to the research question,
and the #link(<ch:summary>)[Summary] presents the main conclusions.
