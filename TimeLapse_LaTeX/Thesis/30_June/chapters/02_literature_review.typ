#import "../template.typ": *

= Review of the Literature <ch:litreview>

#draftnote[this chapter is a skeleton built only from the references already
collected in `references.bib` (all currently placeholder entries that
must be verified against the published record before submission). It is
deliberately scoped to the material needed to justify the unifying
hypothesis of @ch:introduction and the four specific hypotheses tested
in @ch:hyp1, @ch:hyp2, @ch:hyp3 --- add literature on (1) time-lapse/4D
geophysical monitoring, (2) phase-based displacement sensing in other
modalities (e.g. InSAR), and (3) weighted least-squares plane-fitting
precedents, none of which are yet represented in the bibliography.]

== EM Theory and GPR Physics

Ground-penetrating radar images subsurface dielectric contrasts by recording
the reflections of a transmitted electromagnetic pulse; the underlying
electromagnetic principles, antenna behaviour, and propagation losses that
govern the achievable resolution and penetration depth of a survey are
covered in @annan2009 and @daniels2004. Both texts establish
the half-wavelength rule of thumb for resolving two closely spaced
reflectors, which motivates the resolution-limit validation of
@sec:meth-resolution.

== Resolution Problems (Lateral and Vertical) and the Fresnel Zone <sec:lit-resolution>

Two distinct resolution criteria limit what a GPR survey can image, and both
are tested empirically in this thesis.

#para-head[Lateral resolution.] The lateral resolving power of an unmigrated
survey is governed by the _Fresnel zone_ --- the region of the
subsurface from which reflected energy arrives within half a wavelength of
the direct reflection and therefore interferes constructively at the
receiver. Its radius at depth $d$ for a wavelength $lambda$ is
approximately $a_F approx sqrt(lambda d \/ 2)$ @annan2009; two
reflectors closer together than this radius cannot be distinguished without
migration, since their Fresnel zones overlap. Migration (@ch:theory)
collapses this zone back towards a point and is, in principle, capable of
restoring lateral resolution to a fraction of a wavelength --- the limit
quantified empirically in @sec:meth-resolution and tested as a
time-lapse detection problem in @ch:hyp1.

#para-head[Vertical resolution.] Vertical (depth) resolution is instead
governed by the _temporal_ bandwidth of the source pulse: two
reflectors are distinguishable only once their two-way travel-time
separation exceeds roughly one quarter of the dominant period, equivalently
about half a wavelength in depth for the broadband Ricker pulse used
throughout this thesis @daniels2004. Unlike lateral resolution, this
limit is not improved by migration, since migration redistributes energy
spatially but does not compress the pulse in time --- a point that motivates
why the vertical and lateral displacement-detection problems of
@ch:hyp1 are treated as genuinely distinct experiments rather than as
trivial rotations of one another (@sec:th-duality).

== GPR Processing

#para-head[Migration algorithms.] Three families of migration algorithm are
used throughout this thesis, and each has an established lineage in the
geophysical literature. @claerbout1985 introduces the
exploding-reflector model that underlies the zero-offset migration
convention adopted in @ch:theory. @schneider1978 formulates
Kirchhoff (integral, delay-and-sum) migration as the adjoint of a forward
modelling operator, the formulation implemented here via the PyLops
Kirchhoff operator. @gazdag1978 develops the
frequency--wavenumber phase-shift migration algorithm used as an
independent wave-equation-based cross-check; @stolt1978 is the
related constant-velocity Fourier-domain migration that the phase-shift
method generalises to depth-variable continuation. Together these three
algorithms make the resolution and displacement-detection results of
@ch:hyp1, @ch:hyp2, and @ch:hyp3 independent of any single migration
assumption.

#para-head[Simulation tooling.] All synthetic data in this thesis are
generated with gprMax @gprmax, an open-source finite-difference
time-domain electromagnetic simulator purpose-built for GPR. Using a
full-wave simulator rather than a ray-based or convolutional forward model
lets the back-propagation (time-reversal) migration of @ch:theory be
tested directly against the same finite-difference physics used to generate
the data, removing modelling mismatch as a confound when comparing migration
algorithms.

== Phase-Based Interpretation

The use of signal _phase_, rather than amplitude, as the primary
observable is the central methodological choice of this thesis.
@castagna2016 demonstrates phase decomposition as a tool for
seismic interpretation, separating amplitude and phase information that is
otherwise conflated in a conventional seismic or GPR trace; the
phase-decomposition routines referenced in this thesis's processing pipeline
(@ch:methodology) build on this idea. This thesis extends phase-based
analysis from single-survey interpretation to a _cross_-survey,
time-lapse setting, where the quantity of interest is not the phase of a
single image but the phase difference between two migrated images of the
same subsurface region.

The following chapter derives the migration and phase-plane theory
introduced qualitatively above in full, starting from the exploding-reflector
model and the three migration algorithms reviewed here, and building up to
the weighted least-squares phase-plane estimator validated in
@sec:hyp1-phaseplane.
