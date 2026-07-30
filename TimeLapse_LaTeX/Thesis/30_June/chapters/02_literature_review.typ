#import "../template.typ": *

= Review of the Literature <ch:litreview>

// #draftnote[this chapter is a skeleton built only from the references already
// collected in `references.bib` (all currently placeholder entries that
// must be verified against the published record before submission). It is
// deliberately scoped to the material needed to justify the unifying
// hypothesis of @ch:introduction and the four specific hypotheses tested
// in @ch:hyp1, @ch:hyp2, @ch:hyp3 --- add literature on (1) time-lapse/4D
// geophysical monitoring, (2) phase-based displacement sensing in other
// modalities (e.g. InSAR), and (3) weighted least-squares plane-fitting
// precedents, none of which are yet represented in the bibliography.]

// == EM Theory and GPR Physics

// Ground-penetrating radar is governed by the same four Maxwell equations that
// describe all classical electromagnetism, specialised to a low-loss dielectric
// medium through the constitutive relations $bold(J) = sigma bold(E)$,
// $bold(D) = epsilon bold(E)$, and $bold(B) = mu bold(H)$ @annan2009. Combining
// Faraday's and Ampère's laws for a source-free, non-magnetic medium eliminates
// the magnetic field and yields a transverse wave equation for $bold(E)$ whose
// loss term is set by conductivity $sigma$ and whose storage term is set by
// permittivity $epsilon$ @annan2009 @slob2024. When displacement currents
// dominate over conduction currents --- true for GPR frequencies of 10 MHz to a
// few GHz in most earth materials --- the equation supports a propagating,
// weakly attenuated plane wave rather than the purely diffusive response that
// underlies DC resistivity and low-frequency EM sounding methods
// @milsom2011electric. In this propagating regime the wave speed and
// attenuation reduce to $v = c\/sqrt(kappa)$ and
// $alpha = (omega\/2)sqrt(mu epsilon_0 kappa) tan delta$, where
// $kappa = epsilon\/epsilon_0$ is the relative permittivity, or "dielectric
// constant," of the medium and $tan delta = sigma\/(omega epsilon)$ is the loss
// tangent @annan2009 @milsom2011gpr. Because $kappa$ alone sets the
// propagation velocity, and an impedance (and hence reflection-coefficient)
// contrast $Delta Z prop Delta(1\/sqrt(kappa))$ at any interface governs the
// amplitude of a returned signal, GPR is fundamentally a survey of subsurface
// permittivity structure rather than of the conductivity structure probed by
// galvanic and inductive electrical methods @milsom2011electric
// @milsom2011gpr; this distinction motivates treating $kappa$, not $sigma$, as
// the primary state variable whose changes this thesis seeks to detect
// time-lapse.

// The dielectric constant of a natural material is rarely a fixed rock or soil
// property; it is dominated by the volume fraction and phase of water present
// in the pore space, since liquid water's permittivity ($kappa approx 80$) is
// an order of magnitude larger than that of the dry mineral matrix
// ($kappa approx 3$–$8$) or of air ($kappa = 1$) @annan2009. Empirical mixing
// relationships such as the Topp equation exploit this contrast to invert GPR
// velocity for volumetric water content, but the same sensitivity means that
// any process which redistributes water --- infiltration, drainage, or a
// change of phase between liquid and ice --- produces a measurable shift in
// $kappa$ and hence in two-way travel time and reflection amplitude, even when
// the bulk geometry of the target is unchanged. Chen et al. @chen2023
// demonstrate this directly for the water-ice-snow system relevant to
// freeze-thaw monitoring, showing that the relative permittivities of ice and
// water are frequency-dependent and separable, which is precisely the physical
// mechanism by which a subsurface interface can appear, sharpen, or fade
// between repeat GPR surveys without any change in the scatterer's shape or
// position. It is this permittivity-driven, rather than purely geometric,
// origin of GPR reflectivity that justifies treating time-lapse amplitude and
// phase changes as evidence of a genuine subsurface state change, and that
// underlies the reflection and resolution formalism developed in
// @sec:lit-resolution.

// == Resolution Problems (Lateral and Vertical) and the Fresnel Zone <sec:lit-resolution>

// Ground-penetrating radar resolution splits into two independent limits, one
// lateral (or horizontal) and one vertical (or range/depth), and both are
// active research questions rather than settled constants, since the
// literature disagrees on the correct closed-form expression for each even
// though the underlying physical mechanisms are well established
// @milsom2011gpr @rial2007.

// #para-head[Lateral resolution and the Fresnel zone.] An unmigrated GPR
// antenna does not illuminate a single point on a reflecting interface; it
// illuminates a finite footprint whose extent is set by the *first Fresnel
// zone* --- the region of the interface from which returning energy arrives
// within half a wavelength ($lambda\/2$) of the shortest, normal-incidence
// path, and therefore sums constructively at the receiver
// @perezgracia2008 @rial2007. Two point reflectors buried at the same depth
// cannot be resolved as separate anomalies once their separation is smaller
// than the diameter of this footprint. The classical closed-form radius,
// derived from the phase-difference geometry of @perezgracia2008 and used
// throughout the antenna-footprint literature, is
// $ r_F = sqrt(h_1 h_2 v\/(f L)) = sqrt(h v \/ (2f)) = (v\/2) sqrt(t\/f) $
// for a target at depth $h$, two-way travel time $t$, velocity $v$, and
// dominant frequency $f$ @perezgracia2008. Noon et al. treat the same
// footprint area as the effective radar cross-section of a "rough planar"
// target in the radar range equation, linking the Fresnel-zone concept
// directly to detectability rather than only to resolvability
// @noon1998. A competing, elliptical formulation that explicitly incorporates
// the host medium's relative permittivity $epsilon_r$,
// $ A = lambda\/4 + h \/ sqrt(epsilon_r - 1), quad B = A\/2, $
// is preferred by some authors because it accounts for the narrowing of the
// illumination pattern with increasing permittivity @rial2007
// @perezgracia2008. Laboratory comparisons in water and sand show that none
// of the competing closed forms match experiment exactly and that the
// achieved lateral resolution depends strongly on the acceptable level of
// inter-target interference, with the fully separated, interference-free
// distance $D_2$ running roughly 1.5--3$times$ larger than the
// first-hyperbola-visible distance $D_1$ for the same target pair
// @perezgracia2008. Since migration (@ch:theory) is, in principle, capable of
// collapsing the Fresnel zone back toward a diffraction-limited point, this
// scatter in the unmigrated literature values motivates treating the
// achievable post-migration lateral resolution as an empirical quantity to be
// measured directly, rather than assumed from a single formula, in
// @sec:meth-resolution and @ch:hyp1.

// #para-head[Vertical resolution and the Rayleigh criterion.] Vertical
// resolution is governed instead by the temporal bandwidth of the received
// pulse. The standard criterion, adapted from Rayleigh's radar-range
// resolution, states that two equal-strength reflectors are separated in
// range once their travel-time difference exceeds $Delta R = v\/(2B)$, where
// $B$ is the bandwidth of the *received* (not transmitted) signal
// @noon1998. Because attenuation is frequency-dependent, $B$ narrows with
// depth, so a fixed-bandwidth rule of thumb systematically overstates
// resolution at greater penetration depths @noon1998. Experimental
// calibration with bow-tie antennas confirms that vertical resolution in
// practice tracks the effective pulse duration of the source wavelet rather
// than the nominal antenna centre frequency, and additionally depends on the
// electrical contrast of the shallower reflector: a more conductive upper
// interface attenuates more of the signal reaching the lower one, degrading
// resolvability independently of bandwidth @rial2007. Unlike lateral
// resolution, this vertical limit is not improved by migration, since
// migration relocates energy spatially but does not compress the pulse in
// time --- which is why the vertical and lateral displacement-detection
// problems of @ch:hyp1 are treated as genuinely distinct experiments rather
// than as trivial rotations of one another (@sec:th-duality).

// #para-head[Sub-wavelength information in the phase.] The Fresnel-zone and
// Rayleigh limits above both describe the resolution of *discrete,
// separately identifiable* targets from amplitude data, and it is this
// amplitude-based, geometric notion of resolution that migration can, at
// best, push to a fraction of a wavelength. Tsoflias and Hoch show that a
// sub-wavelength layer --- far too thin to be resolved as two separate
// interfaces by the Rayleigh criterion above --- still imprints a
// systematic, quantifiable signature on the *phase* of a transmitted or
// reflected GPR wave, with the sign and magnitude of the polarisation-dependent
// phase shift varying continuously with layer thickness and fluid content
// even when no amplitude criterion would flag the layer as resolvable
// @tsoflias2006. This distinction between amplitude-limited geometric
// resolution and phase-encoded sub-wavelength information is the physical
// basis for treating phase, rather than amplitude, as the primary observable
// of this thesis (@ch:hyp1, @ch:hyp2, @ch:hyp3): the Fresnel and Rayleigh limits
// reviewed here bound what a conventional, amplitude-based reading of a
// migrated image can distinguish, not what the underlying wavefield actually
// encodes about sub-wavelength subsurface change.

// == GPR Processing

// #draftnote[Placeholder]

// // #para-head[Migration algorithms.] Three families of migration algorithm are
// // used throughout this thesis, and each has an established lineage in the
// // geophysical literature. @claerbout1985 introduces the
// // exploding-reflector model that underlies the zero-offset migration
// // convention adopted in @ch:theory. @schneider1978 formulates
// // Kirchhoff (integral, delay-and-sum) migration as the adjoint of a forward
// // modelling operator, the formulation implemented here via the PyLops
// // Kirchhoff operator. @gazdag1978 develops the
// // frequency--wavenumber phase-shift migration algorithm used as an
// // independent wave-equation-based cross-check; @stolt1978 is the
// // related constant-velocity Fourier-domain migration that the phase-shift
// // method generalises to depth-variable continuation. Together these three
// // algorithms make the resolution and displacement-detection results of
// // @ch:hyp1, @ch:hyp2, and @ch:hyp3 independent of any single migration
// // assumption.

// // #para-head[Simulation tooling.] All synthetic data in this thesis are
// // generated with gprMax @gprmax, an open-source finite-difference
// // time-domain electromagnetic simulator purpose-built for GPR. Using a
// // full-wave simulator rather than a ray-based or convolutional forward model
// // lets the back-propagation (time-reversal) migration of @ch:theory be
// // tested directly against the same finite-difference physics used to generate
// // the data, removing modelling mismatch as a confound when comparing migration
// // algorithms.

This review builds, step by step, the argument that motivates the phase-based
time-lapse method of this thesis. It establishes four things in turn: that GPR
reflectivity is a sensitive probe of subsurface water and permittivity
(@sec:lit-fundamentals); that the resolution of a single GPR image is
nonetheless bounded by a wavelength-scale floor, so sub-wavelength change is
invisible in one survey (@sec:lit-resolution); that migration and its
amplitude-based refinements relocate and sharpen energy but never escape that
floor (@sec:lit-amplitude); and that the _phase_ of the wavefield — which every
one of those methods computes and then discards — still carries the
sub-wavelength information amplitude cannot (@sec:lit-phase). The last point is
the pivot on which the rest of the thesis turns.

== Fundamentals of GPR <sec:lit-fundamentals>

Ground-penetrating radar is a _wave-field_ method: unlike the diffusive electrical techniques, it propagates a true wave whose velocity (c = (εμ)⁻¹/²) and reflection strength are set by the subsurface electrical properties (Slob, 2024; Annan, 2009). The property that matters most here is permittivity, and what most controls permittivity is water: liquid water has a relative permittivity of about 80 against 3–8 for most dry minerals, so water content dominates both the wave velocity and the reflection coefficients generated at interfaces (Milsom, 2011). This is what makes GPR a sensitive probe of fluid-related change — fluid entering a fracture, say — and it is the physical basis for reading a time-lapse change in the reflected wavefield as evidence of a subsurface change in water content. Underlying this, GPR obeys Maxwell's equations through the constitutive parameters σ, ε, and μ (Annan, 2009), and is effective only in low-loss media where displacement currents dominate over ohmic dissipation, across the HF–UHF bands (roughly 10 MHz–4 GHz).

== Resolution problems and their bearing on sub-wavelength movement <sec:lit-resolution>

The spatial resolution of a single GPR image is bounded by a wavelength-scale floor: static features smaller than the dominant wavelength—such as the separation between two scatterers or the width of a fluid-filled fracture—cannot be resolved as distinct entities. They remain embedded within a single Fresnel zone or range cell. This static resolution limit is the central obstacle this thesis addresses, motivating the shift to _time-lapse acquisition_. Rather than attempting to resolve sub-wavelength structures in a single baseline survey, repeated monitoring surveys exploit the small but detectable wavefield changes produced by moving scatterers or fluid advance—changes heavily amplified by the strong permittivity contrast of water (Milsom, 2011).
This resolution floor consists of two independent components. Vertical (range) resolution ($Delta r$) is bounded by the pulse width at half amplitude ($W$) and the material velocity ($v$), such that $Delta r >= (W v) / 4$. Lateral resolution ($Delta l$) is governed by the Fresnel zone and degrades with target depth ($r$), following relations such as $Delta l >= sqrt((v r W) / 2)$ (Annan, 2009). By relating $W$ to the bandwidth ($B$) and centre frequency ($f_c$) as $W approx 1/B approx 1/f_c$, and defining the centre wavelength as $lambda_c = v/f_c$, the lateral limit can be expressed as $Delta l = sqrt((r lambda_c) / 2)$. Consequently, two point reflectors buried at the same depth cannot be resolved as separate anomalies if their lateral separation is smaller than the diameter of this antenna footprint (Rial et al., 2007; Pérez-Gracia et al., 2008).
The exact closed form is not settled: Rial (2007) found a Fresnel-zone formulation best matched the observed horizontal resolution of bow-tie antennas, and Pérez-Gracia (2008) distinguished the absolute resolution of a system from the broader antenna footprint needed for clearly separated, interference-free images — both of which remain coarse relative to the scale of an individual scatterer.

== Amplitude-Based Migration: The Traditional Answer to Resolution <sec:lit-amplitude>

Migration is the classical response to the resolution problem posed above: rather than reading reflection amplitudes directly off an unfocused B-scan, it repositions recorded energy into its true subsurface location, collapsing the diffraction hyperbola generated by a point scatterer back toward a focused point wherever the source and receiver wavefields coincide in space and time (Jones, 2014). This thesis implements three such algorithms — Kirchhoff delay-and-sum migration, Gazdag frequency-wavenumber phase-shift migration, and back-propagation (time-reversal) migration — as PyLops-based operators, and each is derived in full in @ch:theory. What matters for the argument here is what all three have in common: although they rely on phase shifts and traveltimes internally to focus the wavefield—whether by weighted backprojection onto isochron surfaces (Kirchhoff), downward-continuation with a phase-shift operator (Gazdag), or time-reversed re-propagation (back-propagation) (Li, 2020; Özdemir, 2014)—they produce a single structural snapshot. Because a standard static image relies strictly on the focused amplitude to define boundaries, the resolution it can deliver is bounded by the same Fresnel-zone / Rayleigh floor identified in @sec:lit-resolution: migration relocates the footprint, but does not by itself shrink it (Jones, 2014; Li, 2020).

A large body of work tries to push this floor lower. Some methods stay within the amplitude imaging condition: deconvolution compresses the source wavelet to separate interfering arrivals (Schmelzbach, 2015; Moghaddam, 2019); reverse-time migration and its least-squares refinement replace ray-based assumptions with a full wave-equation imaging condition (Jones, 2014; Liu, 2024); and full-waveform inversion inverts the entire recorded waveform for quantitative permittivity and conductivity, in principle exceeding the Fresnel-zone limit (Meles, 2011; Wang, 2025). Others abandon the amplitude imaging condition altogether for subspace decomposition — MUSIC and related sparsity methods (Li, 2023; Gao, 2022) — the only class here that demonstrably beats the classical limit rather than merely approaching it. But that gain is bought with structure a routine survey does not supply for free: a multistatic data matrix, an engineered focusing geometry, or an assumed wavelet and background velocity model that a single zero-offset time-lapse survey cannot guarantee.

/* ── Related-work detail deferred (RESTRUCTURE_PLAN Fix #2) ─────────────────
   The deconvolution / RTM / LSRTM / FWI and MUSIC / SREMI comparison below is
   condensed into the single paragraph above. The full version is retained here
   for a possible related-work appendix or the supplementary document; reinstate
   it as live text there rather than in this chapter.

=== Pushing the Amplitude Floor Further: Deconvolution and Iterative Wave-Equation Inversion

Two further lines of work attempt to push this floor lower without abandoning the amplitude-based imaging condition itself. The first attacks the vertical (temporal) component directly: deconvolution inverse-filters a trace to compress the embedded source wavelet and recover the underlying reflectivity, in principle separating interfering arrivals that a broadened pulse would otherwise blur together (Schmelzbach, 2015; Moghaddam, 2019). Standard least-squares deconvolution assumes a minimum-phase wavelet and favours smooth rather than spiky solutions, which suits GPR data poorly; sparse (ℓ1-regularised) variants instead handle the mixed-phase GPR wavelet correctly and restore amplitude at anomalies otherwise hidden in a conventional section (Moghaddam, 2019; Schmelzbach, 2015).

The second line refines the migration operator itself. Reverse Time Migration (RTM) replaces the ray-based hyperbola of Kirchhoff migration with a full wave-equation, cross-correlation imaging condition, removing the ray-based assumptions that limit resolution under complex velocity structure (Jones, 2014; Geng, 2022; Liu, 2023). Because RTM is only an adjoint, not an inverse, operator, it still blurs and mis-scales amplitude; Least-Squares RTM (LSRTM) corrects this by explicitly inverting for the Hessian, at a computational cost that efficient one-step formulations make practical for repeat, time-lapse imaging (Liu, 2024). Full Waveform Inversion (FWI) goes further still, inverting the entire recorded waveform — not just picked reflectivity — for quantitative permittivity and conductivity, and can in principle exceed the Fresnel-zone limit altogether, provided the resulting non-linear inversion is kept from cycle-skipping (Meles, 2011; Wang, 2025). Deconvolution, RTM, LSRTM, and FWI each ask more of the same underlying amplitude and waveform information, at correspondingly higher modelling and computational cost, but none of them step outside the amplitude-based imaging condition itself.

=== Subspace and Sparsity Methods: Genuine Super-Resolution, at a Price

A distinct family of methods abandons the amplitude imaging condition entirely in favour of subspace decomposition, and is the only class reviewed here that demonstrably exceeds the classical Fresnel/Rayleigh limit rather than merely approaching it. MUSIC-based (Multiple Signal Classification) methods split the recorded data into orthogonal signal and noise subspaces via singular value decomposition and localise reflectors from the orthogonality between the noise subspace and the medium's Green's function, achieving sub-meter or finer resolution in tunnel-array and time-reversal GPR settings (Li, 2023; Karami, 2026). A related dual-sparsity approach, SREMI, instead reframes imaging as an inverse problem regularised for both sparse reflectivity and spatial continuity, resolving thin layers that one-pass RTM blurs together (Gao, 2022). These gains are real, but they are bought with additional structure that a routine field survey does not supply for free: MUSIC-type methods generally require a multistatic data matrix or an engineered focusing geometry, and sparsity-based inversion is only as good as its assumed source wavelet and background velocity model (Li, 2023; Karami, 2026; Gao, 2022) — assumptions a single zero-offset time-lapse survey, of the kind used throughout this thesis, cannot always guarantee.
─────────────────────────────────────────────────────────────────────────── */

Across this entire progression — from Kirchhoff's delay-and-sum stack to MUSIC's subspace decomposition — resolution is treated as something to be won from amplitude: by repositioning it (migration), by sharpening it (deconvolution, LSRTM), by inverting for it directly (FWI), or by decomposing it into orthogonal subspaces (MUSIC, SREMI). Yet @sec:lit-resolution already showed that while the structural boundaries of a sub-wavelength layer cannot be resolved by amplitude alone (as wave interference yields a complex, non-unique amplitude response), it still imprints a measurable, continuously varying shift on *phase* (Tsoflias and Hoch, 2006) — information every algorithm above computes internally and then discards. This suggests a different question is worth asking of a single time-lapse pair, without invoking full waveform inversion or a multistatic array: not whether amplitude can be pushed past the classical limit, but whether the phase information already present in a routine migrated image is sufficient, on its own, to detect a sub-wavelength change. That question motivates the shift to phase-based interpretation taken up next, and is the guiding hypothesis carried forward into @ch:theory and @ch:hyp1.



== Moving from Amplitude- to Phase-Based Interpretation <sec:lit-phase>

Using signal _phase_ rather than amplitude as the primary observable is the central methodological choice of this thesis. Phase decomposition is already established as a single-survey interpretation tool: Castagna (2016) separates amplitude and phase information that a conventional trace conflates, and the phase-decomposition routines in this thesis's processing pipeline build on that idea. What is new here is extending phase analysis from single-survey interpretation to a _cross_-survey, time-lapse setting, where the quantity of interest is not the phase of one image but the phase _difference_ between two migrated images of the same region. @ch:theory derives that extension in full — from the migration algorithms reviewed above, through the Fourier shift theorem, to the weighted least-squares phase-plane estimator that recovers sub-wavelength displacement from a baseline/monitor pair.

// The use of signal _phase_, rather than amplitude, as the primary
// observable is the central methodological choice of this thesis.
// @castagna2016 demonstrates phase decomposition as a tool for
// seismic interpretation, separating amplitude and phase information that is
// otherwise conflated in a conventional seismic or GPR trace; the
// phase-decomposition routines referenced in this thesis's processing pipeline
// (@ch:methodology) build on this idea. This thesis extends phase-based
// analysis from single-survey interpretation to a _cross_-survey,
// time-lapse setting, where the quantity of interest is not the phase of a
// single image but the phase difference between two migrated images of the
// same subsurface region.

// The following chapter derives the migration and phase-plane theory
// introduced qualitatively above in full, starting from the exploding-reflector
// model and the three migration algorithms reviewed here, and building up to
// the weighted least-squares phase-plane estimator validated in
// @sec:hyp1-phaseplane.
