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

== Fundamentals of GPR

Ground penetrating radar is governed by Maxwell's equations combined with the constitutive parameters of the subsurface — electrical conductivity (σ), dielectric permittivity (ε), and magnetic permeability (μ) — which together describe how a material responds to an applied electromagnetic field (Annan, 2009). Slob (2024) frames GPR explicitly as a wave-field method, distinct from diffusive geophysical techniques, in which the propagation velocity is dictated by the medium's electrical properties through c = (εμ)⁻¹/², and surveys are most commonly conducted in common-offset mode, with transmitter and receiver held at fixed separation and moved jointly along the profile. Operating across the HF, VHF, and UHF bands (roughly 10 MHz–4 GHz), GPR is only effective in low-loss media, where energy storage through displacement currents dominates over ohmic dissipation (Annan, 2009; Milsom, 2011). Because water has a relative permittivity of around 80 compared with 3–8 for most dry minerals, water content is the dominant control on both wave velocity and the reflection coefficients generated at subsurface interfaces (Milsom, 2011) — a sensitivity that underlies the use of GPR to detect fluid-related changes such as those occurring within fractures.

== Resolution problems and their bearing on sub-wavelength movement

A central limitation of GPR for tracking fine-scale changes is that its resolving power is fundamentally tied to the signal wavelength and depth of investigation, rather than being a fixed, arbitrarily small quantity. Vertical (range) resolution is bounded by the pulse width and material velocity (Δr ≥ Wv/4), while lateral resolution is governed by the Fresnel zone and therefore degrades with increasing target depth, following relations such as Δl ≥ √(vrW/2) (Annan, 2009). Experimental work by Rial (2007) confirmed that a Fresnel-zone formulation best matches observed horizontal resolution for bow-tie antennas, and Pérez-Gracia (2008) further distinguished between the absolute resolution capability of a system and the broader antenna footprint required for clearly separated, interference-free images — both of which remain coarse relative to the scale of individual scatterers. These constraints mean that any physical displacement of a scatterer, or migration of fluid within a fracture, that is smaller than the dominant wavelength cannot be resolved as a distinct spatial feature in a single GPR image; the target simply remains embedded within one Fresnel zone or range cell. This is precisely the motivation for timelapse approaches: rather than attempting to spatially resolve sub-wavelength motion directly, repeated acquisitions exploit small but detectable changes in the reflected wavefield — arising from the strong permittivity contrast associated with water (Milsom, 2011) — to infer movement that lies below the classical resolution limit.

== Traditional Amplitude-Based Migration and Imaging

Traditional migration algorithms reposition recorded amplitude data into their true subsurface geological locations, working from the fundamental imaging condition that a reflector exists wherever the source and receiver wavefields coincide in space and time (Jones, 2014). The most widely used amplitude-based approach is Kirchhoff migration, which performs a weighted backprojection of the recorded wavefield onto isochron surfaces of equal traveltime, a process closely related to diffraction (hyperbolic) stacking, in which waveforms are summed along a theoretical moveout curve to collapse a hyperbolic arrival into a single focused point (Li, 2020; Özdemir, 2014). In GPR imaging specifically, this is implemented under the exploding-source model, where a fictitious source is assumed to radiate at t = 0 and the wave velocity is halved to compensate for the one-way travel assumption; Kirchhoff migration improves on simple hyperbolic summation by incorporating correction factors for spherical spreading, directivity, and wave phase and amplitude, giving it better performance on weakly reflective targets at higher computational cost (Özdemir, 2014). An alternative amplitude-based scheme is Gazdag (phase-shift) migration, which iteratively applies a phase-shift operator in the frequency-wavenumber domain to downward-extrapolate the wavefield to zero time, though it is highly sensitive to the frequency sampling interval (Özdemir, 2014). Back-projection tomography, adapted from medical imaging, offers a further amplitude-based alternative that spreads one-dimensional range profiles across the image space and supports real-time processing, at the cost of increased vulnerability to numerical noise (Özdemir, 2014). Because these ray-based and stacking methods rely on the kinematic timing and amplitude of arrivals rather than exploiting signal-noise subspace separation, their achievable resolution remains inherently bounded by the size of the Fresnel zone (Jones, 2014; Li, 2020).

=== Deconvolution as an Approach to Improving Resolution

Deconvolution addresses the vertical (temporal) resolution limits of GPR by inverse-filtering the recorded trace to remove the embedded source wavelet and recover the underlying reflectivity series, correcting for distortions introduced by the transmitter pulse, antenna-ground coupling, and the Earth's inherent low-pass filtering, which broadens the wavelet in time and absorbs its higher frequencies (Schmelzbach, 2015; Moghaddam, 2019). Because wavelet compression in the time domain requires widening the signal's frequency bandwidth, successful deconvolution allows closely spaced subsurface layers and interfering events to be separated, yielding a section that more faithfully represents the true stratigraphic structure (Moghaddam, 2019; Schmelzbach, 2015). However, standard least-squares deconvolution is often poorly suited to GPR data for two reasons: it assumes a minimum-phase wavelet with front-loaded energy, whereas GPR antennas typically generate mixed-phase wavelets with energy concentrated in the center of the pulse, and its underlying ℓ2-norm formulation promotes smoothness rather than spikiness, tending to produce attenuated, noisy results rather than genuine resolution gains (Schmelzbach, 2015; Moghaddam, 2019). Sparse deconvolution methods, which instead regularize the inversion with an ℓ1-type norm to favor sparse, "spiky" reflectivity solutions, have been shown to handle mixed-phase wavelets more effectively and to restore signal amplitudes associated with anomalies that remain blurred or hidden in conventional least-squares sections (Moghaddam, 2019).

=== Wave-Based Inversion Methods: From RTM to FWI

Among wave-equation migration approaches, Reverse Time Migration (RTM) forms the baseline: it forward-propagates a synthetic source wavefield while simultaneously back-propagating the recorded data, cross-correlating the two at each time step to locate reflectors wherever the wavefields coincide, an imaging condition that allows RTM to handle steep dips and laterally varying velocities without the Fresnel-zone-limited resolution of ray-based methods (Jones, 2014; Geng, 2022; Liu, 2023). However, because standard RTM is only an adjoint — rather than an inverse — operator, it introduces "Hessian effects" such as amplitude imbalance and blurring, which Least-Squares RTM (LSRTM) addresses by explicitly inverting for the Hessian to recover true-amplitude, deblurred reflectivity images; this comes at high iterative cost, though efficient "one-step" formulations using Wiener-filter preconditioning approximate the Hessian inverse cheaply enough to make LSRTM practical for high-resolution, time-lapse monitoring of subtle subsurface changes (Liu, 2024). Full Waveform Inversion (FWI) represents a further increase in complexity and resolving power: rather than migrating amplitude picks or reflectivity alone, it inverts the entire waveform — including reflected, diffracted, and scattered energy — for quantitative maps of permittivity and conductivity through iterative nonlinear least-squares minimisation, which in principle allows it to overcome the Fresnel-zone resolution limit entirely and achieve sub-wavelength resolution, provided that non-linearity and cycle-skipping are managed through strategies such as progressive bandwidth expansion (Meles, 2011; Wang, 2025). Taken together, these three methods trace a clear progression in how resolution is treated: RTM improves geometric accuracy of imaging under complex velocity structure, LSRTM restores amplitude fidelity and deblurs the image by inverting for the imaging operator itself, and FWI pushes past the classical resolution limit altogether by exploiting the full information content of the waveform.

=== Attempts at Obtaining Super-Resolution: Subspace and Signal Classification Methods

A distinct class of methods departs from conventional amplitude-based migration by exploiting subspace decomposition rather than kinematic timing alone, allowing resolution to exceed the classical Rayleigh/Fresnel limit. MUSIC-based (Multiple Signal Classification) approaches use singular value decomposition to split the acquisition data into orthogonal signal and noise subspaces, then exploit the orthogonality between the noise subspace and the medium's Green's function to localise reflectors with a sharpness far exceeding standard migration; in tunnel-array acoustic prospecting this has enabled sub-meter detection of small geological bodies, while a confocal imaging condition — requiring excitation and reception points to coincide spatially — further sharpens the resulting pseudospectrum and improves efficiency relative to RTM (Li, 2023). Building on this principle, High-Resolution Time-Reversal MUSIC (HRTR) combines time-reversal focusing with MUSIC-based subspace separation to achieve resolutions of λ/10 or finer directly from conventional GPR A- and B-scans, without requiring a full multistatic data matrix or additional hardware; this allows shallow, near-surface targets that are normally masked by overlapping reflections to be separated, and permits low-frequency signals to be used for deeper penetration without sacrificing spatial resolution (Karami, 2026). A related but distinct strategy is pursued by SREMI, which reframes imaging as a dual-sparsity inverse problem, promoting sparsity in the reflectivity model (via an L1 norm) alongside spatial continuity (via total variation), solved iteratively using a dual-sparsity FISTA algorithm; this produces markedly sharper images with a wider vertical wavenumber content than one-pass RTM, resolving thin, closely spaced layers that amplitude-based migration blurs together, though its performance remains sensitive to the accuracy of the assumed source wavelet and background permittivity model (Gao, 2022). Collectively, these subspace and sparsity-driven techniques indicate that genuine sub-wavelength resolution in GPR is attainable, but only by moving beyond the amplitude- and traveltime-based imaging conditions that constrain traditional migration.



== Moving from Amplitude- to Phase-Based Interpretation

#draftnote[Placeholder]

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
