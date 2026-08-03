#import "../template.typ": *

= Theoretical Background <ch:theory>

This chapter derives the theory on which the rest of the thesis rests. It first
establishes the three migration algorithms used throughout (@sec:th-migration),
then develops the core result of this thesis --- that a sub-wavelength shift
between two migrated images becomes a recoverable linear phase ramp: the 2D
Fourier shift theorem (@sec:th-fourier-shift), its isolation via the
cross-spectrum (@sec:th-cross-spectrum), and the weighted least-squares fit
that recovers the displacement from it (@sec:th-wls). A supporting view
then builds spatial intuition for that result: a spatial-domain picture of
the lateral/vertical asymmetry (@sec:th-duality), confirmed independently in
simulation and extended with a complementary time-frequency perspective in
the Supplementary Material, §S5. The shared simulation and processing
pipeline that applies this theory to data is described separately in
@ch:methodology.

== GPR Migration Fundamentals <sec:th-migration>

This section derives the three migration algorithms used throughout this
thesis --- Kirchhoff delay-and-sum, Gazdag phase-shift, and
back-propagation --- each formulated under the shared zero-offset,
exploding-reflector convention introduced below. Although they differ
substantially in computational strategy and underlying approximations, all
three collapse the same raw hyperbolic B-scan into a focused image of the
true subsurface reflectivity.

All data in this thesis are zero-offset (collocated transmitter/receiver)
B-scans, which makes the _exploding-reflector model_ applicable
@claerbout1985. Every reflector in the subsurface is treated as an active source that radiates a pulse upward at $t=0$, which is recorded by
receivers at the surface. Because the real two-way travel time corresponds to
a wave travelling down and back up once, this fictitious one-way exploding
source must propagate at half the true medium velocity,
$ #vmig = v / 2 , $ <eq:vmig>
so that the one-way travel time in the exploding-reflector model exactly
matches the two-way travel time of the real survey. Every migration algorithm
used in this thesis operates under this convention.

=== Kirchhoff (Delay-and-Sum) Migration

Physically, Kirchhoff migration treats every subsurface point as a Huygens
secondary source, collapsing the diffraction hyperbola recorded from that
point back onto its apex, the point's true location @yao2012@ozdemir2014.
What separates true Kirchhoff migration from a plain hyperbolic
(diffraction) stack is that it is derived from the Kirchhoff integral
solution of the scalar wave equation, which additionally weights every
summed contribution by an obliquity factor $cos theta$ and a
spherical-spreading term $1 \/ sqrt(#vmig thin r)$ @ozdemir2014@smitha2016.
@schneider1978 formalises this as a linear inverse problem: Kirchhoff
migration is the mathematical adjoint $K^upright(T)$ of a forward
(demigration) operator $K$ that predicts recorded data $d(t,x)$ from a
reflectivity model $m(z,x)$ --- exactly the correlation-type "imaging
condition" that @jones2014 describes generally as building an image
wherever downgoing (source) and upcoming (receiver) wavefield contributions
coincide in space and time. This thesis implements Kirchhoff migration with
the PyLops Python library, whose `pylops.waveeqprocessing.Kirchhoff` class
follows this same _operator-based_ design: the object defines the forward
demigration operator once, and the delay-and-sum migrated image used
throughout this thesis is obtained simply by applying its adjoint (`.H`) to
the recorded data, so that the identical operator can later be reused,
unmodified, inside an iterative least-squares migration @pylops.

// Kirchhoff migration is formulated as a linear forward operator $K$ mapping a
// reflectivity model $m(z,x)$ to recorded data $d(t,x)$ by summing the
// reflectivity along the travel-time hyperbola of every trace,
// $
//   t(x, x_s, z) = r / #vmig + t_0 , quad r = sqrt(z^2 + (x - x_s)^2) ,
// $ <eq:kirchhoff-traveltime>
// where $x_s$ is the source/receiver position and $t_0$ a static time-zero
// correction. Migration is the adjoint operation: every recorded sample is
// smeared back (_delay-and-sum_) along its hyperbola of possible
// origins, and the image is the sum over all traces,
// $ m_"mig" = K^upright(T) d . $ <eq:kirchhoff-adjoint>
// This thesis uses the PyLops zero-offset Kirchhoff operator with an
// analytic-signal Ricker wavelet matched to the source, evaluated with a finite
// migration aperture @schneider1978.

=== Gazdag Phase-Shift Migration

Gazdag (phase-shift) migration is a one-way wave-equation
wavefield-extrapolation method that processes the *entire* recorded
wavefield at once, rather than summing along per-pixel travel-time curves
as Kirchhoff migration does. The B-scan is Fourier transformed to the
frequency-wavenumber domain $(omega, #kx)$ and downward-continued one
depth step $delta z$ at a time by repeatedly multiplying its spectrum by
the phase-shift term $e^(-j #kz delta z)$ derived from the factorised
one-way wave equation @schuster2017; at each step the exploding-reflector
imaging condition reads off the $delta t = 0$ component and adds it into
the image at that depth, the same downward-continuation-then-imaging-condition
procedure @jones2014 describes generally for wavefield-extrapolation
migration. This matches the GPR-specific phase-shift factor
$K = e^(j #kz delta z)$ applied step by step in @ozdemir2014@smitha2016.
Because the whole wavefield is extrapolated via FFT rather than by tracing
individual rays, Gazdag migration avoids Kirchhoff's high-frequency,
single-arrival approximation. It correctly handles multi-path energy, but
each step is restricted to a laterally invariant $v(z)$ velocity model,
whereas Kirchhoff's per-pixel hyperbolas tolerate arbitrary $v(x,z)$
@jones2014@schneider1978@ozdemir2014. Frequency-wavenumber pairs for which
the argument of #kz is negative are evanescent and numerically unstable to
continue, so they are excluded from the sum rather than left to grow
exponentially.

This thesis implements Gazdag migration with PyLops's
`pylops.waveeqprocessing.PhaseShift` operator: as with the Kirchhoff
operator above, it builds the depth-stepped continuation directly from the
phase-shift term for a given velocity and step size, and the migrated
image is obtained by applying it recursively and reading off the
$delta t = 0$ imaging condition at each depth @pylops.

// Gazdag migration works entirely in the frequency--wavenumber ($f$-$#kx$)
// domain @gazdag1978. The recorded wavefield is downward-continued one
// depth step $delta z$ at a time by multiplying its 2D temporal-frequency /
// horizontal-wavenumber spectrum by a phase-shift operator,
// $
//   U(z + delta z, #kx, omega) = U(z, #kx, omega) e^(j #kz delta z) , quad
//   #kz = sqrt((omega / #vmig)^2 - #kx^2) ,
// $ <eq:gazdag>
// and the image is built by applying the imaging condition --- extracting the
// $t=0$ component of the continued field --- at every depth step. Bins for
// which the argument of the square root in @eq:gazdag is negative
// correspond to evanescent energy and are set to zero before continuation;
// otherwise the unstable exponential growth of an imaginary $#kz$ produces
// migration "smile" artefacts. This thesis pads each B-scan with a $5\%$
// cosine taper and $100\%$ zero-padding in $x$ before transforming, to
// suppress wrap-around.

=== Back-Propagation (Time-Reversal) Migration

Unlike Kirchhoff and Gazdag migration, back-propagation (time-reversal)
migration uses no approximate operator at all: each recorded trace is
reversed in time, normalised, and re-injected as a source at its original
receiver position into a full electromagnetic simulation of the medium at
the migration velocity $#vmig$. By the time-reversal symmetry of the wave
equation, the back-propagated field refocuses at the true scatterer
location at the focusing time $t_"focus" = T - t_0$ ($T$ the trace length
and $t_0$ the static time-zero correction that shifts each trace so its
recorded onset coincides with the wave's true departure from the antenna),
and the migrated image is simply the field snapshot read off at that
instant. This thesis performs the re-injection and back-propagation
numerically with gprMax, the same open-source finite-difference
time-domain (FDTD) electromagnetic solver used for all forward modelling
in this thesis @gprmax --- using gprMax as the wave-equation solver that
carries out the extrapolation itself, the same role it plays in
@geng2022's gprMax-based reverse-time migration (RTM) of GPR data.

Geng and Ye’s RTM additionally forward-models the source-side wavefield with gprMax and cross-correlates it against the reverse-time-extrapolated receiver-side wavefield 
at every time step to build the image @geng2022, 
a genuinely two-wavefield imaging condition. While this full cross-correlation approach offers distinct advantages in noise suppression and phase accuracy—particularly when a high-quality wavelet is used for the forward extrapolation—the added computational expense is not strictly required here. Because every survey in this thesis is zero-offset and already treated under the exploding-reflector convention of @sec:th-migration
@claerbout1985, the recorded data effectively represent the source-side wavefield. Consequently, reading off a single back-propagated wavefield at
$t_"focus"$ provides sufficient structural fidelity for this analysis; this thesis therefore implements back-propagation rather than full cross-correlation RTM. Because it makes no high-frequency, single-arrival, or $v(z)$-only approximation beyond the exploding-reflector velocity halving itself, this method serves as a useful independent numerical check on the other two.

// As an independent, purely numerical cross-check of the two analytic methods
// above, every B-scan is also migrated by literal time-reversal: each trace is
// reversed in time, normalised, and re-injected as a source at its original
// receiver position into a finite-difference time-domain (gprMax) model of a
// homogeneous medium at $#vmig$. By the time-reversal symmetry of the wave
// equation, the back-propagated field refocuses at the true scatterer location
// at the focusing time
// $ t_"focus" = T - t_0 , $ <eq:backprop-focus>
// where $T$ is the trace length. The migrated image is read off as the field
// snapshot at $t_"focus"$, either as the full electric-field magnitude
// $||bold(E)||$ or as the single polarised component $E_z$. Unlike
// Kirchhoff and Gazdag migration, this method makes no high-frequency or
// zero-offset approximation beyond the exploding-reflector velocity halving
// itself, which makes it a useful independent check on the other two.

== The 2D Fourier Shift Theorem <sec:th-fourier-shift>

Let the baseline migrated image be $b(z,x)$. If a point scatterer
translates by a vertical distance $#Dz$ and a lateral distance $#Dx$ between
the baseline and monitor survey, the monitor image is, to the extent that
migration is linear. The two surveys are migrated with the same velocity
model, a perfectly translated copy of the baseline image,
$ m(z,x) = b(z - #Dz, space x - #Dx) . $ <eq:spatial-shift>
Taking the 2D continuous Fourier transform of the baseline image,
$
  B(#kz, #kx) = integral_(-oo)^(oo) integral_(-oo)^(oo)
    b(z,x) e^(-j(#kz z + #kx x)) dif z dif x ,
$ <eq:2d-fourier>
the Fourier shift theorem states that the spatial translation in
@eq:spatial-shift becomes a linear phase rotation of the spectrum,
$ M(#kz, #kx) = B(#kz, #kx) dot e^(-j(#kz #Dz + #kx #Dx)) . $ <eq:fourier-shift-theorem>
This is the same 2D shift theorem that underlies classical Fourier-domain
image-registration methods, where it is used to recover a rigid translation
between two images from the phase of their cross-power spectrum
@decastro1987.
The amplitude spectrum is unchanged, $|M|=|B|$, but every frequency
coordinate $(#kz, #kx)$ acquires a phase shift proportional to that
coordinate. This is the central fact the rest of the chapter exploits: a
displacement that is invisible in the migrated amplitude image
(@sec:meth-resolution) is, in principle, exactly recoverable from the phase of
the spectrum.

== Isolating the Phase Plane: The Cross-Spectrum <sec:th-cross-spectrum>

The absolute phase of a single migrated image is not directly useful: it
depends on the shape of the source wavelet and on the (arbitrary) position of
the target within the image, and is, in general, a chaotic, wrapped function
of $(#kz, #kx)$. The displacement information in @eq:fourier-shift-theorem
is obtained by forming the complex _cross-spectrum_ between the baseline
and monitor spectra @decastro1987,
$ #XS (#kz, #kx) = B(#kz, #kx) dot M^*(#kz, #kx) . $ <eq:cross-spectrum-def>
Recovering a relative shift or delay from the phase of a cross-spectrum
between two otherwise-similar signals is an established technique used in for example cross-spectral analysis of ambient seismic noise @clarke2011, and the cross-spectrum
itself is a standard tool for isolating a correlated signal shared between
two measurements in precision metrology @nelson2014.

Substituting @eq:fourier-shift-theorem,
$
  #XS (#kz, #kx)
  = B(#kz, #kx) dot [B(#kz, #kx) e^(-j(#kz #Dz + #kx #Dx))]^*
  = |B(#kz, #kx)|^2 e^(j(#kz #Dz + #kx #Dx)) .
$ <eq:cross-spectrum-expand>
Because $|B|^2$ is real and positive, extracting the phase of $#XS$ via
$Phi = angle #XS$ makes it vanish completely, leaving
$
  Phi(#kz, #kx) = #kz #Dz + #kx #Dx ,
$ <eq:phase-plane>
a perfectly flat plane through the origin: its slope along $#kz$ is exactly
the vertical displacement $#Dz$, and its slope along $#kx$ is exactly the
lateral displacement $#Dx$.

#para-head[Units of the cross-phase spectrum.] The value of $Phi$ is an
angle (rad or deg); its domain, the wavenumbers $(#kz, #kx)$, has units of
$"rad" dot "m"^(-1)$. The slope of the plane in @eq:phase-plane,
$partial Phi \/ partial #kz$, therefore has units of
$"rad" \/ ("rad" dot "m"^(-1)) = "m"$: the radians cancel, and
what remains is a physical distance. This is the precise sense in which an
angular rotation in the frequency domain converts directly into a
millimetre-scale displacement.

== Weighted Least-Squares Plane Fitting <sec:th-wls>

A discretised GPR image yields not one but $N$ frequency-bin observations of
@eq:phase-plane inside a chosen passband (@sec:th-mask-weight),
each of the form
$ Phi_i = k_(z,i) thin #Dz + k_(x,i) thin #Dx + c , quad i = 1 , dots , N , $ <eq:phase-plane-discrete>
where $c$ is a small constant absorbing calibration bias. Collecting
all $N$ equations into the overdetermined linear system $A bold(u) = bold(Phi)$,
$
  underbrace(
    mat(delim: "[",
      k_(z,1), k_(x,1), 1;
      k_(z,2), k_(x,2), 1;
      dots.v, dots.v, dots.v;
      k_(z,N), k_(x,N), 1),
    A
  )
  underbrace(
    mat(delim: "[", #Dz, #Dx, c)^upright(T),
    bold(u)
  )
  =
  underbrace(
    mat(delim: "[", Phi_1; Phi_2; dots.v; Phi_N),
    bold(Phi)
  ) .
$ <eq:design-matrix>

=== Why ordinary least squares fails

Ordinary least squares (OLS) minimises
$S(bold(u)) = (bold(Phi) - A bold(u))^upright(T) (bold(Phi) - A bold(u))$,
giving the closed-form solution
$bold(u) = (A^upright(T) A)^(-1) A^upright(T) bold(Phi)$.
This treats every frequency bin as equally reliable, but in a GPR spectrum a
bin at the antenna's peak power is far more reliable than one at the edge of
the band, which is dominated by background noise. OLS gives both bins equal
weight, letting noisy bins corrupt the plane fit.

=== The weighted solution <sec:th-mask-weight>

Two safeguards make the fit robust enough for sub-wavelength accuracy:

#linebreak()

#para-head[A. The band-pass mask.] The system in @eq:design-matrix is
restricted to bins inside the coherent envelope of the source wavelet,
$|#kz|, |#kx| < 1.4 k_(z c)$, where $k_(z c)$ is the dominant vertical
wavenumber of the pulse. This keeps the total phase rotation $Phi_i$ inside
$plus.minus pi$, preventing the fit from wrapping. The threshold $1.4$ is an
empirical choice, tuned by trial and error to balance the two failure modes
of the mask: too large, and bins far from the spectral peak accumulate
enough phase to wrap around $plus.minus pi$; too small, and bins that still
carry useful, coherent phase information are excluded from the fit.

#linebreak()

#para-head[B. Amplitude weighting.] A diagonal weight matrix $W$ is built from
the cross-spectrum magnitude, $W_(i i) = |italic("XS")_i|$, so that high-energy bins
dominate the fit and noise-floor bins are suppressed, the same
reliability-weighting principle used to stabilise weighted least-squares
phase-to-displacement inversion in differential InSAR time-series processing
@falabella2020 and coherence-weighted phase-delay regression in ambient-noise
seismology @clarke2011. The cost function becomes the energy-weighted residual,
$
  S_W (bold(u)) = sum_(i=1)^(N) |italic("XS")_i|^2
    (Phi_i - (k_(z,i) thin #Dz + k_(x,i) thin #Dx + c))^2 ,
$ <eq:wls-cost>
which is minimised by
$ bold(u) = (A^upright(T) W^2 A)^(-1) A^upright(T) W^2 bold(Phi) . $ <eq:wls-solution>
In practice, computing $(A^upright(T) W^2 A)^(-1)$ directly squares the
condition number of the system; the implementation instead pre-multiplies
both sides by $W$,
$ (W A) bold(u) = (W bold(Phi)) , $ <eq:wls-preweighted>
and solves the pre-weighted system via
singular value decomposition, which is equivalent to @eq:wls-solution
but numerically far more stable.

Beyond absorbing calibration bias, the constant column of
@eq:design-matrix also lets the same fit isolate a frequency-independent
phase offset caused by a _material_ change at the target (e.g. a
sub-wavelength fracture filling with fluid) in the intercept $c$, cleanly
separated from the geometric shift $(#Dz, #Dx)$ --- the derivation of this
decoupling is given in the Supplementary Material, §S6.

@fig:phaseplane-schematic summarises the full pipeline derived above on a
moving PEC point scatterer example --- Gazdag-migrated Lateral movement data at the
smallest scale tested, $1\/32 lambda$ (@ch:hyp1): the baseline image is the
B-scan of the scatterer at its initial position and the monitor image is
the B-scan after it has been laterally translated by $1\/32 lambda$; both
are Fourier transformed, their
cross-spectrum isolates a linear phase ramp via the shift theorem of
@sec:th-fourier-shift, and the weighted least-squares fit of @sec:th-wls
recovers the sub-wavelength displacement from that ramp's slope. For this
example, the phase-plane fit recovered the true lateral displacement with
no error.

#figure(
  img("H1_038_Phase-Plane_Workflow_--_Gazdag_Lateral_132lambda.png"),
  caption: [Phase-plane shift-estimation pipeline on a synthetic example
    (Gazdag migration, Lateral movement, $1\/32 lambda$, @ch:hyp1). Top row:
    (a) migrated baseline image; (b) migrated monitor image; (c) their
    difference. Bottom row: (d) cross-spectrum phase, opaque inside and
    transparent outside the fitting mask ($|#XS| > 10%$ of its peak and
    $|#kz|, |#kx| < 1.4 k_(z c)$, @sec:th-mask-weight); (e) cross-spectrum
    energy $|#XS|$, with the same mask outlined; (f) the weighted
    least-squares plane fit along $#kx$ (@sec:th-wls), points coloured by
    fit weight $|#XS|$.],
) <fig:phaseplane-schematic>

== Space-Wavenumber Duality: Lateral versus Vertical Asymmetry <sec:th-duality>

Having derived the phase plane formally (@sec:th-fourier-shift through
@sec:th-wls), it is worth building the spatial intuition behind it: _why_ a
sub-wavelength displacement leaves a detectable signature in the phase of a
migrated image, and why the lateral and vertical directions do not behave
alike. This section examines how lateral and vertical displacements appear in
the _spatial_ phase of a migrated scatterer, revealing an asymmetry between the
two directions that the Fourier treatment above inherits.

=== Lateral direction: position as a proxy for wavenumber

Near the apex of a migrated point scatterer, the lateral instantaneous phase
$phi(x)$ is well approximated by a parabola,
$phi(x) approx 1/2 C x^2 + phi_0$,
where $C$ is the spatial curvature of the focused pulse. Its spatial
derivative, the instantaneous lateral wavenumber, is then
$ k_(x,"inst")(x) = (partial phi) / (partial x) approx C thin x , $ <eq:kx-inst>
i.e. _linearly proportional to lateral position_: the migrated image
behaves like a prism, naturally separating horizontal spatial frequencies
across its width. If the scatterer shifts laterally by $#Dx$, a first-order
Taylor expansion of $phi_2(x) = phi_1(x - #Dx)$ gives a phase-difference image
$
  Delta phi(x) = phi_2(x) - phi_1(x) approx -k_(x,"inst")(x) thin #Dx
    = -C thin x thin #Dx ,
$ <eq:dphi-lateral>
i.e. a straight _slope_ in $x$, whose steepness is proportional to $#Dx$.

=== Vertical direction: a constant plateau, not a slope

The vertical trace through a migrated scatterer is, by contrast, a
compressed source pulse oscillating at the dominant vertical carrier
wavenumber $k_(z c) = 4 pi f_c \/ v$. Near the peak, its phase is linear in depth,
$phi(z) approx k_(z c)(z - z_0) + phi_0$, so the instantaneous vertical
wavenumber is
$
  k_(z,"inst")(z) = (partial phi) / (partial z) approx k_(z c)
  = "constant" ,
$ <eq:kz-inst>
independent of $z$: unlike the lateral case, depth position is _not_ a
proxy for vertical wavenumber, because the vertical trace does not spatially
separate its frequency content the way the migrated lateral profile does.
A vertical displacement $#Dz$ then produces
$ Delta phi(z) approx -k_(z c) thin #Dz , $ <eq:dphi-vertical>
which contains no remaining $z$ dependence: $partial [Delta phi(z)] \/ partial z approx 0$.
A vertical shift therefore appears in the spatial domain as a flat
_plateau_ whose constant level is proportional to $#Dz$, not as a slope.

=== Reconciling the two pictures

Combining @eq:dphi-lateral and @eq:dphi-vertical, the full 2D spatial phase
difference near the apex $(z_0, 0)$ of a scatterer displaced by $(#Dz, #Dx)$ is
$
  Delta phi(z, x) approx -(k_(z c) thin #Dz + (C thin x) thin #Dx) ,
$ <eq:dphi-2d>
which has the same structure as the phase plane derived formally in
@sec:th-fourier-shift above, with local approximations
$k_(z,"inst") -> #kz$ and $C thin x -> #kx$. The formal Fourier derivation
is preferred as the primary quantitative tool for three reasons: (i) the
parabolic approximation in @eq:kx-inst only holds near the apex of a difference migrated B-scan, whereas the
Fourier plane is exactly flat everywhere inside the wavenumber passband; (ii) amplitude
weighting by $|#XS|$ has no clean spatial-domain analogue; and (iii) in the
spatial picture a calibration-bias intercept mixes irrecoverably with the
vertical term $-k_(z c) #Dz$ in @eq:dphi-2d, whereas the Fourier-domain fit
keeps them exactly orthogonal (@sec:th-wls).

#supp-note[This asymmetry is independently confirmed in simulation by a local instantaneous-phase diagnostic applied to the same lateral and vertical-displacement datasets. It is further developed with a
complementary time-frequency perspective, in the Supplementary Material,
§S5.]
