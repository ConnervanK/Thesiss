#import "../template.typ": *

= Theoretical Background and Methodology <ch:theory>

This chapter covers both the theoretical foundations and the simulation
and processing methodology used throughout the thesis. The first part
(@sec:th-migration through @sec:th-local) derives the three migration
algorithms and the 2D phase-plane shift-estimation method in full. The second
part (@ch:methodology) describes the shared forward-modelling and processing
pipeline and closes with a validation of the amplitude-based resolution floor
(@sec:meth-resolution) that motivates the phase-based approach.

== GPR Migration Fundamentals <sec:th-migration>

All data in this thesis are zero-offset (collocated transmitter/receiver)
B-scans, which makes the _exploding-reflector model_ applicable
@claerbout1985: every reflector in the subsurface is treated as if it
were an active source that radiates a pulse upward at $t=0$, recorded by
receivers at the surface. Because the real two-way travel time corresponds to
a wave travelling down and back up once, this fictitious one-way exploding
source must propagate at half the true medium velocity,
$ #vmig = v / 2 , $ <eq:vmig>
so that the one-way travel time in the exploding-reflector model exactly
matches the two-way travel time of the real survey. Every migration algorithm
used in this thesis operates under this convention.

=== Kirchhoff (Delay-and-Sum) Migration

Kirchhoff migration is formulated as a linear forward operator $K$ mapping a
reflectivity model $m(z,x)$ to recorded data $d(t,x)$ by summing the
reflectivity along the travel-time hyperbola of every trace,
$
  t(x, x_s, z) = r / #vmig + t_0 , quad r = sqrt(z^2 + (x - x_s)^2) ,
$ <eq:kirchhoff-traveltime>
where $x_s$ is the source/receiver position and $t_0$ a static time-zero
correction. Migration is the adjoint operation: every recorded sample is
smeared back (_delay-and-sum_) along its hyperbola of possible
origins, and the image is the sum over all traces,
$ m_"mig" = K^upright(T) d . $ <eq:kirchhoff-adjoint>
This thesis uses the PyLops zero-offset Kirchhoff operator with an
analytic-signal Ricker wavelet matched to the source, evaluated with a finite
migration aperture @schneider1978.

=== Gazdag Phase-Shift Migration

Gazdag migration works entirely in the frequency--wavenumber ($f$-$#kx$)
domain @gazdag1978. The recorded wavefield is downward-continued one
depth step $delta z$ at a time by multiplying its 2D temporal-frequency /
horizontal-wavenumber spectrum by a phase-shift operator,
$
  U(z + delta z, #kx, omega) = U(z, #kx, omega) e^(j #kz delta z) , quad
  #kz = sqrt((omega / #vmig)^2 - #kx^2) ,
$ <eq:gazdag>
and the image is built by applying the imaging condition --- extracting the
$t=0$ component of the continued field --- at every depth step. Bins for
which the argument of the square root in @eq:gazdag is negative
correspond to evanescent energy and are set to zero before continuation;
otherwise the unstable exponential growth of an imaginary $#kz$ produces
migration "smile" artefacts. This thesis pads each B-scan with a $5\%$
cosine taper and $100\%$ zero-padding in $x$ before transforming, to
suppress wrap-around.

=== Back-Propagation (Time-Reversal) Migration

As an independent, purely numerical cross-check of the two analytic methods
above, every B-scan is also migrated by literal time-reversal: each trace is
reversed in time, normalised, and re-injected as a source at its original
receiver position into a finite-difference time-domain (gprMax) model of a
homogeneous medium at $#vmig$. By the time-reversal symmetry of the wave
equation, the back-propagated field refocuses at the true scatterer location
at the focusing time
$ t_"focus" = T - t_0 , $ <eq:backprop-focus>
where $T$ is the trace length. The migrated image is read off as the field
snapshot at $t_"focus"$, either as the full electric-field magnitude
$||bold(E)||$ or as the single polarised component $E_z$. Unlike
Kirchhoff and Gazdag migration, this method makes no high-frequency or
zero-offset approximation beyond the exploding-reflector velocity halving
itself, which makes it a useful independent check on the other two.

== Space-Wavenumber Duality: Lateral versus Vertical Asymmetry <sec:th-duality>

Before deriving the formal Fourier-domain result, it is useful to build
spatial intuition for _why_ a sub-wavelength displacement leaves a detectable
signature in the phase of a migrated image. This section does so by examining
how lateral and vertical displacements appear in the _spatial_ phase of a
migrated scatterer, and reveals an important asymmetry between the two
directions that persists into the full Fourier treatment.

=== Lateral direction: position as a proxy for wavenumber

Near the apex of a migrated point scatterer, the lateral instantaneous phase
$phi(x)$ is well approximated by a parabola,
$phi(x) approx 1/2 C x^2 + phi_0$,
where $C$ is the spatial curvature of the focused pulse. Its spatial
derivative, the instantaneous lateral wavenumber, is then
$ k_(x,"inst")(x) = partial phi / partial x approx C thin x , $ <eq:kx-inst>
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
  k_(z,"inst")(z) = partial phi / partial z approx k_(z c)
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
@sec:th-fourier-shift below, with local approximations
$k_(z,"inst") -> #kz$ and $C thin x -> #kx$. The formal Fourier derivation
is preferred as the primary quantitative tool for three reasons: (i) the
parabolic approximation in @eq:kx-inst only holds near the apex, whereas the
Fourier plane is exactly flat everywhere inside the passband; (ii) amplitude
weighting by $|#XS|$ has no clean spatial-domain analogue; and (iii) in the
spatial picture a calibration-bias intercept mixes irrecoverably with the
vertical term $-k_(z c) #Dz$ in @eq:dphi-2d, whereas the Fourier-domain fit
keeps them exactly orthogonal (@sec:th-wls).

== The 2D Fourier Shift Theorem <sec:th-fourier-shift>

Let the baseline migrated image be $b(z,x)$. If a point scatterer
translates by a vertical distance $#Dz$ and a lateral distance $#Dx$ between
the baseline and monitor survey, the monitor image is, to the extent that
migration is linear and the two surveys are migrated with the same velocity
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
is isolated by forming the complex _cross-spectrum_ between the baseline
and monitor spectra,
$ #XS (#kz, #kx) = B(#kz, #kx) dot M^*(#kz, #kx) . $ <eq:cross-spectrum-def>
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
    mat(delim: "[", #Dz; #Dx; c),
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
bin at the antenna's peak power is far more reliable than a bin at the edge of
the band, which is dominated by background noise. OLS gives both bins equal
weight, letting noisy bins corrupt the plane fit.

=== The weighted solution <sec:th-mask-weight>

Two safeguards make the fit robust enough for sub-wavelength accuracy.

#para-head[A. The band-pass mask.] The system in @eq:design-matrix is
restricted to bins inside the coherent envelope of the source wavelet,
$|#kz|, |#kx| < 1.4 k_(z c)$, where $k_(z c)$ is the dominant vertical
wavenumber of the pulse. This keeps the total phase rotation $Phi_i$ inside
$plus.minus pi$, preventing the fit from wrapping.

#para-head[B. Amplitude weighting.] A diagonal weight matrix $W$ is built from
the cross-spectrum magnitude, $W_(i i) = |italic("XS")_i|$, so that high-energy bins
dominate the fit and noise-floor bins are suppressed. The cost function
becomes the energy-weighted residual,
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
and solves the pre-weighted system with `numpy.linalg.lstsq` via
singular value decomposition, which is equivalent to @eq:wls-solution
but numerically far more stable.

== Decoupling Geometric Movement from Material Change <sec:th-material-change>

The third column of the design matrix in @eq:design-matrix is the
constant vector $bold(1)$, deliberately included alongside the two
wavenumber columns. This section shows, for theoretical completeness, why that
design choice lets the same fit simultaneously measure displacement _and_
detect a change in the dielectric material at the target, and why the two
never contaminate each other.

=== Why a sub-wavelength fracture produces a frequency-independent phase shift

Consider a fracture of thickness $d$ much smaller than the wavelength,
$d lt.double lambda$, embedded in a host matrix of permittivity $epsilon.alt_m$ and filled
with a material of permittivity $epsilon.alt_f$. The wave reflects once at the top
and once at the bottom of the fracture; for $d lt.double lambda$ these two
reflections overlap into a single composite wavelet, and the total
reflection coefficient can be approximated as
$ R_"total" approx R_"top" + R_"bottom" e^(-j 2 k d) . $ <eq:thinlayer-1>
A first-order expansion of the exponential for $k d lt.double 1$ gives
$
  R_"total" approx -j omega thin (d sqrt(epsilon.alt_m)) / (2 c_0)
    ((epsilon.alt_m - epsilon.alt_f) / epsilon.alt_f) ,
$ <eq:thinlayer-2>
where $c_0$ is the speed of light in vacuum. The leading $-j omega$ factor
means a thin layer reflects the time-derivative of the incident pulse.

When the fracture fluid changes between the baseline and monitor survey
(e.g. air, $epsilon.alt_"air" = 1$, displaced by water, $epsilon.alt_"water" = 81$)
while the thickness $d$ stays fixed, the ratio of the two reflection
coefficients is
$
  Delta R = R_"mon" / R_"base"
  = (epsilon.alt_"air" (epsilon.alt_m - epsilon.alt_"water")) / (epsilon.alt_"water" (epsilon.alt_m - epsilon.alt_"air")) .
$ <eq:thinlayer-ratio>
Both the geometric factor $d$ and the derivative factor $-j omega$ cancel
exactly in this ratio, leaving a complex constant, independent of $omega$
and therefore independent of $#kz$ and $#kx$. In the phase domain this is a
uniform rotation
$ Delta R approx alpha e^(-j #Dtheta) , $ <eq:thinlayer-phase>
applied equally to every frequency in the pulse bandwidth: a flat,
frequency-independent phase offset $#Dtheta$, with an amplitude attenuation
$alpha$.

=== Why the WLS fit puts material change exactly into the intercept

For a fluid-substitution event with no mechanical movement
($#Dz = #Dx = 0$), @eq:thinlayer-phase means every observation in
@eq:phase-plane-discrete is the same constant, $Phi_i = #Dtheta$ for
all $i$. Looking at the three columns of $A$ in @eq:design-matrix:

- the $#kz$ column spans negative to positive wavenumbers --- any
  non-zero $#Dz$ tilts the predicted plane, which can only _increase_
  the residual against a perfectly flat target, so the optimum is $#Dz = 0$;
- the $#kx$ column behaves identically, forcing $#Dx = 0$;
- the constant column is exactly $bold(1)$, so setting $c = #Dtheta$
  matches the flat target with zero residual.

This is possible because the wavenumber columns and the constant column are
linearly independent: the matrix inversion in @eq:wls-solution decouples them
exactly. If a target both moves _and_ changes material at once, the fitted
plane both tilts (giving $#Dz, #Dx$) and shifts vertically (giving
$c = #Dtheta$), and the two effects remain perfectly separable in the same
single fit.

== Time-Frequency Perspective: Local Phase and Spectral-Line Analysis <sec:th-local>

The phase-plane fit treats the migrated image as a whole. @sec:hyp1-phaseplane
(@sec:tlp-spectral-line, @sec:tlp-spectrogram, @sec:tlp-stft) instead
analyses individual unmigrated or migrated _traces_ with a localised
time-frequency transform, which gives access to _when_ (in two-way
time) a phase change occurs, complementing the spatial picture above.

=== The localised Fourier shift theorem

Let $s_1(t)$ be a baseline trace. A sliding-window transform (Gaussian or
Hanning window $w$) gives a localised spectrum
$
  S_1(tau, omega) = integral_(-oo)^(oo) s_1(t) thin w(t - tau) thin e^(-j omega t) thin dif t ,
$ <eq:stft-def>
at window centre $tau$ and angular frequency $omega = 2 pi f$. If the monitor
trace is a delayed, phase-rotated copy, $s_2(t) = s_1(t - #Dt) e^(j #Dtheta)$,
and the delay $#Dt$ is small compared with the window width, the window
itself barely shifts and the _localised Fourier shift theorem_ applies,
$ S_2(tau, omega) approx S_1(tau, omega) e^(-j omega #Dt) e^(j #Dtheta) . $ <eq:local-shift>
Forming the local cross-spectrum exactly as in @eq:cross-spectrum-def
and taking its angle, the baseline phase and amplitude cancel, leaving the
local analogue of @eq:phase-plane,
$ Delta Phi(tau, f) approx -2 pi f thin #Dt + #Dtheta . $ <eq:local-phase-line>
At the two-way time $tau_0$ of the target reflection, a straight-line fit of
$Delta Phi$ against $f$ has slope $-2 pi #Dt$ (the mechanical shift) and
intercept $#Dtheta$ (a calibration or material-change offset) --- the exact
time-domain counterpart of the slope/intercept decomposition in @sec:th-wls.
Converting between the two-way-time and depth pictures uses the standard relation
$ #Dt = (2 #Dz) / v , $ <eq:dt-dz>
so that $omega #Dt equiv #kz #Dz$ with $#kz = 2 omega \/ v$: the temporal-frequency
slope and the vertical-wavenumber slope are the same physical quantity, viewed
in two different but exactly equivalent coordinate systems.

=== Three diagnostic views used in @sec:hyp1-phaseplane

Three derived plots make @eq:local-phase-line directly visible in the
data, and are used repeatedly in the figures of @sec:hyp1-phaseplane:

/ Spectral line, $Delta Phi(f)$ at fixed $tau_0$: a straight line
  through the origin for pure mechanical movement; a flat line offset from
  zero for a pure phase-rotation offset; a sloped line with non-zero intercept
  for a combination of the two.

/ Cross-phase spectrogram, $Delta Phi(tau, f)$: outside the target
  reflection this is incoherent, salt-and-pepper phase noise; at the
  target's two-way time a coherent window appears, showing a vertical
  fringe pattern for movement or a uniform colour block for a pure phase
  offset.

/ Polar vector rotation: the complex coefficient at the dominant
  frequency and peak two-way time, plotted as a vector in the complex
  plane for baseline and monitor; a phase offset rotates this vector with
  negligible length change for a purely geometric shift.

These local, trace-based views and the global 2D wavenumber fit of
@sec:th-wls are not competing techniques: they are Fourier duals of the
same underlying physics, related by a spatial Fourier transform of the
time-frequency decomposition with respect to the lateral coordinate $x$,
which reduces (under a smooth-window approximation) directly to the global
spectrum $U(omega, #kx)$ used throughout this chapter. The wavenumber-domain
fit remains the primary quantitative tool because it linearises the spatial
curvature of @eq:kx-inst exactly and admits amplitude weighting
natively (@sec:th-mask-weight); the local, time-frequency view is used
in @sec:hyp1-phaseplane to localise _where_ (in $tau$) a phase
anomaly originates, which the global fit alone cannot show.

// ─────────────────────────────────────────────────────────────────────────────
//  METHODOLOGY
// ─────────────────────────────────────────────────────────────────────────────

== Methodology <ch:methodology>

This section describes the simulation and processing pipeline shared by the
experiments in @ch:hyp1 and @ch:hyp2. Each experiment varies the scatterer
configuration and the displacement under test, but reuses the same forward
model, signal-conditioning steps, migration implementations, and (from @ch:hyp1
onward) the same phase-plane estimator. The section closes with a short
validation experiment (@sec:meth-resolution) that establishes the
amplitude-based resolution floor of the migration algorithms before they are
used to test any of the hypotheses.

=== Forward Modelling with gprMax

All B-scans are simulated with the open-source finite-difference time-domain
solver gprMax @gprmax. The source is a Ricker wavelet with centre frequency
$f_c = 15 "GHz"$ and time-zero offset $t_0 = 0.0943 "ns"$ (@fig:res-setup
(a)), chosen so that its usable bandwidth defines the dominant wavelength
$lambda$ used to express every displacement scale in this thesis ($2 lambda$
down to $1 \/ 32 lambda$). The computational domain is discretised on a
uniform $1 "mm"$ grid with perfectly-matched-layer (PML) absorbing
boundaries.

#draftnote[the figure titles encode the domain extent as "4010 m"; the same
auto-titling code elsewhere strips decimal points from floats (e.g. a depth of
0.676 m appears as `0676_m`, and a position of 2.0 m appears as `20_m`), so
this almost certainly reads as a domain of ≈4.01 m rather than 4010 m ---
confirm against the notebook before quoting a final value.]

A zero-offset (collocated transmitter and receiver) survey is simulated by
sweeping a single transmitter--receiver pair across the surface.

=== Scatterer and Medium Models

One target geometry is used across the synthetic experiments:

- *Point scatterers* (@sec:meth-resolution, @ch:hyp1): perfect-electric-conductor
  (PEC) cylinders of radius $r = 28 "mm"$, buried at a depth of $0.676 "m"$
  in ice. @sec:meth-resolution places two such cylinders at a swept
  separation; @ch:hyp1 instead holds one cylinder fixed as a baseline and
  displaces a second, in the lateral, vertical, or diagonal direction, by
  the same family of sub-wavelength steps.

=== Signal Conditioning Pipeline <sec:meth-conditioning>

Every raw B-scan is processed identically before migration:

+ *Background subtraction.* A background-only simulation (no scatterer) is
  subtracted trace-by-trace to suppress the direct air/ground wave and
  isolate the scatterer reflection (e.g. @fig:res-bscans).

+ *Tapering and $t_0$ alignment.* An exponential decay taper suppresses
  late-arriving energy, a cosine end-taper removes hyperbola tails at the
  edge of the migration aperture, and a static shift aligns the surface
  reflection to $t = 0$ (e.g. @fig:res-taper).

+ *Noise injection (where stated).* @ch:hyp2 contaminates the conditioned
  B-scan with synthetic Laplace-distributed noise at $10%$ of the signal
  standard deviation, fitted from real field data.

=== Migration Algorithms Implemented

Every conditioned B-scan in @sec:meth-resolution, @ch:hyp1, and @ch:hyp2
is migrated with all three algorithms derived in @sec:th-migration:
Kirchhoff delay-and-sum (PyLops zero-offset operator), Gazdag $f$-$k$
phase-shift migration, and gprMax-based time-reversal back-propagation. All
three share the implementation in `helper_functions/migration.py`
(`PylopsKirchoffMigration`, `gazdag_migration`, and `write_backprop_files`)
and the exploding-reflector convention $#vmig = v \/ 2$ of @eq:vmig, so that
the same velocity model and the same migration aperture are used for a
baseline/monitor pair, which is required for the displacement estimate of
@sec:meth-phaseplane to be valid.

=== The 2D Phase-Plane Shift-Estimation Pipeline <sec:meth-phaseplane>

The theory of @sec:th-fourier-shift and @sec:th-wls is applied to the
migrated images produced above via the function `estimate_shift_2d` defined in
`TimeLapse_Processing.ipynb`. The end-to-end workflow is:

+ *Migrate* the baseline and monitor B-scan with an identical velocity model
  and aperture (@sec:th-migration).

+ *Crop a region of interest (ROI)* tightly around the target (in this
  thesis, a window of $plus.minus 2.5 lambda$ about the scatterer apex,
  located on the Hilbert-envelope peak of the baseline image), so that
  static background structure elsewhere in the image cannot drag the fitted
  plane towards zero.

+ *Taper, then take the 2D FFT* of both cropped images and form the
  cross-spectrum $#XS$ of @eq:cross-spectrum-def.

+ *Mask and weight*: restrict the fit to $|#kz|, |#kx| < 1.4 k_(z c)$ and
  weight every bin by $|#XS|$ (@sec:th-mask-weight).

+ *Solve* the weighted least-squares system of @eq:wls-preweighted for
  $(#Dz, #Dx, c)$.

#para-head[Practical considerations carried over from field application.]
Two points are worth flagging for any future extension of this pipeline to
real field surveys (@ch:discussion), even though the synthetic data in this
thesis sidesteps them by construction: the ROI window should be roughly
3--4 times the target pulse width (too small clips the signal under tapering,
too large admits background noise into the fit); and a genuine change in
soil/ice moisture between surveys changes the true wave velocity $v$, which
the estimator would otherwise misinterpret as a spurious vertical shift $#Dz$
unless the velocity is recalibrated against a known static reflector first.

=== Validation: Amplitude Resolution Floor of the Migration Algorithms <sec:meth-resolution>

Before any of the hypotheses can be tested with confidence, the migration
pipeline above is validated on a simple, well-understood baseline problem: how
closely can two _stationary_ point scatterers be spaced before the migration
algorithms can no longer tell them apart? This validation experiment uses two
PEC cylinders illuminated by a zero-offset GPR B-scan, migrated with all three
algorithms. Its result is the amplitude-based resolution floor against which
every displacement-detection result in @ch:hyp1 and @ch:hyp2 is later compared.

@fig:res-setup shows the forward-model setup: the Ricker source wavelet, the
gprMax domain and grid, and the swept separation between two PEC cylinder
scatterers ($r = 28 "mm"$, depth $0.676 "m"$), from $2 lambda$ down to
$1 \/ 16 lambda$.

#figure(
  grid(columns: (1fr, 1fr, 1fr), gutter: 0.8em,
    img("RES_001_Ricker_Wavelet_f_c__15_GHz_t0__0943_ns.png"),
    img("RES_002_Resolution_Study__Model_Geometry__domain_4010_m_Δx__1_mm_PML.png"),
    img("RES_003_Scatterer_Positions__PEC_Cylinders__r__28_mm_depth__0676_m.png"),
  ),
  caption: [Forward-model setup for the resolution validation: (a) the Ricker
    source wavelet ($f_c = 15 "GHz"$, $t_0 = 0.0943 "ns"$); (b) the gprMax
    domain and grid; (c) the swept separation between the two PEC cylinder
    scatterers, $r = 28 "mm"$, depth $0.676 "m"$.],
) <fig:res-setup>

@fig:res-bscans shows the simulated zero-offset B-scans before and after
background subtraction, and @fig:res-taper the effect of the standard
tapering and $t_0$-shift conditioning of @sec:meth-conditioning on the
$2 lambda$ separation scenario.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("RES_004_GPR_B-Scans__Background_and_Separation_Models.png"),
    img("RES_005_GPR_B-Scans__Background_Subtracted.png"),
  ),
  caption: [Raw zero-offset B-scans for the resolution validation, before (a)
    and after (b) background subtraction.],
) <fig:res-bscans>

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("RES_006_Effect_of_Tapering_and_t0_Shift__2λ_dataset_single_trace.png"),
    img("RES_007_B-scan_effect_of_tapering_and_t0_shift__2λ_dataset.png"),
  ),
  caption: [Effect of tapering and the $t_0$ shift on the $2 lambda$ dataset:
    (a) a single representative trace; (b) the complete B-scan.],
) <fig:res-taper>

@fig:res-kirchhoff, @fig:res-gazdag, and @fig:res-backprop show the migrated
image for every separation scenario, for Kirchhoff, Gazdag, and
back-propagation migration respectively, each with a zoomed view around the
true scatterer depth.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("RES_008_Kirchhoff_Migration__All_Datasets____f_c15_GHz____aperture40.png"),
    img("RES_009_Kirchhoff_Migration_zoomed____f_c15_GHz____aperture40.png"),
  ),
  caption: [Kirchhoff migration of the resolution validation ($f_c = 15 "GHz"$,
    aperture $= 40$ traces): (a) all separation scenarios; (b) zoomed view.],
) <fig:res-kirchhoff>

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("RES_010_Gazdag_Phase-Shift_Migration__All_Datasets____f_c15_GHz____z.png"),
    img("RES_011_Gazdag_Phase-Shift_Migration_zoomed____f_c15_GHz.png"),
  ),
  caption: [Gazdag phase-shift migration of the resolution validation
    ($f_c = 15 "GHz"$): (a) all separation scenarios; (b) zoomed view.],
) <fig:res-gazdag>

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("RES_012_Back-Propagation_E__All_Datasets____focus_at_1906_ns.png"),
    img("RES_013_Back-Propagation_E_zoomed____focus_at_1906_ns.png"),
    img("RES_014_Back-Propagation_Ez__All_Datasets____focus_at_1906_ns.png"),
    img("RES_015_Back-Propagation_Ez_zoomed____focus_at_1906_ns.png"),
  ),
  caption: [Time-reversal back-propagation migration of the resolution
    validation, focused at $t = 1906 "ns"$: field-magnitude image (a, b) and
    the $E_z$ component (c, d), each with a zoomed view around the scatterer
    depth.],
) <fig:res-backprop>

@fig:res-psf (a) overlays the signed migrated amplitude from all three
algorithms at $f_c = 15 "GHz"$, and @fig:res-psf (b) plots the normalised
lateral point-spread function (PSF) extracted at the true scatterer depth as
a function of separation.

#figure(
  grid(columns: (1fr, 1fr), gutter: 0.8em,
    img("RES_016_Migration_Comparison__Signed_Amplitude____f_c15_GHz____apert.png"),
    img("RES_017_Normalised_Lateral_PSF_at_True_Scatterer_Depth.png"),
  ),
  caption: [(a) Signed migrated amplitude for Kirchhoff, Gazdag, and
    back-propagation migration overlaid at $f_c = 15 "GHz"$; (b) the
    normalised lateral point-spread function at the true scatterer depth,
    swept across separations from $2 lambda$ to $1 \/ 16 lambda$.],
) <fig:res-psf>

All three migration algorithms collapse the two scatterer hyperbolae into
distinguishable amplitude peaks for separations down to roughly half a
wavelength, but the two peaks progressively merge into a single lobe as the
separation shrinks further: the amplitude image alone cannot certify two
scatterers, or a sub-wavelength displacement of one scatterer, below this
floor.

#draftnote[state the precise separation at which the methods stop resolving
two distinguishable PSF peaks, read directly off @fig:res-psf (b), and
comment on any difference between the three algorithms.]

This amplitude-based floor is the motivation for the phase-plane approach
developed in this chapter and tested in @ch:hyp1 and @ch:hyp2.
