#import "../template.typ": *

= Decoupling Geometric Movement from Material Change <sec:th-material-change>

The third column of the design matrix in @eq:design-matrix is the
constant vector $bold(1)$, deliberately included alongside the two
wavenumber columns. This appendix shows, for theoretical completeness, why that
design choice lets the same fit simultaneously measure displacement _and_
detect a change in the dielectric material at the target, and why the two
never contaminate each other.

== Why a sub-wavelength fracture produces a frequency-independent phase shift

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

== Why the WLS fit puts material change exactly into the intercept

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
