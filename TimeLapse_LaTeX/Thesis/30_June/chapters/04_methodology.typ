#import "../template.typ": *

= General Materials and Methods <ch:methodology>

This chapter gathers the simulation and processing methodology shared by the
synthetic experiments of @ch:hyp1 and @ch:hyp2. Everything common to the experiments in @ch:hyp1, @ch:hyp2, and @ch:hyp3 is described here once --- the gprMax forward model,
the signal-conditioning steps, the three migration implementations, and the 2D
phase-plane estimator derived in @sec:th-wls --- so that each experimental
chapter only needs to describe the physical behaviour of a moving scatterer or fluid front. The chapter
closes with a validation experiment (@sec:meth-resolution)
that establishes the amplitude-based resolution floor of the migration
algorithms before they are used to test any hypothesis.

== Forward Modelling with gprMax

All B-scans are simulated with the open-source finite-difference time-domain
solver gprMax @gprmax. The source is a Ricker wavelet with centre frequency
$f_c = 1.5 "GHz"$ and time-zero offset $t_0 = 0.943 "ns"$ (Supplementary
Material, §S1.1) --- gprMax delays the Ricker pulse by $t_0$ so that it
starts from an amplitude of effectively zero rather than being switched on
abruptly at $t=0$, avoiding the high-frequency numerical noise such a
discontinuity would otherwise inject into the FDTD update --- chosen so
that its usable bandwidth defines the dominant
wavelength $lambda$ used to express every displacement scale in this thesis
($2 lambda$ down to $1 \/ 32 lambda$). The computational domain is discretised
on a uniform $1 "mm"$ grid --- chosen against $4 "GHz"$, the highest
significant frequency component of the Ricker wavelet, so that even that
shortest wavelength in ice is sampled by more than ten cells, keeping
numerical dispersion of the FDTD grid negligible --- with
perfectly-matched-layer (PML) absorbing boundaries, and measures
$4 times 1 "m"$ in $(x, z)$. A zero-offset
(collocated transmitter and receiver) survey is simulated by sweeping a
single transmitter--receiver pair across the surface in $1 "cm"$ steps,
with the receiver trailing the source by a fixed $10 "cm"$ offset; each
resulting position yields one A-scan, and $380$ such steps are stacked to
form each B-scan. Although gprMax records all six field/flux-density
components at every receiver, only the $E_z$ component (the vertical
electric field, in $"V"\/"m"$) is extracted to build every A-scan and
B-scan in this thesis. Each simulation runs for a $20 "ns"$ time window, stepped
at $Delta t approx 2.36 "ps"$ per iteration ($8481$ iterations in total); gprMax
determines this timestep itself from the $1 "mm"$ spatial discretisation,
choosing the largest value that still satisfies the Courant-Friedrichs-Lewy
(CFL) stability criterion for the FDTD update @gprmax.


== Scatterer and Medium Models

One target geometry is used across the synthetic experiments. *Point
scatterers* (@sec:meth-resolution, @ch:hyp1) are perfect-electric-conductor
(PEC) cylinders of radius $r = 28 "mm"$, buried at a depth of $0.676 "m"$ in
ice. @sec:meth-resolution places two such cylinders at a swept separation.
@ch:hyp1's vertical and diagonal scenarios instead use a single cylinder,
simulated once at its original position (baseline) and again after a
sub-wavelength displacement (monitor); its lateral scenarios additionally
keep a second cylinder fixed alongside the moving one, present in both the
baseline and monitor B-scans. This fixed scatterer serves as a check that
differencing the monitor and baseline images successfully cancels a static
target, leaving only the moving one behind. All three directions sweep the
same family of sub-wavelength steps.

#linebreak()

The host medium is ice with relative permittivity $epsilon_r = 3.15$ and
electric conductivity $sigma = 1 times 10^(-6) "S"\/"m"$, giving a
propagation velocity $v = c \/ sqrt(epsilon_r) approx 0.169 "m"\/"ns"$ and,
at the source centre frequency $f_c = 1.5 "GHz"$, a dominant wavelength
$lambda approx 112.6 "mm"$. The scatterer itself is modelled as an ideal
PEC boundary rather than being assigned a permittivity or conductivity, so
that its response is governed purely by its geometry. Both materials are
non-magnetic: relative permeability $mu_r = 1$ and magnetic loss $sigma^* =
0$. The cylinder radius,
$r = 28 "mm"$, is therefore about a quarter of the dominant wavelength
($approx 0.25 lambda$), and its burial depth of $0.676 "m"$ places it
roughly six wavelengths ($approx 6 lambda$) below the surface --- well
beyond the near-field region, so the scatterer is illuminated by a locally
planar wavefront.

== Signal Conditioning Pipeline <sec:meth-conditioning>

Every raw B-scan is processed identically before migration:

+ *Background subtraction.* A background-only simulation (no scatterer) is
  subtracted trace-by-trace to suppress the direct air/ground wave and isolate
  the scatterer reflection (e.g. @fig:res-bscans).

+ *Tapering and $t_0$ alignment.* An exponential decay taper suppresses
  late-arriving energy, a cosine end-taper removes hyperbola tails at the edge
  of the migration aperture, and a static shift aligns the surface reflection
  to $t = 0$ (e.g. @fig:res-bscans).

+ *Noise injection (where stated).* @ch:hyp2 contaminates the conditioned
  B-scan with synthetic Laplace-distributed noise at $10%$ of the signal
  standard deviation, fitted from real field data.

== Migration Algorithms Implemented

Every conditioned B-scan in @sec:meth-resolution, @ch:hyp1, and @ch:hyp2 is
migrated with all three algorithms derived in @sec:th-migration --- Kirchhoff
delay-and-sum (PyLops zero-offset operator), Gazdag $f$-$k$ phase-shift
migration, and gprMax-based time-reversal back-propagation. All three share the exploding-reflector
convention $#vmig = v \/ 2$ of @eq:vmig, so that the same velocity model and the
same migration aperture are used for a baseline/monitor pair, which is required
for the displacement estimate of @sec:meth-phaseplane to be valid.

== The 2D Phase-Plane Shift-Estimation Pipeline <sec:meth-phaseplane>

The theory of @sec:th-fourier-shift and @sec:th-wls is applied to the
migrated images produced above via the function `estimate_shift_2d`,
reproduced in full in the Supplementary Material, §S1.3. The end-to-end
workflow is:

+ *Migrate* the baseline and monitor B-scan with an identical velocity model
  and aperture (@sec:th-migration).

+ *Crop a region of interest (ROI)* tightly around the target (in this
  thesis, a window of $plus.minus 2.5 lambda$ about the scatterer apex,
  located on the Hilbert-envelope peak of the baseline image). The apex is
  found via the high-amplitude region of the difference (monitor minus
  baseline) B-scan, which highlights the rough area the scatterer has moved
  through, so that the ROI stays clear of static background noise elsewhere
  in the image.

+ *Taper, then take the 2D FFT* of both cropped images and form the
  cross-spectrum $#XS$ of @eq:cross-spectrum-def.

+ *Mask and weight*: mask out frequency bins below an amplitude threshold,
  keeping only the high-energy part of the cross-spectrum phase, and weight
  the surviving bins by $|#XS|$ in the WLS fit; further restrict the fit to
  $|#kz|, |#kx| < 1.4 k_(z c)$, where the $1.4$ passband factor is an
  empirical value that must be tuned to avoid phase wraparound
  (@sec:th-mask-weight).

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

== Validation: Amplitude Resolution Floor of the Migration Algorithms <sec:meth-resolution>

Before any of the hypotheses can be tested with confidence, the signal conditioning
pipeline above is validated on a simple, well-understood baseline problem: how
closely can two _stationary_ point scatterers be spaced before the migration
algorithms can no longer tell them apart? This validation experiment uses two
PEC cylinders illuminated by a zero-offset GPR B-scan, migrated with all three
algorithms. Its result is the amplitude-based resolution floor against which
every displacement-detection result in @ch:hyp1 and @ch:hyp2 is later compared.

@fig:res-setup shows the forward-model setup: the gprMax domain and grid,
and the swept separation between two PEC cylinder scatterers ($r = 28 "mm"$,
depth $0.676 "m"$), from $2 lambda$ down to $1 \/ 16 lambda$. The source is a
Ricker wavelet with centre frequency $f_c = 1.5 "GHz"$ and time-zero offset
$t_0 = 0.943 "ns"$, shown in the Supplementary Material, §S1.1.

#figure(
  subfigs(cols: 1,
    img("RES_002_Resolution_Study__Model_Geometry__domain_4010_m_Δx__1_mm_PML.png"),
    img("RES_003_Scatterer_Positions__PEC_Cylinders__r__28_mm_depth__0676_m.png", width: 70%),
  ),
  caption: [Forward-model setup for the resolution validation: (a) the gprMax
    domain and grid; (b) the swept separation between the two PEC cylinder
    scatterers, $r = 28 "mm"$, depth $0.676 "m"$.],
) <fig:res-setup>

#supp-note[The Ricker source wavelet used throughout this thesis is shown in
the Supplementary Material, §S1.1.]

@fig:res-bscans shows the full signal-conditioning pipeline of
@sec:meth-conditioning on the $2 lambda$ separation scenario: the simulated
zero-offset B-scans before and after background subtraction, and the effect
of tapering and the $t_0$ shift on both a single trace and the complete
B-scan. The direct wave itself is removed by the background subtraction of
@sec:meth-conditioning and so is not visible in the single-trace panels; had
it still been present, it would start exactly at $t = 0$. gprMax's source
wavelet is instead excited with a small built-in delay ($t_0 approx
0.943 "ns"$ for the Ricker wavelet used throughout this thesis) purely to
avoid the numerical onset artefact of starting the simulation abruptly at
the wavelet's peak; the $t_0$ shift undoes exactly this delay so the surface
reflection aligns back to $t = 0$ in the conditioned trace.

#figure(
  subfigs(cols: 1,
    img("RES_004_GPR_B-Scans__Background_and_Separation_Models.png"),
    img("RES_005_GPR_B-Scans__Background_Subtracted.png"),
    img("RES_006_Effect_of_Tapering_and_t0_Shift__2λ_dataset_single_trace.png"),
    img("RES_007_B-scan_effect_of_tapering_and_t0_shift__2λ_dataset.png"),
  ),
  caption: [Signal-conditioning pipeline for the resolution validation's
    $2 lambda$ dataset (@sec:meth-conditioning): (a) raw zero-offset B-scans,
    background and separation models; (b) after background subtraction; (c)
    effect of tapering and the $t_0$ shift on a single representative trace;
    (d) the complete B-scan after tapering and the $t_0$ shift. Dashed lines
    in (a)--(b) mark the true baseline/monitor scatterer positions, the same
    convention used in @fig:res-psf.],
) <fig:res-bscans>


@fig:res-psf (a) overlays the signed migrated amplitude from all three
algorithms at $f_c = 1.5 "GHz"$, and @fig:res-psf (b) zooms on the
Baseline-versus-Monitor point-spread function at the true scatterer depth,
for every separation scenario and method, mirroring the amplitude-test
figures used throughout @ch:hyp1 and @ch:hyp2 (@sec:hyp1-lat-amplitude). The shaded band marks the FWHM (the full width, along $x$, of the migrated
point-spread function's main lobe at half its peak amplitude --- a standard
measure of how broad, and hence how resolvable, that lobe is), measured
from a scatterer's point-spread function. We use one FWHM value based on
the reference ($2 lambda$ case) for all the other separation cases. The FWHM value calculated in the reference case is used as a benchmark as the two scatterers are still clearly separated. The dashed
lines mark the true position of each scatterer.

#page(flipped: true)[
#figure(
  subfigs(cols: 2,
    img("RES_016_Resolution_Study_--_Migration_Comparison_f_c15_GHz_aperture4.png", width: 90%),
    img("RES_017_Resolution_Study_--_PSF_Zoom_Baseline_vs_Monitor.png", width: 87%),
  ),
  caption: [(a) Signed migrated amplitude for Kirchhoff, Gazdag, and
    back-propagation migration overlaid at $f_c = 1.5 "GHz"$; (b) zoomed
    Baseline-versus-Monitor PSF at the true scatterer depth, swept across
    separations from $2 lambda$ to $1 \/ 16 lambda$: shaded = reference
    ($2 lambda$) FWHM at each true scatterer position; dashed = true
    scatterer positions.],
) <fig:res-psf>
]

#supp-note[The individual migrated image for every separation scenario, for
Kirchhoff, Gazdag, and back-propagation migration respectively (zoomed
around the true scatterer depth), is provided in the Supplementary
Material, §S1.2.]

All three migration algorithms collapse the two scatterer hyperbolae into distinguishable amplitude peaks for separations down to roughly half a wavelength, but the two peaks progressively merge into a single lobe as the separation shrinks further. The amplitude image alone cannot certify two scatterers, or a sub-wavelength displacement of one scatterer, below this floor—formally known as the Rayleigh diffraction limit. This fundamental amplitude-based limit is the motivation for the phase-plane approach developed in @ch:theory and tested in @ch:hyp1 and @ch:hyp2.

#linebreak()

Read directly off @fig:res-psf (b), the three algorithms differ noticeably
in exactly where that floor falls: Kirchhoff stops resolving two
distinguishable peaks below a separation of $1 lambda$, Gazdag remains
resolvable down to $1 \/ 2 lambda$, and back-propagation --- the
best-performing of the three --- barely still resolves two peaks at
$1 \/ 4 lambda$.
