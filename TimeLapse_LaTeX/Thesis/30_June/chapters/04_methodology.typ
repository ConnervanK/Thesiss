#import "../template.typ": *

= General Materials and Methods <ch:methodology>

This chapter gathers the simulation and processing methodology shared by the
synthetic experiments of @ch:hyp1 and @ch:hyp2. Everything common to two or
more of those experiments is described here once --- the gprMax forward model,
the signal-conditioning steps, the three migration implementations, and the 2D
phase-plane estimator derived in @sec:th-wls --- so that each experimental
chapter need only add the scatterer configuration and displacement specific to
it. The chapter closes with a validation experiment (@sec:meth-resolution)
that establishes the amplitude-based resolution floor of the migration
algorithms before they are used to test any hypothesis.

== Forward Modelling with gprMax

All B-scans are simulated with the open-source finite-difference time-domain
solver gprMax @gprmax. The source is a Ricker wavelet with centre frequency
$f_c = 1.5 "GHz"$ and time-zero offset $t_0 = 0.943 "ns"$ (Supplementary
Material, §S1.1), chosen so that its usable bandwidth defines the dominant
wavelength $lambda$ used to express every displacement scale in this thesis
($2 lambda$ down to $1 \/ 32 lambda$). The computational domain is discretised
on a uniform $1 "mm"$ grid with perfectly-matched-layer (PML) absorbing
boundaries. A zero-offset (collocated transmitter and receiver) survey is
simulated by sweeping a single transmitter--receiver pair across the surface.

#draftnote[the figure titles encode the domain extent as "4010 m"; the same
auto-titling code elsewhere strips decimal points from floats (e.g. a depth of
0.676 m appears as `0676_m`, and a position of 2.0 m appears as `20_m`), so
this almost certainly reads as a domain of ≈4.01 m rather than 4010 m ---
confirm against the notebook before quoting a final value.]

== Scatterer and Medium Models

One target geometry is used across the synthetic experiments. *Point
scatterers* (@sec:meth-resolution, @ch:hyp1) are perfect-electric-conductor
(PEC) cylinders of radius $r = 28 "mm"$, buried at a depth of $0.676 "m"$ in
ice. @sec:meth-resolution places two such cylinders at a swept separation;
@ch:hyp1 instead holds one cylinder fixed as a baseline and displaces a second,
in the lateral, vertical, or diagonal direction, by the same family of
sub-wavelength steps.

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
migration, and gprMax-based time-reversal back-propagation. All three share the
implementation in `helper_functions/migration.py` (`PylopsKirchoffMigration`,
`gazdag_migration`, and `write_backprop_files`) and the exploding-reflector
convention $#vmig = v \/ 2$ of @eq:vmig, so that the same velocity model and the
same migration aperture are used for a baseline/monitor pair, which is required
for the displacement estimate of @sec:meth-phaseplane to be valid.

== The 2D Phase-Plane Shift-Estimation Pipeline <sec:meth-phaseplane>

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

== Validation: Amplitude Resolution Floor of the Migration Algorithms <sec:meth-resolution>

Before any of the hypotheses can be tested with confidence, the migration
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
B-scan.

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

#supp-note[The individual migrated image for every separation scenario, for
Kirchhoff, Gazdag, and back-propagation migration respectively (zoomed
around the true scatterer depth), is provided in the Supplementary
Material, §S1.2.]

@fig:res-psf (a) overlays the signed migrated amplitude from all three
algorithms at $f_c = 1.5 "GHz"$, and @fig:res-psf (b) zooms on the
Baseline-versus-Monitor point-spread function at the true scatterer depth,
for every separation scenario and method, mirroring the amplitude-test
figures used throughout @ch:hyp1 and @ch:hyp2 (@sec:hyp1-lat-amplitude):
the shaded band marks the FWHM measured once from the widest ($2 lambda$)
separation, where the two scatterers' responses do not yet overlap, and the
dashed lines mark the true position of each scatterer.

#page(flipped: true)[
#figure(
  subfigs(cols: 1,
    img("RES_016_Resolution_Study_--_Migration_Comparison_f_c15_GHz_aperture4.png"),
    img("RES_017_Resolution_Study_--_PSF_Zoom_Baseline_vs_Monitor.png"),
  ),
  caption: [(a) Signed migrated amplitude for Kirchhoff, Gazdag, and
    back-propagation migration overlaid at $f_c = 1.5 "GHz"$; (b) zoomed
    Baseline-versus-Monitor PSF at the true scatterer depth, swept across
    separations from $2 lambda$ to $1 \/ 16 lambda$: shaded = reference
    ($2 lambda$) FWHM at each true scatterer position; dashed = true
    scatterer positions.],
) <fig:res-psf>
]

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
developed in @ch:theory and tested in @ch:hyp1 and @ch:hyp2.
