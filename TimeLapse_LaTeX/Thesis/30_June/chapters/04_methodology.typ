#import "../template.typ": *

= Methodology <ch:methodology>

This chapter describes the simulation and processing pipeline shared by the
experiments in @ch:hyp1, @ch:hyp2, and @ch:hyp3. Each experiment varies the
scatterer/material configuration and the displacement or material change
under test, but reuses the same forward model, signal-conditioning steps,
migration implementations, and (from @ch:hyp1 onward) the same phase-plane
estimator. The chapter closes with a short validation experiment
(@sec:meth-resolution) that establishes the amplitude-based resolution floor
of the migration algorithms before they are used to test any of the three
hypotheses.

== Forward Modelling with gprMax

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

== Scatterer and Medium Models

Two target geometries are used across the experiments:

- *Point scatterers* (@sec:meth-resolution, @ch:hyp1): perfect-electric-conductor
  (PEC) cylinders of radius $r = 28 "mm"$, buried at a depth of $0.676 "m"$
  in ice. @sec:meth-resolution places two such cylinders at a swept
  separation; @ch:hyp1 instead holds one cylinder fixed as a baseline and
  displaces a second, in the lateral, vertical, or diagonal direction, by
  the same family of sub-wavelength steps.

- *A distributed fluid front* (@ch:hyp2): a sub-wavelength fracture whose
  air-filled and water-filled regions are separated by a front that advances
  laterally between baseline and monitor surveys, following the thin-layer
  reflectivity model of @sec:th-material-change.

== Signal Conditioning Pipeline <sec:meth-conditioning>

Every raw B-scan is processed identically before migration:

+ *Background subtraction.* A background-only simulation (no scatterer) is
  subtracted trace-by-trace to suppress the direct air/ground wave and
  isolate the scatterer reflection (e.g. @fig:res-bscans).

+ *Tapering and $t_0$ alignment.* An exponential decay taper suppresses
  late-arriving energy, a cosine end-taper removes hyperbola tails at the
  edge of the migration aperture, and a static shift aligns the surface
  reflection to $t = 0$ (e.g. @fig:res-taper).

+ *Noise injection (where stated).* @ch:hyp3 contaminates the conditioned
  B-scan with synthetic Laplace-distributed noise at $10%$ of the signal
  standard deviation, fitted from real field data.

== Migration Algorithms Implemented

Every conditioned B-scan in @sec:meth-resolution, @ch:hyp1, @ch:hyp2, and
@ch:hyp3 is migrated with all three algorithms derived in @sec:th-migration:
Kirchhoff delay-and-sum (PyLops zero-offset operator), Gazdag $f$-$k$
phase-shift migration, and gprMax-based time-reversal back-propagation. All
three share the implementation in `helper_functions/migration.py`
(`PylopsKirchoffMigration`, `gazdag_migration`, and `write_backprop_files`)
and the exploding-reflector convention $#vmig = v \/ 2$ of @eq:vmig, so that
the same velocity model and the same migration aperture are used for a
baseline/monitor pair, which is required for the displacement estimate of
@sec:meth-phaseplane to be valid.

== The 2D Phase-Plane Shift-Estimation Pipeline <sec:meth-phaseplane>

@sec:th-fourier-shift, @sec:th-wls, and @sec:th-material-change are applied
to the migrated images produced above via the function `estimate_shift_2d`
defined in `TimeLapse_Processing.ipynb`. The end-to-end workflow is:

+ *Migrate* the baseline and monitor B-scan with an identical velocity model
  and aperture (@sec:th-migration).

+ *Crop a region of interest (ROI)* tightly around the target (in this
  thesis, a window of $plus.minus 2.5 lambda$ about the scatterer or fracture
  apex, located on the Hilbert-envelope peak of the baseline image), so that
  static background structure elsewhere in the image cannot drag the fitted
  plane towards zero.

+ *Taper, then take the 2D FFT* of both cropped images and form the
  cross-spectrum $#XS$ of @eq:cross-spectrum-def.

+ *Mask and weight*: restrict the fit to $|#kz|, |#kx| < 1.4 k_(z c)$ and
  weight every bin by $|#XS|$ (@sec:th-mask-weight).

+ *Solve* the weighted least-squares system of @eq:wls-preweighted for
  $(#Dz, #Dx, c)$.

For the distributed fluid front of @ch:hyp2, step 2 instead slides a small
window along the known fracture geometry, and step 5 is read for its intercept
$c$ at each position to build a spatial profile of material change
(@eq:intercept-material), rather than a single displacement estimate.

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

Before any of the three hypotheses can be tested with confidence, the
migration pipeline above is validated on a simple, well-understood baseline
problem: how closely can two _stationary_ point scatterers be spaced before
the migration algorithms above can no longer tell them apart? This validation
experiment uses two of the PEC cylinders described above, illuminated by a
zero-offset GPR B-scan, and migrated with all three algorithms. Its result is
the amplitude-based resolution floor against which every
displacement-detection result in @ch:hyp1, @ch:hyp2, and @ch:hyp3 is later
compared.

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
developed in @ch:theory and tested against Hypotheses 1--3 in @ch:hyp1,
@ch:hyp2, and @ch:hyp3.
