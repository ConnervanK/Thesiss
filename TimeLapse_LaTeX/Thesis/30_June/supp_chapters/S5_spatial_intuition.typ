#import "../template.typ": *

= Spatial and Time-Frequency Intuition for the Lateral/Vertical Asymmetry <supp:spatial-intuition>

Chapter 3's derivation of the lateral/vertical asymmetry (Theoretical
Background, "Space-Wavenumber Duality: Lateral versus Vertical Asymmetry")
is analytic: a parabolic approximation of the migrated point-spread function
shows that lateral position acts locally as a proxy for wavenumber, while
vertical position does not, so a lateral displacement leaves a phase
_slope_ and a vertical displacement leaves a phase _plateau_. This chapter
builds spatial intuition for that result with an independent, local
diagnostic --- the spatial instantaneous phase, extracted trace-by-trace
with a Riesz/Hilbert transform --- applied directly to the same
Gazdag-migrated lateral- and vertical-displacement datasets used throughout
the thesis. Because this local, trace-based estimator has so far only been
validated qualitatively against Chapter 5's global phase-plane (WLS) fit,
and not yet compared to it quantitatively on the same data, it also
constitutes the supporting evidence for "Hypothesis 1.5: Local
Phase-Gradient Methods", discussed in Chapter 5 and revisited in Chapter 8
--- the figures below should be read as a qualitative complement to, not a
replacement for, Hypothesis 1's primary evidence.

@supp:spatial-lateral and @supp:spatial-vertical confirm the lateral/vertical
asymmetry directly in simulation. @supp:spatial-timefreq then places the
underlying time-frequency machinery in its correct context: a complementary,
physically equivalent lens on the same duality result, most of which remains
unexplored future work rather than a demonstrated result.

== Lateral Displacement: Instantaneous-Phase Confirmation <supp:spatial-lateral>

The 2D Riesz transform is applied directly to the Gazdag-migrated baseline
and monitor images, producing an instantaneous-phase difference map and, at
the scatterer depth, a 1D phase cross-section whose slope at the
zero-crossing encodes the lateral displacement $#Dx$ --- the simulation-side
counterpart of the parabolic-phase argument in Chapter 3 ("Lateral
direction: position as a proxy for wavenumber"). @fig:tlp-instphase-horiz-example
shows the representative $1\/4 lambda$ case ($#Dx = 28.0 "mm"$), and
@fig:tlp-instphase-horiz-summary summarises the fitted zero-crossing slope
across all seven simulated scales, from $2 lambda$ down to $1\/32 lambda$:
from $1\/2 lambda$ down to $1\/16 lambda$ the slope closely tracks
proportionality with separation, roughly halving each time the separation
halves, the simulation-side signature of the same linear relationship
Chapter 3 derives analytically. This proportionality degrades at the two
largest scales tested: the $1 lambda$ slope no longer exceeds the
$1\/2 lambda$ slope, and at $2 lambda$ the phase has wrapped enough that no
clean zero-crossing exists at all, marking the edge of validity of the
local near-apex approximation underlying this diagnostic.

#figure(
  subfigs(cols: 2,
    img("TLP_059_Instantaneous_Phase__Gazdag____Baseline_vs_¼λ___Δx__280_mm.png"),
    img("TLP_060_Phase_cross-section__z__676_mm____Gazdag____¼λ__Δx__280_mm.png"),
    img("TLP_061_Phase_cross-section__zoomed__½λ_around_scatterers.png"),
  ),
  caption: [Instantaneous-phase analysis for the representative
    $#Dx = 1\/4 lambda = 28.0 "mm"$ lateral displacement, Gazdag-migrated:
    (a) the 2D wrapped phase-difference map; (b) the 1D phase cross-section
    at the scatterer depth $z = 676 "mm"$; (c) the same cross-section zoomed
    to $plus.minus 1\/2 lambda$ around the scatterer.],
) <fig:tlp-instphase-horiz-example>

#figure(
  img("TLP_071_Gazdag__Δφ_zero-crossing_slope_vs_scatterer_separation.png", width: 70%),
  caption: [Fitted zero-crossing slope of the lateral instantaneous
    phase-difference cross-section, as a function of true scatterer
    displacement, across all seven scales from $2 lambda$ to $1\/32 lambda$
    (Gazdag-migrated); no zero-crossing exists at $2 lambda$. The slope
    scales approximately in proportion to separation from $1\/2 lambda$
    down to $1\/16 lambda$ (halving each time the separation halves),
    consistent with Chapter 3's linear near-apex approximation.],
) <fig:tlp-instphase-horiz-summary>

== Vertical Displacement: Instantaneous-Phase Confirmation <supp:spatial-vertical>

An equivalent instantaneous-phase analysis is applied to the
vertical-displacement dataset, with the cross-section now taken along $z$ at
the fixed lateral position $x = 2.0 "m"$ of the scatterer, across six scales
from $1 lambda$ to $1\/32 lambda$. @fig:tlp-instphase-vert-example shows the
representative $1\/8 lambda$ case ($#Dz = 14.0 "mm"$); unlike the lateral
case, the migrated image does not separate phase spatially along $z$, so
instead of a zero-crossing slope, the diagnostic reports the mean phase
difference $lr(⟨ Delta phi ⟩)$ between the two scatterer depths
(panel (c), dashed purple line). @fig:tlp-instphase-vert-summary summarises
this mean phase difference across five of the six scales simulated (the
$1 lambda$ case is excluded: at that separation the two scatterer echoes no
longer overlap coherently, so no pixel between them clears the amplitude
mask and no mean can be formed).

#figure(
  subfigs(cols: 2,
    img("TLP_081_VerticalTimeLapse__Instantaneous_Phase__Gazdag____Baseline_v.png"),
    img("TLP_082_Phase_cross-section__x__2000_mm____Gazdag____⅛λ__Δz__140_mm.png"),
    img("TLP_083_Phase_cross-section__zoomed__½λ_around_scatterer_depths.png"),
  ),
  caption: [Instantaneous-phase analysis for the representative
    $#Dz = 1\/8 lambda = 14.0 "mm"$ vertical displacement, Gazdag-migrated:
    (a) the 2D wrapped phase-difference map; (b) the 1D phase cross-section
    at $x = 2000 "mm"$; (c) the same cross-section zoomed to
    $plus.minus 1\/2 lambda$ around the scatterer depths, with the mean
    phase difference between the two scatterer depths marked.],
) <fig:tlp-instphase-vert-example>

#figure(
  img("TLP_090_VerticalTimeLapse__Gazdag__Mean_Δφ_vs_scatterer_vertical_shi.png", width: 70%),
  caption: [Mean instantaneous phase difference $lr(⟨ Delta phi ⟩)$
    as a function of true vertical scatterer shift, Gazdag-migrated, across
    five scales from $1\/2 lambda$ to $1\/32 lambda$ ($1 lambda$ excluded,
    see text). Unlike the lateral case (@fig:tlp-instphase-horiz-summary),
    this is a plateau level rather than a cross-section slope, consistent
    with the asymmetry derived in Chapter 3. The raw mean phase difference
    is wrapped to $(-pi, pi]$ by construction; at $1\/4 lambda$ and
    $1\/2 lambda$ the true, physically continuous phase has already
    accumulated more than one full cycle, so the wrapped value folds back
    into $(-pi, pi]$ and appears discontinuous with the trend at smaller
    scales. Both points have $2 pi$ subtracted here to restore the
    physically continuous, monotonically decreasing curve],
) <fig:tlp-instphase-vert-summary>

== Time-Frequency Perspective: A Complementary, Largely Unexplored View <supp:spatial-timefreq>

Chapter 5's phase-plane fit treats the migrated image as a whole. A
localised, trace-based alternative instead analyses individual unmigrated or
migrated _traces_ with a sliding time-frequency transform, which gives
access to _when_ (in two-way time) a phase change occurs, complementing the
spatial picture of @supp:spatial-lateral and @supp:spatial-vertical above.
This view is developed analytically below for completeness, but --- unlike
the instantaneous-phase confirmation above --- it remains mostly future
research. It is included here as physical intuition for, not a
demonstrated confirmation of, the duality result.

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
Forming the local cross-spectrum in the same way as the global cross-spectrum
of Chapter 3 (baseline spectrum times the conjugate of the monitor spectrum)
and taking its angle, the baseline phase and amplitude cancel, leaving the
local, time-domain analogue of the global phase-plane relation derived in
Chapter 3,
$ Delta Phi(tau, f) approx -2 pi f thin #Dt + #Dtheta . $ <eq:local-phase-line>
At the two-way time $tau_0$ of the target reflection, a straight-line fit of
$Delta Phi$ against $f$ has slope $-2 pi #Dt$ (the mechanical shift) and
intercept $#Dtheta$ (a calibration or material-change offset) --- the exact
time-domain counterpart of the slope/intercept decomposition in Chapter 3's
Weighted Least-Squares Plane Fitting. Converting between the two-way-time
and depth pictures uses the standard relation
$ #Dt = (2 #Dz) / v , $ <eq:dt-dz>
so that $omega #Dt equiv #kz #Dz$ with $#kz = 2 omega \/ v$: the temporal-frequency
slope and the vertical-wavenumber slope are the same physical quantity, viewed
in two different but exactly equivalent coordinate systems.

=== Three diagnostic views for a future local-phase estimator

Three derived plots would make the local phase line directly visible in the
data, mirroring the views used for the global fit in Chapter 5:

/ Spectral line, $Delta Phi(f)$ at fixed $tau_0$: evaluating the local
  phase-line relation above at the target's two-way time gives a straight
  line in $f$ whose slope and intercept separate the two causes directly ---
  a line through the origin (slope only, $#Dtheta = 0$) for pure mechanical
  movement; a flat line offset from zero (intercept only, $#Dt = 0$) for a
  pure phase-rotation offset; a sloped line with non-zero intercept for a
  combination of the two.

/ Cross-phase spectrogram, $Delta Phi(tau, f)$: the same relation, now
  shown across every window centre $tau$, not only $tau_0$. Outside the
  target reflection the baseline and monitor traces share no coherent
  phase, so this is incoherent, salt-and-pepper noise; at the target's
  two-way time a coherent window appears, and within it the $-2 pi f #Dt$
  term winds through repeated $2 pi$ cycles as $f$ increases, producing a
  vertical fringe pattern for movement, versus a uniform colour block where
  the phase is constant across $f$ for a pure phase offset.

/ Polar vector rotation: the complex coefficient at the dominant
  frequency and peak two-way time, plotted as a vector in the complex
  plane for baseline and monitor. Both $#Dt$ and $#Dtheta$ enter the
  localised Fourier shift theorem above as unit-magnitude phase factors, so
  either cause --- or their combination --- rotates this vector with
  negligible change in length; a genuine amplitude change (e.g. altered
  reflectivity) would
  instead shrink or grow it. This view therefore distinguishes a purely
  geometric/phase shift from an amplitude-affecting one, though separating
  $#Dt$ from $#Dtheta$ still requires the multi-frequency slope of the
  spectral-line view above.

These local, trace-based views and the global 2D wavenumber fit of Chapter
3's Weighted Least-Squares Plane Fitting are not competing techniques: they
are Fourier duals of the same underlying physics, related by a spatial
Fourier transform of the time-frequency decomposition with respect to the
lateral coordinate $x$, which reduces (under a smooth-window approximation)
directly to the global spectrum $U(omega, #kx)$ used throughout Chapter 3.
The wavenumber-domain fit remains the primary quantitative tool used
throughout this thesis because it linearises the spatial curvature of the
migrated point-spread function exactly and admits amplitude weighting
natively; the local, time-frequency view would, if implemented, let a
future analysis localise _where_ (in $tau$) a phase anomaly originates,
which the global fit alone cannot show.
