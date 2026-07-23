#import "../template.typ": *

= Hypothesis 3: Supplementary Figures <supp:hyp3>

This chapter mirrors the section structure of Chapter 6, "Hypothesis 3:
Generalisation to Complex Scenes and Real Field Data", of the main thesis,
in the same way as @supp:hyp1 and @supp:hyp2. The Region of Influence
Workflow section below now collects the Chase/Wait/Pull rectangular-window
stage figures and the sliding-window scan, moved out of the main chapter's
"Region of Influence Workflow" section, which retains only the Push-stage
example and a #raw("supp-note") pointer to this chapter.

== Field Data Explanation <supp:hyp3-fielddata>

== Processing <supp:hyp3-fd-processing>

=== Cross-Profile Registration and Amplitude Normalisation <supp:hyp3-fd-crossprofile>

== Region of Influence Workflow <supp:hyp3-fd-roi>

=== Rectangular Window <supp:hyp3-fd-roi-rect>

The Chase, Wait, and Pull stages of the fixed rectangular ROI, following the
same construction shown in the main thesis for the Push stage: the
time-lapse difference image and its monogenic envelope with
the ROI overlaid, computed independently for Kirchhoff-BP, Gazdag, and
back-propagation. The coherent patch inside the ROI shrinks progressively
from Chase through Wait to Pull, mirroring the weakening signal expected as
the fluid front decelerates.

#figure(
  cimg("FD_stage_diff_envelope_chasing.png"),
  caption: [Chase stage (profiles 5→9): time-lapse difference and monogenic
    envelope, ROI overlaid, for Kirchhoff-BP, Gazdag and back-propagation.],
) <fig:fd-stage-chase>

#figure(
  cimg("FD_stage_diff_envelope_waiting.png"),
  caption: [Wait stage (profiles 10→20): time-lapse difference and monogenic
    envelope, ROI overlaid, for Kirchhoff-BP, Gazdag and back-propagation.],
) <fig:fd-stage-wait>

#figure(
  cimg("FD_stage_diff_envelope_pulling.png"),
  caption: [Pull stage (profiles 21→38): time-lapse difference and monogenic
    envelope, ROI overlaid, for Kirchhoff-BP, Gazdag and back-propagation.],
) <fig:fd-stage-pull>

=== Sliding Window <supp:hyp3-fd-roi-sliding>

A fixed-size window slid across the entire difference B-scan, with the
phase-plane fit repeated independently in every window position to build up
a 2-D map of $#Dz$ and $#Dx$. The window size is
set as a multiple of the dominant wavelength $lambda = v \/ f_0 = 1.0 "m"$,
with independent multiples for the depth and radial extent: the
sub-horizontal reflector geometry means coherent structure in depth extends
over many wavelengths ($k_z approx 0$), while the two energy lobes that carry
the displacement information sit at $k_x approx plus.minus 8 "rad/m"$ (a
radial wavelength of only $approx 0.8 "m"$), so a window narrower than about
one radial wavelength cannot resolve a stable phase gradient in $#Dx$ and
simply fits noise. The window advances by a fixed number of pixels each
step -- large enough to keep the scan fast, with an option to shrink it once
the window size itself has been validated.

A second gate is needed beyond the ROI-window's own amplitude threshold:
because that threshold is relative to *each window's own* spectral peak, a
window sitting entirely in noise can still return a "confident"-looking but
meaningless fit. An additional absolute-energy gate -- requiring a window's
mean $|Delta"amplitude"|$ to reach a fraction of the pair's global
98^"th"-percentile amplitude -- discards these before they reach the map.

@fig:fd-sliding-qc shows the resulting window grid for the Push-stage pair
(window $3.5 "m"$ deep × $0.96 "m"$ radial, $4424$ overlapping positions),
overlaid on the difference B-scan with one window highlighted at full size
for scale. @fig:fd-sliding-maps shows the resulting $#Dz$ / $#Dx$ maps for
all four pairs. In every pair, the surviving (non-blanked) region collapses
onto a single spatially-coherent patch at approximately $68$--$83 "m"$ depth
and $4$--$8 "m"$ radial distance -- the same location as the hand-picked
rectangular ROI of @supp:hyp3-fd-roi-rect, recovered here with no prior knowledge of where the
reflector was. The Push and Chase maps (1→3, 3→8) show the largest-magnitude,
best-defined patches, matching the strongest expected displacement; the Wait
map (8→20) is both smaller in magnitude (colour range roughly a third of
Push/Chase) and less coherent, consistent with near-zero net movement during
the paused-injection stage. This is an independent cross-check of both the
reflector's *location* and the qualitative *stage-to-stage trend* established
by the rectangular-window results, obtained by a method that never had the
ROI told to it.

#figure(
  cimg("FD_sliding_window_qc.png"),
  caption: [Sliding-window scan grid for the Push-stage pair (1→3), Gazdag
    migration: every scanned window overlaid (thin yellow) on the time-lapse
    difference B-scan, with one window highlighted (green) at full size.],
) <fig:fd-sliding-qc>

#figure(
  subfigs(cols: 1,
    cimg("FD_sliding_window_maps_1_to_3.png"),
    cimg("FD_sliding_window_maps_3_to_8.png"),
    cimg("FD_sliding_window_maps_8_to_20.png"),
    cimg("FD_sliding_window_maps_20_to_38.png"),
  ),
  caption: [Sliding-window $#Dz$ / $#Dx$ maps, Gazdag migration, for the four
    representative pairs: (a) 1→3 (Push), (b) 3→8 (Chase), (c) 8→20 (Wait),
    (d) 20→38 (Pull). Windows failing either the WLS mask-count or the
    absolute-energy gate are left blank.],
) <fig:fd-sliding-maps>

=== Picking Pixels <supp:hyp3-fd-roi-picking>

== Phase-Plane Fit Workflow <supp:hyp3-fd-phaseplane>

=== Kirchhoff-BP Diagnostic <supp:hyp3-fd-disp-kirchhoff>

#figure(
  cimg("FD_displacement_diagnostics_kirchhoff_bp.png"),
  caption: [Kirchhoff-BP: WLS cross-spectrum phase-plane diagnostics for the
    four stage-boundary pairs (rows) -- cross-spectrum phase, cross-spectrum
    energy with WLS threshold contour, fitted plane, 1-D $#kz$ slice, 1-D
    $#kx$ slice (columns). Same panel layout as the Gazdag diagnostic in the
    "Phase-Plane Fit Workflow" section of Chapter 6.],
) <fig:fd-disp-kirchhoff>

=== Back-Propagation Diagnostic <supp:hyp3-fd-disp-backprop>

#figure(
  cimg("FD_displacement_diagnostics_backprop.png"),
  caption: [Back-propagation: WLS cross-spectrum phase-plane diagnostics for
    the four stage-boundary pairs (rows), same panel layout as the Gazdag
    diagnostic in the "Phase-Plane Fit Workflow" section of Chapter 6. This
    is the 14-profile, homogeneous-domain result discussed in Chapter 6's
    Interpretation section as disagreeing with Kirchhoff/Gazdag on the sign
    of $#Dz$ before the corrected back-propagation pipeline's independent
    re-derivation.],
) <fig:fd-disp-backprop>

== Corrected Back-Propagation Displacement Re-Estimation <supp:hyp3-fd-bp-corrected>

== Interpretation <supp:hyp3-fd-interpretation>
