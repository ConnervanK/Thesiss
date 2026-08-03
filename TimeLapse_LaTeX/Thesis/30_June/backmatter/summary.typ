#import "../template.typ": *

#heading(level: 1, numbering: none, outlined: true)[Summary] <ch:summary>

// Draft summary --- confirm the headline figures against the final results
// tables and the General Discussion before submission.

This thesis addressed one unifying hypothesis --- the research question: _can
time-lapse ground-penetrating radar accurately track subwavelength movement, a
form of super-resolution, by analysing phase changes in migrated images rather
than their amplitude?_ A single migrated amplitude image is a static snapshot
that cannot show change at all --- and even differencing two such images,
baseline against monitor, still collapses below a wavelength-scale
resolution floor, such as Rayleigh's limit. The method instead migrates the baseline and monitor
survey, forms their Fourier-domain cross-spectrum, and
fits its phase with a weighted least-squares plane whose two slopes give the
vertical and lateral displacement estimates of the movement of a scatterer or moving fluid front. The hypothesis was tested through
three experiments of increasing realism, on synthetic gprMax data across three
migration algorithms --- Kirchhoff, Gazdag phase-shift, and time-reversal
back-propagation --- and on real borehole field data.

#linebreak()

*Hypothesis 1* established that amplitude differencing of migrated time-lapse
images fails below a resolution floor of roughly half a wavelength, whereas the
phase-plane estimator recovers lateral, vertical, and diagonal displacements
down to $1 \/ 32 lambda$ --- to a few tenths of a millimetre, in the sub-wavelength range --- exactly where amplitude based separation has already
collapsed, and does so for a distributed fluid front as well as for a point
scatterer.

#linebreak()

*Hypothesis 2* showed that, under a realistic Laplace noise model, no single
migration technique is unconditionally the most noise-robust. By raw accuracy
Kirchhoff is the best, but it alone manufactures false coherent structure from
pure noise; back-propagation with sign-bit time-reversal --- which clamps every
noise spike to $plus.minus 1$ while preserving phase --- is the more
conservative choice; Gazdag is the weakest on both counts.

#linebreak()

*Hypothesis 3* applied the full pipeline to real borehole GPR field data from a
controlled fluid-injection experiment (38 profiles, 6 June 2016). It produced
stable displacement estimates: Kirchhoff-BP, Gazdag, and back-propagation ---
three independently-processed migration techniques --- agree on the direction of the recovered movement,
whose dominant vertical component also matches the injection geometry.
Fluid is actively pumped in at $78.7 "m"$ depth, within the same depth band
as the tracked reflector, and Gazdag and back-propagation put the front's
upward motion at roughly $1.6$--$2.1 "m"$ during Push and Chase, reversing
downward during Wait and Pull once pumping stops. A smaller lateral shift toward the borehole
during Push, reversing away from it in later stages, is not directly
explained by the injection geometry and may instead reflect the fracture
network's orientation relative to the borehole.

#linebreak()

The unifying hypothesis is therefore supported: time-lapse GPR can recover
subwavelength movement from the phase of migrated images at scales where the
amplitude difference-image cannot, subject to a migration-choice trade-off
under noise. That three independently-processed migration techniques recover closely matching
displacement estimates from the same real field data is strong evidence
that the phase-plane fit measures a genuine physical signal.
