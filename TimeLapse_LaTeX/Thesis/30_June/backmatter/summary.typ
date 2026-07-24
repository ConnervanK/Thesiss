#import "../template.typ": *

#heading(level: 1, numbering: none, outlined: true)[Summary] <ch:summary>

// Draft summary --- confirm the headline figures against the final results
// tables and the General Discussion before submission.

This thesis addressed one unifying hypothesis --- the research question: _can
time-lapse ground-penetrating radar accurately track subwavelength movement, a
form of super-resolution, by analysing phase changes in migrated images rather
than their amplitude?_ Because a single migrated amplitude image cannot resolve
a change smaller than a wavelength-scale floor, the method instead migrates a
baseline and a monitor survey, forms their Fourier-domain cross-spectrum, and
fits its phase with a weighted least-squares plane whose two slopes give the
vertical and lateral displacement and whose intercept separates geometric
movement from sub-wavelength material change. The hypothesis was tested through
three experiments of increasing realism, on synthetic gprMax data across three
migration algorithms --- Kirchhoff, Gazdag phase-shift, and time-reversal
back-propagation --- and on real borehole field data.

*Hypothesis 1* established that amplitude differencing of migrated time-lapse
images fails below a resolution floor of roughly half a wavelength, whereas the
phase-plane estimator recovers lateral, vertical, and diagonal displacements
down to $1 \/ 32 lambda$ --- to a few tenths of a millimetre, several orders of
magnitude below the wavelength scale --- exactly where amplitude has already
collapsed, and does so for a distributed fluid front as well as for a point
scatterer. An optional local phase-gradient analysis (Hypothesis 1.5) is
consistent with this result wherever both were computed, but remains
exploratory.

*Hypothesis 2* showed that, under a realistic Laplace noise model, no single
migration technique is unconditionally the most noise-robust. By raw accuracy
Kirchhoff is the best, but it alone manufactures false coherent structure from
pure noise; back-propagation with sign-bit time-reversal --- which clamps every
noise spike to $plus.minus 1$ while preserving phase --- is the more
conservative choice; Gazdag is the weakest on both counts.

*Hypothesis 3* applied the full pipeline to real borehole GPR field data from a
controlled fluid-injection experiment (38 profiles, 6 June 2016). It produced
stable, cross-technique-corroborated displacement estimates --- of order
$1.3$--$1.9 "m"$ during the Push stage, near-zero during Wait, and a partial
reversal during Pull --- with an independent back-propagation re-derivation
agreeing on sign and order of magnitude. Whether the inferred direction is
physically expected for this experiment remains open, and generalisation to
complex multi-scatterer scenes is left to future work.

The unifying hypothesis is therefore supported: time-lapse GPR can recover
subwavelength movement from the phase of migrated images at scales where the
amplitude image cannot --- subject to a migration-choice trade-off under noise,
and to independent confirmation of the field-data interpretation.
