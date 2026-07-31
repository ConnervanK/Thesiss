#import "../template.typ": *

#heading(level: 1, numbering: none, outlined: true)[Abstract]

// Draft abstract --- confirm the headline figures against the final results
// tables before submission.

Ground-penetrating radar (GPR) resolves subsurface structure only down to a
wavelength-scale floor, so a displacement or material change smaller than a
fraction of a wavelength cannot be read from a migrated amplitude difference image ---
precisely the regime in which time-lapse monitoring of ground movement or
fluid migration must operate.

#linebreak()

This thesis develops a phase-based alternative:
a baseline and a monitor survey are migrated and compared in the
two-dimensional Fourier domain, where the Fourier shift theorem turns a
sub-wavelength translation into a linear phase ramp across the cross-spectrum.
A weighted least-squares fit of that ramp recovers sub-millimetre lateral and
vertical displacements, or a change in sub-wavelength material properties such as a fracture
filling with fluid.

#linebreak()

The method is tested through three hypotheses:

#linebreak()

// - On clean synthetic gprMax data
// --- across three migration algorithms (Kirchhoff, Gazdag phase-shift, and
// time-reversal back-propagation), displacement scales from $2 lambda$ down to
// $1 \/ 32 lambda$, and lateral, vertical, diagonal, and distributed fluid-front
// targets --- the phase-plane fit recovers displacement to a few tenths of a
// millimetre exactly where amplitude differencing has already collapsed. 
// - Under a
// realistic Laplace noise model the choice of migration algorithm matters:
// Kirchhoff is the most accurate but manufactures false coherent structure from
// noise alone, while back-propagation with sign-bit time-reversal is the more
// conservative choice. 
// - Finally, the full pipeline is applied to real borehole GPR
// field data from a fluid-injection experiment, producing stable,
// cross-technique-corroborated displacement estimates.
- On clean synthetic gprMax data --- across three migration algorithms (Kirchhoff, Gazdag phase-shift, and time-reversal back-propagation), displacement scales from $2 lambda$ down to $1/32 lambda$, and lateral, vertical, diagonal, and distributed fluid-front targets --- the phase-plane fit recovers displacement to a few tenths of a millimetre exactly where amplitude differencing has already collapsed. 

#linebreak()

- Under a realistic Laplace noise model the choice of migration algorithm matters: Kirchhoff is the most accurate but manufactures false coherent structures from noise alone, while back-propagation with sign-bit time-reversal is the more conservative choice. 

#linebreak()

- Finally, the full pipeline is applied to real borehole GPR field data from a fluid-injection experiment, producing stable, cross-technique-corroborated displacement estimates.