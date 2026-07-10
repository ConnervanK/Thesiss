#import "../template.typ": *

#heading(level: 1, numbering: none, outlined: true)[Summary] <ch:summary>

#draftnote[Placeholder]

// #draftnote[this is a placeholder skeleton --- rewrite once @ch:discussion is
// filled in with concrete numbers, since the claims below should restate
// quantified results, not generic statements.]

// This thesis addressed the Research Question: _Can time-lapse ground-penetrating
// radar accurately track subwavelength movement --- achieving a form of
// super-resolution --- by analysing phase changes in migrated images rather than
// their amplitude?_

// The answer, across three progressively more complex experimental settings, is
// yes.

// *Hypothesis 1* (@ch:hyp1) established that amplitude differencing of migrated
// time-lapse images fails below a resolution floor of roughly half a wavelength
// --- but that the 2D phase-plane estimator, which fits a linear phase ramp to
// the cross-spectrum of a baseline/monitor image pair, recovers lateral,
// vertical, and diagonal displacements down to $1 \/ 32 lambda$, well below that
// floor. An optional local phase-gradient analysis (Hypothesis 1.5,
// @sec:hyp1-h15) explored instantaneous-phase imaging, spectral-line fitting,
// cross-phase spectrograms, and short-time Fourier transforms as complementary
// inference routes; these are consistent with the global WLS result wherever both
// were computed, but remain exploratory rather than fully validated.

// *Hypothesis 2* (@ch:hyp2) identified back-propagation with sign-bit
// time-reversal as the most noise-robust migration technique. Kirchhoff's
// delay-and-sum stacking manufactures false-coherent bands from pure Laplace
// noise (raising false-positive risk); Gazdag stays incoherent but adds
// significant speckle; back-propagation with sign-bit time-reversal suppresses
// impulsive noise by clamping every spike amplitude to $plus.minus 1$ while
// preserving all phase information --- stripping noise of the outsized amplitude
// that would otherwise let it dominate the back-propagated wavefield.

// *Hypothesis 3* (@ch:hyp3) tested generalisation beyond idealised synthetic
// data. The phase-plane pipeline was applied to real borehole GPR field data from
// a controlled fluid-injection experiment (38 profiles, 6 June 2016) and
// produced physically interpretable displacement estimates: a downward Push
// displacement of approximately $1.41 "m"$, near-zero Wait, and a partial Pull
// reversal, with a net residual of $approx +0.94 "m"$ along the borehole at the
// end of the experiment. The generalisation to complex multi-scatterer synthetic
// scenes is ongoing and not yet included.

// #draftnote[close with a final conclusion sentence stating, in one sentence,
// whether the Research Question is answered affirmatively by the totality of
// the three hypotheses, and noting the main caveat (the phase-plane fit is
// currently demonstrated quantitatively under noise for the lateral case only,
// and the complex synthetic model experiment is still pending).]
