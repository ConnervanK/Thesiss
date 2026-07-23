#import "../template.typ": *

= Methodology: Supplementary Figures <supp:methodology>

This chapter collects figures moved out of the Theory and Methodology
chapter of the main thesis to keep that chapter's resolution-validation
section concise.

== Forward-Model Source Wavelet <supp:methodology-ricker>

#figure(
  img("RES_001_Ricker_Wavelet_f_c__15_GHz_t0__0943_ns.png"),
  caption: [The Ricker source wavelet used throughout this thesis
    ($f_c = 1.5 "GHz"$, $t_0 = 0.943 "ns"$), common to every gprMax forward
    model in Chapters 4 and 5 ("Hypothesis 1" and "Hypothesis 2") and the
    resolution validation in the Theory and Methodology chapter.],
) <fig:supp-ricker-wavelet>

== Per-Method Migration Images <supp:methodology-migration>

The main chapter shows only the combined signed-amplitude comparison and
PSF zoom across all three methods; this section gives the individual
migrated image for every separation scenario, for each method in turn,
zoomed around the true scatterer depth.

#figure(
  img("RES_009_Kirchhoff_Migration_zoomed____f_c15_GHz____aperture40.png", width: 90%),
  caption: [Kirchhoff migration of the resolution validation ($f_c = 1.5 "GHz"$,
    aperture $= 40$ traces), all separation scenarios, zoomed around the
    scatterer depth.],
) <fig:res-kirchhoff>

#figure(
  img("RES_011_Gazdag_Phase-Shift_Migration_zoomed____f_c15_GHz.png", width: 90%),
  caption: [Gazdag phase-shift migration of the resolution validation
    ($f_c = 1.5 "GHz"$), all separation scenarios, zoomed around the
    scatterer depth.],
) <fig:res-gazdag>

#figure(
  subfigs(cols: 1,
    img("RES_013_Back-Propagation_E_zoomed____focus_at_1906_ns.png"),
    img("RES_015_Back-Propagation_Ez_zoomed____focus_at_1906_ns.png"),
  ),
  caption: [Time-reversal back-propagation migration of the resolution
    validation, focused at $t = 19.06 "ns"$ and zoomed around the scatterer
    depth: (a) field magnitude $||bold(E)||$; (b) the $E_z$ component.],
) <fig:res-backprop>
