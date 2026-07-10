#import "../template.typ": *

= Extended Per-Method Migration Figures for Hypothesis 1 <app:hyp1-methods>

@sec:hyp1-lateral, @sec:hyp1-vertical, and @sec:hyp1-diagonal each migrate
every displacement scenario with Kirchhoff, Gazdag, and back-propagation
migration, but the main chapter shows only the combined-comparison summary
figure per direction (e.g. @fig:tl-summary-amp), to keep the chapter focused on
the detectability result rather than on three near-identical sets of migrated
images. This appendix gives the individual migrated image and time-lapse
difference for each of the three migration algorithms, in each of the three
displacement directions, zoomed around the scatterer.

== Lateral Movement

#figure(
  subfigs(cols: 1,
    img("TL_011_Kirchhoff_Migration__All_8_Datasets____f_c15_GHz____aperture.png"),
    img("TL_013_Kirchhoff_Migration__TimeLapse_Differences_migrated__migrate.png"),
  ),
  caption: [Kirchhoff migration of the lateral time-lapse study ($f_c = 1.5 "GHz"$,
    aperture $= 40$), zoomed around the scatterers: (a) migrated image for all
    scenarios; (b) migrated-monitor-minus-migrated-baseline difference.],
) <fig:tl-kirchhoff>

#figure(
  subfigs(cols: 1,
    img("TL_015_Gazdag_Phase-Shift_Migration__All_7_Datasets____f_c15_GHz.png"),
    img("TL_017_Gazdag_Migration__TimeLapse_Differences_migrated__migrated_b.png"),
  ),
  caption: [Gazdag phase-shift migration of the lateral time-lapse study
    ($f_c = 1.5 "GHz"$), zoomed around the scatterers: (a) migrated image;
    (b) time-lapse difference.],
) <fig:tl-gazdag>

#figure(
  subfigs(cols: 1,
    img("TL_019_Back-Propagation_E__All_7_Datasets____focus_at_1906_ns_zoome.png", width: 90%),
    img("TL_021_Back-Propagation_Ez__All_7_Datasets____focus_at_1906_ns_zoom.png", width: 90%),
    img("TL_023_Back-Propagation__TimeLapse_Differences_Ez_Ez__Ez_baseline.png", width: 90%),
  ),
  caption: [Time-reversal back-propagation migration of the lateral time-lapse
    study, focused at $t = 19.06 "ns"$ and zoomed around the scatterers:
    (a) field magnitude $||bold(E)||$; (b) $E_z$; (c) $E_z$ time-lapse
    difference.],
) <fig:tl-backprop>

== Vertical Movement

#figure(
  subfigs(cols: 1,
    img("VTL_010_Kirchhoff_Migration_zoomed____f_c15_GHz____aperture40.png"),
    img("VTL_012_Kirchhoff_Migration__TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Kirchhoff migration of the vertical time-lapse study
    ($f_c = 1.5 "GHz"$, aperture $= 40$), zoomed around the scatterer:
    (a) migrated image for all scenarios; (b) time-lapse difference.],
) <fig:vtl-kirchhoff>

#figure(
  subfigs(cols: 1,
    img("VTL_014_Gazdag_Phase-Shift_Migration_zoomed____f_cf_c_GHz_GHz.png"),
    img("VTL_016_Gazdag_Migration__TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Gazdag phase-shift migration of the vertical time-lapse study
    ($f_c = 1.5 "GHz"$), zoomed around the scatterer: (a) migrated image;
    (b) time-lapse difference.],
) <fig:vtl-gazdag>

#figure(
  subfigs(cols: 1,
    img("VTL_018_Back-Propagation_E_zoomed____focus_at_1906_ns.png", width: 90%),
    img("VTL_020_Back-Propagation_Ez_zoomed____focus_at_1906_ns.png", width: 90%),
    img("VTL_022_Back-Propagation__TimeLapse_Differences_Ez_zoomed.png", width: 90%),
  ),
  caption: [Time-reversal back-propagation migration of the vertical time-lapse
    study, focused at $t = 19.06 "ns"$ and zoomed around the scatterer:
    (a) field magnitude $||bold(E)||$; (b) $E_z$; (c) $E_z$ time-lapse
    difference.],
) <fig:vtl-backprop>

== Diagonal Movement

#figure(
  subfigs(cols: 1,
    img("DTL_010_Kirchhoff_Migration_zoomed____f_c15_GHz____aperture40.png"),
    img("DTL_012_Kirchhoff_Migration__TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Kirchhoff migration of the diagonal time-lapse study
    ($f_c = 1.5 "GHz"$, aperture $= 40$), zoomed around the scatterer:
    (a) migrated image for all scenarios; (b) time-lapse difference.],
) <fig:dtl-kirchhoff>

#figure(
  subfigs(cols: 1,
    img("DTL_014_Gazdag_Phase-Shift_Migration_zoomed____f_c15_GHz.png"),
    img("DTL_016_Gazdag_Migration__TimeLapse_Differences_zoomed.png"),
  ),
  caption: [Gazdag phase-shift migration of the diagonal time-lapse study
    ($f_c = 1.5 "GHz"$), zoomed around the scatterer: (a) migrated image;
    (b) time-lapse difference.],
) <fig:dtl-gazdag>

#figure(
  subfigs(cols: 1,
    img("DTL_018_Back-Propagation_E_zoomed____focus_at_1906_ns.png", width: 90%),
    img("DTL_020_Back-Propagation_Ez_zoomed____focus_at_1906_ns.png", width: 90%),
    img("DTL_022_Back-Propagation__TimeLapse_Differences_Ez_zoomed.png", width: 90%),
  ),
  caption: [Time-reversal back-propagation migration of the diagonal time-lapse
    study, focused at $t = 19.06 "ns"$ and zoomed around the scatterer:
    (a) field magnitude $||bold(E)||$; (b) $E_z$; (c) $E_z$ time-lapse
    difference.],
) <fig:dtl-backprop>
