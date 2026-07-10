#import "../template.typ": *

= Extended Wavelength-Sweep Figures <app:extended-sweeps>

@sec:hyp1-phaseplane showed one representative scale from each of five
repeated parameter sweeps produced by `TimeLapse_Processing.ipynb`, to keep
the main narrative focused. This appendix gives the complete sweep for each
of those five analyses, across every wavelength scale that was simulated.
Figures here intentionally repeat the scale already highlighted in
@sec:hyp1-phaseplane, so that each full sweep is shown as a single,
self-contained reference figure.

== Full Sweep: Lateral Instantaneous-Phase Analysis <app:sec-instphase-horizontal>

@fig:app-instphase-horiz-full, @fig:app-instphase-horiz-cs-full, and
@fig:app-instphase-horiz-zoom-full give the complete seven-scale lateral sweep
underlying @sec:tlp-instphase-horizontal: the 2D instantaneous
phase-difference map, the phase cross-section at the scatterer depth
$z = 676 "mm"$, and the same cross-section zoomed to
$plus.minus 1\/2 lambda$, for every scale from $2 lambda$ to $1\/32 lambda$.

#figure(
  grid(columns: (1fr, 1fr, 1fr, 1fr), gutter: 0.5em,
    img("TLP_050_Instantaneous_Phase__Gazdag____Baseline_vs_2λ___Δx__2240_mm.png"),
    img("TLP_053_Instantaneous_Phase__Gazdag____Baseline_vs_1λ___Δx__1120_mm.png"),
    img("TLP_056_Instantaneous_Phase__Gazdag____Baseline_vs_½λ___Δx__560_mm.png"),
    img("TLP_059_Instantaneous_Phase__Gazdag____Baseline_vs_¼λ___Δx__280_mm.png"),
    img("TLP_062_Instantaneous_Phase__Gazdag____Baseline_vs_⅛λ___Δx__140_mm.png"),
    img("TLP_065_Instantaneous_Phase__Gazdag____Baseline_vs_¹₁₆λ___Δx__70_mm.png"),
    img("TLP_068_Instantaneous_Phase__Gazdag____Baseline_vs_¹₃₂λ___Δx__35_mm.png"),
  ),
  caption: [Full seven-scale sweep: 2D instantaneous phase-difference map,
    Gazdag-migrated, lateral displacement from $2 lambda$ (top-left) down to
    $1\/32 lambda$ (bottom row): (a) $2 lambda$; (b) $1 lambda$;
    (c) $1\/2 lambda$; (d) $1\/4 lambda$; (e) $1\/8 lambda$;
    (f) $1\/16 lambda$; (g) $1\/32 lambda$.],
) <fig:app-instphase-horiz-full>

#figure(
  grid(columns: (1fr, 1fr, 1fr, 1fr), gutter: 0.5em,
    img("TLP_051_Phase_cross-section__z__676_mm____Gazdag____2λ__Δx__2240_mm.png"),
    img("TLP_054_Phase_cross-section__z__676_mm____Gazdag____1λ__Δx__1120_mm.png"),
    img("TLP_057_Phase_cross-section__z__676_mm____Gazdag____½λ__Δx__560_mm.png"),
    img("TLP_060_Phase_cross-section__z__676_mm____Gazdag____¼λ__Δx__280_mm.png"),
    img("TLP_063_Phase_cross-section__z__676_mm____Gazdag____⅛λ__Δx__140_mm.png"),
    img("TLP_066_Phase_cross-section__z__676_mm____Gazdag____¹₁₆λ__Δx__70_mm.png"),
    img("TLP_069_Phase_cross-section__z__676_mm____Gazdag____¹₃₂λ__Δx__35_mm.png"),
  ),
  caption: [Full seven-scale sweep: phase cross-section at $z = 676 "mm"$,
    Gazdag-migrated, lateral displacement from $2 lambda$ down to
    $1\/32 lambda$: (a) $2 lambda$; (b) $1 lambda$; (c) $1\/2 lambda$;
    (d) $1\/4 lambda$; (e) $1\/8 lambda$; (f) $1\/16 lambda$;
    (g) $1\/32 lambda$.],
) <fig:app-instphase-horiz-cs-full>

#figure(
  grid(columns: (1fr, 1fr, 1fr, 1fr), gutter: 0.5em,
    img("TLP_052_Phase_cross-section__zoomed__½λ_around_scatterers.png"),
    img("TLP_055_Phase_cross-section__zoomed__½λ_around_scatterers.png"),
    img("TLP_058_Phase_cross-section__zoomed__½λ_around_scatterers.png"),
    img("TLP_061_Phase_cross-section__zoomed__½λ_around_scatterers.png"),
    img("TLP_064_Phase_cross-section__zoomed__½λ_around_scatterers.png"),
    img("TLP_067_Phase_cross-section__zoomed__½λ_around_scatterers.png"),
    img("TLP_070_Phase_cross-section__zoomed__½λ_around_scatterers.png"),
  ),
  caption: [Full seven-scale sweep: phase cross-section zoomed to
    $plus.minus 1\/2 lambda$ around the scatterers, Gazdag-migrated, lateral
    displacement from $2 lambda$ down to $1\/32 lambda$: (a) $2 lambda$;
    (b) $1 lambda$; (c) $1\/2 lambda$; (d) $1\/4 lambda$; (e) $1\/8 lambda$;
    (f) $1\/16 lambda$; (g) $1\/32 lambda$.],
) <fig:app-instphase-horiz-zoom-full>

== Full Sweep: Vertical Instantaneous-Phase Analysis <app:sec-instphase-vertical>

@fig:app-instphase-vert-full, @fig:app-instphase-vert-cs-full, and
@fig:app-instphase-vert-zoom-full give the complete six-scale vertical sweep
underlying @sec:tlp-instphase-vertical, from $1 lambda$ to $1\/32 lambda$.

#figure(
  grid(columns: (1fr, 1fr, 1fr), gutter: 0.5em,
    img("TLP_072_VerticalTimeLapse__Instantaneous_Phase__Gazdag____Baseline_v.png"),
    img("TLP_075_VerticalTimeLapse__Instantaneous_Phase__Gazdag____Baseline_v.png"),
    img("TLP_078_VerticalTimeLapse__Instantaneous_Phase__Gazdag____Baseline_v.png"),
    img("TLP_081_VerticalTimeLapse__Instantaneous_Phase__Gazdag____Baseline_v.png"),
    img("TLP_084_VerticalTimeLapse__Instantaneous_Phase__Gazdag____Baseline_v.png"),
    img("TLP_087_VerticalTimeLapse__Instantaneous_Phase__Gazdag____Baseline_v.png"),
  ),
  caption: [Full six-scale sweep: 2D instantaneous phase-difference map,
    Gazdag-migrated, vertical displacement from $1 lambda$ down to
    $1\/32 lambda$: (a) $1 lambda$; (b) $1\/2 lambda$; (c) $1\/4 lambda$;
    (d) $1\/8 lambda$; (e) $1\/16 lambda$; (f) $1\/32 lambda$.],
) <fig:app-instphase-vert-full>

#figure(
  grid(columns: (1fr, 1fr, 1fr), gutter: 0.5em,
    img("TLP_073_Phase_cross-section__x__2000_mm____Gazdag____1λ__Δz__1120_mm.png"),
    img("TLP_076_Phase_cross-section__x__2000_mm____Gazdag____½λ__Δz__560_mm.png"),
    img("TLP_079_Phase_cross-section__x__2000_mm____Gazdag____¼λ__Δz__280_mm.png"),
    img("TLP_082_Phase_cross-section__x__2000_mm____Gazdag____⅛λ__Δz__140_mm.png"),
    img("TLP_085_Phase_cross-section__x__2000_mm____Gazdag____¹₁₆λ__Δz__70_mm.png"),
    img("TLP_088_Phase_cross-section__x__2000_mm____Gazdag____¹₃₂λ__Δz__35_mm.png"),
  ),
  caption: [Full six-scale sweep: phase cross-section at $x = 2000 "mm"$,
    Gazdag-migrated, vertical displacement from $1 lambda$ down to
    $1\/32 lambda$: (a) $1 lambda$; (b) $1\/2 lambda$; (c) $1\/4 lambda$;
    (d) $1\/8 lambda$; (e) $1\/16 lambda$; (f) $1\/32 lambda$.],
) <fig:app-instphase-vert-cs-full>

#figure(
  grid(columns: (1fr, 1fr, 1fr), gutter: 0.5em,
    img("TLP_074_Phase_cross-section__zoomed__½λ_around_scatterer_depths.png"),
    img("TLP_077_Phase_cross-section__zoomed__½λ_around_scatterer_depths.png"),
    img("TLP_080_Phase_cross-section__zoomed__½λ_around_scatterer_depths.png"),
    img("TLP_083_Phase_cross-section__zoomed__½λ_around_scatterer_depths.png"),
    img("TLP_086_Phase_cross-section__zoomed__½λ_around_scatterer_depths.png"),
    img("TLP_089_Phase_cross-section__zoomed__½λ_around_scatterer_depths.png"),
  ),
  caption: [Full six-scale sweep: phase cross-section zoomed to
    $plus.minus 1\/2 lambda$ around the scatterer depths, Gazdag-migrated,
    vertical displacement from $1 lambda$ down to $1\/32 lambda$:
    (a) $1 lambda$; (b) $1\/2 lambda$; (c) $1\/4 lambda$; (d) $1\/8 lambda$;
    (e) $1\/16 lambda$; (f) $1\/32 lambda$.],
) <fig:app-instphase-vert-zoom-full>

== Full Sweep: Spectral-Line Phase Analysis <app:sec-spectral-line>

#draftnote[*Gap:* @sec:tlp-spectral-line cites
a full seven-scale CLSSA spectral-line sweep here, but the figures no longer
exist in `TimeLapse_Figures/` under any filename -- `TimeLapse_Processing.ipynb`
Section 7 and 7b (the cells that produce the Gazdag-migrated and raw-unmigrated
CLSSA decompositions respectively) are currently commented out in their
entirety, so no run of the notebook since they were disabled has regenerated
this output. `_reorg_mapping.txt` confirms the files existed as of the
2026-07-08 figure-store reorganisation (as `TLP_057`--`070`), so the analysis
itself is not lost, only its cells need re-enabling and the notebook
re-running from that point onward (which will renumber every subsequent
figure in the store). Restore this section once that is done.]

== Full Sweep: Cross-Phase Spectrograms <app:sec-spectrogram>

#draftnote[*Gap:* same issue as @app:sec-spectral-line above --
@sec:tlp-spectrogram cites a full seven-scale
cross-phase spectrogram sweep here, but `TimeLapse_Processing.ipynb` Section 8
(the cell that produces it) is currently commented out and the figures no
longer exist under any filename (previously `TLP_071`--`077`, confirmed via
`_reorg_mapping.txt`). Re-enable and re-run that section, then restore this
figure.]

== Full Sweep: Localised Fourier Shift via STFT <app:sec-stft>

#draftnote[*Gap:* same issue again -- @sec:tlp-stft
cites a full seven-output localised-STFT sweep here, but
`TimeLapse_Processing.ipynb` Section 9 (the cell that produces it) is
currently commented out and the figures no longer exist under any filename
(previously `TLP_078`--`084`, confirmed via `_reorg_mapping.txt`). Note this
was already the most uncertain sweep in this appendix even before the figures
disappeared -- the seven outputs shared an identical auto-generated filename
stem and did not encode which displacement scale each panel corresponded to,
so when Section 9 is re-enabled, confirm the per-panel scale labelling
against the notebook's cell order before trusting the captions. Re-enable and
re-run that section, then restore this figure.]
