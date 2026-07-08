# Project Wiki — Subwavelength Imaging in GPR Using Time-Lapse Data

Master directory for the Jupyter-notebook wiki. Each page below documents one notebook in `TimeLapse_Notebooks/` and maps it to the 30 June thesis draft (`TimeLapse_LaTeX/Thesis/30_June`, Typst).

## Notebook Directory

- [Resolution_Playground](Resolution_Playground.md) — Simulates two static subwavelength scatterers at decreasing separations (2λ→¹⁄₁₆λ) and migrates them with Kirchhoff, Gazdag, back-propagation, and least-squares migration to find the classical lateral resolution limit. Establishes the ≈λ/2 amplitude resolution floor of Chapter 3 §3.7.6 (`sec:meth-resolution`) that motivates the thesis's entire phase-based research question.

- [TimeLapse_Playground](TimeLapse_Playground.md) — Generates and migrates the lateral time-lapse study: a point scatterer displaced by 2λ→¹⁄₃₂λ, clean and with Laplace noise, across all three migration methods. Supplies all data and figures for Chapter 5 §5.1 (lateral movement and its amplitude detectability floor) and the noisy lateral migrations of Chapter 6 §6.4.

- [Vertical_TimeLapse_Playground](Vertical_TimeLapse_Playground.md) — The vertical counterpart: a single scatterer sinking by 1λ→¹⁄₃₂λ, migrated clean and noisy, including the only completed noisy sign-bit back-propagation sweep. Provides Chapter 5 §5.2 (vertical movement) and Chapter 6 §6.5, the strongest direct evidence for Hypothesis 2's sign-bit back-propagation claim.

- [Diagonal_TimeLapse_Playground](Diagonal_TimeLapse_Playground.md) — Simulates combined lateral+vertical movement along a fixed 2:1 diagonal (5 scenarios) and analyses the time-lapse difference along the motion direction. Provides Chapter 5 §5.3 (diagonal movement) and Chapter 6 §6.6, and holds the migrated data needed to close the thesis's flagged gap of §5.4.4 (the missing diagonal phase-plane fit).

- [Noise_Playground](Noise_Playground.md) — Fits Laplace vs Gaussian distributions to real field noise at each pipeline stage, saves the project-wide Laplace noise model, and migrates a pure-noise B-scan to expose each algorithm's false-positive behaviour. Supplies the noise-model table and the "Kirchhoff manufactures coherence, Gazdag stays speckled" finding of Chapter 6 §§6.1–6.2, plus the pending sign-bit pure-noise back-propagation test of §6.7.

- [TimeLapse_Cleaning](TimeLapse_Cleaning.md) — Benchmarks denoising strategies and robust estimators (median/wavelet/SVD/F-X cleaning, Huber/RANSAC/GCC fits, median-stacked apex-finding, background removal) against the noisy lateral dataset. Its winning GCC + PSR-fallback + shared-apex design became the noise-robust estimator behind Chapter 6 §6.4's noisy phase-plane results and substantiates the WLS weighting rationale of Chapter 3 §3.5.

- [TimeLapse_Processing](TimeLapse_Processing.md) — The central analysis notebook: applies the 2D WLS phase-plane fit (and GCC fallback) to every migrated dataset, compiles clean-vs-noisy error tables, ranks the three migration methods, and runs the local phase analyses (instantaneous phase, spectral lines, spectrograms, STFT). Generates Chapter 5 §5.4 (phase-plane validation) and §5.5 (Hypothesis 1.5), Chapter 6's quantitative noisy evidence and conclusion, all of Appendix A, and the pipeline described in Chapter 3 §3.7.5.

- [FluidFlow_Playground](FluidFlow_Playground.md) — Simulates a fluid front with a graded permittivity wetting zone (εr 80→1) advancing through a subwavelength fracture by 2λ→¹⁄₃₂λ, clean and noisy, across all three migrations. Realises the Abstract's "distributed fluid-front scenario" and Chapter 3 §3.6's material-change physics, and feeds the fluid-flow columns of Chapter 6's migration comparison as the synthetic bridge to Chapter 7's field experiment.

- [FieldData_Playground](FieldData_Playground.md) — Applies the full pipeline to 38 real borehole GPR profiles from the 2016 Ploemeur fluid-injection experiment: 11-step pre-processing, three migrations, time-lapse differencing, and cubed-weight WLS phase fitting under three tracking strategies. Provides everything in Chapter 7 §7.2 (Hypothesis 3's field-data evidence), including the stage-anchored displacement table (Push +1.41 m down, net +0.94 m), and sources the noise samples behind Chapter 6's Laplace model.

## Literature

- [LITERATURE.md](LITERATURE.md) — Categorized index of all research papers in `TimeLapse_LaTeX/Literature/`, organized by thesis-chapter topic folders, with a title, filename, and short description for each PDF (audited and reorganized 2026-07-08; confirmed duplicate downloads are quarantined in `_Duplicates_Review/`).

## Protocols & Standards

- [FIGURES_PROTOCOL.md](FIGURES_PROTOCOL.md) — The standardized figure pipeline: how notebooks auto-save 300 dpi PNGs into `TimeLapse_Figures/<Study>/<Technique>/` (Kirchhoff / Gazdag / Back-Propagation / Noisy / General, auto-classified from the figure title) via `helper_functions/figures.py`, and how the Typst thesis imports them through the two-stage `img()` router.
