# Literature folder structure

Organized by thesis chapter → subject → topic, matching the section numbers in
`Thesis/27_June/main.tex` (and the current `30_June` Typst chapters). A paper
lives in exactly **one** folder even when it's relevant to more than one
chapter — cite it from wherever you need in the `.bib` file regardless of
where the PDF physically sits.

## Chapter mapping

- `1_Introduction/` — motivation / framing papers.
- `2_Review_of_Literature/` — the literature-review chapter, split into its
  four sections (2-1 EM Theory, 2-2 Resolution/Fresnel Zone, 2-3 GPR
  Processing, 2-4 Phase-Based Interpretation). This is the largest and most
  subdivided folder since it's where nearly all background reading lives.
- `3_Theoretical_Background/` — only holds papers specific to sections that
  aren't already covered by a 2-3/2-4 topic (Fourier Shift Theorem,
  Cross-Spectrum, Space-Wavenumber Duality, Time-Frequency Perspective). See
  **Cross-references** below for 3-1, 3-4, 3-5.
- `4_Methodology/` — gprMax modelling and validation references specific to
  your own pipeline (as opposed to `2-3/Simulation_Tooling_gprMax`, which
  holds gprMax literature in general).
- `5-7_Hypotheses_Field_and_TimeLapse_Studies/` — precedent field/time-lapse
  GPR studies that motivate or compare against Hypotheses 1-3.
- `8_Discussion/` — currently empty; discussion mostly re-cites earlier
  chapters rather than needing new sources.
- `Thesis_Writing_Resources/` — writing/stats guides, not thesis content.

## Cross-references (same folder serves multiple sections)

- **Section 3-1** (Kirchhoff / Gazdag / Back-Propagation fundamentals) →
  papers are filed under
  `2_Review_of_Literature/2-3_GPR_Processing/Migration_Algorithms/{Kirchhoff, Gazdag_Phase-Shift, Back-Propagation_Time-Reversal}`.
- **Section 3-4** (WLS Plane Fitting) and **3-5** (Geometric vs Material
  Decoupling) → foundational math is in
  `2_Review_of_Literature/2-4_Phase-Based_Interpretation/WLS_Phase-Difference_Fitting`.
- **Section 4-1** (gprMax forward modelling) → general gprMax
  literature/theory is in
  `2_Review_of_Literature/2-1_EM_Theory_and_GPR_Physics/Forward_Modelling_Principles`
  and `2-3_GPR_Processing/Simulation_Tooling_gprMax`; `4_Methodology/4-1_gprMax_Models`
  is reserved for anything specific to your own model setup only.

## `_Duplicates_Review/`

19 PDFs that turned out to be duplicate downloads of a paper already filed
elsewhere (same title/content, different filename — e.g. a DOI-named
ScienceDirect download next to a hand-renamed copy). Nothing was deleted;
skim this folder and delete once you've confirmed you don't need the second
copy (e.g. for annotations).
