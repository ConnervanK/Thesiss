#import "../template.typ": *

= Hypothesis 3: Supplementary Figures <supp:hyp3>

This chapter collects supplementary detail for Chapter 7, "Hypothesis 3:
Generalisation to Real Field Data", of the main thesis, in the same way as
@supp:hyp1 and @supp:hyp2 -- though, unlike those two chapters, it no longer
mirrors Chapter 7's full section structure one-for-one. Chapter 7 was
substantially condensed: automated ROI-selection (rectangular-window and
sliding-window scanning) and the consecutive/fixed-baseline
trajectory-chaining strategies were retired from the thesis entirely, in
favour of exclusively manual napari-based ROI picking and stage-anchored
chaining, so the figures that used to document those retired methods
(per-stage rectangular-window envelopes, the sliding-window scan grid and
displacement maps, and the per-technique five-panel phase-plane diagnostics)
have been removed rather than relocated here, along with the now-empty
section headings that used to hold them. What remains below -- the
per-pick WLS-versus-RANSAC diagnostics -- is the only genuinely
supplementary detail left for methods still used in Chapter 7; every other
aspect of the field-data workflow is fully covered in the main chapter
text with nothing held back here.

== Interpretation: WLS vs. RANSAC <supp:hyp3-fd-interpretation>

Per-pick inlier/outlier phase panels for the WLS-versus-RANSAC comparison in
the "Interpretation (WLS vs. RANSAC)" section of Chapter 7, for every stage
and migration technique, in both the wavenumber-domain and amplitude-domain
napari picking variants. Each panel overlays the painted or gated cell
selection on the cross-spectrum phase and marks the RANSAC inliers and the
rejected outliers.

#supp-note[These per-pick panels are produced by the field-data notebook
(`fielddata/output/roi_phase/`); the two aggregate WLS-versus-RANSAC summaries
are shown in the main text. #strong[Not yet embedded here] -- see the
accompanying chat message for why (unresolved filename/stage-label
mismatches found in the candidate panel files) and what's needed to finish
this section.]
