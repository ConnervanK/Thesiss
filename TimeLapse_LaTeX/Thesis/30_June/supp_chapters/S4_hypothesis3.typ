#import "../template.typ": *

= Hypothesis 3: Supplementary Figures <supp:hyp3>

This chapter mirrors the section structure of Chapter 7, "Hypothesis 3:
Generalisation to Real Field Data", of the main thesis, in the same way as
@supp:hyp1 and @supp:hyp2. Chapter 7 was substantially condensed: automated
ROI-selection (rectangular-window and sliding-window scanning) and the
consecutive/fixed-baseline trajectory-chaining strategies were retired from
the thesis entirely, in favour of exclusively manual napari-based ROI
picking and stage-anchored chaining, so the figures that used to document
those retired methods (per-stage rectangular-window envelopes, the
sliding-window scan grid and displacement maps, and the per-technique
five-panel phase-plane diagnostics) have been removed rather than relocated
here. What remains below is genuinely supplementary detail for methods still
used in Chapter 7.

== Field Data Explanation <supp:hyp3-fielddata>

== Processing and Migration <supp:hyp3-fd-processing>

=== Cross-Profile Registration and Amplitude Normalisation <supp:hyp3-fd-crossprofile>

== Timelapse Differencing and Region of Influence Selection <supp:hyp3-fd-roi>

== WLS Cross-Spectrum Phase Plane Fit <supp:hyp3-fd-phaseplane>

=== Back-Propagation Cross-Check and Stage Displacements <supp:hyp3-fd-bp-corrected>

== Interpretation (WLS vs. RANSAC) <supp:hyp3-fd-interpretation>

=== WLS vs RANSAC <supp:hyp3-fd-ransac>

Per-pick inlier/outlier phase panels for the WLS-versus-RANSAC comparison in
the "Interpretation (WLS vs. RANSAC)" section of Chapter 7, for every stage
and migration technique, in both the wavenumber-domain and amplitude-domain
napari picking variants. Each panel overlays the painted or gated cell
selection on the cross-spectrum phase and marks the RANSAC inliers and the
rejected outliers.

#supp-note[These per-pick panels are produced by the field-data notebook
(`fielddata/output/roi_phase/`); the two aggregate WLS-versus-RANSAC summaries
are shown in the main text.]
