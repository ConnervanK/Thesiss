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
the "Interpretation (WLS vs. RANSAC)" section of Chapter 7, for a
representative selection of stage/technique/domain combinations -- the two
aggregate summaries covering every combination are shown in the main text,
§7.5; the five picks below are chosen to span the full
range of behaviour seen there, from near-perfect WLS/RANSAC agreement to the
most unstable amplitude-domain picks. Each panel overlays the painted or
gated cell selection on the cross-spectrum phase and amplitude, with the
$x$-/$z$-direction 1-D fits alongside (WLS top row, RANSAC bottom row,
rejected outliers marked with a red cross).

#figure(
  cimg("FD_ransac_vs_wls_phase_amp_fit_kirchhoff_bp_chasing_3_to_8.png"),
  caption: [Kirchhoff-BP, Chase-stage pair (3→8), amplitude-domain picking:
    the worst-case amplitude-domain pick referenced in the main text's
    "Interpretation (WLS vs. RANSAC)" section, §7.5, where RANSAC shifts the
    estimate the furthest of any pick tested ($#Dz$: $-1.74$ to
    $-1.78 "m"$; $#Dx$: $+0.02$ to $+0.03 "m"$).],
) <fig:supp-h3-ransac-kirchhoffbp-chase-bscan>

#figure(
  cimg("FD_ransac_vs_wls_phase_amp_fit_kirchhoff_bp_pulling_20_to_38.png"),
  caption: [Kirchhoff-BP, Pull-stage pair (20→38), amplitude-domain picking:
    an even more extreme case than the Chase-stage pick above -- WLS and
    RANSAC disagree on the *sign* of $#Dz$ ($+0.21 "m"$ versus
    $-0.17 "m"$), visible in the scattered, non-linear phase points that
    neither fit tracks cleanly.],
) <fig:supp-h3-ransac-kirchhoffbp-pull-bscan>

#figure(
  cimg("FD_ransac_vs_wls_phase_amp_fit_kspace_gazdag_chasing_3_to_8.png"),
  caption: [Gazdag, Chase-stage pair (3→8), wavenumber-domain picking: the
    direct k-space-picking counterpart to the amplitude-domain example
    figure in the main text, §7.5 (same stage and technique, opposite
    picking domain). WLS and
    RANSAC agree to four decimal places ($#Dz = -1.495 "m"$,
    $#Dx = +0.064 "m"$ for both) -- a hand-curated k-space pick needs
    essentially no further robustification.],
) <fig:supp-h3-ransac-gazdag-chase-kspace>

#figure(
  cimg("FD_ransac_vs_wls_phase_amp_fit_backprop_chasing_3_to_8.png"),
  caption: [Back-propagation, Chase-stage pair (3→8), amplitude-domain
    picking: the same instability seen for Kirchhoff-BP above, here for the
    third migration technique ($#Dx$: $+0.04$ to $+0.08 "m"$), showing the
    amplitude-domain gate's sensitivity is not specific to one migration
    algorithm.],
) <fig:supp-h3-ransac-backprop-chase-bscan>

#figure(
  cimg("FD_ransac_vs_wls_phase_amp_fit_kspace_gazdag_waiting_8_to_20.png"),
  caption: [Gazdag, Wait-stage pair (8→20), wavenumber-domain picking:
    another near-perfect-agreement case ($#Dz = +0.628 "m"$,
    $#Dx = +0.012 "m"$ for both WLS and RANSAC), confirming the tight
    k-space/amplitude-domain contrast above holds beyond the Chase stage.],
) <fig:supp-h3-ransac-gazdag-wait-kspace>
