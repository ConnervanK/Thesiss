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

== Phase-Plane Shift-Estimation Function <supp:methodology-code>

The weighted least-squares phase-plane shift estimator applied throughout
this thesis (§4.5, "The 2D Phase-Plane Shift-Estimation Pipeline") is
implemented as the function `estimate_shift_2d` in
`helper_functions/WLS.py`, reproduced below in full for completeness. The
listing below additionally exposes the passband half-width and amplitude
threshold as named parameters (`band_factor`, `amp_thr`) and adds a
`kx_pos_only` option, generalising the fixed values used elsewhere in this
thesis and the one-sided-in-$kx$ masking needed for the field-data
processing pipeline; the defaults reproduce exactly the fixed $1.4$
passband and $10%$ threshold used throughout the rest of this document.

#block(
  fill: rgb("#F5F5F5"),
  stroke: 0.4pt + rgb("#CCCCCC"),
  inset: 0.8em,
  radius: 2pt,
  width: 100%,
)[
#text(size: 8pt)[
```python
def estimate_shift_2d(
    base, mon, dz_g, dx_g, kz_cent,
    kz_pos_only=False, kx_pos_only=False, force_dz_zero=False,
    weighted=True, band_factor=1.4, amp_thr=0.10,
    return_fit_points=False,
):
    """
    Estimate sub-pixel 2-D shifts between two images via weighted least-squares
    phase-plane fitting of their cross-spectrum.

    The cross-spectrum of a rigidly shifted pair satisfies:
        phi(kz, kx) = kz*dz + kx*dx + phi_0
    where phi_0 absorbs any constant phase offset (e.g. from a polarity flip).
    This function fits that plane to the measured cross-spectrum phase using the
    cross-spectrum magnitude as weights, so noisy / low-energy frequency bins
    contribute less to the solution.

    Input preparation (caller's responsibility):
        base and mon should be cropped to a sub-window centred on the region of
        interest before calling this function. The recommended approach is to
        compute the envelope (e.g. |hilbert(timelapse_diff)|) of the time-lapse
        difference image and use its peak location to define the crop centre.
        Fitting the phase plane over the full image mixes contributions from
        regions with no time-lapse signal, diluting the WLS weights and
        potentially biasing the shift estimate.

    Pre-processing:
        - A 2-D Tukey window (alpha=0.15) is applied before FFT to suppress
          spectral leakage from non-periodic edges.
        - Only bins inside a band |kz|, |kx| < band_factor * kz_cent AND with
          magnitude above amp_thr of the peak are included in the fit, limiting
          the WLS to the informative part of the spectrum and excluding the DC
          bin.

    Args:
        base (ndarray): Reference (baseline) image of shape (Nz, Nx).
        mon (ndarray): Monitor image of shape (Nz, Nx); same grid as base.
        dz_g (float): Grid spacing along the z (depth/time) axis [same unit as
            the desired dz_est output, e.g. m or ns].
        dx_g (float): Grid spacing along the x (lateral) axis [m].
        kz_cent (float): Central wavenumber of the dominant signal band [rad/m
            or rad/ns]. Used to define the spectral fitting window; should match
            the approximate centre of the wavelet's kz spectrum.
        kz_pos_only (bool): If True, restrict WLS to bins with KZ > 0. Use when
            base and mon are analytic signals (one-sided spectrum), so that the
            negative-frequency mirror image does not double-count phase.
        kx_pos_only (bool): If True, restrict WLS to bins with KX > 0. Useful when
            the +kx and -kx lobes carry opposing phase offsets (e.g. a hand-picked
            k-space region straddling both mirror lobes): fitting through both at
            once forces a single shared phi_0 onto data that actually needs phi_0
            and -phi_0, biasing the fitted plane.
        force_dz_zero (bool): If True, drop the KZ column and fit only
            (dx, phi_0), returning dz_est = 0.0. Use when vertical movement is
            known to be absent: avoids ill-conditioning that arises with a
            narrowband wavelet where the kz*dz and phi_0 terms become nearly
            indistinguishable on the one-sided spectrum.
        weighted (bool): If True (default), weight each masked bin by its
            cross-spectrum magnitude (WLS). If False, every masked bin gets
            equal weight (ordinary least squares) -- the mask/band selection
            is unchanged, only the fit itself stops discounting low-energy
            bins. Exists for the OLS-vs-WLS comparison in Chapter 5.4.
        band_factor (float): Half-width of the rectangular |kz|, |kx| passband, in
            multiples of kz_cent; bins with |kz| or |kx| >= band_factor * kz_cent
            are excluded from the fit regardless of amplitude. Default 1.4.
        amp_thr (float): Minimum cross-spectrum magnitude admitted to the fit, as
            a fraction of the peak magnitude |XS|.max(); bins below this threshold
            are treated as noise floor and excluded. Default 0.10.
        return_fit_points (bool): If True, also return a dict with the
            per-bin (kz, kx, phi, weight) arrays actually used in the fit,
            for diagnostic scatter plots.

    Returns:
        dz_est (float): Estimated shift along the z axis [units of dz_g].
            Always 0.0 when force_dz_zero=True.
        dx_est (float): Estimated shift along the x axis [units of dx_g].
        phi_0 (float): Fitted constant phase offset [rad].
        XS (ndarray): Complex cross-spectrum array of shape (Nz, Nx), useful
            for diagnostic plots.
        kz_ax (ndarray): 1-D kz wavenumber axis [rad / unit of dz_g].
        kx_ax (ndarray): 1-D kx wavenumber axis [rad/m].
        fit_points (dict, only if return_fit_points=True): {'kz', 'kx',
            'phi', 'weight'} -- the masked, per-bin arrays the fit was run
            on (weight is always the cross-spectrum magnitude, regardless of
            the `weighted` flag, so it can be used to visualise what WLS
            would/does emphasise).
    """
    Nz, Nx = base.shape

    kz_ax = np.fft.fftfreq(Nz, d=dz_g) * 2 * np.pi
    kx_ax = np.fft.fftfreq(Nx, d=dx_g) * 2 * np.pi
    KZ, KX = np.meshgrid(kz_ax, kx_ax, indexing='ij')

    taper = np.outer(tukey(Nz, alpha=0.15), tukey(Nx, alpha=0.15))

    base_fft = np.fft.fft2(base * taper)
    mon_fft  = np.fft.fft2(mon  * taper)
    XS       = base_fft * np.conj(mon_fft)

    w   = np.abs(XS)
    phi = np.angle(XS)

    band = ((np.abs(KZ) < band_factor * kz_cent) &
            (np.abs(KX) < band_factor * kz_cent))
    mask = (w > amp_thr * w.max()) & band & ((np.abs(KZ) + np.abs(KX)) > 0)
    if kz_pos_only:
        mask &= (KZ > 0)
    if kx_pos_only:
        mask &= (KX > 0)

    W = w[mask]
    Wfit = W if weighted else np.ones_like(W)
    if force_dz_zero:
        # 2-parameter fit: phi = kx*dx + phi_0  (dz = 0 by physics)
        A = np.column_stack([KX[mask], np.ones(mask.sum())])
        c = np.linalg.lstsq(A * Wfit[:, None], phi[mask] * Wfit, rcond=None)[0]
        dz_est, dx_est, phi_0 = 0.0, c[0], c[1]
    else:
        # 3-parameter fit: phi = kz*dz + kx*dx + phi_0
        A = np.column_stack([KZ[mask], KX[mask], np.ones(mask.sum())])
        c = np.linalg.lstsq(A * Wfit[:, None], phi[mask] * Wfit, rcond=None)[0]
        dz_est, dx_est, phi_0 = c[0], c[1], c[2]

    if return_fit_points:
        fit_points = dict(kz=KZ[mask], kx=KX[mask], phi=phi[mask], weight=W)
        return dz_est, dx_est, phi_0, XS, kz_ax, kx_ax, fit_points
    return dz_est, dx_est, phi_0, XS, kz_ax, kx_ax
```
]
]
