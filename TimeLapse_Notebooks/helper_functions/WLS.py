import numpy as np
import matplotlib.pyplot as plt
from scipy.signal         import hilbert as sp_hilbert
from scipy.signal.windows import tukey

def estimate_shift_2d(base, mon, dz_g, dx_g, kz_cent, kz_pos_only=False, force_dz_zero=False):
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
        - Only bins inside a band |kz|, |kx| < 1.4 * kz_cent AND with magnitude
          above 10 % of the peak are included in the fit, limiting the WLS to
          the informative part of the spectrum and excluding the DC bin.

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
        force_dz_zero (bool): If True, drop the KZ column and fit only
            (dx, phi_0), returning dz_est = 0.0. Use when vertical movement is
            known to be absent: avoids ill-conditioning that arises with a
            narrowband wavelet where the kz*dz and phi_0 terms become nearly
            indistinguishable on the one-sided spectrum.

    Returns:
        dz_est (float): Estimated shift along the z axis [units of dz_g].
            Always 0.0 when force_dz_zero=True.
        dx_est (float): Estimated shift along the x axis [units of dx_g].
        phi_0 (float): Fitted constant phase offset [rad].
        XS (ndarray): Complex cross-spectrum array of shape (Nz, Nx), useful
            for diagnostic plots.
        kz_ax (ndarray): 1-D kz wavenumber axis [rad / unit of dz_g].
        kx_ax (ndarray): 1-D kx wavenumber axis [rad/m].
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

    band = (np.abs(KZ) < 1.4 * kz_cent) & (np.abs(KX) < 1.4 * kz_cent)
    mask = (w > 0.10 * w.max()) & band & ((np.abs(KZ) + np.abs(KX)) > 0)
    if kz_pos_only:
        mask &= (KZ > 0)

    W = w[mask]
    if force_dz_zero:
        # 2-parameter fit: phi = kx*dx + phi_0  (dz = 0 by physics)
        A = np.column_stack([KX[mask], np.ones(mask.sum())])
        c = np.linalg.lstsq(A * W[:, None], phi[mask] * W, rcond=None)[0]
        return 0.0, c[0], c[1], XS, kz_ax, kx_ax
    else:
        # 3-parameter fit: phi = kz*dz + kx*dx + phi_0
        A = np.column_stack([KZ[mask], KX[mask], np.ones(mask.sum())])
        c = np.linalg.lstsq(A * W[:, None], phi[mask] * W, rcond=None)[0]
        return c[0], c[1], c[2], XS, kz_ax, kx_ax
