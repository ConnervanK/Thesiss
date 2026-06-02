"""
phase_decomposition.py
======================
Spectral and phase decomposition of a single GPR trace, following the
methodology of Castagna et al. (2016) "Phase decomposition", Interpretation.

Key equations
-------------
CWT coefficients (complex Morlet):
    W(f, t) — complex array, shape (n_freqs, n_time)
    A(f, t) = |W(f, t)|          instantaneous amplitude
    θ(f, t) = angle(W(f, t))     instantaneous phase

Phase gather (Castagna 2016, eq. 3):
    S'(θ, t) = ∫ A(f,t) cos(θ(f,t) − θ) df

Expanded (cos subtraction identity) for efficient vectorised computation:
    S'(θ, t) = cos(θ)·R(t) + sin(θ)·I(t)
    where  R(t) = ∫ Re[W] df,   I(t) = ∫ Im[W] df   (trapezoid over freq)

Edge-effect handling
--------------------
At low frequencies the complex Morlet wavelet has a large time-support.
Near the trace edges this produces the "cone of influence" (COI) — spurious
oscillations in the amplitude and phase spectra.  Two mitigations are applied:

1.  Reflect-padding: the trace is mirrored at both ends by `n_pad` samples
    before the CWT, then trimmed back to the original length.  This pushes
    the COI boundary outside the signal, eliminating edge noise.
    n_pad is set to 2 × σ_wavelet(f_min), the e-folding half-width of the
    Gaussian envelope at the lowest requested frequency.

2.  COI boundary overlay: dashed white lines mark the time from each edge
    at which wavelet support falls to 1/e.  Data inside the COI should be
    interpreted with caution.

Usage
-----
    from phase_decomposition import phase_decomposition, plot_decomposition
    A, theta_deg, pg, freqs, t_ns, theta_ang = phase_decomposition(trace, dt)
    plot_decomposition(A, theta_deg, pg, freqs, t_ns, theta_ang)

    # or run standalone for a demo:
    python phase_decomposition.py
"""

import numpy as np
import matplotlib.pyplot as plt
import pywt


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def phase_decomposition(
    trace,
    dt,
    f_min_GHz: float = 0.5,
    f_max_GHz: float = 5.0,
    n_scales: int = 60,
    n_theta: int = 181,
    wavelet: str = "cmor1.5-1.0",
    reflect_pad: bool = True,
):
    """
    Spectral and phase decomposition of a 1-D GPR trace.

    Parameters
    ----------
    trace : 1-D ndarray
        Amplitude values of the trace (arbitrary units).
    dt : float
        Sampling interval in **nanoseconds**.
    f_min_GHz, f_max_GHz : float
        Frequency band for the CWT and phase gather [GHz].
        Default lower limit 0.5 GHz; a 1.5 GHz Ricker has negligible energy
        below ~0.5 GHz, so lower frequencies would be noise-only.
    n_scales : int
        Number of CWT scales (log-spaced across the frequency band).
    n_theta : int
        Number of phase angles in the gather (default 181 → −180° … +180°).
    wavelet : str
        PyWavelets complex wavelet name.  'cmor{B}-{C}' where B is the
        bandwidth and C the centre frequency (cycles/sample).
    reflect_pad : bool
        If True (default), reflect-pad the trace before CWT to suppress
        edge artefacts (cone-of-influence noise) at low frequencies.

    Returns
    -------
    A : ndarray (n_scales, n_time)
        Instantaneous amplitude A(f, t).
    theta_2d_deg : ndarray (n_scales, n_time)
        Instantaneous phase θ(f, t) in degrees.
    phase_gather : ndarray (n_theta, n_time)
        Phase gather S'(θ, t).
    freqs_GHz : ndarray (n_scales,)
        Frequency axis in GHz.
    time_ns : ndarray (n_time,)
        Two-way time axis in nanoseconds.
    theta_deg : ndarray (n_theta,)
        Phase axis in degrees for the phase gather.
    """
    trace = np.asarray(trace, dtype=float)
    n_time = len(trace)
    time_ns = np.arange(n_time) * dt

    # --- wavelet parameters -----------------------------------------------
    cwavelet = pywt.ContinuousWavelet(wavelet)
    f_c = cwavelet.center_frequency          # cycles/sample (typically 1.0)
    B   = float(wavelet.split("-")[0].replace("cmor", ""))  # bandwidth

    # --- scales: f_GHz = f_c / (scale * dt)  →  scale = f_c / (f_GHz * dt)
    freq_axis = np.geomspace(f_min_GHz, f_max_GHz, n_scales)  # ascending
    scales = f_c / (freq_axis * dt)          # descending (large → small)

    # --- reflect-pad to suppress COI edge artefacts ----------------------
    if reflect_pad:
        # e-folding half-width of Gaussian envelope at the lowest frequency:
        #   σ(s) = sqrt(B/2) * s  [samples]
        n_pad = int(np.ceil(2.0 * np.sqrt(B / 2.0) * scales.max()))
        n_pad = min(n_pad, n_time)            # cap at trace length
        trace_work = np.pad(trace, n_pad, mode="reflect")
    else:
        n_pad = 0
        trace_work = trace

    # --- CWT (complex Morlet → complex coefficients) ---------------------
    coeffs, freqs_GHz = pywt.cwt(
        trace_work, scales, wavelet, sampling_period=dt
    )
    # coeffs: complex (n_scales, n_time_padded); freqs_GHz in GHz

    # trim back to original signal length
    if n_pad > 0:
        coeffs = coeffs[:, n_pad : n_pad + n_time]

    # --- Amplitude and phase ----------------------------------------------
    A             = np.abs(coeffs)                    # (n_scales, n_time)
    theta_2d_rad  = np.angle(coeffs)                  # radians
    theta_2d_deg  = np.rad2deg(theta_2d_rad)          # degrees

    # --- Phase gather (vectorised) ----------------------------------------
    # Pre-integrate Re and Im parts over frequency using the trapezoidal rule.
    R = np.trapezoid(np.real(coeffs), freqs_GHz, axis=0)   # (n_time,)
    I = np.trapezoid(np.imag(coeffs), freqs_GHz, axis=0)   # (n_time,)

    theta_deg = np.linspace(-180.0, 180.0, n_theta)
    theta_rad = np.deg2rad(theta_deg)

    # S'(θ, t) = cos(θ)·R(t) + sin(θ)·I(t)  →  shape (n_theta, n_time)
    phase_gather = (
        np.outer(np.cos(theta_rad), R)
        + np.outer(np.sin(theta_rad), I)
    )

    return A, theta_2d_deg, phase_gather, freqs_GHz, time_ns, theta_deg


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

def _coi_boundary(freqs_GHz, time_ns, B, f_c=1.0):
    """
    Return the COI boundary time [ns] from each edge for each frequency.
    σ(f) = sqrt(B/2) * f_c / f_GHz  [ns]
    """
    coi = np.sqrt(B / 2.0) * f_c / freqs_GHz     # ns from each edge
    # Clip so lines stay inside the time axis
    t_max = time_ns[-1]
    return np.clip(coi, 0, t_max / 2.0)


def plot_decomposition(
    A,
    theta_2d_deg,
    phase_gather,
    freqs_GHz,
    time_ns,
    theta_deg,
    title: str = "Spectral & Phase Decomposition",
    db_clip: float = 40.0,
    show_coi: bool = True,
    wavelet: str = "cmor1.5-1.0",
):
    """
    Three-panel figure:
      1. Amplitude Spectrum  A(f, t)  [jet, dB-normalised]
      2. Phase Spectrum      θ(f, t)  [RdBu, ±180°]
      3. Phase Gather        S'(θ,t)  [RdBu, symmetric clip]

    Parameters
    ----------
    show_coi : bool
        Overlay dashed white lines marking the cone-of-influence boundary.
    db_clip : float
        Dynamic range in dB for the amplitude display.
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 6), sharey=True)
    fig.suptitle(title, fontsize=13, fontweight="bold")

    t0, t1 = time_ns[0], time_ns[-1]
    # freqs_GHz[0]=f_min, freqs_GHz[-1]=f_max  (ascending output from pywt)
    f_left, f_right = freqs_GHz[0], freqs_GHz[-1]

    # imshow: A.T has shape (n_time, n_scales) → time on y-axis, freq on x-axis
    extent_spec = [f_left, f_right, t1, t0]   # [x_left, x_right, y_bottom, y_top]

    # ── Plot 1: Amplitude Spectrum (dB) ──────────────────────────────────
    ax = axes[0]
    A_norm = A / (A.max() + 1e-30)
    A_dB = 20.0 * np.log10(A_norm + 1e-30)
    A_dB = np.clip(A_dB, -db_clip, 0.0)

    im0 = ax.imshow(
        A_dB.T,
        aspect="auto",
        origin="upper",
        extent=extent_spec,
        cmap="jet",
        vmin=-db_clip,
        vmax=0.0,
        interpolation="bilinear",
    )
    cb0 = fig.colorbar(im0, ax=ax, fraction=0.046, pad=0.04)
    cb0.set_label("Amplitude [dB]", fontsize=9)
    ax.set_xlabel("Frequency [GHz]", fontsize=10)
    ax.set_ylabel("TWT [ns]", fontsize=10)
    ax.set_title("Amplitude Spectrum", fontsize=11)

    # ── Plot 2: Phase Spectrum (degrees) ─────────────────────────────────
    ax = axes[1]
    im1 = ax.imshow(
        theta_2d_deg.T,
        aspect="auto",
        origin="upper",
        extent=extent_spec,
        cmap="RdBu",
        vmin=-180.0,
        vmax=180.0,
        interpolation="nearest",
    )
    cb1 = fig.colorbar(im1, ax=ax, fraction=0.046, pad=0.04)
    cb1.set_label("Phase [°]", fontsize=9)
    ax.set_xlabel("Frequency [GHz]", fontsize=10)
    ax.set_title("Phase Spectrum", fontsize=11)

    # ── COI overlay on amplitude and phase panels ─────────────────────────
    if show_coi:
        B = float(wavelet.split("-")[0].replace("cmor", ""))
        f_c = 1.0   # centre freq = 1.0 cycles/sample for 'cmor?-1.0'
        coi_t = _coi_boundary(freqs_GHz, time_ns, B, f_c)
        for ax_sp in (axes[0], axes[1]):
            ax_sp.plot(freqs_GHz, coi_t,              "w--", lw=1.0, alpha=0.75,
                       label="COI boundary")
            ax_sp.plot(freqs_GHz, time_ns[-1] - coi_t, "w--", lw=1.0, alpha=0.75)

    # ── Plot 3: Phase Gather ─────────────────────────────────────────────
    ax = axes[2]
    pg_clip = np.percentile(np.abs(phase_gather), 99)
    im2 = ax.imshow(
        phase_gather.T,    # (n_time, n_theta)
        aspect="auto",
        origin="upper",
        extent=[theta_deg[0], theta_deg[-1], t1, t0],
        cmap="RdBu",
        vmin=-pg_clip,
        vmax=pg_clip,
        interpolation="bilinear",
    )
    cb2 = fig.colorbar(im2, ax=ax, fraction=0.046, pad=0.04)
    cb2.set_label("Amplitude", fontsize=9)
    ax.set_xlabel("Phase [°]", fontsize=10)
    ax.set_title("Phase Gather  S'(θ, t)", fontsize=11)

    # Mark ±90° — most sensitive to thin-bed impedance variations
    for ph in (-90, 90):
        ax.axvline(ph, color="lime", lw=1.0, ls="--", alpha=0.85)
    ax.axvline(0, color="white", lw=0.6, ls=":", alpha=0.5)
    ax.set_xticks([-180, -90, 0, 90, 180])

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# CLSSA — Constrained Least-Squares Spectral Analysis
# ---------------------------------------------------------------------------

def clssa_phase_decomposition(
    trace,
    dt,
    f_min_GHz: float = 0.1,
    f_max_GHz: float = 4.0,
    n_freqs: int = 60,
    n_theta: int = 181,
    win_ns: float = 5.0,
    alpha: float = 1e-2,
):
    """
    CLSSA + phase decomposition of a 1-D GPR trace.

    Unlike CWT, the tapering window is applied to the BASIS FUNCTIONS
    (sinusoids inside G), not to the data.  The inversion therefore
    recovers the spectrum of the RAW trace, not the windowed trace.

    Forward model for each time centre t_c:
        d  =  G m  +  noise
        G[j,k]  =  h[j] * exp(i 2pi f_k (t_j - t_c))   [tapered sinusoid]
        h[j]    =  Hanning window

    Tikhonov solution (normalised columns so alpha is scale-invariant):
        m  =  (G_n^H G_n  +  alpha I)^{-1} G_n^H d  /  col_norms
    equivalently the Woodbury form from the problem statement:
        m  =  G_n^H (G_n G_n^H  +  alpha I)^{-1} d  /  col_norms

    Parameters
    ----------
    trace     : 1-D ndarray
    dt        : sampling interval [ns]
    f_min_GHz, f_max_GHz : frequency band [GHz]
    n_freqs   : number of discrete frequency components
    n_theta   : number of phase angles in the phase gather
    win_ns    : full sliding-window length [ns]
    alpha     : Tikhonov regularisation (fraction of normalised Gram diagonal)

    Returns
    -------
    A, theta_2d_deg, phase_gather, freqs_GHz, time_ns, theta_deg
        Same shapes/units as phase_decomposition().
    """
    trace = np.asarray(trace, dtype=float)
    n = len(trace)
    time_ns = np.arange(n) * dt

    freqs_GHz = np.linspace(f_min_GHz, f_max_GHz, n_freqs)

    half = int(round(win_ns / dt / 2))
    n_win = 2 * half + 1

    # Reflect-pad the trace so every output step uses a full window.
    # This eliminates slow per-sample edge solves and pushes boundary
    # artefacts outside the signal — same principle as CWT reflect_pad.
    trace_pad = np.pad(trace, half, mode="reflect")   # length n + 2*half

    # Hanning window — applied to kernel sinusoids (NOT to the data)
    h = np.hanning(n_win)                              # (n_win,)
    t_loc = (np.arange(n_win) - half) * dt             # centred at 0

    # Central kernel G_c: (n_win, n_freqs)
    G_c = h[:, None] * np.exp(
        1j * 2 * np.pi * freqs_GHz[None, :] * t_loc[:, None]
    )
    # Normalise columns so alpha has consistent meaning across window lengths
    col_norms = np.linalg.norm(G_c, axis=0)            # (n_freqs,)
    G_n = G_c / col_norms[None, :]

    # Regularised Gram matrix (precomputed once — same for all steps)
    A_c   = G_n.conj().T @ G_n + alpha * np.eye(n_freqs)
    A_inv = np.linalg.inv(A_c)
    W     = G_n.conj()                                 # (n_win, n_freqs)

    # ── Vectorised batch over all n steps ────────────────────────────────────
    # sliding_window_view over padded trace: shape (n, n_win)
    D   = np.lib.stride_tricks.sliding_window_view(trace_pad, n_win)  # (n, n_win)
    RHS = D @ W                                        # (n, n_freqs)
    M   = (A_inv @ RHS.T)                             # (n_freqs, n)
    C   = M / col_norms[:, None]

    # ── Amplitude, phase, phase gather ───────────────────────────────────────
    A_out        = np.abs(C)
    theta_2d_deg = np.rad2deg(np.angle(C))

    R = np.trapezoid(np.real(C), freqs_GHz, axis=0)
    I = np.trapezoid(np.imag(C), freqs_GHz, axis=0)

    theta_deg = np.linspace(-180.0, 180.0, n_theta)
    theta_rad = np.deg2rad(theta_deg)
    phase_gather = (
        np.outer(np.cos(theta_rad), R)
        + np.outer(np.sin(theta_rad), I)
    )

    return A_out, theta_2d_deg, phase_gather, freqs_GHz, time_ns, theta_deg


def plot_clssa_decomposition(
    A,
    theta_2d_deg,
    phase_gather,
    freqs_GHz,
    time_ns,
    theta_deg,
    title: str = "CLSSA Spectral & Phase Decomposition",
    db_clip: float = 40.0,
    highlight_neg90: bool = True,
):
    """
    Three-panel figure matching the CWT plot_decomposition layout but labelled
    CLSSA, with optional -90° anomaly highlight in the phase gather.

    Parameters
    ----------
    highlight_neg90 : bool
        If True, add a magenta band at theta = -90 deg in the phase gather.
        A dominant -90 deg phase is diagnostic of a sub-resolution low-impedance
        thin layer (see Castagna et al. 2016).
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 6), sharey=True)
    fig.suptitle(title, fontsize=13, fontweight="bold")

    t0, t1 = time_ns[0], time_ns[-1]
    f0, f1 = freqs_GHz[0], freqs_GHz[-1]
    extent_spec = [f0, f1, t1, t0]

    # ── Plot 1: Amplitude Spectrum (dB) ──────────────────────────────────────
    ax = axes[0]
    A_norm = A / (A.max() + 1e-30)
    A_dB   = np.clip(20.0 * np.log10(A_norm + 1e-30), -db_clip, 0.0)
    im0 = ax.imshow(
        A_dB.T, aspect="auto", origin="upper", extent=extent_spec,
        cmap="jet", vmin=-db_clip, vmax=0.0, interpolation="bilinear",
    )
    fig.colorbar(im0, ax=ax, fraction=0.046, pad=0.04).set_label("Amplitude [dB]", fontsize=9)
    ax.set_xlabel("Frequency [GHz]", fontsize=10)
    ax.set_ylabel("TWT [ns]", fontsize=10)
    ax.set_title("CLSSA Amplitude Spectrum", fontsize=11)

    # ── Plot 2: Phase Spectrum ────────────────────────────────────────────────
    ax = axes[1]
    im1 = ax.imshow(
        theta_2d_deg.T, aspect="auto", origin="upper", extent=extent_spec,
        cmap="RdBu", vmin=-180.0, vmax=180.0, interpolation="nearest",
    )
    fig.colorbar(im1, ax=ax, fraction=0.046, pad=0.04).set_label("Phase [deg]", fontsize=9)
    ax.set_xlabel("Frequency [GHz]", fontsize=10)
    ax.set_title("CLSSA Phase Spectrum", fontsize=11)

    # ── Plot 3: Phase Gather ─────────────────────────────────────────────────
    ax = axes[2]
    pg_clip = np.percentile(np.abs(phase_gather), 99)
    im2 = ax.imshow(
        phase_gather.T, aspect="auto", origin="upper",
        extent=[theta_deg[0], theta_deg[-1], t1, t0],
        cmap="RdBu", vmin=-pg_clip, vmax=pg_clip, interpolation="bilinear",
    )
    fig.colorbar(im2, ax=ax, fraction=0.046, pad=0.04).set_label("Amplitude", fontsize=9)
    ax.set_xlabel("Phase [deg]", fontsize=10)
    ax.set_title("Phase Gather  S'(theta, t)", fontsize=11)

    # Standard ±90° markers
    for ph in (-90, 90):
        ax.axvline(ph, color="lime", lw=1.0, ls="--", alpha=0.85)
    ax.axvline(0, color="white", lw=0.6, ls=":", alpha=0.5)
    ax.set_xticks([-180, -90, 0, 90, 180])

    # -90° anomaly highlight (thin-bed low-impedance indicator)
    if highlight_neg90:
        ax.axvspan(-105, -75, color="magenta", alpha=0.18, zorder=0)
        ax.text(-90, t0 + (t1 - t0) * 0.02, "-90°\nthin-bed",
                color="magenta", fontsize=7, ha="center", va="top", fontweight="bold")

    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Demo / standalone entry point
# ---------------------------------------------------------------------------

def _synthetic_ricker(n, dt, f0_GHz=1.5, t_peak_ns=None):
    """Ricker wavelet with added thin-layer interference."""
    if t_peak_ns is None:
        t_peak_ns = n * dt / 3
    t = np.arange(n) * dt
    u = np.pi * f0_GHz * (t - t_peak_ns)
    ricker = (1.0 - 2.0 * u**2) * np.exp(-(u**2))
    t_thin = t_peak_ns + 0.4 / f0_GHz
    u2 = np.pi * f0_GHz * (t - t_thin)
    thin_layer = 0.4 * (1.0 - 2.0 * u2**2) * np.exp(-(u2**2))
    noise = 0.05 * np.random.default_rng(42).standard_normal(n)
    return ricker + thin_layer + noise


if __name__ == "__main__":
    dt  = 0.004717        # ns  (from notebook)
    n_t = 4241
    f_c = 1.5             # GHz

    print("Phase decomposition demo — synthetic Ricker + thin-layer trace")
    print(f"  dt = {dt:.6f} ns,  n_t = {n_t},  f_c = {f_c} GHz")

    trace = _synthetic_ricker(n_t, dt, f0_GHz=f_c)

    print("Running CWT and building phase gather …", end=" ", flush=True)
    A, theta_2d_deg, phase_gather, freqs_GHz, time_ns, theta_deg = (
        phase_decomposition(
            trace, dt,
            f_min_GHz=0.5,
            f_max_GHz=4.5,
            n_scales=60,
            n_theta=181,
            reflect_pad=True,
        )
    )
    print("done.")
    print(f"  A shape      : {A.shape}  (n_freqs × n_time)")
    print(f"  phase_gather : {phase_gather.shape}  (n_theta × n_time)")
    print(f"  freq range   : {freqs_GHz[0]:.3f} – {freqs_GHz[-1]:.3f} GHz")

    fig = plot_decomposition(
        A, theta_2d_deg, phase_gather, freqs_GHz, time_ns, theta_deg,
        title=f"Phase Decomposition — synthetic GPR trace  (f_c = {f_c} GHz)",
        show_coi=True,
    )
    plt.show()
