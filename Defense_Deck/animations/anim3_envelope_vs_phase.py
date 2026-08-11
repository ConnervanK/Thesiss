"""
anim3_envelope_vs_phase.py   --   SLIDE 7
=========================================
"The envelope lies. The wave underneath doesn't."

The key physical insight, on one scatterer, in three panels:

  (a) A depth trace through the migrated scatterer. Baseline is a grey ghost,
      monitor is orange. As the scatterer creeps down, the ENVELOPE stays
      almost exactly on top of itself while the ZERO-CROSSINGS visibly slide.

  (b) A phasor. The displacement is literally an angle: Δφ = k_c · Δz.
      This is the panel that makes it click for a non-specialist.

  (c) Estimated vs true displacement, accumulating live:
        - reading the ENVELOPE PEAK gives a staircase, quantised by the
          image sample spacing (this is exactly the finite-grid artefact the
          thesis reports for back-propagation at 1/32 λ, where the baseline
          and monitor peaks land on the identical sample and the measured
          separation is exactly zero);
        - reading the PHASE lies on the 1:1 line, continuously, with no floor.

Usage
-----
    python anim3_envelope_vs_phase.py                       # 18 s of motion + 3 s freeze
    python anim3_envelope_vs_phase.py --duration 25         # slower still
    python anim3_envelope_vs_phase.py --frames 30           # fast layout preview
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Arc, FancyArrow

from gpr_common import (LAMBDA, ImageGrid, apply_house_style, readout,
                        slide_title, footnote, save_animation, save_poster,
                        common_cli, frame_counts, C_DARK, C_ACCENT, C_WARM, C_GREY,
                        C_PANEL, FIGSIZE)

# Stays well inside λ/2 so the phasor never wraps AND the envelope genuinely
# still looks unmoved — pushing further would undercut the slide's own claim.
SHIFT_MAX = 0.30 * LAMBDA


def analytic_signal(x: np.ndarray) -> np.ndarray:
    """Hilbert analytic signal via FFT (avoids a scipy dependency)."""
    n = x.size
    X = np.fft.fft(x)
    h = np.zeros(n)
    if n % 2 == 0:
        h[0] = h[n // 2] = 1.0
        h[1:n // 2] = 2.0
    else:
        h[0] = 1.0
        h[1:(n + 1) // 2] = 2.0
    return np.fft.ifft(X * h)


def main():
    args = common_cli(__doc__.splitlines()[1]).parse_args()
    n_frames, n_hold = frame_counts(args)

    apply_house_style(base=14.5)
    os.makedirs(args.outdir, exist_ok=True)

    shifts = np.linspace(0.0, SHIFT_MAX, n_frames)
    shifts = np.concatenate([shifts, np.repeat(shifts[-1], n_hold)])

    grid = ImageGrid(n=256, lam=LAMBDA, samples_per_lambda=16)
    S = grid.psf_spectrum(aperture_deg=50.0, bw_frac=0.85)
    ix0 = grid.n // 2                      # trace at x = 0

    # Effective vertical carrier.
    #
    # It is tempting to use k_c = 2*pi/lambda, but that is subtly wrong: an
    # aperture-limited migrated image has k_z = k*cos(theta), so its *vertical*
    # carrier is the energy-weighted mean of k_z over the aperture, which sits
    # a few percent below k_c. Using k_c instead would bias every phase
    # estimate low by ~8% -- visible as the teal line drifting off the 1:1
    # line in panel (c).
    #
    # This is precisely why the thesis fits a PLANE across the whole measured
    # band rather than assuming a single carrier wavenumber. Worth mentioning
    # if someone asks why the global Fourier fit is preferred over reading a
    # local phase gradient.
    _w = S ** 2
    _pos = grid.KZ > 0
    kc = float(np.sum(grid.KZ[_pos] * _w[_pos]) / np.sum(_w[_pos]))

    # ---- baseline trace and its reference sample -------------------------
    img_b = grid.psf_image(S, 0.0, 0.0)
    NORM = np.abs(img_b[:, ix0]).max()     # normalise so the y-axis reads 0..1
    tr_b = img_b[:, ix0] / NORM
    A_b = analytic_signal(tr_b)
    env_b = np.abs(A_b)
    iref = int(np.argmax(env_b))           # fixed measurement point
    z = grid.z
    z_peak_b = z[iref]

    zoom = 0.90 * LAMBDA
    amp = 1.0

    # ---------------- figure ----------------
    fig = plt.figure(figsize=FIGSIZE)
    # NOTE ON WORDING: the migrated envelope of a single scatterer really does
    # translate by Δz -- claiming it "does not move" would be overselling, and
    # a sharp committee will catch it. What is actually true, and what this
    # animation shows, is that the envelope can only be READ on a discrete
    # grid (a staircase), while the phase beneath it is continuous and exact.
    # The title is worded to match what panel (c) actually demonstrates.
    slide_title(fig,
                "Amplitude reads in steps. Phase reads continuously.",
                f"Depth trace through one migrated scatterer  ·  λ = {LAMBDA:.0f} m  ·  "
                f"image sample spacing = λ/{grid.samples_per_lambda:.0f}")
    footnote(fig, r"$\Delta\varphi = k_c\,\Delta z$   —   an angle, not an amplitude")

    gs = fig.add_gridspec(2, 2, left=0.055, right=0.965, top=0.755, bottom=0.115,
                          wspace=0.26, hspace=0.50,
                          width_ratios=[1.55, 1.0], height_ratios=[1.0, 1.0])
    ax_tr = fig.add_subplot(gs[:, 0])
    ax_ph = fig.add_subplot(gs[0, 1])
    ax_es = fig.add_subplot(gs[1, 1])

    # ---- (a) the trace ---------------------------------------------------
    ax_tr.set_title("(a)  The wave under the envelope", color=C_DARK, pad=24)
    ax_tr.set_xlabel("depth offset  [m]")
    ax_tr.set_ylabel("normalised amplitude")
    ax_tr.set_xlim(z_peak_b - zoom, z_peak_b + zoom)
    ax_tr.set_ylim(-1.30 * amp, 1.95 * amp)
    ax_tr.axhline(0, color=C_GREY, lw=0.8, zorder=0)

    ax_tr.plot(z, tr_b, color=C_GREY, lw=2.2, ls="--", label="baseline", zorder=3)
    ax_tr.plot(z, env_b, color=C_GREY, lw=1.4, alpha=0.75, zorder=2)
    ax_tr.plot(z, -env_b, color=C_GREY, lw=1.4, alpha=0.75, zorder=2)

    (ln_mon,) = ax_tr.plot([], [], color=C_WARM, lw=3.0, label="monitor", zorder=6)
    (ln_env,) = ax_tr.plot([], [], color=C_WARM, lw=1.5, alpha=0.8, zorder=5)
    (ln_env2,) = ax_tr.plot([], [], color=C_WARM, lw=1.5, alpha=0.8, zorder=5)
    (mk_zc_b,) = ax_tr.plot([], [], "v", color=C_DARK, ms=11, zorder=8)
    (mk_zc_m,) = ax_tr.plot([], [], "v", color=C_WARM, ms=11, zorder=8)
    ax_tr.legend(loc="upper right", fontsize=13)

    txt_tr = readout(ax_tr, "", loc="upper left", size=15)

    # ---- (b) the phasor --------------------------------------------------
    ax_ph.set_title("(b)  Displacement is an angle", color=C_DARK, pad=20)
    ax_ph.set_aspect("equal")
    ax_ph.set_xlim(-1.45, 1.45)
    ax_ph.set_ylim(-1.45, 1.45)
    ax_ph.axis("off")
    th = np.linspace(0, 2 * np.pi, 400)
    ax_ph.plot(np.cos(th), np.sin(th), color=C_GREY, lw=1.6)
    ax_ph.plot([-1.25, 1.25], [0, 0], color=C_GREY, lw=0.8)
    ax_ph.plot([0, 0], [-1.25, 1.25], color=C_GREY, lw=0.8)
    ax_ph.annotate("", xy=(1.0, 0.0), xytext=(0, 0),
                   arrowprops=dict(arrowstyle="-|>", lw=2.6, color=C_GREY))
    ax_ph.text(1.10, -0.16, "baseline", color=C_GREY, fontsize=12.5)
    arrow_mon = ax_ph.annotate("", xy=(1.0, 0.0), xytext=(0, 0),
                               arrowprops=dict(arrowstyle="-|>", lw=3.4,
                                               color=C_WARM))
    arc = Arc((0, 0), 1.0, 1.0, theta1=0, theta2=0, lw=3.0, color=C_ACCENT)
    ax_ph.add_patch(arc)
    txt_ang = ax_ph.text(0, -1.36, "", ha="center", va="top", fontsize=16,
                         fontweight="bold", color=C_ACCENT)

    # ---- (c) estimate vs truth ------------------------------------------
    ax_es.set_title("(c)  Recovered displacement", color=C_DARK, pad=20)
    ax_es.set_xlabel("true Δz  [λ]")
    ax_es.set_ylabel("estimate  [λ]")
    ax_es.set_xlim(0, SHIFT_MAX / LAMBDA)
    ax_es.set_ylim(0, SHIFT_MAX / LAMBDA * 1.12)
    ax_es.set_facecolor(C_PANEL)
    lim = SHIFT_MAX / LAMBDA
    ax_es.plot([0, lim], [0, lim], color=C_DARK, ls=":", lw=1.6,
               label="perfect (1:1)", zorder=2)
    (ln_env_est,) = ax_es.plot([], [], color=C_GREY, lw=3.0, drawstyle="steps-post",
                               label="from envelope peak", zorder=4)
    (ln_ph_est,) = ax_es.plot([], [], color=C_ACCENT, lw=3.2,
                              label="from phase", zorder=5)
    ax_es.legend(loc="upper left", fontsize=12)

    hist_true, hist_env, hist_ph = [], [], []

    # ---------------- frame update ----------------
    def update(i):
        d = shifts[i]
        if i == 0:
            hist_true.clear(); hist_env.clear(); hist_ph.clear()

        img_m = grid.psf_image(S, 0.0, d)
        tr_m = img_m[:, ix0] / NORM
        A_m = analytic_signal(tr_m)
        env_m = np.abs(A_m)

        ln_mon.set_data(z, tr_m)
        ln_env.set_data(z, env_m)
        ln_env2.set_data(z, -env_m)

        # zero-crossing markers, just left of each envelope peak
        def first_zero_after(tr, i_start):
            for j in range(i_start, tr.size - 1):
                if tr[j] * tr[j + 1] <= 0:
                    f = tr[j] / (tr[j] - tr[j + 1] + 1e-30)
                    return z[j] + f * (z[j + 1] - z[j])
            return np.nan

        zc_b = first_zero_after(tr_b, iref)
        zc_m = first_zero_after(tr_m, iref)
        mk_zc_b.set_data([zc_b], [1.16 * amp])
        mk_zc_m.set_data([zc_m], [1.16 * amp])

        # -- estimate 1: envelope peak, quantised by the sample grid
        env_est = z[int(np.argmax(env_m))] - z_peak_b

        # -- estimate 2: phase at the fixed reference sample
        #    m(z) = b(z-d)  =>  Δφ = -k_c·d   =>   d = -Δφ / k_c
        dphi = np.angle(A_m[iref] * np.conj(A_b[iref]))
        ph_est = -dphi / kc

        hist_true.append(d / LAMBDA)
        hist_env.append(env_est / LAMBDA)
        hist_ph.append(ph_est / LAMBDA)
        ln_env_est.set_data(hist_true, hist_env)
        ln_ph_est.set_data(hist_true, hist_ph)

        # phasor
        ang = -dphi                     # positive as the scatterer descends
        arrow_mon.xy = (np.cos(ang), np.sin(ang))
        arc.theta2 = np.degrees(ang)
        txt_ang.set_text(f"Δφ = {np.degrees(ang):5.1f}°")

        env_err = abs(env_est - d) * 1000
        txt_tr.set_text(
            f"true Δz          {d*1000:6.0f} mm   ({d/LAMBDA:.3f} λ)\n"
            f"envelope peak    {env_est*1000:6.0f} mm   err {env_err:5.0f} mm\n"
            f"phase            {ph_est*1000:6.0f} mm   err "
            f"{abs(ph_est-d)*1000:5.1f} mm")
        txt_tr.set_family("DejaVu Sans Mono")

        return [ln_mon, ln_env, ln_env2, mk_zc_b, mk_zc_m,
                ln_env_est, ln_ph_est, arrow_mon, arc, txt_ang, txt_tr]

    anim = FuncAnimation(fig, update, frames=len(shifts),
                         interval=1000 // args.fps, blit=False, repeat=False)

    stem = os.path.join(args.outdir, "anim3_envelope_vs_phase")
    print("Rendering anim3  (slide 7)  — envelope vs phase")
    save_animation(anim, stem, fps=args.fps, dpi=args.dpi, prefer=args.format)
    if not args.no_poster:
        for k in range(len(shifts)):     # replay so the history curves are full
            update(k)
        save_poster(fig, stem, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
