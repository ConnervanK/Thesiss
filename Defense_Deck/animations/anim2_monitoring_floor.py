"""
anim2_monitoring_floor.py   --   SLIDE 6
========================================
"Subtracting two images doesn't help."

Now it IS time-lapse. One scatterer, two surveys. The scatterer moves by a
displacement that shrinks from 2 lambda down to lambda/32, and we watch the
obvious approach fail:

  (a) Baseline migrated image        static reference
  (b) Monitor migrated image         visually identical below ~lambda/2
  (c) Difference, FIXED colour scale  fades to nothing
  (d) Profiles overlaid + difference  the two lobes stop separating

The critical honesty of this animation is the FIXED colour scale on panel (c).
Auto-scaling a difference image always shows *something*, which is exactly how
people fool themselves. Pinned to the baseline's own peak, the difference
visibly dies -- and the running "peak diff = X% of baseline" readout says how
dead it is.

This is the slide that creates the need for the phase method.

Usage
-----
    python anim2_monitoring_floor.py                       # 18 s of motion + 3 s freeze
    python anim2_monitoring_floor.py --duration 25         # slower still
    python anim2_monitoring_floor.py --frames 30           # fast layout preview
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from gpr_common import (LAMBDA, ImageGrid, apply_house_style, lateral_profile,
                        fwhm, readout, verdict_box, set_verdict,
                        slide_title, footnote,
                        save_animation, save_poster, common_cli, frame_counts,
                        C_DARK, C_ACCENT, C_WARM, C_GREY, C_GOOD, C_BAD,
                        CMAP_IMG, FIGSIZE)

SHIFT_START = 2.00 * LAMBDA
SHIFT_END = LAMBDA / 32.0

# Scale-marker rungs annotated on the sweep bar
RUNGS = [(2.0, "2λ"), (1.0, "1λ"), (0.5, "½λ"), (0.25, "¼λ"),
         (0.125, "⅛λ"), (0.0625, "¹⁄₁₆λ"), (0.03125, "¹⁄₃₂λ")]


def main():
    args = common_cli(__doc__.splitlines()[1]).parse_args()
    n_frames, n_hold = frame_counts(args)

    apply_house_style(base=14.5)
    os.makedirs(args.outdir, exist_ok=True)

    shifts = np.geomspace(SHIFT_START, SHIFT_END, n_frames)
    shifts = np.concatenate([shifts, np.repeat(shifts[-1], n_hold)])

    grid = ImageGrid(n=256, lam=LAMBDA, samples_per_lambda=16)
    S = grid.psf_spectrum(aperture_deg=50.0, bw_frac=0.85)

    img_base = grid.psf_image(S, 0.0, 0.0)
    vmax = np.abs(img_base).max()
    xb, prof_base, _ = lateral_profile(img_base, grid)
    prof_base_n = prof_base / np.abs(prof_base).max()
    FWHM_REF = fwhm(xb, prof_base)

    zoom = 1.5 * LAMBDA

    # ---------------- figure ----------------
    fig = plt.figure(figsize=FIGSIZE)
    slide_title(fig,
                "Differencing two images inherits the same floor",
                f"One scatterer, two surveys  ·  λ = {LAMBDA:.0f} m  ·  "
                "panel (c) is pinned to the baseline's own colour scale")
    footnote(fig, "Time-lapse monitoring — amplitude only")

    gs = fig.add_gridspec(2, 3, left=0.055, right=0.965, top=0.755, bottom=0.135,
                          wspace=0.26, hspace=0.68, height_ratios=[1.0, 0.72])
    ax_b = fig.add_subplot(gs[0, 0])
    ax_m = fig.add_subplot(gs[0, 1])
    ax_d = fig.add_subplot(gs[0, 2])
    ax_p = fig.add_subplot(gs[1, :])

    def style_img(ax, title):
        ax.set_xlim(-zoom, zoom)
        ax.set_ylim(zoom, -zoom)
        ax.set_title(title, color=C_DARK, pad=22)
        ax.set_xlabel("x  [m]")
        ax.set_xticks([-1, 0, 1])
        ax.set_yticks([-1, 0, 1])

    im_b = ax_b.imshow(img_base, aspect="auto", cmap=CMAP_IMG,
                       extent=grid.extent, vmin=-vmax, vmax=vmax,
                       interpolation="bilinear")
    style_img(ax_b, "(a)  Baseline")
    ax_b.set_ylabel("depth offset  [m]")

    im_m = ax_m.imshow(img_base, aspect="auto", cmap=CMAP_IMG,
                       extent=grid.extent, vmin=-vmax, vmax=vmax,
                       interpolation="bilinear")
    style_img(ax_m, "(b)  Monitor")

    im_d = ax_d.imshow(np.zeros_like(img_base), aspect="auto", cmap=CMAP_IMG,
                       extent=grid.extent, vmin=-vmax, vmax=vmax,
                       interpolation="bilinear")
    style_img(ax_d, "(c)  Difference")

    # (d) profiles
    ax_p.set_xlim(-zoom, zoom)
    ax_p.set_ylim(-1.18, 1.60)
    ax_p.axhline(0.0, color=C_GREY, lw=0.8, zorder=0)
    ax_p.set_xlabel("x  [m]")
    ax_p.set_ylabel("normalised amplitude")
    ax_p.plot(xb, prof_base_n, color=C_GREY, lw=2.4, ls="--",
              label="baseline", zorder=3)
    (ln_mon,) = ax_p.plot([], [], color=C_WARM, lw=3.0,
                          label="monitor", zorder=5)
    (ln_dif,) = ax_p.plot([], [], color=C_DARK, lw=2.4,
                          label="difference (same scale)", zorder=4)
    # Legend upper-right and verdict upper-left, both INSIDE the wide profile
    # panel: the curve only occupies the centre, so neither can ever collide.
    ax_p.legend(loc="upper right", ncol=1, fontsize=12)

    txt_shift = readout(ax_b, "", loc="upper left", size=14)
    txt_diff = readout(ax_d, "", loc="lower left", size=14)
    verdict = verdict_box(ax_p, loc="upper left", size=14)

    def update(i):
        d = shifts[i]
        img_mon = grid.psf_image(S, d, 0.0)
        img_dif = img_mon - img_base

        im_m.set_data(img_mon)
        im_d.set_data(img_dif)

        _, prof_mon, _ = lateral_profile(img_mon, grid)
        scale = np.abs(prof_base).max()
        ln_mon.set_data(xb, prof_mon / scale)
        ln_dif.set_data(xb, (prof_mon - prof_base) / scale)

        peak_pct = 100.0 * np.abs(img_dif).max() / vmax
        ratio = d / FWHM_REF
        ok = ratio >= 1.0

        txt_shift.set_text(f"true shift Δx = {d/LAMBDA:.3f} λ  ({d*1000:.0f} mm)")
        txt_diff.set_text(f"peak |diff| = {peak_pct:.0f}% of baseline")

        label = "RESOLVED" if ok else "UNRESOLVABLE from amplitude"
        set_verdict(verdict, f"{label}\nseparation / FWHM = {ratio:.2f}", ok)

        return [im_m, im_d, ln_mon, ln_dif, txt_shift, txt_diff, verdict]

    anim = FuncAnimation(fig, update, frames=len(shifts),
                         interval=1000 // args.fps, blit=False, repeat=False)

    stem = os.path.join(args.outdir, "anim2_monitoring_floor")
    print("Rendering anim2  (slide 6)  — amplitude differencing fails")
    save_animation(anim, stem, fps=args.fps, dpi=args.dpi, prefer=args.format)
    if not args.no_poster:
        update(len(shifts) - 1)
        save_poster(fig, stem, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
