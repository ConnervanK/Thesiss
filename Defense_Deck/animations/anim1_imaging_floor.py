"""
anim1_imaging_floor.py   --   SLIDES 3 & 4
==========================================
"Where the image gives up."

Two STATIC point scatterers walk toward each other. Three panels tell the
whole imaging-resolution story in one continuous motion:

  (a) Raw B-scan          two hyperbolas merging into one
  (b) Migrated image      two focused lobes merging into one
  (c) Lateral profile     two peaks merging into one, with the live
                          Rayleigh ratio and a RESOLVED / MERGED verdict

The verdict badge flips from green to red as the separation crosses the
half-wavelength floor. That flip is the moment the slide exists for.

Nothing here is time-lapse. This is the *static* problem, deliberately, so
that slide 5 can then pivot to monitoring.

Usage
-----
    python anim1_imaging_floor.py                       # 18 s of motion + 3 s freeze
    python anim1_imaging_floor.py --duration 25         # slower still
    python anim1_imaging_floor.py --frames 30           # fast layout preview
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from gpr_common import (LAMBDA, V, FC, DEPTH, ImageGrid, apply_house_style,
                        synth_bscan, lateral_profile, fwhm, two_peaks_resolved,
                        readout, verdict_box, set_verdict, slide_title, footnote,
                        save_animation, save_poster, common_cli, frame_counts,
                        C_DARK, C_ACCENT, C_WARM, C_GREY, C_GOOD, C_BAD,
                        CMAP_IMG, FIGSIZE)

# --------------------------------------------------------------------------
# Sweep: separation from 2 lambda down to 1/8 lambda, geometric so the
# interesting sub-wavelength end gets plenty of frames.
# --------------------------------------------------------------------------
SEP_START = 2.00 * LAMBDA
SEP_END = 0.125 * LAMBDA
# Frame count now comes from --duration / --fps (see gpr_common.frame_counts).

# Raw B-scan geometry
NX_TRACES = 260
X_HALF = 6.5                       # +/- extent of the B-scan panel [m]
T_MAX = 160.0                      # [ns]
NT = 700


def build_sweep(n_frames):
    return np.geomspace(SEP_START, SEP_END, n_frames)


def main():
    args = common_cli(__doc__.splitlines()[1]).parse_args()
    n_frames, n_hold = frame_counts(args)

    apply_house_style(base=15)
    os.makedirs(args.outdir, exist_ok=True)

    seps = build_sweep(n_frames)
    seps = np.concatenate([seps, np.repeat(seps[-1], n_hold)])

    # ---------------- static setup ----------------
    grid = ImageGrid(n=256, lam=LAMBDA, samples_per_lambda=16)
    S = grid.psf_spectrum(aperture_deg=50.0, bw_frac=0.85)

    x_tr = np.linspace(-X_HALF, X_HALF, NX_TRACES)
    t = np.linspace(0.0, T_MAX, NT)

    # Reference FWHM: measured ONCE from a single isolated scatterer, exactly
    # as the thesis does (the 2-lambda reference case).
    ref_img = grid.psf_image(S, 0.0, 0.0)
    xp, ref_prof, _ = lateral_profile(ref_img, grid)
    FWHM_REF = fwhm(xp, ref_prof)

    # ---------------- figure ----------------
    fig = plt.figure(figsize=FIGSIZE)
    slide_title(fig,
                "Two reflectors, closing in",
                f"Zero-offset GPR  ·  illustrative scale: λ = {LAMBDA:.0f} m  "
                f"(v = {V:.2f} m/ns, fᴄ = {FC*1000:.0f} MHz)")
    footnote(fig, "Static imaging — one survey, one snapshot")

    gs = fig.add_gridspec(1, 3, left=0.055, right=0.965, top=0.76, bottom=0.13,
                          wspace=0.30, width_ratios=[1.0, 1.0, 1.15])
    ax_b = fig.add_subplot(gs[0, 0])
    ax_m = fig.add_subplot(gs[0, 1])
    ax_p = fig.add_subplot(gs[0, 2])

    # (a) raw B-scan
    B0 = synth_bscan(x_tr, t, [(-SEP_START / 2, DEPTH, 1.0),
                               (+SEP_START / 2, DEPTH, 1.0)])
    vmaxB = np.abs(B0).max()
    im_b = ax_b.imshow(B0, aspect="auto", cmap=CMAP_IMG,
                       extent=(x_tr[0], x_tr[-1], t[-1], t[0]),
                       vmin=-vmaxB, vmax=vmaxB, interpolation="bilinear")
    ax_b.set_title("(a)  Raw B-scan", color=C_DARK, pad=26)
    ax_b.set_xlabel("antenna position  x  [m]")
    ax_b.set_ylabel("two-way time  [ns]")
    ax_b.set_ylim(150, 60)

    # (b) migrated image  (zoom to +/- 1.6 lambda about the targets)
    zoom = 1.6 * LAMBDA
    img0 = grid.psf_image(S, -SEP_START / 2, 0) + grid.psf_image(S, SEP_START / 2, 0)
    vmaxM = np.abs(img0).max()
    # aspect="auto" (not "equal") so all three panel titles sit on one line
    im_m = ax_m.imshow(img0, aspect="auto", cmap=CMAP_IMG,
                       extent=grid.extent, vmin=-vmaxM, vmax=vmaxM,
                       interpolation="bilinear")
    ax_m.set_xlim(-zoom, zoom)
    ax_m.set_ylim(zoom, -zoom)
    ax_m.set_title("(b)  After migration", color=C_DARK, pad=26)
    ax_m.set_xlabel("x  [m]")
    ax_m.set_ylabel("depth offset  [m]")
    ax_m.set_xticks([-1, 0, 1])
    ax_m.set_yticks([-1, 0, 1])   # keeps the top tick clear of (a)'s readout

    # (c) lateral profile
    ax_p.set_title("(c)  Amplitude across the targets", color=C_DARK, pad=26)
    ax_p.set_xlabel("x  [m]")
    ax_p.set_ylabel("normalised amplitude")
    ax_p.set_xlim(-zoom, zoom)
    ax_p.set_ylim(-0.55, 1.68)   # headroom so the verdict box clears the peak
    ax_p.axhline(0.0, color=C_GREY, lw=0.8, zorder=0)

    (ln_prof,) = ax_p.plot([], [], color=C_ACCENT, lw=3.0, zorder=5)
    fill_l = ax_p.axvspan(0, 0, color=C_GREY, alpha=0.16, lw=0)
    fill_r = ax_p.axvspan(0, 0, color=C_GREY, alpha=0.16, lw=0)
    vl_l = ax_p.axvline(0, color=C_WARM, ls="--", lw=1.6)
    vl_r = ax_p.axvline(0, color=C_WARM, ls="--", lw=1.6)
    ax_p.plot([], [], color=C_WARM, ls="--", lw=1.6, label="true positions")
    ax_p.plot([], [], color=C_GREY, lw=7, alpha=0.35,
              label=f"reference FWHM = {FWHM_REF/LAMBDA:.2f} λ")
    ax_p.legend(loc="lower left", fontsize=11.5)

    txt_sep = readout(ax_b, "", loc="upper left", size=15)
    # verdict + ratio in ONE box, pinned inside the axes (legend is lower-left,
    # the profile peaks at centre) so it can never collide with anything
    verdict = verdict_box(ax_p, loc="upper right", size=14)

    # ---------------- frame update ----------------
    def update(i):
        sep = seps[i]
        h = sep / 2.0

        # (a)
        B = synth_bscan(x_tr, t, [(-h, DEPTH, 1.0), (+h, DEPTH, 1.0)])
        im_b.set_data(B)
        im_b.set_clim(-vmaxB, vmaxB)

        # (b)
        img = grid.psf_image(S, -h, 0.0) + grid.psf_image(S, +h, 0.0)
        im_m.set_data(img)
        im_m.set_clim(-vmaxM, vmaxM)

        # (c)
        xx, prof, _ = lateral_profile(img, grid)
        prof_n = prof / (np.abs(prof).max() + 1e-30)
        ln_prof.set_data(xx, prof_n)

        vl_l.set_xdata([-h, -h])
        vl_r.set_xdata([+h, +h])
        # axvspan returns a Rectangle: move it in data coords via x / width
        for patch, centre in ((fill_l, -h), (fill_r, +h)):
            patch.set_x(centre - FWHM_REF / 2)
            patch.set_width(FWHM_REF)

        ratio = sep / FWHM_REF
        ok = two_peaks_resolved(prof_n) and ratio >= 1.0

        txt_sep.set_text(f"sep = {sep/LAMBDA:.3f} λ  ({sep*1000:.0f} mm)")

        label = "RESOLVED" if ok else "MERGED — one lobe"
        set_verdict(verdict, f"{label}\nseparation / FWHM = {ratio:.2f}", ok)
        ln_prof.set_color(C_ACCENT if ok else C_BAD)

        return [im_b, im_m, ln_prof, vl_l, vl_r, fill_l, fill_r,
                txt_sep, verdict]

    anim = FuncAnimation(fig, update, frames=len(seps),
                         interval=1000 // args.fps, blit=False, repeat=False)

    stem = os.path.join(args.outdir, "anim1_imaging_floor")
    print("Rendering anim1  (slides 3 & 4)  — imaging resolution floor")
    save_animation(anim, stem, fps=args.fps, dpi=args.dpi, prefer=args.format)
    if not args.no_poster:
        update(len(seps) - 1)
        save_poster(fig, stem, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
