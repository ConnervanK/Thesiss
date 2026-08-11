"""
anim4_kdomain_phase_ramp.py   --   SLIDES 8 & 10
================================================
"A shift in space is a tilt in wavenumber."

This is the transition the whole method rests on: a moving scatterer in the
image domain, and the SAME event rendered as the cross-spectrum phase, live,
side by side. Flat colour means no shift. A colour ramp means a shift, and
the ramp's slope IS the displacement.

Three renders come out of this one script:

  --mode lateral   SLIDE 8   image | phase | 1-D fit along k_x | error sweep
  --mode vertical            same, but the ramp runs along k_z
  --mode both      SLIDE 10  lateral vs vertical, sharing one error sweep

The displacement is recovered live by the same masked, amplitude-weighted
least-squares fit used in the thesis, so the number printed on screen is the
method's own output -- not an annotation.

Everything is expressed in WAVELENGTHS, not millimetres. The whole argument is
about resolving a fraction of a wavelength, so lambda is the unit that carries
the meaning; a reader should never have to divide by 1000 in their head to see
that 0.0156 lambda is one sixty-fourth of a wavelength.

The sweep and the error panel
-----------------------------
The displacement sweeps GEOMETRICALLY from 1/64 lambda to 1/2 lambda, so the
x-axis is the same ladder of scenario scales the thesis uses, and the error
panel accumulates live underneath it. Measured behaviour of this estimator:

    shift      lateral err    vertical err
    1/64 l        0.29 %         0.05 %
    1/4  l        0.29 %         0.05 %
    3/8  l        0.29 %         4.0  %
    1/2  l        0.29 %        79    %

Flat, then a cliff. That cliff is the estimator's unambiguous range, and it is
worth showing rather than hiding: the aperture limits |k_x| to about 0.77 k_c
so a lateral shift stays unwrapped to roughly 0.65 lambda, whereas k_z reaches
k_c so a vertical shift wraps at 0.5 lambda. Same plane, same fit -- the two
axes simply have different unambiguous ranges. If someone asks why the field
chapter restricts its k-space ROI, this panel is the answer.

Usage
-----
    python anim4_kdomain_phase_ramp.py --mode lateral
    python anim4_kdomain_phase_ramp.py --mode both
    python anim4_kdomain_phase_ramp.py --mode all        # renders all three
    python anim4_kdomain_phase_ramp.py --duration 25     # slower
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from gpr_common import (LAMBDA, ImageGrid, apply_house_style, cross_spectrum,
                        fit_phase_plane, readout, slide_title, footnote,
                        save_animation, save_poster, common_cli, frame_counts,
                        C_DARK, C_ACCENT, C_WARM, C_GREY, C_PANEL,
                        CMAP_IMG, CMAP_PHASE, FIGSIZE)

# Geometric sweep over the thesis's own ladder of displacement scales.
SHIFT_MIN = LAMBDA / 64.0
SHIFT_MAX = LAMBDA / 2.0
K_VIEW = 2.2                   # k-axis limit, in units of k_c
SUBSAMPLE = 9                  # thin the scatter so the fit panel stays legible
ERR_THRESHOLD = 5.0            # %, the same threshold as the thesis's map

# Tick ladder for the error panel's x-axis, in units of lambda
LAM_TICKS = [1 / 64, 1 / 32, 1 / 16, 1 / 8, 1 / 4, 1 / 2]
LAM_LABELS = ["¹⁄₆₄λ", "¹⁄₃₂λ", "¹⁄₁₆λ", "⅛λ", "¼λ", "½λ"]


def make_error_panel(ax, title="(d)  Estimate error across the sweep"):
    """Log-log error panel: |error| in % of the true shift, vs shift in λ."""
    ax.set_title(title, color=C_DARK, pad=22)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(SHIFT_MIN / LAMBDA * 0.85, SHIFT_MAX / LAMBDA * 1.18)
    ax.set_ylim(1e-3, 3e2)
    ax.set_xlabel("true displacement   [λ]")
    ax.set_ylabel("|error|   [% of true]")
    ax.set_facecolor(C_PANEL)
    ax.set_xticks(LAM_TICKS)
    ax.set_xticklabels(LAM_LABELS, fontsize=10.5)
    ax.tick_params(axis="x", which="minor", bottom=False)
    ax.grid(True, which="major", color="white", lw=1.0, alpha=0.9)
    ax.axhline(ERR_THRESHOLD, color=C_GREY, ls="--", lw=1.6, zorder=2)
    ax.text(SHIFT_MAX / LAMBDA * 1.10, ERR_THRESHOLD * 1.4,
            f"{ERR_THRESHOLD:.0f}% threshold", fontsize=10.5, color=C_GREY,
            va="bottom", ha="right")
    return ax


# ==========================================================================
#  shared per-frame physics
# ==========================================================================

class Engine:
    def __init__(self):
        self.grid = ImageGrid(n=256, lam=LAMBDA, samples_per_lambda=16)
        self.S = self.grid.psf_spectrum(aperture_deg=50.0, bw_frac=0.85)
        self.kc = 2.0 * np.pi / LAMBDA
        self.img_base = self.grid.psf_image(self.S, 0.0, 0.0)
        self.vmax = np.abs(self.img_base).max()

    def frame(self, dx, dz):
        g = self.grid
        img_mon = g.psf_image(self.S, dx, dz)
        XS = cross_spectrum(self.img_base, img_mon, g)
        est_dz, est_dx, _, mask = fit_phase_plane(XS, g)
        return img_mon, XS, mask, est_dx, est_dz


# ==========================================================================
#  panel builders
# ==========================================================================

def make_image_panel(ax, eng, title, arrow_dir):
    zoom = 1.35 * LAMBDA
    im = ax.imshow(eng.img_base, aspect="auto", cmap=CMAP_IMG,
                   extent=eng.grid.extent, vmin=-eng.vmax, vmax=eng.vmax,
                   interpolation="bilinear")
    ax.set_xlim(-zoom, zoom)
    ax.set_ylim(zoom, -zoom)
    ax.set_title(title, color=C_DARK, pad=22)
    ax.set_xlabel("x  [m]")
    ax.set_ylabel("depth offset  [m]")
    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([-1, 0, 1])
    # motion arrow
    if arrow_dir == "lateral":
        ax.annotate("", xy=(0.72, 0.0), xytext=(0.06, 0.0),
                    arrowprops=dict(arrowstyle="-|>", lw=3.0, color=C_WARM))
    else:
        ax.annotate("", xy=(0.0, 0.72), xytext=(0.0, 0.06),
                    arrowprops=dict(arrowstyle="-|>", lw=3.0, color=C_WARM))
    return im


def make_phase_panel(ax, eng, title):
    g = eng.grid
    kc = eng.kc
    # origin="lower" so k_z increases upward, as a wavenumber axis should
    ext = (g.kx[0] / kc, g.kx[-1] / kc, g.kz[0] / kc, g.kz[-1] / kc)
    rgba = np.zeros((g.n, g.n, 4))
    im = ax.imshow(rgba, aspect="auto", extent=ext, origin="lower",
                   interpolation="nearest", zorder=3)
    ax.set_facecolor(C_PANEL)
    ax.set_xlim(-K_VIEW, K_VIEW)
    ax.set_ylim(-K_VIEW, K_VIEW)
    ax.set_title(title, color=C_DARK, pad=22)
    ax.set_xlabel("$k_x$   [$k_c$]")
    ax.set_ylabel("$k_z$   [$k_c$]")
    ax.axhline(0, color=C_GREY, lw=0.7, zorder=4)
    ax.axvline(0, color=C_GREY, lw=0.7, zorder=4)
    return im


def phase_rgba(XS, mask):
    """Phase as RGBA: opaque inside the fitting mask, transparent outside."""
    cmap = plt.get_cmap(CMAP_PHASE)
    phi = np.angle(XS)
    rgba = cmap((phi + np.pi) / (2 * np.pi))
    rgba[..., 3] = np.where(mask, 1.0, 0.06)
    return rgba


# ==========================================================================
#  render: single-direction, three panels (slide 8)
# ==========================================================================

def render_single(mode, args):
    eng = Engine()
    g, kc = eng.grid, eng.kc
    lateral = (mode == "lateral")

    n_frames, n_hold = frame_counts(args)
    shifts = np.geomspace(SHIFT_MIN, SHIFT_MAX, n_frames)
    shifts = np.concatenate([shifts, np.repeat(shifts[-1], n_hold)])

    fig = plt.figure(figsize=FIGSIZE)
    slide_title(
        fig,
        "Move in space → tilt in wavenumber" if lateral
        else "Vertical motion tilts the other axis",
        ("Lateral displacement  ·  the ramp runs along $k_x$" if lateral
         else "Vertical displacement  ·  the ramp runs along $k_z$")
        + "  ·  swept from ¹⁄₆₄λ to ½λ")
    footnote(fig, "Fourier shift theorem:  a translation rotates phase, "
                  "leaves amplitude untouched")

    gs = fig.add_gridspec(2, 2, left=0.062, right=0.945, top=0.755, bottom=0.135,
                          wspace=0.40, hspace=0.66)
    ax_i = fig.add_subplot(gs[0, 0])
    ax_k = fig.add_subplot(gs[0, 1])
    ax_f = fig.add_subplot(gs[1, 0])
    ax_e = fig.add_subplot(gs[1, 1])

    im_i = make_image_panel(ax_i, eng, "(a)  Migrated image", mode)
    im_k = make_phase_panel(ax_k, eng, "(b)  Cross-spectrum phase")

    # (c) 1-D fit
    axis_lbl = "$k_x$" if lateral else "$k_z$"
    ax_f.set_title(f"(c)  Fit the slope along {axis_lbl}", color=C_DARK, pad=22)
    ax_f.set_xlabel(f"{axis_lbl}   [$k_c$]")
    ax_f.set_ylabel("cross-spectrum phase  [rad]")
    ax_f.set_xlim(-K_VIEW, K_VIEW)
    ax_f.set_ylim(-np.pi * 1.08, np.pi * 1.08)
    ax_f.set_facecolor(C_PANEL)
    ax_f.axhline(0, color=C_GREY, lw=0.8)
    ax_f.set_yticks([-np.pi, 0, np.pi])
    ax_f.set_yticklabels([r"$-\pi$", "0", r"$\pi$"])
    # fit line UNDER the scatter, so the measurements stay visible
    (ln_fit,) = ax_f.plot([], [], color=C_WARM, lw=3.4, zorder=3,
                          label="weighted least-squares fit")
    sc = ax_f.scatter([], [], s=26, c=[], cmap="viridis", vmin=0, vmax=1,
                      alpha=0.9, zorder=5, edgecolors="none")
    ax_f.legend(loc="lower right", fontsize=11.5)
    cb = fig.colorbar(sc, ax=ax_f, fraction=0.05, pad=0.03)
    cb.set_label("fit weight  |XS|", fontsize=11.5)

    # (d) error sweep, accumulating live
    make_error_panel(ax_e)
    col = C_ACCENT if lateral else C_WARM
    (ln_err,) = ax_e.plot([], [], color=col, lw=3.2, zorder=6,
                          label="lateral" if lateral else "vertical")
    (mk_err,) = ax_e.plot([], [], "o", color=col, ms=9, zorder=7)
    ax_e.legend(loc="lower right", fontsize=11.5)

    txt_i = readout(ax_i, "", loc="upper left", size=15)
    # The numeric readout belongs beside the error curve, not on top of the
    # fit scatter — panel (d)'s top-left is empty by construction.
    txt_f = readout(ax_e, "", loc="upper left", size=14, color=C_ACCENT)
    txt_f.set_family("DejaVu Sans Mono")

    hist_d, hist_e = [], []

    def update(idx):
        d = shifts[idx]
        if idx == 0:
            hist_d.clear(); hist_e.clear()
        dx, dz = (d, 0.0) if lateral else (0.0, d)
        img_mon, XS, mask, est_dx, est_dz = eng.frame(dx, dz)

        im_i.set_data(img_mon)
        im_k.set_data(phase_rgba(XS, mask))

        # 1-D projection: remove the already-fitted orthogonal contribution
        phi = np.angle(XS)
        w = np.abs(XS)
        if lateral:
            resid = phi[mask] - g.KZ[mask] * est_dz
            kax = g.KX[mask] / kc
            slope = est_dx
        else:
            resid = phi[mask] - g.KX[mask] * est_dx
            kax = g.KZ[mask] / kc
            slope = est_dz
        ww = w[mask]
        ww = ww / (ww.max() + 1e-30)

        sel = slice(None, None, SUBSAMPLE)
        sc.set_offsets(np.column_stack([kax[sel], resid[sel]]))
        sc.set_array(ww[sel])

        # draw the fitted line only across the k-range the mask actually spans
        if kax.size:
            kk = np.linspace(kax.min(), kax.max(), 60)
        else:
            kk = np.linspace(-1, 1, 60)
        ln_fit.set_data(kk, slope * kk * kc)

        # Everything in wavelengths: the argument is about fractions of λ,
        # so λ is the unit that carries the meaning.
        true_l = d / LAMBDA
        est_l = (est_dx if lateral else est_dz) / LAMBDA
        err_pct = abs(est_l - true_l) / true_l * 100.0

        hist_d.append(true_l)
        hist_e.append(max(err_pct, 1e-3))       # log axis: clamp off zero
        ln_err.set_data(hist_d, hist_e)
        mk_err.set_data([hist_d[-1]], [hist_e[-1]])

        txt_i.set_text(f"true shift = {true_l:.4f} λ")
        txt_f.set_text(f"true       {true_l:8.4f} λ\n"
                       f"recovered  {est_l:8.4f} λ\n"
                       f"error      {err_pct:8.2f} %")
        return [im_i, im_k, sc, ln_fit, ln_err, mk_err, txt_i, txt_f]

    anim = FuncAnimation(fig, update, frames=len(shifts),
                         interval=1000 // args.fps, blit=False, repeat=False)

    stem = os.path.join(args.outdir, f"anim4_kdomain_{mode}")
    print(f"Rendering anim4 [{mode}]  (slide 8)  — shift → phase ramp")
    save_animation(anim, stem, fps=args.fps, dpi=args.dpi, prefer=args.format)
    if not args.no_poster:
        update(len(shifts) - 1)
        save_poster(fig, stem, dpi=300)
    plt.close(fig)


# ==========================================================================
#  render: side-by-side comparison (slide 10)
# ==========================================================================

def render_both(args):
    eng = Engine()
    g, kc = eng.grid, eng.kc

    n_frames, n_hold = frame_counts(args)
    shifts = np.geomspace(SHIFT_MIN, SHIFT_MAX, n_frames)
    shifts = np.concatenate([shifts, np.repeat(shifts[-1], n_hold)])

    fig = plt.figure(figsize=FIGSIZE)
    slide_title(fig,
                "Lateral tilts one way. Vertical tilts the other.",
                "Two independent slopes of the same phase plane  ·  "
                "swept from ¹⁄₆₄λ to ½λ")
    footnote(fig, "Φ(k_z, k_x) = k_z·Δz + k_x·Δx")

    # Left two columns: the two directions. Right column: ONE shared error
    # panel spanning both rows, so the two curves can be compared directly.
    gs = fig.add_gridspec(2, 3, left=0.058, right=0.965, top=0.755, bottom=0.135,
                          wspace=0.46, hspace=0.66,
                          width_ratios=[1.0, 1.0, 1.15])
    axes = {
        ("lateral", "img"): fig.add_subplot(gs[0, 0]),
        ("lateral", "k"): fig.add_subplot(gs[0, 1]),
        ("vertical", "img"): fig.add_subplot(gs[1, 0]),
        ("vertical", "k"): fig.add_subplot(gs[1, 1]),
    }
    ax_e = fig.add_subplot(gs[:, 2])

    im, txt = {}, {}
    for row, mode in enumerate(("lateral", "vertical")):
        ai = axes[(mode, "img")]
        ak = axes[(mode, "k")]
        im[(mode, "img")] = make_image_panel(
            ai, eng, f"{'LATERAL' if row == 0 else 'VERTICAL'} motion", mode)
        im[(mode, "k")] = make_phase_panel(ak, eng, "cross-spectrum phase")
        txt[mode] = readout(ak, "", loc="upper left", size=13, color=C_ACCENT)
        txt[mode].set_family("DejaVu Sans Mono")

    # shared error panel
    make_error_panel(ax_e, "Estimate error across the sweep")
    lines, marks, hist = {}, {}, {"lateral": ([], []), "vertical": ([], [])}
    for mode, col in (("lateral", C_ACCENT), ("vertical", C_WARM)):
        (lines[mode],) = ax_e.plot([], [], color=col, lw=3.2, zorder=6,
                                   label=mode)
        (marks[mode],) = ax_e.plot([], [], "o", color=col, ms=9, zorder=7)
    ax_e.legend(loc="lower right", fontsize=12)

    # Shared colourbar for phase, placed from the real axes geometry so it
    # lands in the gap between the k column and the error panel no matter how
    # the gridspec is later tweaked.
    sm = plt.cm.ScalarMappable(cmap=CMAP_PHASE,
                               norm=plt.Normalize(-np.pi, np.pi))
    kbox_top = axes[("lateral", "k")].get_position()
    kbox_bot = axes[("vertical", "k")].get_position()
    ebox = ax_e.get_position()
    gap_x = kbox_top.x1 + 0.10 * (ebox.x0 - kbox_top.x1)
    cax = fig.add_axes([gap_x, kbox_bot.y0, 0.011,
                        kbox_top.y1 - kbox_bot.y0])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_ticks([-np.pi, 0, np.pi])
    cb.set_ticklabels([r"$-\pi$", "0", r"$\pi$"])
    cb.ax.tick_params(labelsize=10.5)
    # No set_label(): a right-hand colourbar label reaches straight into the
    # error panel's y-label. A caption above the bar stays in its own lane.
    fig.text(cax.get_position().x0 + 0.005, kbox_top.y1 + 0.018,
             "phase", fontsize=10.5, color=C_GREY, ha="center", va="bottom")

    def update(idx):
        d = shifts[idx]
        true_l = d / LAMBDA
        if idx == 0:
            for k in hist:
                hist[k][0].clear(); hist[k][1].clear()

        artists = []
        for mode in ("lateral", "vertical"):
            dx, dz = (d, 0.0) if mode == "lateral" else (0.0, d)
            img_mon, XS, mask, est_dx, est_dz = eng.frame(dx, dz)
            im[(mode, "img")].set_data(img_mon)
            im[(mode, "k")].set_data(phase_rgba(XS, mask))
            txt[mode].set_text(f"Δx = {est_dx/LAMBDA:+8.4f} λ\n"
                               f"Δz = {est_dz/LAMBDA:+8.4f} λ")

            est_l = (est_dx if mode == "lateral" else est_dz) / LAMBDA
            err = max(abs(est_l - true_l) / true_l * 100.0, 1e-3)
            hist[mode][0].append(true_l)
            hist[mode][1].append(err)
            lines[mode].set_data(*hist[mode])
            marks[mode].set_data([true_l], [err])
            artists += [im[(mode, "img")], im[(mode, "k")], txt[mode],
                        lines[mode], marks[mode]]
        return artists

    anim = FuncAnimation(fig, update, frames=len(shifts),
                         interval=1000 // args.fps, blit=False, repeat=False)

    stem = os.path.join(args.outdir, "anim4_kdomain_both")
    print("Rendering anim4 [both]  (slide 10)  — lateral vs vertical")
    save_animation(anim, stem, fps=args.fps, dpi=args.dpi, prefer=args.format)
    if not args.no_poster:
        update(len(shifts) - 1)
        save_poster(fig, stem, dpi=300)
    plt.close(fig)


# ==========================================================================

def main():
    p = common_cli(__doc__.splitlines()[1])
    p.add_argument("--mode", default="all",
                   choices=["lateral", "vertical", "both", "all"])
    args = p.parse_args()

    apply_house_style(base=14.5)
    os.makedirs(args.outdir, exist_ok=True)

    modes = ["lateral", "vertical", "both"] if args.mode == "all" else [args.mode]
    for m in modes:
        if m == "both":
            render_both(args)
        else:
            render_single(m, args)


if __name__ == "__main__":
    main()
