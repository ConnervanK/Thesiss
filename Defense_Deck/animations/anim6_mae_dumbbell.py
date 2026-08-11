"""
anim6_mae_dumbbell.py   --   SLIDES 12 / 13 / 15
================================================
Two plots built from the thesis's own headline MAE tables.

  --which clean   a static dot plot of the CLEAN synthetic results
  --which noisy   an ANIMATION that starts from exactly that clean plot and
                  lets every marker travel to its noisy value, drawing the
                  connecting bar as it goes. You watch noise degrade the
                  method rather than being shown a before/after pair.

The second grows out of the first on purpose: show the clean plot, then play
the animation, and the audience sees the damage happen.

Where the numbers come from
---------------------------
Verbatim from the thesis:

    clean   Table "Mean absolute phase-plane displacement error [mm]"   (tab:h1-mae)
    noisy   its noisy counterpart                                       (tab:h2-mae)

Both are the sub-half-wavelength regime only.

A NOTE ON THE UNITS -- please read before defending this figure
---------------------------------------------------------------
The thesis's own dumbbell (H2_030) plots MAE as a percentage OF THE TRUE
DISPLACEMENT. That number is NOT reproducible from the tables: recomputing it
from the per-scenario supplementary tables does not return the reported mm MAE
either. For lateral/clean/back-propagation the per-scenario values give
18.62 mm over five scenarios or 23.27 mm over four, against 23.645 mm reported
-- so the exact scenario set and rounding behind the headline table are not
recoverable from the document.

Rather than reverse-engineer it and put an invented number on a defense slide,
this script does something that cannot be wrong: it plots the thesis's mm
values verbatim and expresses them as a percentage OF THE WAVELENGTH,

    error [% of λ]  =  MAE [mm] / 112.6 [mm] × 100

which is a pure unit conversion. It is also arguably the more useful framing
for this talk: "the error is 0.02% of a wavelength" speaks directly to the
resolution-floor argument the whole deck is built on.

If you want the thesis's own "% of true displacement" metric instead, fill in
PCT_TRUE below from your results notebook and run with --metric pct_true. The
script refuses to plot that metric until the table is complete -- it will not
guess.

Usage
-----
    python anim6_mae_dumbbell.py                 # both plots
    python anim6_mae_dumbbell.py --which clean
    python anim6_mae_dumbbell.py --metric mm     # raw millimetres instead
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from gpr_common import (apply_house_style, slide_title, footnote,
                        save_animation, save_poster, common_cli,
                        C_DARK, C_GREY, C_PANEL, C_GOOD, C_BAD, FIGSIZE)

LAMBDA_MM = 112.6            # thesis wavelength, 1.5 GHz in ice
FLOOR_MM = 5e-4              # log-axis floor for values reported as 0.000

MOVEMENTS = ["Lateral", "Vertical", "Diagonal", "FluidFlow"]
METHODS = ["Kirchhoff", "Gazdag", "Back-prop"]

# Colours chosen to match the thesis's own detectability / dumbbell figures,
# so the audience reads them as the same three methods.
MCOL = {"Kirchhoff": "#D9531E", "Gazdag": "#0E9C6E", "Back-prop": "#7B6FD0"}

# -- thesis tab:h1-mae  (clean, mm) ----------------------------------------
CLEAN_MM = {
    "Lateral":   {"Back-prop": 23.645, "Gazdag": 0.000,  "Kirchhoff": 0.027},
    "Vertical":  {"Back-prop": 11.950, "Gazdag": 16.079, "Kirchhoff": 16.799},
    "Diagonal":  {"Back-prop": 0.683,  "Gazdag": 0.003,  "Kirchhoff": 0.016},
    "FluidFlow": {"Back-prop": 0.967,  "Gazdag": 0.172,  "Kirchhoff": 0.335},
}
# -- thesis tab:h2-mae  (noisy, mm) ----------------------------------------
NOISY_MM = {
    "Lateral":   {"Back-prop": 14.381, "Gazdag": 22.254, "Kirchhoff": 0.236},
    "Vertical":  {"Back-prop": 13.410, "Gazdag": 34.476, "Kirchhoff": 17.485},
    "Diagonal":  {"Back-prop": 0.685,  "Gazdag": 23.308, "Kirchhoff": 0.065},
    "FluidFlow": {"Back-prop": 15.758, "Gazdag": 0.980,  "Kirchhoff": 0.234},
}
# Thesis reports these means for the noisy case only.
NOISY_MEAN_MM = {"Back-prop": 11.059, "Gazdag": 20.254, "Kirchhoff": 4.505}

# Fill this in ONLY from your own results if you want the thesis's own metric.
# Leave as None and --metric pct_true will refuse to run.
PCT_TRUE = None


def clean_mean_mm(method):
    """Derived, not quoted: the thesis prints means for the noisy case only."""
    return float(np.mean([CLEAN_MM[mv][method] for mv in MOVEMENTS]))


def convert(mm, metric):
    mm = max(mm, FLOOR_MM)
    if metric == "mm":
        return mm
    return mm / LAMBDA_MM * 100.0        # pct_lambda


AXIS_LABEL = {
    "mm": "mean absolute displacement error   [mm]",
    "pct_lambda": "mean absolute displacement error   [% of λ]",
}


# ==========================================================================
#  Row layout
# ==========================================================================

def build_rows(include_mean=True):
    """Returns list of (label, movement, method, y) top-to-bottom."""
    rows, y = [], 0.0
    groups = MOVEMENTS + (["MEAN"] if include_mean else [])
    for gi, mv in enumerate(groups):
        for method in METHODS:
            rows.append(dict(movement=mv, method=method, y=y))
            y += 1.0
        y += 0.7                                   # gap between groups
    return rows


def value_for(row, which, metric):
    mv, me = row["movement"], row["method"]
    if mv == "MEAN":
        mm = clean_mean_mm(me) if which == "clean" else NOISY_MEAN_MM[me]
    else:
        mm = (CLEAN_MM if which == "clean" else NOISY_MM)[mv][me]
    return convert(mm, metric), mm


# ==========================================================================
#  Shared canvas
# ==========================================================================

def label_side(ax, v, v_other=None):
    """
    Place the value label clear of BOTH ends of a dumbbell.

    Taking only one endpoint into account is not enough: when the clean and
    noisy markers are close together, a label offset from one of them lands on
    the other. Offset from the outer edge of the pair instead.
    """
    v_other = v if v_other is None else v_other
    lo, hi = np.log10(ax.get_xlim())
    right_end = max(v, v_other)
    left_end = min(v, v_other)
    if (np.log10(right_end) - lo) / (hi - lo) > 0.68:
        return left_end / 1.9, "right"
    return right_end * 1.28, "left"


def make_canvas(metric, title, subtitle, note):
    fig = plt.figure(figsize=FIGSIZE)
    slide_title(fig, title, subtitle)
    footnote(fig, note)

    ax = fig.add_axes([0.205, 0.125, 0.735, 0.585])
    ax.set_xscale("log")
    ax.set_facecolor(C_PANEL)
    ax.set_xlabel(AXIS_LABEL[metric])
    if metric == "pct_lambda":
        ax.set_xlim(3e-4, 1.2e2)
    else:
        ax.set_xlim(3e-4, 1.4e2)
    ax.grid(True, which="major", axis="x", color="white", lw=1.1)
    ax.grid(True, which="minor", axis="x", color="white", lw=0.5, alpha=0.6)

    rows = build_rows()
    ax.set_ylim(rows[-1]["y"] + 0.9, rows[0]["y"] - 0.9)
    ax.set_yticks([r["y"] for r in rows])
    ax.set_yticklabels([r["method"] for r in rows], fontsize=11)

    # group bands + labels
    for gi, mv in enumerate(MOVEMENTS + ["MEAN"]):
        sel = [r for r in rows if r["movement"] == mv]
        y0, y1 = sel[0]["y"] - 0.5, sel[-1]["y"] + 0.5
        if gi % 2 == 0:
            ax.axhspan(y0, y1, color="white", alpha=0.55, zorder=0)
        lbl = "MEAN\n(all four)" if mv == "MEAN" else mv
        ax.text(-0.235, (y0 + y1) / 2, lbl, transform=ax.get_yaxis_transform(),
                ha="left", va="center", fontsize=13, fontweight="bold",
                color=C_DARK)
    return fig, ax, rows


def annotate_zero_floor(ax, rows, which, metric):
    for r in rows:
        _, mm = value_for(r, which, metric)
        if mm < FLOOR_MM * 2:
            ax.text(convert(FLOOR_MM, metric) * 1.35, r["y"],
                    "reported 0.000 mm", fontsize=9.5, color=C_GREY,
                    va="center", ha="left", style="italic")


# ==========================================================================
#  Plot 1 — clean, static
# ==========================================================================

def render_clean(args):
    metric = args.metric
    fig, ax, rows = make_canvas(
        metric,
        "Clean synthetic data: where the phase fit lands",
        "Mean absolute error, sub-½λ regime  ·  thesis Table 5.4 (tab:h1-mae)",
        "λ = 112.6 mm  ·  MEAN row for clean data is derived from the four "
        "movement values")

    for r in rows:
        v, mm = value_for(r, "clean", metric)
        c = MCOL[r["method"]]
        ax.plot([ax.get_xlim()[0], v], [r["y"], r["y"]], color=c, lw=2.0,
                alpha=0.30, zorder=2)
        ax.plot([v], [r["y"]], "o", color=c, ms=13, zorder=5,
                markeredgecolor="white", markeredgewidth=1.4)
        if mm >= FLOOR_MM * 2:          # the floor annotation says it instead
            lx, ha = label_side(ax, v)
            ax.text(lx, r["y"], f"{mm:.3f} mm", fontsize=10.5,
                    color=C_DARK, va="center", ha=ha)
    annotate_zero_floor(ax, rows, "clean", metric)

    for m in METHODS:
        ax.plot([], [], "o", color=MCOL[m], ms=11, label=m)
    # legend above the panel: inside it would sit on the MEAN group
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.015), fontsize=12,
              ncol=3, borderaxespad=0.0)

    stem = os.path.join(args.outdir, "anim6_dumbbell_clean")
    out = stem + ".png"
    fig.savefig(out, dpi=300, facecolor="white")
    print(f"  -> {out}  (300 dpi still)")
    plt.close(fig)


# ==========================================================================
#  Plot 2 — noisy, animated out of the clean state
# ==========================================================================

def render_noisy(args):
    metric = args.metric
    fig, ax, rows = make_canvas(
        metric,
        "Adding noise: watch the accuracy degrade",
        "Laplace noise at 10% of each B-scan's own signal std  ·  "
        "thesis Tables 5.4 → 6.4",
        "λ = 112.6 mm  ·  bar colour: red = worse under noise, green = better")

    fps = args.fps
    n_start = int(2.2 * fps)                 # hold on the clean state
    n_move = int(9.0 * fps)
    n_end = int(args.hold * fps)

    art = {}
    for r in rows:
        c = MCOL[r["method"]]
        (bar,) = ax.plot([], [], lw=6.5, solid_capstyle="round", alpha=0.55,
                         zorder=3)
        (mk_clean,) = ax.plot([], [], "o", color="white", ms=13, zorder=5,
                              markeredgecolor=c, markeredgewidth=2.4)
        (mk_now,) = ax.plot([], [], "o", color=c, ms=13, zorder=6,
                            markeredgecolor="white", markeredgewidth=1.4)
        txt = ax.text(0, r["y"], "", fontsize=10.5, color=C_DARK,
                      va="center", ha="left")
        art[(r["movement"], r["method"])] = (bar, mk_clean, mk_now, txt)

    for m in METHODS:
        ax.plot([], [], "o", color=MCOL[m], ms=11, label=m)
    ax.plot([], [], "o", color="white", ms=11, markeredgecolor=C_GREY,
            markeredgewidth=2.0, label="clean")
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.015), fontsize=11.5,
              ncol=4, borderaxespad=0.0)

    # banner right-aligned so it shares the strip above the panel with the
    # legend instead of landing on top of it
    banner = ax.text(1.0, 1.02, "", transform=ax.transAxes, ha="right",
                     va="bottom", fontsize=16, fontweight="bold", color=C_DARK)

    def update(i):
        if i < n_start:
            u = 0.0
        elif i < n_start + n_move:
            p = (i - n_start) / n_move
            u = 0.5 - 0.5 * np.cos(np.pi * p)          # ease in/out
        else:
            u = 1.0

        for r in rows:
            v0, mm0 = value_for(r, "clean", metric)
            v1, mm1 = value_for(r, "noisy", metric)
            # travel geometrically: the axis is log, so interpolate in log
            v = 10 ** (np.log10(v0) + u * (np.log10(v1) - np.log10(v0)))
            mm = mm0 + u * (mm1 - mm0)
            worse = v1 > v0
            bar, mk_clean, mk_now, txt = art[(r["movement"], r["method"])]

            mk_clean.set_data([v0], [r["y"]])
            mk_now.set_data([v], [r["y"]])
            if u > 0.002:
                bar.set_data([v0, v], [r["y"], r["y"]])
                bar.set_color(C_BAD if worse else C_GOOD)
            else:
                bar.set_data([], [])
            lx, ha = label_side(ax, v, v0)
            txt.set_position((lx, r["y"]))
            txt.set_ha(ha)
            txt.set_text(f"{mm:.3f} mm")

        if u < 0.01:
            banner.set_text("CLEAN")
            banner.set_color(C_DARK)
        elif u < 0.99:
            banner.set_text(f"noise  →  {u*10:.1f}%  of signal std")
            banner.set_color(C_GREY)
        else:
            banner.set_text("NOISY   —   10% Laplace")
            banner.set_color(C_BAD)
        return []

    frames = n_start + n_move + n_end
    if args.frames:
        frames = args.frames
        n_start, n_move, n_end = frames // 5, frames // 2, frames // 5

    anim = FuncAnimation(fig, update, frames=frames,
                         interval=1000 // fps, blit=False, repeat=False)
    stem = os.path.join(args.outdir, "anim6_dumbbell_noisy")
    save_animation(anim, stem, fps=fps, dpi=args.dpi, prefer=args.format)
    if not args.no_poster:
        update(frames - 1)
        save_poster(fig, stem, dpi=300)
    plt.close(fig)


# ==========================================================================

def main():
    p = common_cli(__doc__.splitlines()[1])
    p.add_argument("--which", default="both",
                   choices=["clean", "noisy", "both"])
    p.add_argument("--metric", default="pct_lambda",
                   choices=["pct_lambda", "mm", "pct_true"])
    args = p.parse_args()

    if args.metric == "pct_true" and PCT_TRUE is None:
        raise SystemExit(
            "--metric pct_true needs the PCT_TRUE table filled in.\n"
            "It cannot be derived from the thesis document: recomputing MAE\n"
            "from the per-scenario supplementary tables does not reproduce\n"
            "the reported mm values, so the scenario set behind the headline\n"
            "table is not recoverable. Fill PCT_TRUE from your results\n"
            "notebook, or use the default --metric pct_lambda.")

    apply_house_style(base=14.5)
    os.makedirs(args.outdir, exist_ok=True)

    if args.which in ("clean", "both"):
        print("Rendering anim6 clean dumbbell  (static)")
        render_clean(args)
    if args.which in ("noisy", "both"):
        print("Rendering anim6 noisy dumbbell  (animated degradation)")
        render_noisy(args)


if __name__ == "__main__":
    main()
