"""
gpr_style.py
============
Shared visual style and export helpers for GPR animations across the project.

This module is deliberately physics-free: it only owns colours, matplotlib
rcParams, slide furniture (titles, readouts, verdict badges) and the
non-looping GIF/MP4 export machinery. `gpr_common.py` (Defense_Deck's
illustrative Gaussian-blob demos) and `gif_maker.ipynb` (the real-physics
Ricker/Kirchhoff/Gazdag/back-propagation animations) both import from here so
the two families of animation look like one system despite very different
underlying computations.

Author: shared style extracted for Conner's MSc defense deck + gif_maker.ipynb.
"""

from __future__ import annotations

import os
import sys

import matplotlib

matplotlib.use("Agg")  # headless-safe; no window is ever opened
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter

# --------------------------------------------------------------------------
# 1. House style  (clean minimal / white, per the deck design)
# --------------------------------------------------------------------------

C_DARK = "#16223A"    # near-black navy: text, axes
C_ACCENT = "#00748C"  # teal: baseline / "the method works" colour
C_WARM = "#C1440E"    # burnt orange: monitor / "the thing that moved"
C_GREY = "#8892A6"    # muted grey: ghosts, de-emphasised elements
C_GOOD = "#1B7F4B"    # green: RESOLVED badge
C_BAD = "#B3261E"     # red: MERGED / FAILED badge
C_PANEL = "#F4F6F9"   # very light grey: readout panel fills

CMAP_IMG = "RdBu_r"    # signed amplitude, matches the thesis figures
CMAP_PHASE = "RdBu_r"  # cross-spectrum phase ramp reads as a clean tilt

# 16:9 at PowerPoint's native widescreen size. Multi-panel wide layouts should
# default to this; a genuinely square/single-panel figure may deviate, but
# stay wide rather than tall wherever the panel count allows it.
FIGSIZE = (13.333, 7.5)


def apply_house_style(base: float = 15.0) -> None:
    """Set matplotlib rcParams for projector-legible, minimal figures."""
    plt.rcParams.update({
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.facecolor": "white",
        "font.family": "DejaVu Sans",
        "font.size": base,
        "axes.titlesize": base * 1.15,
        "axes.titleweight": "bold",
        "axes.labelsize": base * 0.95,
        "axes.edgecolor": C_DARK,
        "axes.labelcolor": C_DARK,
        "axes.linewidth": 1.1,
        "xtick.color": C_DARK,
        "ytick.color": C_DARK,
        "xtick.labelsize": base * 0.82,
        "ytick.labelsize": base * 0.82,
        "text.color": C_DARK,
        "legend.frameon": False,
        "legend.fontsize": base * 0.85,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.autolayout": False,
    })


# --------------------------------------------------------------------------
# 2. Slide furniture
# --------------------------------------------------------------------------

def readout(ax, text: str, loc: str = "upper left",
            color: str = C_DARK, size: float = 15, weight: str = "bold"):
    """A boxed numeric readout pinned inside an axes."""
    pos = {
        "upper left": (0.025, 0.975, "left", "top"),
        "upper right": (0.975, 0.975, "right", "top"),
        "lower left": (0.025, 0.03, "left", "bottom"),
        "lower right": (0.975, 0.03, "right", "bottom"),
    }[loc]
    return ax.text(pos[0], pos[1], text, transform=ax.transAxes,
                   ha=pos[2], va=pos[3], fontsize=size, color=color,
                   fontweight=weight, family="DejaVu Sans",
                   bbox=dict(boxstyle="round,pad=0.42", fc="white",
                             ec=color, lw=1.4, alpha=0.94))


def verdict_box(ax, loc: str = "upper left", size: float = 15):
    """
    A filled verdict box pinned INSIDE an axes.

    Placing it inside (rather than floating above the axes, where matplotlib
    will happily let it collide with a neighbouring panel's title or run off
    the figure edge) means it can never overlap anything, at any frame.

    Update it with  set_verdict(box, "MERGED", ok=False).
    """
    t = readout(ax, "", loc=loc, size=size, color=C_BAD)
    t.set_color("white")
    return t


def set_verdict(box, text: str, ok: bool):
    c = C_GOOD if ok else C_BAD
    box.set_text(text)
    patch = box.get_bbox_patch()
    patch.set_facecolor(c)
    patch.set_edgecolor(c)
    patch.set_alpha(1.0)


def badge(ax, text: str, ok: bool, loc=(0.5, 1.06), size: float = 17):
    """A RESOLVED / MERGED style verdict badge above an axes."""
    c = C_GOOD if ok else C_BAD
    return ax.text(loc[0], loc[1], text, transform=ax.transAxes,
                   ha="center", va="bottom", fontsize=size,
                   color="white", fontweight="bold",
                   bbox=dict(boxstyle="round,pad=0.40", fc=c, ec=c))


def slide_title(fig, title: str, subtitle: str = "", y: float = 0.965):
    """Title block matching the PowerPoint deck's typography."""
    fig.text(0.045, y, title, fontsize=27, fontweight="bold",
             color=C_DARK, ha="left", va="top")
    if subtitle:
        fig.text(0.045, y - 0.058, subtitle, fontsize=15.0,
                 color=C_GREY, ha="left", va="top")
    # thin accent rule, clear of the subtitle's descenders
    fig.lines.append(plt.Line2D([0.045, 0.955], [y - 0.108, y - 0.108],
                                transform=fig.transFigure,
                                color=C_ACCENT, lw=2.2))


def footnote(fig, text: str):
    fig.text(0.955, 0.022, text, fontsize=11, color=C_GREY,
             ha="right", va="bottom", style="italic")


# --------------------------------------------------------------------------
# 3. Export  --  non-looping by construction
# --------------------------------------------------------------------------

class PlayOncePillowWriter(PillowWriter):
    """
    PillowWriter that omits the GIF Netscape loop extension entirely, so the
    animation plays exactly once and freezes on the final frame.

    (matplotlib's stock PillowWriter hardcodes loop=0 = loop forever.)
    """

    def finish(self):
        if not self._frames:
            raise RuntimeError(
                "No frames were captured — the update() function raised. "
                "Scroll up: the real traceback is printed above this one.")
        self._frames[0].save(
            self.outfile,
            save_all=True,
            append_images=self._frames[1:],
            duration=int(1000.0 / self.fps),
            # no `loop=` kwarg -> play once
        )


def ffmpeg_available() -> bool:
    try:
        return animation.writers.is_available("ffmpeg")
    except Exception:
        return False


def save_animation(anim, outstem: str, fps: int = 25, dpi: int = 200,
                   prefer: str = "auto") -> str:
    """
    Write the animation to MP4 (preferred, PowerPoint-native) or GIF.

    Both formats are non-looping: MP4 because PowerPoint does not loop video
    by default, GIF because of PlayOncePillowWriter above.
    """
    os.makedirs(os.path.dirname(os.path.abspath(outstem)) or ".", exist_ok=True)

    use_mp4 = (prefer == "mp4") or (prefer == "auto" and ffmpeg_available())
    if use_mp4:
        out = outstem + ".mp4"
        try:
            writer = animation.FFMpegWriter(
                fps=fps, bitrate=-1,
                extra_args=["-pix_fmt", "yuv420p", "-vcodec", "libx264",
                            "-preset", "slow", "-crf", "18"])
            anim.save(out, writer=writer, dpi=dpi)
            print(f"  -> {out}  ({fps} fps, {dpi} dpi)")
            return out
        except Exception as exc:                       # pragma: no cover
            print(f"  !! ffmpeg failed ({exc}); falling back to GIF",
                  file=sys.stderr)

    out = outstem + ".gif"
    anim.save(out, writer=PlayOncePillowWriter(fps=fps), dpi=min(dpi, 120))
    print(f"  -> {out}  ({fps} fps, plays once)")
    return out


def save_poster(fig, outstem: str, dpi: int = 300) -> str:
    """Write the final frame as a 300-dpi PNG, for printing or a static slide."""
    out = outstem + "_final.png"
    fig.savefig(out, dpi=dpi, facecolor="white", bbox_inches=None)
    print(f"  -> {out}  ({dpi} dpi still)")
    return out


def hold(frames, n_hold: int):
    """Append n_hold repeats of the last value so the animation freezes on it."""
    import numpy as np
    if n_hold <= 0:
        return frames
    return np.concatenate([frames, np.repeat(frames[-1], n_hold)])
