"""
anim5_crossspectrum_derivation.py   --   SLIDE 9
================================================
The cross-spectrum derivation, as a four-act build.

Two panels, side by side, that change what they REPRESENT four times while
staying in the same place on screen. That is the whole trick: the audience
never loses track of "left = baseline, right = monitor", so each new domain
lands as a re-reading of the same two objects rather than a new figure.

  ACT 1   model        left: static scatterer at the baseline position
                       right: the same scatterer moving to its monitor position
                       →  m(x,z) = b(x-Δx, z-Δz)

  ACT 2   migrated     both panels dissolve into their MIGRATED counterparts.
                       The sharp point becomes a broad lobe -- and the relation
                       above still holds, because migration is linear and both
                       surveys use the same operator.

  ACT 3   spectrum     both panels dissolve into their 2-D wavenumber spectra.
                       →  B(kx,kz) = ∫∫ b e^{-j(kx x + kz z)} dx dz
                       →  M(kx,kz) = B(kx,kz) · e^{-j(kx Δx + kz Δz)}
                       The payoff: |M| and |B| are IDENTICAL. A displacement is
                       completely invisible in the amplitude spectrum.

  ACT 4   cross        the pair is replaced by the cross-spectrum's amplitude
                       and phase.
                       →  XS = B · M*
                       →  XS = |B|² · e^{+j(kx Δx + kz Δz)}
                       →  ∠XS = Φ(kx,kz) = kx Δx + kz Δz
                       |XS| still says nothing. The phase is a clean plane, and
                       its two slopes ARE the displacement.

Working example: Δx = ¼λ, Δz = ⅛λ. Diagonal on purpose, so the phase plane
tilts along both axes at once and Φ = kx Δx + kz Δz has two live terms rather
than one. Total phase excursion stays under π, so nothing wraps.

Output
------
By default this writes FIVE files:

    anim5_act1_model.mp4        show your formula, then play this
    anim5_act2_migrated.mp4     ... and so on
    anim5_act3_spectrum.mp4
    anim5_act4_cross.mp4
    anim5_full.mp4              all four, continuous, for rehearsal

Each standalone act opens by holding the PREVIOUS act's final state for a
moment, so clicking from one clip to the next is seamless -- no jump cut.

Usage
-----
    python anim5_crossspectrum_derivation.py                 # all five
    python anim5_crossspectrum_derivation.py --acts 3        # just act 3
    python anim5_crossspectrum_derivation.py --no-equations  # panels only
    python anim5_crossspectrum_derivation.py --frames 40     # layout preview

--no-equations is the one to use if you are revealing the formulas as
PowerPoint text yourself and do not want them duplicated inside the video.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from gpr_common import (LAMBDA, ImageGrid, apply_house_style, cross_spectrum,
                        slide_title, save_animation, save_poster, common_cli,
                        C_DARK, C_ACCENT, C_WARM, C_GREY, C_PANEL, CMAP_IMG,
                        FIGSIZE)

# --------------------------------------------------------------------------
# The worked example
# --------------------------------------------------------------------------
DX = 0.250 * LAMBDA          # lateral displacement
DZ = 0.125 * LAMBDA          # vertical displacement
K_VIEW = 2.2                 # k-axis limit, in units of k_c
CMAP_MAG = "magma"           # sequential, for |spectrum| and |XS|
CMAP_PHI = "RdBu_r"          # diverging, for the phase plane

# Magnitude panels are gamma-compressed for display. |B| already spans several
# decades, and |XS| = |B|^2 squares that: on a linear scale the panel is a
# black rectangle with one bright dot, which teaches nothing. The colourbar is
# labelled "normalised" so the compression is stated rather than smuggled in.
GAMMA = 0.45


def _mag_display(a):
    return (a / (a.max() + 1e-30)) ** GAMMA

# --------------------------------------------------------------------------
# Act structure, in seconds (scaled by --duration / 26)
# --------------------------------------------------------------------------
ACTS = [
    dict(key="model",    dissolve=0.0, motion=6.0, hold=2.0),
    dict(key="migrated", dissolve=1.6, motion=0.0, hold=3.4),
    dict(key="spectrum", dissolve=1.6, motion=0.0, hold=4.4),
    dict(key="cross",    dissolve=1.6, motion=0.0, hold=5.4),
]
BASE_TOTAL = sum(a["dissolve"] + a["motion"] + a["hold"] for a in ACTS)  # 26 s
LEAD_IN = 0.9                # seconds of previous state at the top of a clip

# --------------------------------------------------------------------------
# Equations, revealed progressively within their act
# --------------------------------------------------------------------------
EQ = {
    "model": [
        r"$m(x,z)\;=\;b\,(\,x-\Delta x,\;\; z-\Delta z\,)$",
    ],
    "migrated": [
        r"$m(x,z)\;=\;b\,(\,x-\Delta x,\;\; z-\Delta z\,)$"
        r"$\quad$ — still true after migration",
    ],
    "spectrum": [
        r"$B(k_x,k_z)\;=\;\int_{-\infty}^{\infty}\!\!\int_{-\infty}^{\infty}"
        r"b(x,z)\;e^{-j(k_x x\,+\,k_z z)}\;dx\,dz$",
        r"$M(k_x,k_z)\;=\;B(k_x,k_z)\;\cdot\;"
        r"e^{-j(k_x \Delta x\,+\,k_z \Delta z)}$",
    ],
    "cross": [
        r"$XS(k_x,k_z)\;=\;B(k_x,k_z)\;\cdot\;M^{*}(k_x,k_z)$",
        r"$XS(k_x,k_z)\;=\;|B(k_x,k_z)|^{2}\;\cdot\;"
        r"e^{+j(k_x \Delta x\,+\,k_z \Delta z)}$",
        r"$\angle\,XS(k_x,k_z)\;=\;\Phi(k_x,k_z)\;=\;"
        r"k_x \Delta x\,+\,k_z \Delta z$",
    ],
}

PANEL_TITLES = {
    "model":    ("BASELINE   $b(x,z)$", "MONITOR   $m(x,z)$"),
    "migrated": ("migrated baseline", "migrated monitor"),
    "spectrum": ("$|B(k_x,k_z)|$", "$|M(k_x,k_z)|$"),
    "cross":    ("$|XS|$   —   amplitude", "$\\angle XS$   —   phase"),
}

# The one-line takeaway stamped between the panels, per act
VERDICTS = {
    "model":    ("", ""),
    "migrated": ("", ""),
    "spectrum": ("|M| = |B|", "the shift is invisible in amplitude"),
    "cross":    ("|XS| = |B|²", "still nothing …  but ∠XS is a plane"),
}


# ==========================================================================
#  Physics: every representation, precomputed once
# ==========================================================================

class Derivation:
    def __init__(self):
        self.grid = ImageGrid(n=256, lam=LAMBDA, samples_per_lambda=16)
        g = self.grid
        self.S = g.psf_spectrum(aperture_deg=50.0, bw_frac=0.85)
        self.kc = 2.0 * np.pi / LAMBDA

        # --- true reflectivity model: a compact blob at each position
        self.sigma = 0.045 * LAMBDA
        self.model_b = self._blob(0.0, 0.0)

        # --- migrated images
        self.mig_b = g.psf_image(self.S, 0.0, 0.0)
        self.mig_m = g.psf_image(self.S, DX, DZ)
        self.mig_vmax = np.abs(self.mig_b).max()

        # --- spectra
        B = g.spectrum_of(self.mig_b)
        M = g.spectrum_of(self.mig_m)
        self.absB = np.abs(B)
        self.absM = np.abs(M)
        self.dispB = _mag_display(self.absB)
        self.dispM = _mag_display(self.absM)

        # --- cross-spectrum (tapered, exactly as the estimator does it)
        XS = cross_spectrum(self.mig_b, self.mig_m, g)
        self.absXS = np.abs(XS)
        self.phiXS = np.angle(XS)
        # Blank the phase where there is no coherent energy: outside the
        # wavelet's own band the "phase" is numerical noise and showing it
        # would bury the ramp we are trying to point at.
        self.xs_mask = self.absXS > 0.03 * self.absXS.max()
        self.dispXS = _mag_display(self.absXS)

        self.extent_space = g.extent
        self.extent_k = (g.kx[0] / self.kc, g.kx[-1] / self.kc,
                         g.kz[0] / self.kc, g.kz[-1] / self.kc)

    def _blob(self, dx, dz):
        g = self.grid
        X, Z = np.meshgrid(g.x, g.z)
        return np.exp(-(((X - dx) ** 2 + (Z - dz) ** 2) / (2 * self.sigma ** 2)))

    def model_monitor(self, u: float):
        """Monitor model with the scatterer u of the way to (DX, DZ)."""
        return self._blob(u * DX, u * DZ)

    # -- what each panel shows, per act ------------------------------------

    def panel_data(self, state, side, u=1.0):
        """
        Returns (data, cmap, vmin, vmax, extent, cbar_label, cbar_ticks,
                 xlabel, ylabel, alpha_mask)
        """
        if state == "model":
            d = self.model_b if side == 0 else self.model_monitor(u)
            return dict(data=d, cmap="Greys", vmin=0, vmax=1,
                        extent=self.extent_space, clabel="reflectivity",
                        cticks=None, xlabel="x   [m]", ylabel="z   [m]",
                        origin="upper", mask=None)

        if state == "migrated":
            d = self.mig_b if side == 0 else self.mig_m
            v = self.mig_vmax
            return dict(data=d, cmap=CMAP_IMG, vmin=-v, vmax=v,
                        extent=self.extent_space, clabel="amplitude",
                        cticks=None, xlabel="x   [m]", ylabel="z   [m]",
                        origin="upper", mask=None)

        if state == "spectrum":
            d = self.dispB if side == 0 else self.dispM
            return dict(data=d, cmap=CMAP_MAG, vmin=0, vmax=1,
                        extent=self.extent_k, clabel="magnitude (normalised)",
                        cticks=None, xlabel="$k_x$   [$k_c$]",
                        ylabel="$k_z$   [$k_c$]", origin="lower", mask=None)

        # cross
        if side == 0:
            return dict(data=self.dispXS, cmap=CMAP_MAG, vmin=0,
                        vmax=1, extent=self.extent_k,
                        clabel="$|XS|$  (normalised)", cticks=None,
                        xlabel="$k_x$   [$k_c$]", ylabel="$k_z$   [$k_c$]",
                        origin="lower", mask=None)
        return dict(data=self.phiXS, cmap=CMAP_PHI, vmin=-np.pi, vmax=np.pi,
                    extent=self.extent_k, clabel="phase   [rad]",
                    cticks=([-np.pi, 0, np.pi], [r"$-\pi$", "0", r"$\pi$"]),
                    xlabel="$k_x$   [$k_c$]", ylabel="$k_z$   [$k_c$]",
                    origin="lower", mask=self.xs_mask)


# ==========================================================================
#  Timeline
# ==========================================================================

def build_timeline(sel, fps, scale):
    """
    sel : list of act indices (0-based) to include, in order.

    Each frame is dict(state, prev, alpha, u, prog) where `prog` is progress
    through the act's hold phase, used to reveal equations one at a time.
    """
    frames = []

    # A clip that does not start at act 1 opens by holding the previous
    # act's final state, so clicking between clips never jump-cuts.
    if sel[0] > 0:
        prev = ACTS[sel[0] - 1]["key"]
        for _ in range(max(1, int(LEAD_IN * fps * scale))):
            frames.append(dict(state=prev, alpha=1.0, u=1.0, prog=1.0))

    for i in sel:
        a = ACTS[i]
        prev = ACTS[i - 1]["key"] if i > 0 else None

        n_dis = int(a["dissolve"] * fps * scale)
        for f in range(n_dis):
            p = (f + 1) / n_dis
            if p < 0.5:                       # fade the old representation out
                frames.append(dict(state=prev, alpha=1.0 - 2 * p,
                                   u=1.0, prog=0.0))
            else:                             # fade the new one in
                frames.append(dict(state=a["key"], alpha=2 * p - 1.0,
                                   u=1.0, prog=0.0))

        n_mot = int(a["motion"] * fps * scale)
        for f in range(n_mot):
            u = (f + 1) / n_mot
            # ease-in-out so the scatterer starts and stops gently
            u = 0.5 - 0.5 * np.cos(np.pi * u)
            frames.append(dict(state=a["key"], alpha=1.0, u=u, prog=0.0))

        n_hold = int(a["hold"] * fps * scale)
        for f in range(n_hold):
            frames.append(dict(state=a["key"], alpha=1.0, u=1.0,
                               prog=(f + 1) / n_hold))

    return frames


def n_equations_visible(state, prog):
    """Reveal an act's equations one at a time across its hold phase."""
    eqs = EQ.get(state, [])
    if len(eqs) <= 1:
        return len(eqs)
    # first equation immediately, the rest spread over the first 70% of hold
    step = 0.70 / (len(eqs) - 1)
    return min(len(eqs), 1 + int(prog / step + 1e-9))


# ==========================================================================
#  Render
# ==========================================================================

def render(sel, stem, args, deriv):
    fps = args.fps
    scale = args.duration / BASE_TOTAL if args.duration else 1.0
    frames = build_timeline(sel, fps, scale)
    if args.frames:                       # preview: uniformly thin the timeline
        idx = np.linspace(0, len(frames) - 1, args.frames).astype(int)
        frames = [frames[i] for i in idx]
    # freeze on the last frame
    frames += [dict(frames[-1])] * max(1, int(args.hold * fps))

    fig = plt.figure(figsize=FIGSIZE)
    slide_title(fig,
                "From a moved scatterer to a phase plane",
                f"Worked example:  Δx = ¼λ,  Δz = ⅛λ   ·   λ = {LAMBDA:.0f} m")

    # Vertical budget, from the top: title block, then an equation band deep
    # enough for act 4's three stacked lines, then the panels, then the verdict
    # strip. Act 4 is the constraint -- everything is sized so IT fits.
    show_eq = not args.no_equations
    eq_top = 0.790 if show_eq else 0.80
    eq_step = 0.055
    panel_top = 0.575 if show_eq else 0.730

    gs = fig.add_gridspec(1, 2, left=0.075, right=0.885, top=panel_top,
                          bottom=0.170, wspace=0.42)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]

    ims, cbs, caxes, titles = [], [], [], []
    for k, ax in enumerate(axes):
        spec = deriv.panel_data("model", k)
        im = ax.imshow(spec["data"], cmap=spec["cmap"], vmin=spec["vmin"],
                       vmax=spec["vmax"], extent=spec["extent"],
                       origin=spec["origin"], aspect="auto",
                       interpolation="bilinear")
        ax.set_facecolor("white")
        ims.append(im)
        titles.append(ax.set_title("", color=C_DARK, pad=13, fontsize=16))

        box = ax.get_position()
        cax = fig.add_axes([box.x1 + 0.014, box.y0, 0.013, box.height])
        caxes.append(cax)
        cbs.append(fig.colorbar(im, cax=cax))
        cbs[-1].ax.tick_params(labelsize=10.5)

    # zoom windows
    zoom_space = 1.15 * LAMBDA

    # ghost trail of the positions the scatterer has passed through
    trail = axes[1].scatter([], [], s=26, facecolors="none",
                            edgecolors=C_WARM, linewidths=1.4, zorder=6)
    (trail_head,) = axes[1].plot([], [], "o", color=C_WARM, ms=9, zorder=7)
    arrow = axes[1].annotate("", xy=(0, 0), xytext=(0, 0),
                             arrowprops=dict(arrowstyle="-|>", lw=2.4,
                                             color=C_WARM), zorder=6)

    # equation band
    eq_texts = []
    if show_eq:
        for i in range(3):
            eq_texts.append(fig.text(0.5, eq_top - i * eq_step, "", ha="center",
                                     va="top", fontsize=15, color=C_DARK))

    # verdict stamp, centred between the panels
    verdict_a = fig.text(0.5, 0.070, "", ha="center", va="bottom",
                         fontsize=19, fontweight="bold", color=C_ACCENT)
    verdict_b = fig.text(0.5, 0.026, "", ha="center", va="bottom",
                         fontsize=13.5, color=C_GREY, style="italic")

    current = {"state": None}

    def apply_state(state, u):
        for k, ax in enumerate(axes):
            spec = deriv.panel_data(state, k, u)
            im = ims[k]
            im.set_data(spec["data"])
            im.set_cmap(spec["cmap"])
            im.set_clim(spec["vmin"], spec["vmax"])
            im.set_extent(spec["extent"])
            im.set_norm(im.norm)                # keep the colorbar in step

            if spec["mask"] is not None:
                im.set_alpha(np.where(spec["mask"], 1.0, 0.03))
            else:
                im.set_alpha(None)

            if state in ("model", "migrated"):
                ax.set_xlim(-zoom_space, zoom_space)
                ax.set_ylim(zoom_space, -zoom_space)
                ax.set_xticks([-1, 0, 1])
                ax.set_yticks([-1, 0, 1])
            else:
                ax.set_xlim(-K_VIEW, K_VIEW)
                ax.set_ylim(-K_VIEW, K_VIEW)
                ax.set_xticks([-2, -1, 0, 1, 2])
                ax.set_yticks([-2, -1, 0, 1, 2])

            ax.set_xlabel(spec["xlabel"])
            ax.set_ylabel(spec["ylabel"])
            titles[k].set_text(PANEL_TITLES[state][k])

            cbs[k].update_normal(im)
            cbs[k].set_label(spec["clabel"], fontsize=11.5)
            if spec["cticks"]:
                cbs[k].set_ticks(spec["cticks"][0])
                cbs[k].set_ticklabels(spec["cticks"][1])
            cbs[k].ax.tick_params(labelsize=10.5)

        current["state"] = state

    def update(i):
        fr = frames[i]
        state, alpha, u, prog = fr["state"], fr["alpha"], fr["u"], fr["prog"]

        if state != current["state"] or state == "model":
            apply_state(state, u)

        # dissolve: scale the whole panel's opacity
        for k, im in enumerate(ims):
            spec = deriv.panel_data(state, k, u)
            if spec["mask"] is not None:
                im.set_alpha(np.where(spec["mask"], alpha, 0.03 * alpha))
            else:
                im.set_alpha(alpha)
            caxes[k].set_alpha(alpha)

        # act 1 only: the moving scatterer, its trail and its arrow
        showing_model = (state == "model")
        if showing_model and u > 0.02:
            n = max(2, int(u * 7))
            us = np.linspace(0, u, n)[:-1]
            trail.set_offsets(np.column_stack([us * DX, us * DZ]))
            trail_head.set_data([u * DX], [u * DZ])
            arrow.xy = (u * DX, u * DZ)
            arrow.set_visible(alpha > 0.4)
        else:
            trail.set_offsets(np.empty((0, 2)))
            trail_head.set_data([], [])
            arrow.set_visible(False)
        trail.set_alpha(alpha if showing_model else 0.0)
        trail_head.set_alpha(alpha if showing_model else 0.0)

        # equations
        if show_eq:
            eqs = EQ.get(state, [])
            n_vis = n_equations_visible(state, prog)
            for j, t in enumerate(eq_texts):
                if j < n_vis and j < len(eqs):
                    t.set_text(eqs[j])
                    # newest equation in accent, settled ones in dark
                    t.set_color(C_ACCENT if j == n_vis - 1 else C_DARK)
                    t.set_alpha(alpha)
                else:
                    t.set_text("")

        va, vb = VERDICTS[state]
        verdict_a.set_text(va if prog > 0.25 else "")
        verdict_b.set_text(vb if prog > 0.25 else "")

        return ims + [trail, trail_head, verdict_a, verdict_b] + eq_texts

    anim = FuncAnimation(fig, update, frames=len(frames),
                         interval=1000 // fps, blit=False, repeat=False)

    save_animation(anim, stem, fps=fps, dpi=args.dpi, prefer=args.format)
    if not args.no_poster:
        update(len(frames) - 1)
        save_poster(fig, stem, dpi=300)
    plt.close(fig)


# ==========================================================================

def main():
    p = common_cli(__doc__.splitlines()[1])
    p.add_argument("--acts", default="all",
                   help="'all' (five files), or a comma list like '1,3'")
    p.add_argument("--no-equations", action="store_true",
                   help="omit the formulas from the video — use this if you "
                        "are revealing them as PowerPoint text yourself")
    p.set_defaults(duration=BASE_TOTAL)
    args = p.parse_args()

    apply_house_style(base=14.5)
    os.makedirs(args.outdir, exist_ok=True)
    deriv = Derivation()

    names = ["act1_model", "act2_migrated", "act3_spectrum", "act4_cross"]

    if args.acts == "all":
        wanted = [0, 1, 2, 3]
        make_full = True
    else:
        wanted = [int(x) - 1 for x in args.acts.split(",")]
        make_full = False

    for i in wanted:
        stem = os.path.join(args.outdir, f"anim5_{names[i]}")
        print(f"Rendering anim5 act {i+1}  ({ACTS[i]['key']})")
        render([i], stem, args, deriv)

    if make_full:
        stem = os.path.join(args.outdir, "anim5_full")
        print("Rendering anim5 full  (all four acts, continuous)")
        render([0, 1, 2, 3], stem, args, deriv)


if __name__ == "__main__":
    main()
