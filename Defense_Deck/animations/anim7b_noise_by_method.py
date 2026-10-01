"""
anim7b_noise_by_method.py   --   SLIDE 14b  (companion to anim7)
==================================================================
anim7 asks "does noise break the phase fit?" for one migration method.
This asks the question the thesis actually answers differently per method:
does noise break the fit the SAME WAY under Kirchhoff, Gazdag, and
back-propagation?

Same B-scan geometry, same true displacement, same single noise
realisation, same weighted phase-plane fit as anim7 -- the only thing that
changes across the three columns is the migration algorithm standing
between the noisy data and the fit.

  (a)/(b)/(c)  migrated image, one column per method, noise ramping
               together on all three
  (d)          recovered-displacement error vs. noise level, one line per
               method, sharing the axes so the comparison is direct

Honest construction / what this does NOT show
-----------------------------------------------
Kirchhoff (delay-and-sum) and Gazdag (f-k phase-shift) are both LINEAR in
the B-scan, so -- exactly as in anim7 -- the clean and noise images are
migrated once each and blended per frame.

Back-propagation here uses DIRECT, unnormalised injection of the B-scan
into the reverse-time wave equation (matching gif_maker.ipynb's own
back-propagation cell). That keeps it linear too, so it can share the same
precompute-once trick -- four migrations total per method, not one per
frame, which is what makes hundreds of frames of three real migrations
affordable.

The thesis's own back-propagation is NOISE-ROBUST specifically because it
uses SIGN-BIT excitation, a nonlinear clip that this animation does not
reproduce (a sign-bit ramp would need a full reverse-time simulation at
every noise level -- expensive, and a separate, narrower claim). So do not
read this animation as "here is why back-propagation wins" -- for that,
anim6_mae_dumbbell.py plots the thesis's actual sign-bit numbers. Read this
one as "here is how three migration algorithms, run identically, respond
to the same growing noise on the same input" -- which is the fair,
like-for-like comparison the algorithms alone support.

Usage
-----
    python anim7b_noise_by_method.py
    python anim7b_noise_by_method.py --max-noise 0.10   # strict thesis level
    python anim7b_noise_by_method.py --frames 30         # layout preview
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.interpolate import RegularGridInterpolator

from gpr_common import (LAMBDA, V, apply_house_style, fit_phase_plane,
                        readout, slide_title, footnote, save_animation,
                        save_poster, common_cli, frame_counts, C_DARK,
                        C_GREY, C_PANEL, FIGSIZE)

from anim7_noise_propagation import (Pipeline as KirchhoffPipeline, DEPTH,
                                     DX_TRUE, THESIS_LEVEL,
                                     MAX_NOISE_DEFAULT)

METHODS = ["Kirchhoff", "Gazdag", "Back-prop"]
# same hues as anim6_mae_dumbbell.py so the audience reads them as the same
# three methods across the whole deck
MCOL = {"Kirchhoff": "#D9531E", "Gazdag": "#0E9C6E", "Back-prop": "#7B6FD0"}

V_MIG = V / 2.0                    # exploding-reflector velocity (Gazdag + back-prop)


# ==========================================================================
#  Pipeline3: Kirchhoff (inherited) + Gazdag + back-propagation, all on the
#  SAME output grid, all evaluated by the SAME cross-spectrum + fit.
# ==========================================================================

class Pipeline3(KirchhoffPipeline):
    def __init__(self, max_noise):
        super().__init__(max_noise)             # bs_b, bs_m, nz_b, nz_m, Kirchhoff op
        self.dx_tr = self.x_tr[1] - self.x_tr[0]

        print("  migrating (Gazdag)...")
        self.mig_b_clean_gz = self._migrate_gazdag(self.bs_b)
        self.mig_m_clean_gz = self._migrate_gazdag(self.bs_m)
        self.mig_b_noise_gz = self._migrate_gazdag(self.nz_b)
        self.mig_m_noise_gz = self._migrate_gazdag(self.nz_m)

        print("  migrating (back-propagation)...")
        self._setup_backprop_grid()
        self.mig_b_clean_bp = self._migrate_backprop(self.bs_b)
        self.mig_m_clean_bp = self._migrate_backprop(self.bs_m)
        self.mig_b_noise_bp = self._migrate_backprop(self.nz_b)
        self.mig_m_noise_bp = self._migrate_backprop(self.nz_m)

        self.vmax_gz = np.abs(self.mig_b_clean_gz).max()
        self.vmax_bp = np.abs(self.mig_b_clean_bp).max()

    # -- Gazdag: f-k phase-shift, adapted from gif_maker.ipynb's verified
    #    cell (anti-wrap padding, exploding-reflector kz, t=0 imaging
    #    condition). Depths are stepped directly at self.grid's own z
    #    positions and the native-spacing output row is 1-D interpolated
    #    onto self.grid.x, so the result lands pixel-for-pixel on the same
    #    grid Kirchhoff already uses -- no need for a separate resample step
    #    downstream. --------------------------------------------------------

    def _migrate_gazdag(self, bscan):
        nt, ntr = bscan.shape
        pad_t = nt
        pad_x = ntr // 2
        b_pad = np.pad(bscan, ((0, pad_t), (pad_x, pad_x)), mode="constant")

        P_fk = np.fft.fft2(b_pad)
        freqs = np.fft.fftfreq(nt + pad_t, d=self.dt)
        kxs = np.fft.fftfreq(ntr + 2 * pad_x, d=self.dx_tr)
        KX, F = np.meshgrid(kxs, freqs)
        omega = 2 * np.pi * F
        kx = 2 * np.pi * KX

        kz_sq = (omega / V_MIG) ** 2 - kx ** 2
        kz = np.where(kz_sq >= 0, np.sqrt(np.maximum(kz_sq, 0)), 0.0)
        kz = np.sign(omega) * kz

        z_targets = DEPTH + self.grid.z
        img = np.zeros((self.grid.n, self.grid.n))
        for j, z in enumerate(z_targets):
            shifted = P_fk * np.exp(1j * kz * z)
            xt_pad = np.real(np.fft.ifft2(shifted))
            row_native = xt_pad[0, pad_x:pad_x + ntr]         # t=0 slice, native x
            img[j, :] = np.interp(self.grid.x, self.x_tr, row_native)
        return img

    # -- Back-propagation: raw (unnormalised) FDTD time-reversal, adapted
    #    from gif_maker.ipynb's verified cell. Direct injection keeps the
    #    wave equation linear in the B-scan, so it shares Kirchhoff/Gazdag's
    #    precompute-once trick -- the honest tradeoff being that this is
    #    NOT the thesis's sign-bit excitation (see module docstring).
    #    The FDTD grid runs on the antenna's own x-spacing; only the small
    #    window around the target is cropped and bilinear-interpolated onto
    #    self.grid for a like-for-like comparison with the other two. -------

    def _setup_backprop_grid(self):
        self._bp_dx = self.dx_tr
        z_max = DEPTH + 1.5 * LAMBDA * 1.3
        self._bp_nz = int(z_max / self._bp_dx) + 1
        self._bp_z = np.arange(self._bp_nz) * self._bp_dx
        self._bp_x = self.x_tr

        C = V_MIG * self.dt / self._bp_dx
        self._bp_C2 = C ** 2

    def _migrate_backprop(self, bscan):
        nz, nx = self._bp_nz, len(self._bp_x)
        nt = bscan.shape[0]
        P_past = np.zeros((nz, nx))
        P_now = np.zeros((nz, nx))
        C2 = self._bp_C2

        for step in range(nt - 1, -1, -1):
            lap = (P_now[:-2, 1:-1] + P_now[2:, 1:-1] +
                  P_now[1:-1, :-2] + P_now[1:-1, 2:] -
                  4 * P_now[1:-1, 1:-1])
            P_next = np.zeros_like(P_now)
            P_next[1:-1, 1:-1] = 2 * P_now[1:-1, 1:-1] - P_past[1:-1, 1:-1] + C2 * lap
            P_next[0, :] = bscan[step, :]
            P_past, P_now = P_now, P_next

        interp = RegularGridInterpolator((self._bp_z, self._bp_x), P_now,
                                         bounds_error=False, fill_value=0.0)
        Zt, Xt = np.meshgrid(DEPTH + self.grid.z, self.grid.x, indexing="ij")
        return interp((Zt, Xt))

    # -- per-frame, one blended image + one fit per method ------------------

    def frame(self, frac):
        out = {}
        specs = {
            "Kirchhoff": (self.mig_b_clean, self.mig_b_noise,
                         self.mig_m_clean, self.mig_m_noise),
            "Gazdag": (self.mig_b_clean_gz, self.mig_b_noise_gz,
                      self.mig_m_clean_gz, self.mig_m_noise_gz),
            "Back-prop": (self.mig_b_clean_bp, self.mig_b_noise_bp,
                         self.mig_m_clean_bp, self.mig_m_noise_bp),
        }
        for name, (cb, nb, cm, nm) in specs.items():
            mig_b = cb + frac * nb
            mig_m = cm + frac * nm
            # zero-padded spectra (self.gk/self._pad, inherited from
            # KirchhoffPipeline) -- without this the 128-sample image alone
            # gives dk too coarse for a good phase-plane fit, exactly the
            # problem anim7's own docstring flags and works around
            B = self.gk.spectrum_of(self._pad(mig_b))
            M = self.gk.spectrum_of(self._pad(mig_m))
            XS = B * np.conj(M)
            est_dz, est_dx, _, mask = fit_phase_plane(XS, self.gk)
            out[name] = (mig_b, est_dx)
        return out


# ==========================================================================

def main():
    p = common_cli(__doc__.splitlines()[1])
    p.add_argument("--max-noise", type=float, default=MAX_NOISE_DEFAULT,
                   help="peak noise level, as a fraction of signal std")
    args = p.parse_args()

    apply_house_style(base=14)
    os.makedirs(args.outdir, exist_ok=True)

    print("Building forward model and all three migration operators ...")
    pipe = Pipeline3(args.max_noise)
    g = pipe.grid

    n_frames, n_hold = frame_counts(args)
    n_lead = max(1, n_frames // 8)
    fracs = np.concatenate([
        np.zeros(n_lead),
        np.linspace(0.0, args.max_noise, n_frames),
        np.full(n_hold, args.max_noise)])

    # ---------------- figure ----------------
    fig = plt.figure(figsize=FIGSIZE)
    slide_title(fig,
               "Three migration methods, the same growing noise",
               f"Noise ramped 0 → {args.max_noise*100:.0f}% of signal std  ·  "
               f"true Δx = {DX_TRUE/LAMBDA:.2f} λ  ·  same noise draw, same fit")
    footnote(fig, "back-prop uses direct (non-sign-bit) injection here -- "
                  "see the module docstring before citing thesis noise-robustness")

    gs = fig.add_gridspec(2, 3, left=0.055, right=0.965, top=0.72,
                          bottom=0.10, wspace=0.32, hspace=0.62,
                          height_ratios=[1.0, 0.85])
    axes_img = [fig.add_subplot(gs[0, i]) for i in range(3)]
    ax_err = fig.add_subplot(gs[1, :])

    ims = {}
    txts = {}
    vmax_map = {"Kirchhoff": pipe.mig_vmax, "Gazdag": pipe.vmax_gz,
               "Back-prop": pipe.vmax_bp}
    for ax, name in zip(axes_img, METHODS):
        vmax = vmax_map[name]
        im = ax.imshow(np.zeros((g.n, g.n)), aspect="auto", extent=g.extent,
                       cmap="seismic", vmin=-0.45 * vmax, vmax=0.45 * vmax,
                       interpolation="bilinear")
        ax.set_title(name, color=MCOL[name], fontweight="bold", pad=12)
        ax.set_xlabel("x  [m]")
        ax.set_ylabel("depth offset  [m]")
        ax.set_xlim(-1.0, 1.0)
        ax.set_ylim(1.0, -1.0)
        txt = readout(ax, "", loc="upper left", size=11, color=C_DARK)
        txt.set_family("DejaVu Sans Mono")
        ims[name] = im
        txts[name] = txt

    # Shared banner above all three panels -- pinning this inside any one
    # narrow panel (as "upper right") let a long string span the whole
    # panel width and collide with that panel's own recovered/error
    # readout at "upper left".
    txt_noise = fig.text(0.5, 0.78, "", ha="center", va="bottom",
                         fontsize=15, fontweight="bold", color=C_DARK)

    # -- error-vs-noise panel, one line per method ---------------------------
    ax_err.set_title("(d)  Recovered-displacement error vs. noise level", color=C_DARK, pad=12)
    ax_err.set_xlabel("noise level  [% of signal std]")
    ax_err.set_ylabel("error in recovered  Δx  [%]")
    ax_err.set_facecolor(C_PANEL)
    ax_err.set_xlim(0, args.max_noise * 100)
    ax_err.set_yscale("log")
    ax_err.axvline(THESIS_LEVEL * 100, color=C_GREY, ls=":", lw=1.6)
    ax_err.text(THESIS_LEVEL * 100, 1.0, "  thesis level (10%)", color=C_GREY,
               fontsize=10, rotation=90, va="bottom", ha="left",
               transform=ax_err.get_xaxis_transform())

    lines = {}
    for name in METHODS:
        (ln,) = ax_err.plot([], [], color=MCOL[name], lw=2.6, label=name)
        lines[name] = ln
    ax_err.legend(loc="upper left", fontsize=11)

    # Indexed by frame position, not appended -- update(i) must be
    # idempotent, since save_poster() re-calls update() on the final frame
    # AFTER the animation writer already rendered it once.
    noise_pct = fracs * 100.0
    err_hist = {name: np.full(len(fracs), np.nan) for name in METHODS}

    def update(i):
        frac = fracs[i]
        result = pipe.frame(frac)

        for name in METHODS:
            mig_b, est_dx = result[name]
            ims[name].set_data(mig_b)
            true_l = DX_TRUE / LAMBDA
            est_l = est_dx / LAMBDA if np.isfinite(est_dx) else np.nan
            err = abs(est_l - true_l) / true_l * 100.0 if np.isfinite(est_l) else np.nan
            txts[name].set_text(f"recovered {est_l:6.3f} λ\nerror     {err:6.2f} %")

            err_hist[name][i] = err
            lines[name].set_data(noise_pct[:i + 1], err_hist[name][:i + 1])

        tag = "  <- thesis level" if abs(frac - THESIS_LEVEL) < 0.012 else ""
        txt_noise.set_text(f"noise = {frac*100:.0f}% of signal std{tag}")

        # autoscale the error axis so early near-zero errors don't collapse it
        all_err = [err_hist[name][j] for name in METHODS for j in range(i + 1)
                  if np.isfinite(err_hist[name][j]) and err_hist[name][j] > 0]
        if all_err:
            ax_err.set_ylim(max(1e-3, min(all_err) * 0.5), max(all_err) * 2.0)

        return (list(ims.values()) + list(txts.values()) + list(lines.values())
               + [txt_noise])

    anim = FuncAnimation(fig, update, frames=len(fracs),
                         interval=1000 // args.fps, blit=False, repeat=False)

    stem = os.path.join(args.outdir, "anim7b_noise_by_method")
    print("Rendering anim7b  (slide 14b)  -- noise across three migration methods")
    save_animation(anim, stem, fps=args.fps, dpi=args.dpi, prefer=args.format)
    if not args.no_poster:
        update(len(fracs) - 1)
        save_poster(fig, stem, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
