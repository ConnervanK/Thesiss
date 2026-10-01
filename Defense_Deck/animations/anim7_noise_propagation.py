"""
anim7_noise_propagation.py   --   SLIDE 14
==========================================
Laplace noise travelling the whole pipeline, from B-scan to phase-plane fit.

Five panels, all updating together as the noise level ramps:

  (a) B-scan          the clean hyperbola, then noise growing on top of it
  (b) migrated        what that noise becomes after migration
  (c) cross-spectrum  the phase ramp acquiring speckle
      phase
  (d) the fit         phase vs k_x with the weighted least-squares line
  (e) residual        the same points with the fitted plane subtracted, on a
                      ZOOMED axis -- this is where the point cloud is actually
                      visible

READ THIS BEFORE USING THE ANIMATION
------------------------------------
The obvious version of this slide -- "add the thesis's 10% Laplace noise and
watch the straight line become a point cloud" -- does not survive contact with
the numbers. Measured on this pipeline:

    noise (x signal std)   image SNR   fit error   residual rms [rad]
             0.00               inf      0.17 %          0.013
             0.10             1019       0.24 %          0.014
             0.50              204       0.48 %          0.019
             1.00              102       0.79 %          0.028
             4.00               26       2.43 %          0.094

At the thesis's own 10% level the fit barely notices. The reason is that
migration is itself a powerful noise suppressor: roughly 150 traces stack
coherently for the signal and incoherently for the noise, so a 10% B-scan
perturbation arrives at the image as about 0.8% -- and then the |XS| amplitude
mask keeps only the highest-SNR bins on top of that.

So this animation does NOT pretend the fit collapses at 10%. It ramps further
(default 100% of signal std, with the thesis's 10% marked as a waypoint) and
puts the residual on a zoomed axis in panel (e), where the cloud opening up is
genuinely visible rather than imagined.

That is better defense material than the dramatic version would have been. It
lets you say: the phase fit does not fall apart under additive noise -- the
degradation the thesis reports comes from specific migration failures (Gazdag's
cross-axis leakage, wrapping near lambda/2, ROI mis-localisation), not from the
cross-spectrum turning to mush. If someone asks "surely noise destroys the
phase?", this is the answer, with a number attached.

Honest construction
-------------------
  * the B-scan is a Ricker wavelet on the true travel-time hyperbola
  * Laplace noise is added to the B-SCAN, before migration
  * both surveys are migrated by an actual diffraction stack
  * the cross-spectrum and the masked, amplitude-weighted fit are the same
    functions every other animation here uses

One noise realisation is generated up front and scaled as the level ramps, so
growth is smooth rather than flickering. Migration is linear, so the migrated
noise is computed once and scaled too -- which is what makes hundreds of frames
of real migration cheap.

Usage
-----
    python anim7_noise_propagation.py
    python anim7_noise_propagation.py --max-noise 0.10   # strict thesis level
    python anim7_noise_propagation.py --max-noise 4.0    # push until it breaks
    python anim7_noise_propagation.py --frames 30        # layout preview
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from gpr_common import (LAMBDA, V, FC, ImageGrid, apply_house_style, ricker,
                        cross_spectrum, fit_phase_plane, readout, slide_title,
                        footnote, save_animation, save_poster, common_cli,
                        frame_counts, C_DARK, C_ACCENT, C_WARM, C_GREY,
                        C_PANEL, C_BAD, CMAP_IMG, FIGSIZE)

# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------
DEPTH = 1.5 * LAMBDA           # scatterer depth [m]
DX_TRUE = 0.25 * LAMBDA        # lateral displacement between surveys
NTR = 201                      # traces
X_HALF = 4.0                   # antenna sweep half-width [m]
T_MAX = 110.0                  # [ns]
NT = 801
APERTURE = 3.0                 # migration aperture half-width [m]

# Wavenumber resolution is set by the image's spatial EXTENT, not by how many
# samples it has. A 2.4 m migrated window gives dk = 2*pi/2.4 = 2.6 rad/m,
# against k_c = 6.28 -- barely two samples per k_c, which turns the phase panel
# into blocks. So the tapered image is zero-padded out to NPAD before the FFT:
# same physics, same sample spacing, just a finely interpolated spectrum.
NIMG = 128                     # migrated image samples per side
NPAD = 640                     # padded size for the FFT  (~12 samples / k_c)
K_VIEW = 2.0
SUBSAMPLE = 23
MAX_NOISE_DEFAULT = 1.00       # ramp to 1x signal std; 10% is a waypoint

SEED = 20260811
THESIS_LEVEL = 0.10            # the level used throughout the thesis


def _tukey(n, alpha):
    """Tukey window: cosine tapers at the edges, flat in the middle."""
    w = np.ones(n)
    e = max(1, int(alpha * n / 2))
    t = 0.5 * (1.0 - np.cos(np.pi * np.arange(e) / e))
    w[:e] = t
    w[-e:] = t[::-1]
    return w


# ==========================================================================
#  Real forward model + real diffraction-stack migration
# ==========================================================================

class Pipeline:
    def __init__(self, max_noise):
        self.max_noise = max_noise

        # image grid: NIMG samples spanning ~2.4 m, centred on the scatterer
        self.grid = ImageGrid(n=NIMG, lam=LAMBDA, samples_per_lambda=53)
        # padded grid, same sample spacing -- this is the one the FFT and the
        # plane fit actually live on
        self.gk = ImageGrid(n=NPAD, lam=LAMBDA, samples_per_lambda=53)
        g = self.grid
        self.kc = 2.0 * np.pi / LAMBDA
        # Tukey, not Hanning. A Hanning taper is still falling steeply at the
        # monitor lobe's off-centre position, so it weights the two surveys
        # differently and breaks the exact shift relation -- measured as a
        # 10% bias in the recovered displacement at ZERO noise. A Tukey window
        # with a flat centre leaves both lobes untouched and drops that bias
        # to 0.2%.
        wx = _tukey(NIMG, 0.40)
        self._taper = np.outer(wx, wx)

        self.x_tr = np.linspace(-X_HALF, X_HALF, NTR)
        self.t = np.linspace(0.0, T_MAX, NT)
        self.dt = self.t[1] - self.t[0]

        # ---- clean B-scans ------------------------------------------------
        self.bs_b = self._bscan(0.0)
        self.bs_m = self._bscan(DX_TRUE)

        # ---- one fixed noise realisation per survey, scaled later ---------
        rng = np.random.default_rng(SEED)
        sig = self.bs_b.std()
        nb = rng.laplace(0.0, 1.0, self.bs_b.shape)
        nm = rng.laplace(0.0, 1.0, self.bs_m.shape)
        self.nz_b = nb / nb.std() * sig          # std == signal std
        self.nz_m = nm / nm.std() * sig

        # ---- migration operator, precomputed --------------------------------
        self._build_migration_operator()

        # Migration is LINEAR, so migrate the clean data and the noise
        # separately, once. Every frame is then a weighted sum of two
        # precomputed images instead of a fresh migration.
        self.mig_b_clean = self._migrate(self.bs_b)
        self.mig_m_clean = self._migrate(self.bs_m)
        self.mig_b_noise = self._migrate(self.nz_b)
        self.mig_m_noise = self._migrate(self.nz_m)

        self.bs_vmax = np.abs(self.bs_b).max()
        self.mig_vmax = np.abs(self.mig_b_clean).max()

    # -- forward -----------------------------------------------------------

    def _bscan(self, xs):
        r = np.hypot(self.x_tr - xs, DEPTH)
        t0 = 2.0 * r / V
        gain = 1.0 / np.sqrt(r)
        return gain[None, :] * ricker(self.t[:, None] - t0[None, :], FC)

    # -- migration ---------------------------------------------------------

    def _build_migration_operator(self):
        g = self.grid
        X0 = g.x[None, :, None]                     # image x
        Z0 = (DEPTH + g.z)[:, None, None]           # image depth
        XT = self.x_tr[None, None, :]               # trace x

        r = np.sqrt((XT - X0) ** 2 + Z0 ** 2)
        tt = 2.0 * r / V
        idx = np.clip(np.round(tt / self.dt).astype(np.int32), 0, NT - 1)

        # aperture taper + obliquity, exactly the weights a Kirchhoff-style
        # delay-and-sum applies
        off = np.abs(XT - X0)
        w = np.cos(np.arctan2(off, Z0))
        w = np.where(off <= APERTURE, w, 0.0)
        w *= 0.5 * (1.0 + np.cos(np.pi * np.clip(off / APERTURE, 0, 1)))

        self._flat = (idx * NTR + np.arange(NTR)[None, None, :]).astype(np.int64)
        self._w = (w / np.sqrt(r)).astype(np.float32)

    def _migrate(self, bscan):
        return (bscan.ravel()[self._flat] * self._w).sum(axis=-1)

    # -- per-frame ---------------------------------------------------------

    def _pad(self, a):
        out = np.zeros((NPAD, NPAD))
        o = (NPAD - NIMG) // 2
        out[o:o + NIMG, o:o + NIMG] = a * self._taper
        return out

    def frame(self, frac):
        bs_b = self.bs_b + frac * self.nz_b
        mig_b = self.mig_b_clean + frac * self.mig_b_noise
        mig_m = self.mig_m_clean + frac * self.mig_m_noise
        B = self.gk.spectrum_of(self._pad(mig_b))
        M = self.gk.spectrum_of(self._pad(mig_m))
        XS = B * np.conj(M)
        est_dz, est_dx, _, mask = fit_phase_plane(XS, self.gk)
        return bs_b, mig_b, XS, mask, est_dx, est_dz


def phase_rgba(XS, mask):
    cmap = plt.get_cmap("RdBu_r")
    rgba = cmap((np.angle(XS) + np.pi) / (2 * np.pi))
    rgba[..., 3] = np.where(mask, 1.0, 0.05)
    return rgba


# ==========================================================================

def main():
    p = common_cli(__doc__.splitlines()[1])
    p.add_argument("--max-noise", type=float, default=MAX_NOISE_DEFAULT,
                   help="peak noise level, as a fraction of signal std")
    args = p.parse_args()

    apply_house_style(base=14)
    os.makedirs(args.outdir, exist_ok=True)

    print("Building forward model and migration operator ...")
    pipe = Pipeline(args.max_noise)
    g = pipe.grid

    n_frames, n_hold = frame_counts(args)
    n_lead = max(1, n_frames // 8)                 # dwell on the clean state
    fracs = np.concatenate([
        np.zeros(n_lead),
        np.linspace(0.0, args.max_noise, n_frames),
        np.full(n_hold, args.max_noise)])

    # ---------------- figure ----------------
    fig = plt.figure(figsize=FIGSIZE)
    slide_title(fig,
                "Noise propagates all the way to the fit",
                f"Laplace noise added to the B-scan, ramped 0 → "
                f"{args.max_noise*100:.0f}% of signal std  ·  "
                f"true Δx = {DX_TRUE/LAMBDA:.2f} λ")
    footnote(fig, "one noise realisation, scaled — so the growth is smooth "
                  "rather than flickering")

    gs = fig.add_gridspec(2, 3, left=0.062, right=0.905, top=0.775,
                          bottom=0.105, wspace=0.38, hspace=0.60,
                          width_ratios=[1.0, 1.0, 1.05])
    ax_bs = fig.add_subplot(gs[0, 0])
    ax_mig = fig.add_subplot(gs[0, 1])
    ax_k = fig.add_subplot(gs[0, 2])
    ax_fit = fig.add_subplot(gs[1, 0:2])
    ax_res = fig.add_subplot(gs[1, 2])

    # (a) B-scan
    im_bs = ax_bs.imshow(pipe.bs_b, aspect="auto", cmap=CMAP_IMG,
                         extent=(pipe.x_tr[0], pipe.x_tr[-1],
                                 pipe.t[-1], pipe.t[0]),
                         vmin=-0.28 * pipe.bs_vmax, vmax=0.28 * pipe.bs_vmax,
                         interpolation="bilinear")
    ax_bs.set_title("(a)  B-scan  (clipped)", color=C_DARK, pad=16)
    ax_bs.set_xlabel("antenna position  x  [m]")
    ax_bs.set_ylabel("two-way time  [ns]")
    ax_bs.set_ylim(95, 20)

    # (b) migrated image
    im_mig = ax_mig.imshow(pipe.mig_b_clean, aspect="auto", cmap=CMAP_IMG,
                           extent=g.extent, vmin=-0.45 * pipe.mig_vmax,
                           vmax=0.45 * pipe.mig_vmax, interpolation="bilinear")
    ax_mig.set_title("(b)  Migrated  (clipped)", color=C_DARK, pad=16)
    ax_mig.set_xlabel("x  [m]")
    ax_mig.set_ylabel("depth offset  [m]")
    ax_mig.set_xlim(-1.0, 1.0)
    ax_mig.set_ylim(1.0, -1.0)

    # (c) cross-spectrum phase
    gk = pipe.gk
    ext_k = (gk.kx[0] / pipe.kc, gk.kx[-1] / pipe.kc,
             gk.kz[0] / pipe.kc, gk.kz[-1] / pipe.kc)
    im_k = ax_k.imshow(np.zeros((gk.n, gk.n, 4)), aspect="auto", extent=ext_k,
                       origin="lower", interpolation="nearest", zorder=3)
    ax_k.set_facecolor(C_PANEL)
    ax_k.set_title("(c)  Cross-spectrum phase", color=C_DARK, pad=16)
    ax_k.set_xlabel("$k_x$   [$k_c$]")
    ax_k.set_ylabel("$k_z$   [$k_c$]")
    ax_k.set_xlim(-K_VIEW, K_VIEW)
    ax_k.set_ylim(-K_VIEW, K_VIEW)

    # (d) the fit
    ax_fit.set_title("(d)  Phase-plane fit along $k_x$", color=C_DARK, pad=16)
    ax_fit.set_xlabel("$k_x$   [$k_c$]")
    ax_fit.set_ylabel("cross-spectrum phase  [rad]")
    ax_fit.set_xlim(-K_VIEW, K_VIEW)
    ax_fit.set_ylim(-np.pi * 1.08, np.pi * 1.08)
    ax_fit.set_facecolor(C_PANEL)
    ax_fit.set_yticks([-np.pi, 0, np.pi])
    ax_fit.set_yticklabels([r"$-\pi$", "0", r"$\pi$"])
    ax_fit.axhline(0, color=C_GREY, lw=0.8)
    (ln_fit,) = ax_fit.plot([], [], color=C_WARM, lw=3.2, zorder=3,
                            label="weighted least-squares fit")
    sc = ax_fit.scatter([], [], s=22, c=[], cmap="viridis", vmin=0, vmax=1,
                        alpha=0.9, zorder=5, edgecolors="none")
    ax_fit.legend(loc="lower right", fontsize=11)

    # (e) residual, zoomed -- the panel where the cloud is actually visible
    ax_res.set_title("(e)  Residual  (zoomed)", color=C_DARK, pad=16)
    ax_res.set_xlabel("$k_x$   [$k_c$]")
    ax_res.set_ylabel("phase $-$ fitted plane  [rad]")
    ax_res.set_xlim(-K_VIEW, K_VIEW)
    ax_res.set_ylim(-0.16, 0.16)
    ax_res.set_facecolor(C_PANEL)
    ax_res.axhline(0, color=C_WARM, lw=2.2, zorder=3)
    sc_r = ax_res.scatter([], [], s=20, c=[], cmap="viridis", vmin=0, vmax=1,
                          alpha=0.9, zorder=5, edgecolors="none")
    txt_rms = readout(ax_res, "", loc="upper left", size=12, color=C_DARK)
    txt_rms.set_family("DejaVu Sans Mono")

    # ONE shared colourbar for the weight colouring, which both (d) and (e)
    # use. Attached to the right of (e) so it cannot reach into a neighbour.
    rb = ax_res.get_position()
    cax = fig.add_axes([rb.x1 + 0.012, rb.y0, 0.011, rb.height])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label("fit weight  |XS|", fontsize=10.5)
    cb.ax.tick_params(labelsize=10)

    txt_noise = readout(ax_bs, "", loc="upper left", size=14, color=C_BAD)
    txt_fit = readout(ax_fit, "", loc="upper left", size=13, color=C_ACCENT)
    txt_fit.set_family("DejaVu Sans Mono")

    def update(i):
        frac = fracs[i]
        bs_b, mig_b, XS, mask, est_dx, est_dz = pipe.frame(frac)

        im_bs.set_data(bs_b)
        im_mig.set_data(mig_b)
        im_k.set_data(phase_rgba(XS, mask))

        phi = np.angle(XS)
        w = np.abs(XS)
        resid = phi[mask] - gk.KZ[mask] * est_dz
        kax = gk.KX[mask] / pipe.kc
        g_kx_masked = gk.KX[mask]
        ww = w[mask]
        ww = ww / (ww.max() + 1e-30)

        sel = slice(None, None, SUBSAMPLE)
        sc.set_offsets(np.column_stack([kax[sel], resid[sel]]))
        sc.set_array(ww[sel])

        if kax.size:
            kk = np.linspace(kax.min(), kax.max(), 60)
            ln_fit.set_data(kk, est_dx * kk * pipe.kc)

        true_l = DX_TRUE / LAMBDA
        est_l = est_dx / LAMBDA
        err = abs(est_l - true_l) / true_l * 100.0

        # residual panel
        resid_c = resid - est_dx * g_kx_masked
        sc_r.set_offsets(np.column_stack([kax[sel], resid_c[sel]]))
        sc_r.set_array(ww[sel])
        rms = float(np.sqrt(np.mean(resid_c ** 2)))
        txt_rms.set_text(f"residual rms  {rms:6.3f} rad")

        tag = "  ← thesis level" if abs(frac - THESIS_LEVEL) < 0.012 else ""
        txt_noise.set_text(f"noise = {frac*100:.0f}% of signal std{tag}")
        txt_fit.set_text(f"true       {true_l:7.4f} λ\n"
                         f"recovered  {est_l:7.4f} λ\n"
                         f"error      {err:7.2f} %")
        return [im_bs, im_mig, im_k, sc, sc_r, ln_fit, txt_noise, txt_fit,
                txt_rms]

    anim = FuncAnimation(fig, update, frames=len(fracs),
                         interval=1000 // args.fps, blit=False, repeat=False)

    stem = os.path.join(args.outdir, "anim7_noise_propagation")
    print("Rendering anim7  (slide 14)  — noise through the pipeline")
    save_animation(anim, stem, fps=args.fps, dpi=args.dpi, prefer=args.format)
    if not args.no_poster:
        update(len(fracs) - 1)
        save_poster(fig, stem, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
