"""
anim1_extended_resolution_floor.py   --   EXTENDED IMAGING RESOLUTION
=====================================================================
Five panels tell the complete story: subsurface sources, raw data, migration,
lateral profile, and real-time resolution-metric evolution.

Two point scatterers walk toward each other from 2λ to ⅛λ separation:

  (a) Subsurface layout    the two point sources closing in the ground
  (b) Raw B-scan          hyperbolas converging as sources merge
  (c) Migrated image       real Kirchhoff migration (not illustrative PSF)
  (d) Lateral profile      peaked lobes sharpening & merging with verdict
  (e) Resolution metric    separation / FWHM ratio plotted in real time,
                           showing the crossing through the Rayleigh floor

Usage
-----
    python anim1_extended_resolution_floor.py              # 18 s + 3 s freeze
    python anim1_extended_resolution_floor.py --duration 25  # slower
    python anim1_extended_resolution_floor.py --frames 30    # preview
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gpr_style as gs

# Physical constants -- matched to gpr_common.py so this stays on the same
# scale as the illustrative anim1 and the real-physics gif_maker cells
# (v=0.1 m/ns, lambda=1 m -> fc=0.1 GHz; DEPTH=5 m keeps two-way time inside
# a sane T_MAX without the reflection falling outside the plotted window).
LAMBDA = 1.0                # wavelength [m]
V = 0.10                    # velocity [m/ns]
FC = V / LAMBDA              # centre frequency [GHz]
DEPTH = 5.0                 # depth of targets [m]

# Sweep: separation from 2 lambda down to 1/8 lambda (geometric)
SEP_START = 2.00 * LAMBDA
SEP_END = 0.125 * LAMBDA

# Raw B-scan geometry
NX_TRACES = 260
X_HALF = 6.5               # +/- extent [m]
T_MAX = 160.0              # [ns]
NT = 700

# Migration grid (for real Kirchhoff) -- centred on the TRUE target depth,
# not on zero, since z is an absolute depth axis here.
NX_MIG = 256
NZ_MIG = 256
X_MIG_HALF = 1.8 * LAMBDA
Z_MIG_MIN, Z_MIG_MAX = DEPTH - 1.6 * LAMBDA, DEPTH + 1.6 * LAMBDA


def ricker(t, f):
    """Ricker wavelet: (1 - 2(πft)²) exp(-(πft)²)"""
    return (1 - 2 * (np.pi * f * t)**2) * np.exp(-(np.pi * f * t)**2)


def synth_bscan(x_tr, t, sources):
    """
    Synthesize B-scan from point sources.
    sources: list of (x_pos, depth, amplitude) tuples
    """
    B = np.zeros((len(t), len(x_tr)))
    dt_dt = t[1] - t[0] if len(t) > 1 else 1.0

    for x_src, z_src, amp in sources:
        for i, x_receiver in enumerate(x_tr):
            dist = np.sqrt((x_receiver - x_src)**2 + z_src**2)
            B[:, i] += amp * ricker(t - (2 * dist / V), FC) * (z_src / dist)**2
    return B


def kirchhoff_migrate(bscan, x_tr, t, x_mig, z_mig, velocity):
    """
    Fast Kirchhoff migration: stack along diffraction curves.
    Returns image on (z_mig, x_mig) grid.
    """
    dt = t[1] - t[0] if len(t) > 1 else 1.0
    nt = len(t)

    img = np.zeros((len(z_mig), len(x_mig)))
    X_mig_grid, Z_mig_grid = np.meshgrid(x_mig, z_mig)

    for i, x_receiver in enumerate(x_tr):
        # Distance from each image point to this receiver
        r = np.sqrt((X_mig_grid - x_receiver)**2 + Z_mig_grid**2)
        # Two-way travel time
        t_grid = 2 * r / velocity
        # Interpolate B-scan at this time for all image points
        t_idx = np.clip(t_grid / dt, 0, nt - 1.5).astype(int)
        frac = np.clip(t_grid / dt, 0, nt - 1) - t_idx

        b_contrib = (1 - frac) * bscan[t_idx, i] + frac * bscan[np.minimum(t_idx + 1, nt - 1), i]
        # Obliquity correction: (z / r)
        b_contrib *= Z_mig_grid / (r + 1e-30)
        img += b_contrib

    return img


def lateral_profile(img, z_idx=None):
    """Extract lateral profile at a given depth (or max amplitude depth)."""
    if z_idx is None:
        # Find depth of maximum amplitude
        max_per_z = np.abs(img).max(axis=1)
        z_idx = np.argmax(max_per_z)
    return img[z_idx, :]


def fwhm(profile):
    """Estimate FWHM from a 1D profile (assume it's centered at max)."""
    prof_abs = np.abs(profile)
    max_val = prof_abs.max()
    half_max = max_val / 2.0
    indices_above = np.where(prof_abs >= half_max)[0]
    if len(indices_above) < 2:
        return 0.0
    return indices_above[-1] - indices_above[0]


def two_peaks_resolved(profile, threshold=0.5):
    """Check if two peaks are resolvable in the profile."""
    prof_abs = np.abs(profile)
    # Find local minima; if profile is not well-separated, central dip is shallow
    mid = len(profile) // 2
    max_val = prof_abs.max()
    if prof_abs[mid] < threshold * max_val:
        return True
    return False


def build_sweep(n_frames):
    """Geometric spacing from SEP_START to SEP_END."""
    return np.geomspace(SEP_START, SEP_END, n_frames)


def common_cli(title):
    """Minimal CLI parser."""
    import argparse
    p = argparse.ArgumentParser(description=title, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default="../assets/animations")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--fps", type=int, default=25)
    p.add_argument("--duration", type=float, default=18.0, help="seconds of motion")
    p.add_argument("--hold", type=float, default=3.0, help="seconds frozen on final frame")
    p.add_argument("--format", choices=["auto", "mp4", "gif"], default="auto")
    p.add_argument("--frames", type=int, default=0, help="override frame count for preview")
    p.add_argument("--no-poster", action="store_true", help="skip 300 dpi poster")
    return p


def frame_counts(args):
    """Compute n_frames and n_hold from CLI args."""
    if args.frames:
        return args.frames, max(1, int(args.fps * args.hold / 3.0))
    n_frames = max(10, int(args.fps * args.duration))
    n_hold = max(1, int(args.fps * args.hold))
    return n_frames, n_hold


def main():
    args = common_cli(__doc__.splitlines()[1]).parse_args()
    n_frames, n_hold = frame_counts(args)

    gs.apply_house_style(base=14)
    os.makedirs(args.outdir, exist_ok=True)

    seps = build_sweep(n_frames)
    seps = np.concatenate([seps, np.repeat(seps[-1], n_hold)])

    # Static geometry
    x_tr = np.linspace(-X_HALF, X_HALF, NX_TRACES)
    t = np.linspace(0.0, T_MAX, NT)
    x_mig = np.linspace(-X_MIG_HALF, X_MIG_HALF, NX_MIG)
    z_mig = np.linspace(Z_MIG_MIN, Z_MIG_MAX, NZ_MIG)

    # Reference FWHM: measured from a SINGLE isolated scatterer's migrated
    # response, exactly as the thesis defines the resolution kernel. Using the
    # two-source B0 here instead would pick up the span between the two
    # lobes' half-max crossings once they're resolved, not one lobe's width.
    dx_mig = x_mig[1] - x_mig[0]
    B_ref = synth_bscan(x_tr, t, [(0.0, DEPTH, 1.0)])
    img_ref = kirchhoff_migrate(B_ref, x_tr, t, x_mig, z_mig, V)
    prof_ref = lateral_profile(img_ref)
    FWHM_REF = fwhm(prof_ref) * dx_mig  # physical units [m]

    # Figure with 5 columns
    fig = plt.figure(figsize=(16, 5.2))
    gs_main = fig.add_gridspec(1, 5, left=0.05, right=0.97, top=0.82, bottom=0.12,
                               wspace=0.38, width_ratios=[0.9, 1.0, 1.0, 1.1, 1.2])

    ax_src = fig.add_subplot(gs_main[0, 0])  # Subsurface sources
    ax_b = fig.add_subplot(gs_main[0, 1])    # Raw B-scan
    ax_m = fig.add_subplot(gs_main[0, 2])    # Migrated
    ax_p = fig.add_subplot(gs_main[0, 3])    # Lateral profile
    ax_m_plot = fig.add_subplot(gs_main[0, 4])  # Metric vs separation

    gs.slide_title(fig,
                   "Two reflectors closing in: real migration edition",
                   f"Kirchhoff diffraction-stack migration  ·  "
                   f"λ = {LAMBDA:.0f} m (v = {V:.2f} m/ns, fᴄ = {FC*1000:.1f} MHz)")
    gs.footnote(fig, "Left to right: subsurface geometry | raw B-scan | migrated image | profile | resolution metric")

    # ========== Panel (a): Subsurface sources ==========
    ax_src.set_xlim(-1.2 * LAMBDA, 1.2 * LAMBDA)
    ax_src.set_ylim(DEPTH + 5, -5)
    ax_src.set_title("(a)  Subsurface sources", color=gs.C_DARK, pad=12)
    ax_src.set_xlabel("x  [m]")
    ax_src.set_ylabel("depth  z  [m]")
    ax_src.grid(True, alpha=0.2, color=gs.C_GREY)

    (src_l,) = ax_src.plot([], [], 'o', color=gs.C_ACCENT, ms=10, markeredgecolor="white", markeredgewidth=1.5)
    (src_r,) = ax_src.plot([], [], 'o', color=gs.C_ACCENT, ms=10, markeredgecolor="white", markeredgewidth=1.5)
    (src_sep_line,) = ax_src.plot([], [], '--', color=gs.C_WARM, lw=1.8, alpha=0.6)
    txt_src_sep = gs.readout(ax_src, "", loc="upper left", size=12)

    # ========== Panel (b): Raw B-scan ==========
    B0 = synth_bscan(x_tr, t, [(-SEP_START / 2, DEPTH, 1.0), (+SEP_START / 2, DEPTH, 1.0)])
    vmaxB = np.abs(B0).max()
    im_b = ax_b.imshow(B0, aspect="auto", cmap=gs.CMAP_IMG,
                       extent=(x_tr[0], x_tr[-1], t[-1], t[0]),
                       vmin=-vmaxB, vmax=vmaxB, interpolation="bilinear")
    ax_b.set_title("(b)  Raw B-scan", color=gs.C_DARK, pad=12)
    ax_b.set_xlabel("antenna x  [m]")
    ax_b.set_ylabel("two-way time  [ns]")
    ax_b.set_ylim(150, 60)

    # ========== Panel (c): Migrated image ==========
    zoom = 1.6 * LAMBDA
    img0 = kirchhoff_migrate(B0, x_tr, t, x_mig, z_mig, V)
    vmaxM = np.abs(img0).max()
    im_m = ax_m.imshow(img0, aspect="auto", cmap=gs.CMAP_IMG,
                       extent=(x_mig[0], x_mig[-1], z_mig[-1], z_mig[0]),
                       vmin=-vmaxM, vmax=vmaxM, interpolation="bilinear")
    ax_m.set_xlim(-zoom, zoom)
    ax_m.set_ylim(DEPTH + zoom, DEPTH - zoom)
    ax_m.set_title("(c)  After Kirchhoff migration", color=gs.C_DARK, pad=12)
    ax_m.set_xlabel("x  [m]")
    ax_m.set_ylabel("depth  z  [m]")

    # ========== Panel (d): Lateral profile ==========
    ax_p.set_title("(d)  Amplitude profile", color=gs.C_DARK, pad=12)
    ax_p.set_xlabel("x  [m]")
    ax_p.set_ylabel("normalised amplitude")
    ax_p.set_xlim(-zoom, zoom)
    ax_p.set_ylim(-0.55, 1.68)
    ax_p.axhline(0.0, color=gs.C_GREY, lw=0.8, zorder=0)

    (ln_prof,) = ax_p.plot([], [], color=gs.C_ACCENT, lw=3.0, zorder=5)
    fill_l = ax_p.axvspan(0, 0, color=gs.C_GREY, alpha=0.16, lw=0)
    fill_r = ax_p.axvspan(0, 0, color=gs.C_GREY, alpha=0.16, lw=0)
    vl_l = ax_p.axvline(0, color=gs.C_WARM, ls="--", lw=1.6)
    vl_r = ax_p.axvline(0, color=gs.C_WARM, ls="--", lw=1.6)
    ax_p.plot([], [], color=gs.C_WARM, ls="--", lw=1.6, label="true positions")
    ax_p.plot([], [], color=gs.C_GREY, lw=7, alpha=0.35,
              label=f"FWHM ≈ {FWHM_REF/LAMBDA:.2f} λ")
    ax_p.legend(loc="lower left", fontsize=10.5)

    verdict = gs.verdict_box(ax_p, loc="upper right", size=12)

    # ========== Panel (e): Resolution metric plot ==========
    ax_m_plot.set_title("(e)  Resolution metric", color=gs.C_DARK, pad=12)
    ax_m_plot.set_xlabel("true separation  [λ]")
    ax_m_plot.set_ylabel("separation / FWHM")
    ax_m_plot.grid(True, alpha=0.3, color=gs.C_GREY, which="both")

    # Ratio range spans SEP_END/FWHM_REF .. SEP_START/FWHM_REF; give it headroom
    # on both ends rather than a fixed bound (FWHM_REF depends on aperture/
    # bandwidth of the real migration, not a value chosen to fit an axis).
    ratio_lo = 0.6 * (SEP_END / FWHM_REF)
    ratio_hi = 1.3 * (SEP_START / FWHM_REF)
    x_lo, x_hi = 0.08, 2.5

    # Draw Rayleigh criterion line at y=1
    ax_m_plot.axhline(1.0, color=gs.C_BAD, ls="--", lw=2.0, alpha=0.7, label="Rayleigh floor (1.0)")
    ax_m_plot.fill_between([x_lo, x_hi], [ratio_lo, ratio_lo], [1.0, 1.0],
                           color=gs.C_BAD, alpha=0.08, label="merged")
    ax_m_plot.fill_between([x_lo, x_hi], [1.0, 1.0], [ratio_hi, ratio_hi],
                           color=gs.C_GOOD, alpha=0.08, label="resolved")

    (ln_metric,) = ax_m_plot.plot([], [], '-o', color=gs.C_ACCENT, lw=2.5, ms=5, zorder=5, label="swept separation")
    ax_m_plot.set_xlim(x_lo, x_hi)
    ax_m_plot.set_ylim(ratio_lo, ratio_hi)
    ax_m_plot.set_xscale("log")
    ax_m_plot.set_yscale("log")
    ax_m_plot.legend(loc="upper left", fontsize=10)

    # Track all points for the metric plot
    metric_seps = []
    metric_ratios = []

    # ========== Animation update ==========
    def update(i):
        sep = seps[i]
        h = sep / 2.0

        # Panel (a): Subsurface sources
        src_l.set_data([-h], [DEPTH])
        src_r.set_data([+h], [DEPTH])
        src_sep_line.set_data([-h, +h], [DEPTH, DEPTH])
        txt_src_sep.set_text(f"sep = {sep/LAMBDA:.3f} λ")

        # Panel (b): Raw B-scan
        B = synth_bscan(x_tr, t, [(-h, DEPTH, 1.0), (+h, DEPTH, 1.0)])
        im_b.set_data(B)

        # Panel (c): Migrated image (real Kirchhoff)
        img = kirchhoff_migrate(B, x_tr, t, x_mig, z_mig, V)
        im_m.set_data(img)

        # Panel (d): Lateral profile
        prof = lateral_profile(img)
        prof_n = prof / (np.abs(prof).max() + 1e-30)

        # Convert grid indices to physical positions
        x_prof = np.linspace(-X_MIG_HALF, X_MIG_HALF, len(prof))
        ln_prof.set_data(x_prof, prof_n)

        vl_l.set_xdata([-h, -h])
        vl_r.set_xdata([+h, +h])
        for patch, centre in ((fill_l, -h), (fill_r, +h)):
            patch.set_x(centre - FWHM_REF / 2)
            patch.set_width(FWHM_REF)

        # Ratio uses the FIXED single-scatterer FWHM_REF as the resolution
        # kernel (the Rayleigh-criterion definition) -- not a live refit of
        # this frame's (possibly two-lobed) profile, which would measure the
        # lobe-to-lobe span instead of a single lobe's width once resolved.
        ratio = sep / FWHM_REF
        ok = two_peaks_resolved(prof_n) and ratio >= 1.0

        label = "RESOLVED" if ok else "MERGED"
        gs.set_verdict(verdict, f"{label}\nratio = {ratio:.2f}", ok)
        ln_prof.set_color(gs.C_ACCENT if ok else gs.C_BAD)

        # Panel (e): Update metric plot
        metric_seps.append(sep / LAMBDA)
        metric_ratios.append(ratio)

        if len(metric_seps) > 0:
            ln_metric.set_data(metric_seps, metric_ratios)

        return [im_b, im_m, ln_prof, vl_l, vl_r, fill_l, fill_r,
                src_l, src_r, src_sep_line, txt_src_sep, ln_metric, verdict]

    anim = FuncAnimation(fig, update, frames=len(seps),
                         interval=1000 // args.fps, blit=False, repeat=False)

    stem = os.path.join(args.outdir, "anim1_extended_resolution_floor")
    print("Rendering anim1_extended  —  real Kirchhoff migration with metric tracking")
    gs.save_animation(anim, stem, fps=args.fps, dpi=args.dpi, prefer=args.format)
    if not args.no_poster:
        update(len(seps) - 1)
        gs.save_poster(fig, stem, dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    main()
