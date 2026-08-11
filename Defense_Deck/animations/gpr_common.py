"""
gpr_common.py
=============
Shared physics, styling and export helpers for the MSc defense animations.

Everything here is *illustrative*, not a reproduction of the thesis simulations.
The dominant wavelength is deliberately set to lambda = 1.0 m so that every
number an audience reads off a slide is trivially interpretable:

    lambda      = 1.00 m
    lambda/2    = 0.50 m   <-- classical imaging resolution floor
    lambda/4    = 0.25 m
    lambda/32   = 0.031 m  = 31 mm  <-- smallest scale tested in the thesis

The thesis itself used f_c = 1.5 GHz in ice (lambda = 112.6 mm). The physics is
scale-invariant, so the demonstration is faithful; only the units are friendlier.
Say this out loud on the slide if anyone asks.

Physics model
-------------
A migrated zero-offset image of a point scatterer is modelled directly in the
wavenumber domain, which is both physically correct and exactly what makes the
phase argument visible:

  * Radial support: a Gaussian band centred on k_c = 2*pi/lambda. This is the
    pulse bandwidth -> it sets vertical resolution.
  * Angular support: a smooth cutoff at the migration aperture angle. This is
    the finite aperture -> it sets lateral resolution.

The product is the migrated point-spread function's spectrum. Because a spatial
translation is an exact phase rotation in this domain (the Fourier shift
theorem), sub-wavelength displacement is applied *exactly*, with no
interpolation error to muddy the demonstration.

Author: generated for Conner's MSc defense deck.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass

import numpy as np
import matplotlib

matplotlib.use("Agg")  # headless-safe; no window is ever opened
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter

# --------------------------------------------------------------------------
# 1. Illustrative survey parameters
# --------------------------------------------------------------------------

LAMBDA = 1.0        # dominant wavelength of the migrated image [m]
V = 0.10            # medium velocity [m/ns]  (-> f_c = 0.1 GHz for lambda = 1 m)
FC = V / LAMBDA     # centre frequency [GHz]
DEPTH = 5.0         # scatterer burial depth for the raw B-scan panels [m]

# --------------------------------------------------------------------------
# 2. House style  (clean minimal / white, per the deck design)
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

# 16:9 at PowerPoint's native widescreen size
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
# 3. Wavelets and raw B-scan synthesis
# --------------------------------------------------------------------------

def ricker(t: np.ndarray, fc: float) -> np.ndarray:
    """Ricker (Mexican-hat) wavelet, peak-normalised. t in ns, fc in GHz."""
    a = (np.pi * fc * t) ** 2
    return (1.0 - 2.0 * a) * np.exp(-a)


def synth_bscan(x_traces: np.ndarray,
                t: np.ndarray,
                scatterers,
                v: float = V,
                fc: float = FC,
                spreading: bool = True) -> np.ndarray:
    """
    Zero-offset B-scan: a Ricker wavelet placed along each scatterer's
    two-way travel-time hyperbola. This is the textbook 'convolve a Ricker
    with a point scatterer' construction.

    Parameters
    ----------
    x_traces : (nx,) antenna positions [m]
    t        : (nt,) two-way time axis [ns]
    scatterers : iterable of (x_s, z_s, amplitude)

    Returns
    -------
    (nt, nx) B-scan.
    """
    B = np.zeros((t.size, x_traces.size), dtype=float)
    for xs, zs, amp in scatterers:
        r = np.hypot(x_traces - xs, zs)              # (nx,)
        t0 = 2.0 * r / v                             # two-way time [ns]
        gain = amp / np.sqrt(r) if spreading else np.full_like(r, amp)
        B += gain[None, :] * ricker(t[:, None] - t0[None, :], fc)
    return B


# --------------------------------------------------------------------------
# 4. Migrated image: the wavenumber-domain point-spread function
# --------------------------------------------------------------------------

@dataclass
class ImageGrid:
    """A centred (x, z) image grid plus its matching centred wavenumber grid."""
    n: int = 256
    lam: float = LAMBDA
    samples_per_lambda: int = 16

    def __post_init__(self):
        self.d = self.lam / self.samples_per_lambda           # sample spacing [m]
        idx = np.arange(self.n) - self.n // 2
        self.x = idx * self.d
        self.z = idx * self.d
        k = 2.0 * np.pi * np.fft.fftshift(np.fft.fftfreq(self.n, self.d))
        # array is indexed [iz, ix] -> KX varies along axis 1, KZ along axis 0
        self.KX, self.KZ = np.meshgrid(k, k)
        self.kx = k
        self.kz = k
        self.extent = (self.x[0], self.x[-1], self.z[-1], self.z[0])  # origin upper

    # -- spectrum -----------------------------------------------------------

    def psf_spectrum(self, aperture_deg: float = 50.0,
                     bw_frac: float = 0.85) -> np.ndarray:
        """
        Migrated point-scatterer spectrum: a Gaussian band around k_c
        (pulse bandwidth) truncated by a smooth aperture cutoff.

        aperture_deg : half-aperture of the migration operator
        bw_frac      : fractional bandwidth (FWHM / k_c)
        """
        kc = 2.0 * np.pi / self.lam
        K = np.hypot(self.KX, self.KZ)
        theta = np.arctan2(np.abs(self.KX), np.abs(self.KZ))

        sigma = bw_frac * kc / 2.3548                     # FWHM -> sigma
        radial = np.exp(-0.5 * ((K - kc) / sigma) ** 2)

        th0 = np.radians(aperture_deg)
        soft = np.radians(13.0)   # generous taper -> low PSF sidelobes
        aperture = 0.5 * (1.0 - np.tanh((theta - th0) / soft))

        S = radial * aperture
        S[K < 1e-9] = 0.0
        return S

    # -- image --------------------------------------------------------------

    def psf_image(self, S: np.ndarray, dx: float = 0.0, dz: float = 0.0) -> np.ndarray:
        """
        Migrated image of one point scatterer displaced by (dx, dz).

        The displacement is applied as an exact phase rotation
        exp(-i(kx*dx + kz*dz)) -- the Fourier shift theorem -- so there is no
        interpolation error even at lambda/1000.
        """
        F = S * np.exp(-1j * (self.KX * dx + self.KZ * dz))
        img = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(F)))
        return img.real

    def spectrum_of(self, img: np.ndarray) -> np.ndarray:
        """Centred 2D spectrum of a real image."""
        return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(img)))


# --------------------------------------------------------------------------
# 5. The method itself: cross-spectrum -> masked WLS plane fit
# --------------------------------------------------------------------------

def cross_spectrum(img_base: np.ndarray, img_mon: np.ndarray,
                   grid: ImageGrid, taper: bool = True):
    """
    Cross-spectrum XS = B * conj(M).

    With M(k) = B(k)*exp(-i k.d) this gives XS = |B|^2 * exp(+i k.d), so
    angle(XS) = kz*dz + kx*dx  --  a flat plane through the origin.
    """
    if taper:
        wx = np.hanning(grid.n)
        w = np.outer(wx, wx)
        img_base = img_base * w
        img_mon = img_mon * w
    B = grid.spectrum_of(img_base)
    M = grid.spectrum_of(img_mon)
    return B * np.conj(M)


def fit_phase_plane(XS: np.ndarray, grid: ImageGrid,
                    amp_frac: float = 0.10,
                    band_factor: float = 1.4):
    """
    Weighted least-squares fit of Phi = kz*dz + kx*dx + c.

    Mirrors the thesis estimator:
      * amplitude mask  |XS| > amp_frac * max|XS|
      * band-pass mask  |kz|,|kx| < band_factor * k_c   (prevents wrapping)
      * weights W_ii = |XS_i|, solved via the pre-weighted system for stability

    Returns (dz, dx, c, mask).
    """
    kc = 2.0 * np.pi / grid.lam
    mag = np.abs(XS)
    phi = np.angle(XS)

    mask = (mag > amp_frac * mag.max())
    mask &= (np.abs(grid.KZ) < band_factor * kc)
    mask &= (np.abs(grid.KX) < band_factor * kc)

    if mask.sum() < 10:
        return np.nan, np.nan, np.nan, mask

    kz = grid.KZ[mask]
    kx = grid.KX[mask]
    w = mag[mask]
    A = np.column_stack([kz, kx, np.ones_like(kz)])
    # pre-weight both sides, then solve by SVD (never form A^T W^2 A directly)
    sol, *_ = np.linalg.lstsq(A * w[:, None], phi[mask] * w, rcond=None)
    return sol[0], sol[1], sol[2], mask


# --------------------------------------------------------------------------
# 6. Resolution diagnostics (Rayleigh criterion)
# --------------------------------------------------------------------------

def fwhm(x: np.ndarray, y: np.ndarray) -> float:
    """Full width at half maximum of the main lobe of |y|, by linear interp."""
    a = np.abs(y)
    i = int(np.argmax(a))
    half = 0.5 * a[i]
    if a[i] <= 0:
        return np.nan

    def cross(idx_range):
        prev = i
        for j in idx_range:
            if a[j] < half:
                f = (half - a[j]) / (a[prev] - a[j] + 1e-30)
                return x[j] + f * (x[prev] - x[j])
            prev = j
        return np.nan

    left = cross(range(i - 1, -1, -1))
    right = cross(range(i + 1, a.size))
    if np.isnan(left) or np.isnan(right):
        return np.nan
    return right - left


def two_peaks_resolved(y: np.ndarray, dip_frac: float = 0.26) -> bool:
    """
    Classical Rayleigh test on a profile: are there two local maxima with a
    dip between them deeper than dip_frac of the lower peak?
    """
    a = np.abs(y)
    interior = np.arange(1, a.size - 1)
    peaks = interior[(a[1:-1] > a[:-2]) & (a[1:-1] > a[2:])]
    peaks = peaks[a[peaks] > 0.4 * a.max()]
    if peaks.size < 2:
        return False
    p1, p2 = peaks[0], peaks[-1]
    valley = a[p1:p2 + 1].min()
    return valley < (1.0 - dip_frac) * min(a[p1], a[p2])


def lateral_profile(img: np.ndarray, grid: ImageGrid):
    """Lateral cut through the row of maximum absolute amplitude."""
    row = int(np.unravel_index(np.argmax(np.abs(img)), img.shape)[0])
    return grid.x, img[row, :], row


# --------------------------------------------------------------------------
# 7. Slide furniture
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
# 8. Export  --  non-looping by construction
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

    A note on DPI: PowerPoint plays video at the projector's pixel resolution,
    so 200 dpi (2667x1500, above 1080p) is already oversampled for playback.
    The 300-dpi requirement matters for the *still* poster frame, which
    save_poster() always writes at 300.
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


def hold(frames: np.ndarray, n_hold: int) -> np.ndarray:
    """Append n_hold repeats of the last value so the animation freezes on it."""
    if n_hold <= 0:
        return frames
    return np.concatenate([frames, np.repeat(frames[-1], n_hold)])


def common_cli(description: str):
    """Shared argument parser so every script behaves identically."""
    import argparse
    p = argparse.ArgumentParser(description=description)
    p.add_argument("--outdir", default="output", help="output directory")
    p.add_argument("--dpi", type=int, default=200,
                   help="animation dpi (200 = 2667x1500; use 300 for 4000x2250)")
    p.add_argument("--fps", type=int, default=25)
    p.add_argument("--duration", type=float, default=DEFAULT_DURATION,
                   help=f"seconds of motion (default {DEFAULT_DURATION:.0f}). "
                        "This is the knob for pacing — raise it to slow down.")
    p.add_argument("--hold", type=float, default=DEFAULT_HOLD,
                   help=f"seconds frozen on the final frame "
                        f"(default {DEFAULT_HOLD:.0f})")
    p.add_argument("--frames", type=int, default=None,
                   help="override the frame count outright (use ~30 for a "
                        "fast layout preview); ignores --duration")
    p.add_argument("--format", choices=["auto", "mp4", "gif"], default="auto")
    p.add_argument("--no-poster", action="store_true",
                   help="skip the 300-dpi final-frame PNG")
    return p


# Pacing defaults.
#
# These animations are narrated, not watched in silence: you speak over them for
# 50-75 s per slide. A sweep that finishes in 5 s is over before anyone has read
# the axes. 18 s of motion lets the audience actually track what is changing,
# and the 3 s freeze gives you a stable final frame to keep talking against.
DEFAULT_DURATION = 18.0     # seconds of motion
DEFAULT_HOLD = 3.0          # seconds frozen at the end


def frame_counts(args):
    """
    Resolve (n_motion_frames, n_hold_frames) from the CLI.

    --frames overrides everything (preview mode) and scales the hold down with
    it, so a 30-frame preview does not spend a quarter of its length frozen.
    """
    if args.frames:
        n = max(2, args.frames)
        return n, max(2, n // 6)
    n = max(2, int(round(args.duration * args.fps)))
    hold = max(1, int(round(args.hold * args.fps)))
    return n, hold
