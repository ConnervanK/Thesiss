"""
render_all.py
=============
Render every defense animation in one go.

    python render_all.py                    # 18 s per animation  (~12 min)
    python render_all.py --duration 25      # slower still
    python render_all.py --preview          # 30 frames, low dpi  (~30 s)
    python render_all.py --format gif       # force GIF instead of MP4

Pacing
------
These are narrated animations: you talk over them for 50-75 s per slide, so a
sweep that finishes in 5 s is over before anyone has read the axes. The default
is 18 s of motion followed by a 3 s freeze on the final frame, giving you a
stable image to keep talking against. Raise --duration if you want more.

Outputs land in  ../assets/animations/  which is exactly where build_deck.py
tells PowerPoint to look for them.

What gets produced
------------------
    anim1_imaging_floor          slides 3 & 4   two scatterers merging (illustrative PSF)
    anim1_extended_resolution    backup/detail  real Kirchhoff migration + metric
    anim2_monitoring_floor       slide 6        amplitude differencing fails
    anim3_envelope_vs_phase      slide 7        why phase is the observable
    anim4_kdomain_lateral        slide 8        shift -> phase ramp
    anim4_kdomain_vertical       (backup)       vertical equivalent
    anim4_kdomain_both           slide 10       lateral vs vertical
    anim5_act1..act4             slide 9        the derivation, one clip per
                                                formula reveal
    anim5_full                   (rehearsal)    all four acts, continuous
    anim6_dumbbell_clean.png     slide 12       clean MAE, thesis numbers
    anim6_dumbbell_noisy         slide 15       the same plot degrading under
                                                noise, animated
    anim7_noise_propagation      slide 14       noise through the whole chain
    anim7b_noise_by_method       backup/detail  same noise ramp, Kirchhoff vs
                                                Gazdag vs back-prop side by side

Each also writes a 300-dpi *_final.png so you have a static fallback if the
projector refuses to play video. Put those on a USB stick too.
"""

import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUT = os.path.normpath(os.path.join(HERE, "..", "assets", "animations"))

JOBS = [
    ("anim1_imaging_floor.py", [], "slides 3 & 4"),
    ("anim1_extended_resolution_floor.py", [], "backup/detail — real migration"),
    ("anim2_monitoring_floor.py", [], "slide 6"),
    ("anim3_envelope_vs_phase.py", [], "slide 7"),
    ("anim4_kdomain_phase_ramp.py", ["--mode", "all"], "slides 8 & 10"),
    ("anim5_crossspectrum_derivation.py", [], "slide 9"),
    ("anim6_mae_dumbbell.py", [], "slides 12 & 15"),
    ("anim7_noise_propagation.py", [], "slide 14"),
    ("anim7b_noise_by_method.py", [], "backup/detail — 3-method comparison"),
]


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=DEFAULT_OUT)
    p.add_argument("--preview", action="store_true",
                   help="fast, low-quality pass for checking layout")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--fps", type=int, default=25)
    p.add_argument("--duration", type=float, default=18.0,
                   help="seconds of motion per animation (pacing knob)")
    p.add_argument("--hold", type=float, default=3.0,
                   help="seconds frozen on the final frame")
    p.add_argument("--format", choices=["auto", "mp4", "gif"], default="auto")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    common = ["--outdir", args.outdir, "--format", args.format]
    if args.preview:
        common += ["--frames", "30", "--dpi", "90", "--fps", "12"]
    else:
        common += ["--dpi", str(args.dpi), "--fps", str(args.fps),
                   "--hold", str(args.hold)]

    print(f"Output directory: {args.outdir}")
    if args.preview:
        print("Mode: PREVIEW")
    else:
        print(f"Mode: FULL QUALITY  —  {args.duration:.0f}s motion "
              f"+ {args.hold:.0f}s freeze, {args.fps} fps")
    print()

    t0 = time.time()
    failed = []
    for script, extra, where in JOBS:
        print(f"[{where}]  {script}")
        # anim5 derives its length from its own four-act structure, so a
        # blanket --duration would squash the acts out of proportion.
        per = list(common)
        if not args.preview and not any(k in script for k in ("anim5", "anim6", "anim7")):
            per += ["--duration", str(args.duration)]
        r = subprocess.run([sys.executable, os.path.join(HERE, script)]
                           + extra + per, cwd=HERE)
        if r.returncode != 0:
            failed.append(script)
        print()

    dt = time.time() - t0
    if failed:
        print(f"FAILED: {', '.join(failed)}")
        sys.exit(1)
    print(f"All animations rendered in {dt:.0f} s -> {args.outdir}")


if __name__ == "__main__":
    main()
