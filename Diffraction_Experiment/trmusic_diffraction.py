import argparse
import glob
import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


def _rx_sort_key(name: str) -> int:
    match = re.search(r"(\d+)$", name)
    return int(match.group(1)) if match else 0


def _load_out_file(out_path: Path):
    with h5py.File(out_path, "r") as f:
        dt = float(f.attrs["dt"])
        rx_names = sorted(list(f["rxs"].keys()), key=_rx_sort_key)

        traces = []
        rx_pos = []
        for name in rx_names:
            rx_group = f["rxs"][name]
            traces.append(rx_group["Ez"][:])
            rx_pos.append(rx_group.attrs["Position"])

        tx_pos = f["srcs"]["src1"].attrs["Position"]

    traces = np.asarray(traces, dtype=float)
    rx_pos = np.asarray(rx_pos, dtype=float)
    tx_pos = np.asarray(tx_pos, dtype=float)
    return dt, traces, rx_pos, tx_pos


def _discover_files(glob_pattern: str, must_not_contain: str = ""):
    paths = sorted(Path(p) for p in glob.glob(glob_pattern))
    if must_not_contain:
        paths = [p for p in paths if must_not_contain not in p.name]
    return paths


def _pair_baselines(tx_paths, baseline_paths):
    if len(baseline_paths) == 0:
        return [None] * len(tx_paths)
    if len(tx_paths) == 1:
        preferred = [
            "horizontal_scattering_homogeneous.out",
            "horizontal_scattering_homogeneous_merged.out",
        ]
        by_name = {p.name: p for p in baseline_paths}
        for name in preferred:
            if name in by_name:
                return [by_name[name]]
        return [baseline_paths[0]]
    if len(baseline_paths) == 1:
        return [baseline_paths[0]] * len(tx_paths)
    if len(baseline_paths) == len(tx_paths):
        return baseline_paths
    print(
        "Warning: Baseline count does not match Tx count. "
        "Proceeding without baseline subtraction."
    )
    return [None] * len(tx_paths)


def _preprocess_traces(traces: np.ndarray, dt: float, mute_time_s: float, eps: float) -> np.ndarray:
    """Suppress direct-wave clutter and balance amplitudes before MUSIC."""
    proc = traces.copy()

    # Remove per-trace DC component.
    proc -= np.mean(proc, axis=-1, keepdims=True)

    # Mute early-time direct arrivals with cosine taper.
    t = np.arange(proc.shape[-1]) * dt
    mute_mask = t < mute_time_s
    if np.any(mute_mask):
        taper = 0.5 * (1.0 - np.cos(np.pi * t[mute_mask] / mute_time_s))
        proc[:, :, mute_mask] *= taper

    # Mild time gain control boosts later weak diffractions.
    t_norm = t / (t[-1] + eps)
    gain = 1.0 + 2.0 * t_norm
    proc *= gain[None, None, :]

    # Normalize each array so no trace dominates, but keep relative amp intact.
    rms = np.sqrt(np.mean(proc**2, axis=(1, 2), keepdims=True) + eps)
    proc /= rms
    return proc


def _box_blur_axis(arr: np.ndarray, axis: int, window: int) -> np.ndarray:
    if window <= 1:
        return arr.copy()
    pad = window // 2
    kernel = np.ones(window, dtype=float) / float(window)

    def _conv_1d(vec: np.ndarray) -> np.ndarray:
        vec_pad = np.pad(vec, (pad, pad), mode="edge")
        return np.convolve(vec_pad, kernel, mode="valid")

    return np.apply_along_axis(_conv_1d, axis, arr)


def _local_contrast_ratio(image: np.ndarray, wx: int, wy: int, eps: float) -> np.ndarray:
    smooth = _box_blur_axis(image, axis=0, window=max(3, int(wy) | 1))
    smooth = _box_blur_axis(smooth, axis=1, window=max(3, int(wx) | 1))
    return image / (smooth + eps)


def _estimate_signal_rank(S: np.ndarray, target_rank: int, eps: float = 1e-12) -> int:
    """Estimate signal subspace rank from singular value decay.
    
    Aggressively seeks target_rank; only stops early if very strong gap (>0.3) detected.
    """
    if len(S) <= 1:
        return max(1, target_rank)
    
    S_norm = S / (np.max(np.abs(S)) + eps)
    max_possible = len(S) - 2
    
    # Try to reach target rank; only short-circuit on very strong gaps
    if len(S) >= 3:
        for r in range(1, min(target_rank + 1, max_possible)):
            if r < len(S_norm):
                gap = S_norm[r - 1] - S_norm[r]
                # Only early-exit on very strong gap (not at 0.1)
                if gap > 0.3:  # Increased from 0.1 to 0.3
                    return max(target_rank, min(r + 4, max_possible))
    
    # Default: use target_rank capped by available dimensions
    return max(1, min(target_rank, max_possible))


def _single_tx_coherence_migration(
    traces_tx: np.ndarray,
    rx_x: np.ndarray,
    tx_x: float,
    tx_y: float,
    y_rx: float,
    grid_x: np.ndarray,
    grid_y: np.ndarray,
    v: float,
    dt: float,
    eps: float,
    max_angle_deg: float = 32.0,
    min_depth_m: float = 0.05,
) -> np.ndarray:
    """Coherence-factor weighted diffraction stacking for single-Tx data."""
    num_rx, num_t = traces_tx.shape
    img = np.zeros((len(grid_y), len(grid_x)), dtype=float)
    illum = np.zeros((len(grid_y), len(grid_x)), dtype=float)

    # Envelope-like trace transform to suppress phase oscillation smiles.
    trace_env = np.abs(traces_tx)
    trace_env = _box_blur_axis(trace_env, axis=1, window=9)
    tan_max = np.tan(np.deg2rad(max_angle_deg))

    for iy, y in enumerate(grid_y):
        depth = np.maximum(y_rx - y, eps)
        depth_taper = 1.0 / (1.0 + np.exp(-(depth - min_depth_m) / max(0.005, 0.2 * min_depth_m + eps)))

        d_rx = np.sqrt((grid_x[:, None] - rx_x[None, :]) ** 2 + (y - y_rx) ** 2)
        d_tx = np.sqrt((grid_x - tx_x) ** 2 + (y - tx_y) ** 2)
        t_total = (d_rx + d_tx[:, None]) / v
        idx_f = t_total / dt

        aperture_half = np.maximum(depth * tan_max, 0.02)
        tx_ap = np.abs(grid_x - tx_x) / aperture_half
        w_ap_tx = np.exp(-(tx_ap**4))

        stack = np.zeros(len(grid_x), dtype=float)
        abs_stack = np.zeros(len(grid_x), dtype=float)

        for ir in range(num_rx):
            idx = idx_f[:, ir]
            i0 = np.floor(idx).astype(int)
            frac = idx - i0
            valid = (i0 >= 0) & (i0 < (num_t - 1))
            rx_ap = np.abs(grid_x - rx_x[ir]) / aperture_half
            w_ap_rx = np.exp(-(rx_ap**4))

            vals = np.zeros(len(grid_x), dtype=float)
            i0v = i0[valid]
            vals[valid] = (1.0 - frac[valid]) * trace_env[ir, i0v] + frac[valid] * trace_env[ir, i0v + 1]

            # Obliquity weighting attenuates steep-angle migration smiles.
            ob_rx = depth / np.maximum(d_rx[:, ir], eps)
            ob_tx = depth / np.maximum(d_tx, eps)
            w = np.clip(ob_rx * ob_tx, 0.0, 1.0)
            vals *= (w * w_ap_rx * w_ap_tx * depth_taper)

            stack += vals
            abs_stack += np.abs(vals)
            illum[iy, :] += (w * w_ap_rx * w_ap_tx * depth_taper)

        coherence = (stack**2) / (num_rx * (abs_stack**2) + eps)
        img[iy, :] = np.abs(stack) * np.clip(coherence, 0.0, 1.0)

    illum_norm = illum / (np.max(illum) + eps)
    illum_floor = 0.15
    img /= np.maximum(illum_norm, illum_floor)
    return img


def apply_trmusic_diffraction(
    tx_glob: str = "Diffraction_Experiment/horizontal_scattering_0p5lambda*.out",
    baseline_glob: str = "Diffraction_Experiment/horizontal_scattering_homogeneous*.out",
    out_png: str = "Diffraction_Experiment/trmusic_diffraction_multi_tx.png",
    num_scatterers: int = 2,
    top_freq_fraction: float = 0.6,
    use_local_contrast: bool = True,
    contrast_window_x: int = 31,
    contrast_window_y: int = 31,
    single_tx_mode: str = "auto",
    migration_max_angle_deg: float = 32.0,
    migration_min_depth_m: float = 0.05,
):
    tx_paths = _discover_files(tx_glob, must_not_contain="homogeneous")
    if len(tx_paths) == 0:
        raise FileNotFoundError(f"No Tx files found for pattern: {tx_glob}")

    baseline_paths = _discover_files(baseline_glob)
    paired_baselines = _pair_baselines(tx_paths, baseline_paths)

    traces_list = []
    tx_pos_list = []
    dt = None
    rx_pos_ref = None
    num_t_ref = None

    for tx_file, base_file in zip(tx_paths, paired_baselines):
        dt_i, traces_i, rx_pos_i, tx_pos_i = _load_out_file(tx_file)

        if dt is None:
            dt = dt_i
            rx_pos_ref = rx_pos_i
            num_t_ref = traces_i.shape[1]
        else:
            if abs(dt_i - dt) > 1e-18:
                raise ValueError(f"dt mismatch in file: {tx_file}")
            if traces_i.shape[1] != num_t_ref:
                raise ValueError(f"time-length mismatch in file: {tx_file}")
            if not np.allclose(rx_pos_i, rx_pos_ref):
                raise ValueError(f"receiver geometry mismatch in file: {tx_file}")

        if base_file is not None:
            _, traces_b, rx_pos_b, _ = _load_out_file(base_file)
            if traces_b.shape != traces_i.shape:
                raise ValueError(f"baseline shape mismatch: {base_file}")
            if not np.allclose(rx_pos_b, rx_pos_i):
                raise ValueError(f"baseline receiver geometry mismatch: {base_file}")
            traces_i = traces_i - traces_b

        traces_list.append(traces_i)
        tx_pos_list.append(tx_pos_i)

    traces = np.asarray(traces_list, dtype=float)
    tx_pos = np.asarray(tx_pos_list, dtype=float)

    num_tx, num_rx, num_t = traces.shape

    eps = 1e-12
    traces = _preprocess_traces(traces, dt=dt, mute_time_s=2.0e-9, eps=eps)

    D_f = np.fft.rfft(traces, axis=-1)
    freqs = np.fft.rfftfreq(num_t, d=dt)

    f_min = 0.5e9
    f_max = 3.5e9
    f_idx_valid = np.where((freqs >= f_min) & (freqs <= f_max))[0]
    if len(f_idx_valid) == 0:
        raise RuntimeError("No valid frequency bins in requested band.")

    # Prefer high-SNR bins to avoid flattening from noisy frequencies.
    bin_energy = np.mean(np.abs(D_f[:, :, f_idx_valid]) ** 2, axis=(0, 1))
    keep_n = max(8, int(np.ceil(len(f_idx_valid) * top_freq_fraction)))
    keep_n = min(keep_n, len(f_idx_valid))
    local_idx = np.argsort(bin_energy)[-keep_n:]
    f_idx_use = f_idx_valid[np.sort(local_idx)]
    freq_weights = bin_energy[np.sort(local_idx)]
    freq_weights = freq_weights / (np.sum(freq_weights) + eps)

    c0 = 299792458.0
    eps_r = 6.0
    v = c0 / np.sqrt(eps_r)

    x_min = float(np.min(rx_pos_ref[:, 0])) + 0.02
    x_max = float(np.max(rx_pos_ref[:, 0])) - 0.02
    y_rx = float(np.mean(rx_pos_ref[:, 1]))

    # Search below the antennas in this 2D model (z is constant in gprMax output).
    y_min = 0.02
    y_max = max(0.04, y_rx - 0.03)
    grid_x = np.linspace(x_min, x_max, 180)
    grid_y = np.linspace(y_min, y_max, 160)

    image = np.zeros((len(grid_y), len(grid_x)), dtype=float)

    rx_x = rx_pos_ref[:, 0]
    tx_x_arr = tx_pos[:, 0]
    tx_y_arr = tx_pos[:, 1]

    fallback_mode = single_tx_mode
    if fallback_mode == "auto":
        fallback_mode = "migration"

    print(
        f"Using {num_tx} Tx file(s), {num_rx} Rx, "
        f"{len(f_idx_use)}/{len(f_idx_valid)} selected frequency bins."
    )

    used_mode = "trmusic"
    if num_tx == 1 and fallback_mode == "migration":
        used_mode = "migration"
        image = _single_tx_coherence_migration(
            traces_tx=traces[0],
            rx_x=rx_x,
            tx_x=float(tx_x_arr[0]),
            tx_y=float(tx_y_arr[0]),
            y_rx=y_rx,
            grid_x=grid_x,
            grid_y=grid_y,
            v=v,
            dt=dt,
            eps=eps,
            max_angle_deg=migration_max_angle_deg,
            min_depth_m=migration_min_depth_m,
        )
    else:
        if num_tx == 1:
            used_mode = "music"

        for wi, fi in enumerate(f_idx_use):
            f = freqs[fi]
            k = 2.0 * np.pi * f / v

            # K has shape (Nrx, Ntx) for TR-MUSIC.
            K = D_f[:, :, fi].T

            if num_tx > 1:
                U, S, Vh = np.linalg.svd(K, full_matrices=False)
                signal_dim = _estimate_signal_rank(S, target_rank=num_scatterers, eps=eps)
                signal_dim = max(2, min(signal_dim, len(S) - 2))
                Un = U[:, signal_dim:]
                Vn = Vh.conj().T[:, signal_dim:]
            else:
                # Single-Tx fallback: use receive-only noise subspace.
                R = (K @ K.conj().T) / max(num_tx, 1)
                trace_R = np.trace(R).real
                R += (1e-4 * trace_R / max(num_rx, 1) + eps) * np.eye(num_rx)
                evals, evecs = np.linalg.eigh(R)
                order = np.argsort(evals)[::-1]
                evecs = evecs[:, order]
                signal_dim = 1
                Un = evecs[:, signal_dim:]
                Vn = None

            I_f = np.zeros((len(grid_y), len(grid_x)), dtype=float)

            for iy, y in enumerate(grid_y):
                d_rx = np.sqrt((grid_x[:, None] - rx_x[None, :]) ** 2 + (y - y_rx) ** 2)
                d_tx = np.sqrt((grid_x[:, None] - tx_x_arr[None, :]) ** 2 + (y - tx_y_arr[None, :]) ** 2)

                g_rx = np.exp(-1j * k * d_rx) / np.sqrt(np.maximum(d_rx, eps))
                g_tx = np.exp(-1j * k * d_tx) / np.sqrt(np.maximum(d_tx, eps))

                g_rx /= np.linalg.norm(g_rx, axis=1, keepdims=True) + eps
                g_tx /= np.linalg.norm(g_tx, axis=1, keepdims=True) + eps

                proj_rx = np.sum(np.abs(np.conj(g_rx) @ Un) ** 2, axis=1)
                if Vn is not None:
                    # Using additive norm for stability instead of multiplicative
                    proj_tx = np.sum(np.abs(np.conj(g_tx) @ Vn) ** 2, axis=1)
                    I_f[iy, :] = 1.0 / (proj_rx + proj_tx + eps)
                else:
                    I_f[iy, :] = 1.0 / (proj_rx + eps)

            image += freq_weights[wi] * I_f

    image_to_plot = image
    if use_local_contrast:
        image_to_plot = _local_contrast_ratio(
            image,
            wx=contrast_window_x,
            wy=contrast_window_y,
            eps=eps,
        )

    ref = np.percentile(image_to_plot, 99.7)
    image_db = 10.0 * np.log10((image_to_plot + eps) / (ref + eps))

    plt.figure(figsize=(10, 6))
    plt.imshow(
        image_db,
        extent=[grid_x[0], grid_x[-1], grid_y[0], grid_y[-1]],
        origin="lower",
        aspect="auto",
        cmap="plasma",
        vmin=-35,
        vmax=0,
    )
    plt.colorbar(label="Pseudo-spectrum (dB)")
    if used_mode == "trmusic":
        title = "TR-MUSIC Imaging for Diffraction Experiment (Multi-Tx)"
    elif used_mode == "migration":
        title = "Coherence Migration Fallback (Single-Tx)"
    else:
        title = "MUSIC Imaging Fallback (Single-Tx)"
    if use_local_contrast:
        title += " + Local Contrast"
    plt.title(title)
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.tight_layout()

    out_path = Path(out_png)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved TR-MUSIC image to: {out_path}")


def _build_arg_parser():
    parser = argparse.ArgumentParser(description="Multi-Tx TR-MUSIC for diffraction experiment outputs.")
    parser.add_argument(
        "--tx-glob",
        default="Diffraction_Experiment/horizontal_scattering_0p5lambda*.out",
        help="Glob for diffraction Tx output files.",
    )
    parser.add_argument(
        "--baseline-glob",
        default="Diffraction_Experiment/horizontal_scattering_homogeneous*.out",
        help="Glob for baseline output files. Use empty string to disable subtraction.",
    )
    parser.add_argument(
        "--out",
        default="Diffraction_Experiment/trmusic_diffraction_multi_tx.png",
        help="Output image path.",
    )
    parser.add_argument(
        "--num-scatterers",
        type=int,
        default=2,
        help="Expected number of dominant scatterers for subspace split.",
    )
    parser.add_argument(
        "--top-freq-fraction",
        type=float,
        default=0.6,
        help="Fraction of strongest in-band frequency bins used for imaging (0,1].",
    )
    parser.add_argument(
        "--no-local-contrast",
        action="store_true",
        help="Disable local background-normalization contrast enhancement.",
    )
    parser.add_argument(
        "--contrast-window-x",
        type=int,
        default=31,
        help="Odd window size in x for local contrast background estimation.",
    )
    parser.add_argument(
        "--contrast-window-y",
        type=int,
        default=31,
        help="Odd window size in y for local contrast background estimation.",
    )
    parser.add_argument(
        "--single-tx-mode",
        choices=["auto", "migration", "music"],
        default="auto",
        help="Method used when only one Tx file is available.",
    )
    parser.add_argument(
        "--migration-max-angle-deg",
        type=float,
        default=32.0,
        help="Maximum aperture angle for single-Tx migration fallback.",
    )
    parser.add_argument(
        "--migration-min-depth-m",
        type=float,
        default=0.05,
        help="Minimum depth to image for single-Tx migration fallback.",
    )
    return parser


if __name__ == "__main__":
    parser = _build_arg_parser()
    args = parser.parse_args()

    baseline_glob = args.baseline_glob if args.baseline_glob.strip() else "__none__"
    if baseline_glob == "__none__":
        baseline_glob = ""

    apply_trmusic_diffraction(
        tx_glob=args.tx_glob,
        baseline_glob=baseline_glob,
        out_png=args.out,
        num_scatterers=max(1, int(args.num_scatterers)),
        top_freq_fraction=max(0.1, min(1.0, float(args.top_freq_fraction))),
        use_local_contrast=not bool(args.no_local_contrast),
        contrast_window_x=max(3, int(args.contrast_window_x)),
        contrast_window_y=max(3, int(args.contrast_window_y)),
        single_tx_mode=args.single_tx_mode,
        migration_max_angle_deg=max(10.0, min(80.0, float(args.migration_max_angle_deg))),
        migration_min_depth_m=max(0.0, float(args.migration_min_depth_m)),
    )
