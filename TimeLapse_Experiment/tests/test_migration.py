"""
Tests that Kirchhoff migration recovers a single point scatterer.

Traces are generated analytically: for a scatterer at a known (x, z), the
two-way travel time to each Tx/Rx pair is computed and a Ricker wavelet is
placed at that time in each trace. Migration should then collapse the diffraction
hyperbola back to the true scatterer location.

Two geometries are tested:
  - Common-shot gather: fixed Tx, multiple Rx positions  →  kirchhoff_migration()
  - Common-offset B-scan: moving Tx-Rx pair            →  GPRBScanData.kirchhoff_migration_bscan()

Run with pytest:
    cd TimeLapse_Experiment
    pytest tests/test_migration.py -v

Run standalone (also saves diagnostic plots to tests/outputs/):
    python tests/test_migration.py
"""

import sys
import os

# Make src/ importable regardless of where the script is invoked from
SRC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src")
sys.path.insert(0, SRC_DIR)

import numpy as np

from utils.processing import kirchhoff_migration
from utils.processing_bscan import GPRBScanData


# ---------------------------------------------------------------------------
# Shared physical constants  (match experiment settings)
# ---------------------------------------------------------------------------

F_CENTRAL = 1.5e9                           # Hz
C = 3e8                                     # m/s
PERMITTIVITY_ICE = 3.15
V_ICE = C / np.sqrt(PERMITTIVITY_ICE)       # ~1.69e8 m/s

# True scatterer position
SCATTERER_X = 2.0    # m  (centre of the survey profile)
SCATTERER_Z = 0.5    # m

# Survey geometry
N_TRACES = 40
RX_SPREAD = 2.0                                              # m
RX_X_ARRAY = np.linspace(1.0, 3.0, N_TRACES)                # common-shot receivers / B-scan midpoints
TX_X = 2.0                                                   # common-shot transmitter
TX_RX_OFFSET = 0.1                                           # B-scan Tx–Rx separation (m)

# Time axis  (40 ps step → ~15 samples per Ricker period at 1.5 GHz; 600 samples → 24 ns)
DT = 4e-11
N_TIME = 600
TIME_S = np.arange(N_TIME) * DT        # seconds  (used by B-scan class)
TIME_NS = TIME_S * 1e9                 # nanoseconds  (used by standalone kirchhoff_migration)

# Migration parameters
DZ = 0.005   # m
MAX_DEPTH = 0.8  # m

# Acceptance tolerances for recovered peak location
Z_TOL = 0.05   # m  (10 depth cells at dz=0.005)
X_TOL = 0.20   # m  (~4 receiver spacings at dx≈0.051)


# ---------------------------------------------------------------------------
# Synthetic trace generation
# ---------------------------------------------------------------------------

def ricker_wavelet(t, t0, f):
    """Ricker wavelet centred at t0 (t and t0 in the same units; f in Hz·unit⁻¹)."""
    tau = t - t0
    u = (np.pi * f * tau) ** 2
    return (1.0 - 2.0 * u) * np.exp(-u)


def make_common_shot_traces():
    """
    Common-shot gather: Tx fixed at TX_X, receivers at RX_X_ARRAY.
    Returns traces (N_TRACES × N_TIME) with the scatterer diffraction hyperbola.
    """
    traces = np.zeros((N_TRACES, N_TIME))
    dist_tx = np.sqrt((SCATTERER_X - TX_X) ** 2 + SCATTERER_Z ** 2)
    for i, rx_x in enumerate(RX_X_ARRAY):
        dist_rx = np.sqrt((SCATTERER_X - rx_x) ** 2 + SCATTERER_Z ** 2)
        t_arrival = (dist_tx + dist_rx) / V_ICE   # seconds
        traces[i] = ricker_wavelet(TIME_S, t_arrival, F_CENTRAL)
    return traces


def make_bscan_traces():
    """
    Common-offset B-scan: midpoints at RX_X_ARRAY, fixed Tx–Rx separation.
    Returns traces (N_TRACES × N_TIME) with the scatterer diffraction hyperbola.
    """
    traces = np.zeros((N_TRACES, N_TIME))
    for i, mp in enumerate(RX_X_ARRAY):
        tx_x = mp - TX_RX_OFFSET / 2.0
        rx_x = mp + TX_RX_OFFSET / 2.0
        dist_tx = np.sqrt((SCATTERER_X - tx_x) ** 2 + SCATTERER_Z ** 2)
        dist_rx = np.sqrt((SCATTERER_X - rx_x) ** 2 + SCATTERER_Z ** 2)
        t_arrival = (dist_tx + dist_rx) / V_ICE
        traces[i] = ricker_wavelet(TIME_S, t_arrival, F_CENTRAL)
    return traces


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def peak_location(migrated_img, z_array, x_array):
    """
    Return (x_peak, z_peak) of the maximum absolute amplitude in migrated_img.
    Raw |amplitude| is used for peak finding because the Ricker wavelet peaks
    exactly at the true scatterer location; the Hilbert envelope along depth
    can shift the peak due to asymmetric dT/dz rates across the aperture.
    migrated_img shape: (n_depth, n_x)
    """
    iz, ix = np.unravel_index(np.argmax(np.abs(migrated_img)), migrated_img.shape)
    return float(x_array[ix]), float(z_array[iz])


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_common_shot_migration_recovers_scatterer():
    """
    Common-shot Kirchhoff migration should place the peak at the true scatterer
    position (SCATTERER_X, SCATTERER_Z) within the specified tolerances.
    """
    traces = make_common_shot_traces()

    migrated_img, z_array = kirchhoff_migration(
        traces,
        time_array=TIME_NS,           # nanoseconds — matches internal unit of kirchhoff_migration
        rx_x_array=RX_X_ARRAY,
        tx_x=TX_X,
        velocity=V_ICE,
        max_depth=MAX_DEPTH,
        dz=DZ,
    )

    x_peak, z_peak = peak_location(migrated_img, z_array, RX_X_ARRAY)

    assert abs(z_peak - SCATTERER_Z) <= Z_TOL, (
        f"Common-shot z-peak {z_peak:.4f} m differs from true z={SCATTERER_Z} m "
        f"by more than {Z_TOL} m"
    )
    assert abs(x_peak - SCATTERER_X) <= X_TOL, (
        f"Common-shot x-peak {x_peak:.4f} m differs from true x={SCATTERER_X} m "
        f"by more than {X_TOL} m"
    )

    print(f"  Common-shot: recovered ({x_peak:.3f}, {z_peak:.3f}) m  "
          f"| true ({SCATTERER_X}, {SCATTERER_Z}) m  ✓")


def test_bscan_migration_recovers_scatterer():
    """
    Common-offset B-scan Kirchhoff migration should place the peak at the true
    scatterer position within the specified tolerances.
    """
    traces = make_bscan_traces()
    midpoints = RX_X_ARRAY
    midpoint_spacing = float(midpoints[1] - midpoints[0])

    # Build GPRBScanData without loading a file — inject synthetic time axis directly.
    bscan = GPRBScanData(
        baseline_file="",
        timelapse_file="",
        midpoint_start_x=float(midpoints[0]),
        midpoint_spacing=midpoint_spacing,
        n_traces=N_TRACES,
        tx_rx_offset=TX_RX_OFFSET,
    )
    bscan.time = TIME_S   # seconds — kirchhoff_migration_bscan converts internally to ns
    bscan.dt = DT
    bscan.n_traces = N_TRACES

    migrated_img, z_array = bscan.kirchhoff_migration_bscan(
        traces,
        velocity=V_ICE,
        max_depth=MAX_DEPTH,
        dz=DZ,
    )

    x_peak, z_peak = peak_location(migrated_img, z_array, midpoints)

    assert abs(z_peak - SCATTERER_Z) <= Z_TOL, (
        f"B-scan z-peak {z_peak:.4f} m differs from true z={SCATTERER_Z} m "
        f"by more than {Z_TOL} m"
    )
    assert abs(x_peak - SCATTERER_X) <= X_TOL, (
        f"B-scan x-peak {x_peak:.4f} m differs from true x={SCATTERER_X} m "
        f"by more than {X_TOL} m"
    )

    print(f"  B-scan:      recovered ({x_peak:.3f}, {z_peak:.3f}) m  "
          f"| true ({SCATTERER_X}, {SCATTERER_Z}) m  ✓")


# ---------------------------------------------------------------------------
# Standalone runner — also saves diagnostic plots
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")
    os.makedirs(out_dir, exist_ok=True)

    def _save_migration_plot(migrated_img, z_array, x_array, x_peak, z_peak, title, filename):
        from scipy.signal import hilbert as _hilbert
        env = np.abs(_hilbert(migrated_img, axis=0))
        vmax = np.percentile(env, 99)
        if vmax == 0:
            vmax = 1.0
        fig, ax = plt.subplots(figsize=(8, 5))
        im = ax.imshow(
            env,
            extent=[x_array[0], x_array[-1], z_array[-1], z_array[0]],
            aspect="auto",
            cmap="hot_r",
            vmin=0,
            vmax=vmax,
        )
        ax.axhline(SCATTERER_Z, color="cyan", linewidth=1.2, linestyle="--",
                   label=f"True z = {SCATTERER_Z} m")
        ax.axvline(SCATTERER_X, color="cyan", linewidth=1.2, linestyle=":",
                   label=f"True x = {SCATTERER_X} m")
        ax.plot(x_peak, z_peak, "b*", markersize=12, label=f"Peak ({x_peak:.3f}, {z_peak:.3f}) m")
        ax.set_xlabel("x (m)")
        ax.set_ylabel("Depth (m)")
        ax.set_title(title)
        ax.invert_yaxis()
        ax.legend(fontsize=8)
        plt.colorbar(im, ax=ax, label="Envelope amplitude")
        plt.tight_layout()
        plt.savefig(filename, dpi=120)
        plt.close()
        print(f"  Saved: {filename}")

    def _save_bscan_plot(traces, time_ns, x_array, title, filename):
        vmax = np.max(np.abs(traces))
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.imshow(
            traces.T,
            extent=[x_array[0], x_array[-1], time_ns[-1], time_ns[0]],
            aspect="auto",
            cmap="seismic",
            vmin=-vmax,
            vmax=vmax,
        )
        ax.set_xlabel("Midpoint x (m)")
        ax.set_ylabel("Two-way travel time (ns)")
        ax.set_title(title)
        plt.tight_layout()
        plt.savefig(filename, dpi=120)
        plt.close()
        print(f"  Saved: {filename}")

    print("=== Common-shot gather test ===")
    cs_traces = make_common_shot_traces()
    _save_bscan_plot(cs_traces, TIME_NS, RX_X_ARRAY,
                     "Synthetic common-shot gather — single scatterer",
                     os.path.join(out_dir, "cs_gather.png"))

    cs_migrated, z_array = kirchhoff_migration(
        cs_traces,
        time_array=TIME_NS,
        rx_x_array=RX_X_ARRAY,
        tx_x=TX_X,
        velocity=V_ICE,
        max_depth=MAX_DEPTH,
        dz=DZ,
    )
    cs_x, cs_z = peak_location(cs_migrated, z_array, RX_X_ARRAY)
    _save_migration_plot(cs_migrated, z_array, RX_X_ARRAY, cs_x, cs_z,
                         "Common-shot Kirchhoff migration",
                         os.path.join(out_dir, "cs_migration.png"))
    test_common_shot_migration_recovers_scatterer()

    print("\n=== B-scan test ===")
    bs_traces = make_bscan_traces()
    _save_bscan_plot(bs_traces, TIME_NS, RX_X_ARRAY,
                     "Synthetic B-scan — single scatterer",
                     os.path.join(out_dir, "bscan_gather.png"))

    midpoints = RX_X_ARRAY
    bscan = GPRBScanData(
        baseline_file="",
        timelapse_file="",
        midpoint_start_x=float(midpoints[0]),
        midpoint_spacing=float(midpoints[1] - midpoints[0]),
        n_traces=N_TRACES,
        tx_rx_offset=TX_RX_OFFSET,
    )
    bscan.time = TIME_S
    bscan.dt = DT
    bscan.n_traces = N_TRACES

    bs_migrated, z_array = bscan.kirchhoff_migration_bscan(
        bs_traces,
        velocity=V_ICE,
        max_depth=MAX_DEPTH,
        dz=DZ,
    )
    bs_x, bs_z = peak_location(bs_migrated, z_array, midpoints)
    _save_migration_plot(bs_migrated, z_array, midpoints, bs_x, bs_z,
                         "B-scan Kirchhoff migration",
                         os.path.join(out_dir, "bscan_migration.png"))
    test_bscan_migration_recovers_scatterer()

    print("\nAll tests passed.")
