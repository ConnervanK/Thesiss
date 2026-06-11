"""
GPU runner for TimeLapse_Notebooks gprMax simulations.

Usage — single-trace (Common Shot Gather):
    python run_experiment.py
    python run_experiment.py --prefix fracture_10cm --subdir run_01

Usage — B-scan (Common-Offset, multiple traces):
    python run_experiment.py --bscan --n-traces 350
    python run_experiment.py --bscan --n-traces 350 --prefix run1 --subdir bscan_01

In B-scan mode gprMax writes one .out file per trace; this script merges them
into <stem>_merged.out and removes the numbered files automatically.
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
from pathlib import Path

VCVARS = r"C:\Progra~2\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
NOTEBOOK_DIR = Path(__file__).parent


def _run_gpu(in_file: Path, n: int = 1) -> None:
    py_exe = sys.executable
    ps = f"""
    $envDump = cmd /c '"{VCVARS}" && set'
    $envDump | ForEach-Object {{
        $p = $_ -split '=',2
        if ($p.Length -eq 2) {{ Set-Item -Path env:$($p[0]) -Value $p[1] }}
    }}
    & "{py_exe}" -m gprMax "{in_file}" -n {n} -gpu
    if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
    """
    subprocess.run(["powershell", "-Command", ps], check=True)


def _merge_bscan(in_file: Path) -> None:
    """Merge the per-trace .out files gprMax produces for a B-scan run.

    Replicates the logic of gprMax/tools/outputfiles_merge.py using h5py
    directly so no external tool import is required.
    Output: <stem>_merged.out  (same format, field arrays shaped [timesteps, traces]).
    """
    import glob as _glob
    import re
    import h5py

    base = str(in_file.with_suffix(""))
    files = sorted(
        [f for f in _glob.glob(base + "[0-9]*.out") if "_merged" not in f],
        key=lambda f: int(re.search(r"\d+", f[len(base):]).group()),
    )
    if not files:
        raise FileNotFoundError(f"No numbered .out files found matching {base}*.out")

    n_traces   = len(files)
    merged_out = base + "_merged.out"

    with h5py.File(merged_out, "w") as fm:
        # --- header and pre-allocated datasets from trace 1 ---
        with h5py.File(files[0], "r") as f0:
            for k, v in f0.attrs.items():
                fm.attrs[k] = v
            nrx     = int(f0.attrs["nrx"])
            n_iters = int(f0.attrs["Iterations"])
            for rx in range(1, nrx + 1):
                path = f"/rxs/rx{rx}"
                grp  = fm.create_group(path)
                for comp in f0[path].keys():
                    grp.create_dataset(comp, (n_iters, n_traces),
                                       dtype=f0[f"{path}/{comp}"].dtype)
        # --- fill one column per trace ---
        for col, fpath in enumerate(files):
            with h5py.File(fpath, "r") as fi:
                for rx in range(1, nrx + 1):
                    path = f"/rxs/rx{rx}"
                    for comp in fi[path].keys():
                        fm[f"{path}/{comp}"][:, col] = fi[f"{path}/{comp}"][:]

    for fpath in files:
        os.remove(fpath)

    print(f"  Merged {n_traces} traces -> {in_file.stem}_merged.out")


def run_experiment(
    prefix: str = "",
    subdir: str | None = None,
    bscan: bool = False,
    n_traces: int = 1,
) -> None:
    def fname(suffix: str) -> str:
        return f"{prefix}_{suffix}.in" if prefix else f"{suffix}.in"

    baseline_in  = NOTEBOOK_DIR / fname("baseline")
    timelapse_in = NOTEBOOK_DIR / fname("timelapse")

    mode_label = f"B-scan ({n_traces} traces)" if bscan else "single-trace"
    print(f"=== TimeLapse gprMax GPU Runner  [{mode_label}] ===")
    for f in [baseline_in, timelapse_in]:
        if not f.exists():
            raise FileNotFoundError(f"Input file not found: {f}")

    n = n_traces if bscan else 1

    print(f"[1/2] Baseline:  {baseline_in.name}")
    _run_gpu(baseline_in, n=n)
    if bscan:
        _merge_bscan(baseline_in)

    print(f"[2/2] Timelapse: {timelapse_in.name}")
    _run_gpu(timelapse_in, n=n)
    if bscan:
        _merge_bscan(timelapse_in)

    print("=== Simulations complete ===")

    if subdir:
        dest = NOTEBOOK_DIR / "outputs" / subdir
        dest.mkdir(parents=True, exist_ok=True)

        for in_file in [baseline_in, timelapse_in]:
            # B-scan produces a merged file; single-trace produces a plain .out
            candidates = [
                in_file.parent / f"{in_file.stem}_merged.out",
                in_file.with_suffix(".out"),
            ]
            for out_file in candidates:
                if out_file.exists():
                    shutil.move(str(out_file), str(dest / out_file.name))
                    print(f"Moved: {out_file.name} -> outputs/{subdir}/")
                    break

        for path in glob.glob(str(NOTEBOOK_DIR / "*_snaps*")):
            if os.path.isdir(path):
                shutil.move(path, str(dest / os.path.basename(path)))
                print(f"Moved: {os.path.basename(path)}/ -> outputs/{subdir}/")

        print(f"Outputs in: {dest}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run TimeLapse gprMax experiment on GPU",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--prefix",   default="",    help="File prefix (matches file_prefix in notebook)")
    parser.add_argument("--subdir",   default=None,  help="Output subdirectory under TimeLapse_Notebooks/outputs/")
    parser.add_argument("--bscan",    action="store_true", help="B-scan mode: run all traces and merge outputs")
    parser.add_argument("--n-traces", type=int, default=1,
                        help="Number of traces (printed by BScan_generator.ipynb Cell 2; required with --bscan)")
    args = parser.parse_args()

    if args.bscan and args.n_traces == 1:
        parser.error("--bscan requires --n-traces N  (N > 1)")

    run_experiment(
        prefix=args.prefix,
        subdir=args.subdir,
        bscan=args.bscan,
        n_traces=args.n_traces,
    )
