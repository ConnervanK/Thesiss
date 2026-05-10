"""
GPU runner for TimeLapse_Notebooks gprMax simulations.

Usage:
    python run_experiment.py
    python run_experiment.py --prefix fracture_10cm --subdir run_01
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


def _run_gpu(in_file: Path) -> None:
    py_exe = sys.executable
    ps = f"""
    $envDump = cmd /c '"{VCVARS}" && set'
    $envDump | ForEach-Object {{
        $p = $_ -split '=',2
        if ($p.Length -eq 2) {{ Set-Item -Path env:$($p[0]) -Value $p[1] }}
    }}
    & "{py_exe}" -m gprMax "{in_file}" -n 1 -gpu
    if ($LASTEXITCODE -ne 0) {{ exit $LASTEXITCODE }}
    """
    subprocess.run(["powershell", "-Command", ps], check=True)


def run_experiment(prefix: str = "", subdir: str | None = None) -> None:
    def fname(suffix: str) -> str:
        return f"{prefix}_{suffix}.in" if prefix else f"{suffix}.in"

    baseline_in  = NOTEBOOK_DIR / fname("baseline")
    timelapse_in = NOTEBOOK_DIR / fname("timelapse")

    print("=== TimeLapse gprMax GPU Runner ===")
    for f in [baseline_in, timelapse_in]:
        if not f.exists():
            raise FileNotFoundError(f"Input file not found: {f}")

    print(f"[1/2] Baseline:  {baseline_in.name}")
    _run_gpu(baseline_in)
    print(f"[2/2] Timelapse: {timelapse_in.name}")
    _run_gpu(timelapse_in)
    print("=== Simulations complete ===")

    if subdir:
        dest = NOTEBOOK_DIR / "outputs" / subdir
        dest.mkdir(parents=True, exist_ok=True)

        for in_file in [baseline_in, timelapse_in]:
            out_file = in_file.with_suffix(".out")
            if out_file.exists():
                shutil.move(str(out_file), str(dest / out_file.name))
                print(f"Moved: {out_file.name} -> outputs/{subdir}/")

        for path in glob.glob(str(NOTEBOOK_DIR / "*_snaps*")):
            if os.path.isdir(path):
                shutil.move(path, str(dest / os.path.basename(path)))
                print(f"Moved: {os.path.basename(path)}/ -> outputs/{subdir}/")

        print(f"Outputs in: {dest}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run TimeLapse gprMax experiment on GPU")
    parser.add_argument("--prefix", default="",   help="File prefix (matches file_prefix in notebook)")
    parser.add_argument("--subdir", default=None, help="Output subdirectory under TimeLapse_Notebooks/outputs/")
    args = parser.parse_args()
    run_experiment(prefix=args.prefix, subdir=args.subdir)
