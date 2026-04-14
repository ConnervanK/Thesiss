import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _extract_rx_x_list(in_text: str):
    rx_x = []
    for line in in_text.splitlines():
        if line.startswith("#rx:"):
            parts = line.split()
            if len(parts) >= 4:
                rx_x.append(float(parts[1]))
    if not rx_x:
        raise ValueError("No receiver lines (#rx:) found in template.")
    return rx_x


def _replace_hertzian_x(in_text: str, new_x: float) -> str:
    pattern = r"^#hertzian_dipole:\s+z\s+([0-9eE+\-.]+)\s+([0-9eE+\-.]+)\s+0\s+my_ricker\s*$"
    lines = in_text.splitlines()
    replaced = False
    for i, line in enumerate(lines):
        m = re.match(pattern, line)
        if m:
            y = m.group(2)
            lines[i] = f"#hertzian_dipole: z {new_x:.15g} {y} 0 my_ricker"
            replaced = True
            break
    if not replaced:
        raise ValueError("Could not find #hertzian_dipole line in template.")
    return "\n".join(lines) + "\n"


def _run_cmd(cmd, cwd: Path):
    print(" ".join(str(c) for c in cmd))
    subprocess.run(cmd, cwd=str(cwd), check=True)


def main():
    parser = argparse.ArgumentParser(description="Generate multi-Tx diffraction acquisitions and run true TR-MUSIC.")
    parser.add_argument("--n-tx", type=int, default=8, help="Number of Tx positions.")
    parser.add_argument("--tx-margin", type=float, default=0.08, help="Margin from aperture ends in meters.")
    parser.add_argument("--num-scatterers", type=int, default=8, help="Signal subspace dimension for TR-MUSIC.")
    parser.add_argument("--skip-sim", action="store_true", help="Skip gprMax simulation and only run TR-MUSIC.")
    parser.add_argument("--gpu", action="store_true", help="Use gprMax GPU flag.")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    alt_template = root / "gapped_scattering.in"
    base_template = root / "horizontal_scattering_homogeneous.in"
    trmusic_script = root / "trmusic_diffraction.py"

    if not alt_template.exists() or not base_template.exists():
        raise FileNotFoundError("Template .in files not found in Diffraction_Experiment.")

    alt_text = _read_text(alt_template)
    base_text = _read_text(base_template)
    rx_x = _extract_rx_x_list(alt_text)

    x_min = min(rx_x) + args.tx_margin
    x_max = max(rx_x) - args.tx_margin
    if x_min >= x_max:
        raise ValueError("Invalid Tx range after margin. Reduce --tx-margin.")

    n_tx = max(2, int(args.n_tx))
    tx_positions = [x_min + i * (x_max - x_min) / (n_tx - 1) for i in range(n_tx)]

    work_dir = root / "multi_tx"
    in_alt_dir = work_dir / "inputs_alt"
    in_base_dir = work_dir / "inputs_base"
    out_dir = work_dir / "outputs"

    if not args.skip_sim:
        if work_dir.exists():
            shutil.rmtree(work_dir)
        in_alt_dir.mkdir(parents=True, exist_ok=True)
        in_base_dir.mkdir(parents=True, exist_ok=True)
        out_dir.mkdir(parents=True, exist_ok=True)

        for i, tx_x in enumerate(tx_positions, start=1):
            idx = f"{i:02d}"
            alt_in = in_alt_dir / f"horizontal_scattering_0p5lambda_tx{idx}.in"
            base_in = in_base_dir / f"horizontal_scattering_homogeneous_tx{idx}.in"
            _write_text(alt_in, _replace_hertzian_x(alt_text, tx_x))
            _write_text(base_in, _replace_hertzian_x(base_text, tx_x))

        gpr_flags = ["-n", "1"]
        if args.gpu:
            gpr_flags.append("-gpu")

        for i in range(1, n_tx + 1):
            idx = f"{i:02d}"
            alt_in = in_alt_dir / f"horizontal_scattering_0p5lambda_tx{idx}.in"
            base_in = in_base_dir / f"horizontal_scattering_homogeneous_tx{idx}.in"

            _run_cmd([sys.executable, "-m", "gprMax", str(alt_in)] + gpr_flags, cwd=root)
            _run_cmd([sys.executable, "-m", "gprMax", str(base_in)] + gpr_flags, cwd=root)

            alt_out_src = alt_in.with_suffix(".out")
            base_out_src = base_in.with_suffix(".out")
            alt_out_dst = out_dir / f"horizontal_scattering_0p5lambda_tx{idx}.out"
            base_out_dst = out_dir / f"horizontal_scattering_homogeneous_tx{idx}.out"
            shutil.copy2(alt_out_src, alt_out_dst)
            shutil.copy2(base_out_src, base_out_dst)

    tx_glob = str((out_dir / "horizontal_scattering_0p5lambda_tx*.out").as_posix())
    base_glob = str((out_dir / "horizontal_scattering_homogeneous_tx*.out").as_posix())
    out_png = str((out_dir / "trmusic_true_multi_tx.png").as_posix())

    _run_cmd(
        [
            sys.executable,
            str(trmusic_script),
            "--tx-glob",
            tx_glob,
            "--baseline-glob",
            base_glob,
            "--single-tx-mode",
            "music",
            "--num-scatterers",
            str(max(3, int(args.num_scatterers))),
            "--top-freq-fraction",
            "0.7",
            "--out",
            out_png,
        ],
        cwd=root,
    )

    print(f"Done. True multi-Tx TR-MUSIC output: {out_png}")


if __name__ == "__main__":
    main()
