"""One-off script: merge already-generated B-scan output files and move them to the outputs folder."""
import os
import sys
import shutil

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forward

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
configs_dir = os.path.join(base_dir, "configs")

forward.merge_bscan_outputs(configs_dir, [
    "horizontal_scattering_baseline_bscan",
    "horizontal_scattering_timelapse_bscan",
])

block_label = "1_8Lambda"
dest_dir = os.path.join(base_dir, "data", "outputs", block_label)
os.makedirs(dest_dir, exist_ok=True)

for fname in [
    "horizontal_scattering_baseline_bscan.out",
    "horizontal_scattering_timelapse_bscan.out",
]:
    src = os.path.join(configs_dir, fname)
    if os.path.exists(src):
        shutil.move(src, os.path.join(dest_dir, fname))
        print(f"Moved: {fname} -> {dest_dir}")
    else:
        print(f"Not found after merge: {src}")

print("Done.")
