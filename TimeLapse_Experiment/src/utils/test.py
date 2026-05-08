import sys
sys.path.insert(0, r"c:\Users\Administrator\OneDrive\Thesis\TimeLapse_Experiment\src")
import forward, os

configs_dir = r"c:\Users\Administrator\OneDrive\Thesis\TimeLapse_Experiment\configs"
forward.merge_bscan_outputs(configs_dir, [
    "horizontal_scattering_baseline_bscan",
    "horizontal_scattering_timelapse_bscan",
])

# Then move to the output folder
import shutil
dest = r"c:\Users\Administrator\OneDrive\Thesis\TimeLapse_Experiment\data\outputs\1_8Lambda"
os.makedirs(dest, exist_ok=True)
for f in ["horizontal_scattering_baseline_bscan.out", "horizontal_scattering_timelapse_bscan.out"]:
    src = os.path.join(configs_dir, f)
    if os.path.exists(src):
        shutil.move(src, os.path.join(dest, f))
