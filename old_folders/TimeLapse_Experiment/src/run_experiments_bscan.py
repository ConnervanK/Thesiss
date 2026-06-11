import os
import old_folders.TimeLapse_Experiment.src.forward as forward


class GPRBScanExperiment:
    """Configures and runs a single GPR B-scan forward modeling experiment."""
    def __init__(self, target_dir, model_params):
        self.target_dir = target_dir
        self.model_params = model_params

    def run(self):
        os.makedirs(os.path.join(self.target_dir, "data", "outputs"), exist_ok=True)
        os.chdir(self.target_dir)

        output_context = forward.create_bscan_input_file(**self.model_params)
        block_label = output_context.get('block_label', 'default')
        n_traces = output_context.get('n_traces', 50)

        print(f"Running B-scan simulation for {block_label} ({n_traces} traces)...")
        forward.run_gprmax_bscan_to_subdir(n_traces=n_traces, output_subdir=block_label)
        print(f"Finished B-scan simulation for {block_label}.")


def main():
    target_dir = r"c:\Users\Administrator\OneDrive\Thesis\TimeLapse_Experiment"

    base_params = {
        'f_central': 1.5 * 1e9,
        'c': 3 * 1e8,
        'permittivity_ice': 3.15,
        'permittivity_air': 1,
        'permittivity_fracture': 80,
        'conductivity_ice': 1e-6,
        'conductivity_air': 0,
        'conductivity_fracture': 0.01,
        'fracture_depth': 0.6,
        'depth_below_fracture': 0.1,
        'air_thickness': 0.1,

        # B-scan geometry: midpoints sweep the same positions as the original receivers
        'rx_count': 50,
        'rx_spread': 2.0,
        'tx_rx_offset': 0.1,     # fixed separation between Tx and Rx (m)

        # Block size settings
        'block_size_fraction': 1/8,

        # Time-lapse settings
        'time_lapse_shift': 0.05,
    }

    experiment = GPRBScanExperiment(target_dir, base_params)
    experiment.run()


if __name__ == "__main__":
    main()
