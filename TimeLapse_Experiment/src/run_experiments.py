import os
import forward

class GPRExperiment:
    """Configures and runs a single GPR forward modeling experiment."""
    def __init__(self, target_dir, model_params):
        self.target_dir = target_dir
        self.model_params = model_params
        
    def run(self):
        os.makedirs(os.path.join(self.target_dir, "data", "outputs"), exist_ok=True)
        os.chdir(self.target_dir)

        output_context = forward.create_input_file(**self.model_params)
        block_label = output_context.get('block_label', 'default')

        print(f"Running simulation for {block_label}...")
        forward.run_gprmax_to_subdir(output_subdir=block_label)
        print(f"Finished simulation for {block_label}.")
 
def main():
    target_dir = r"c:\Users\Administrator\OneDrive\Thesis\TimeLapse_Experiment"
    
    # Base parameters matching those available in forward.create_input_file
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
        
        # Receivers positioning
        'rx_count': 20,
        'rx_spread': 2.0,  # Spread uniformly across 2 meters centered at source
        
        # Block size settings
        'block_size_fraction': 1/4, # 1/4 of central wavelength in ice
        
        # Time-lapse settings
        'time_lapse_shift': 0.05    # Shift in meters
    }

    experiment1 = GPRExperiment(target_dir, base_params)
    experiment1.run()
    
    # # Run another one with different block size quickly 
    # params2 = base_params.copy()
    # params2['block_size_fraction'] = 1/2
    # experiment2 = GPRExperiment(target_dir, params2)
    # experiment2.run()

if __name__ == "__main__":
    main()