import os
import numpy as np

def generate_in_file():
    in_file = "RTM/subwavelength_scatters.in"
    
    # Typical GprMax setup for 1.5 GHz in ice
    # Subwavelength at 1.5 GHz in ice (er=6) is ~ 8.16 cm. We use 1 cm radius spheres
    with open(in_file, "w") as f:
        f.write("#title: 2D Sub-wavelength scatterers in ice\n")
        f.write("#domain: 1.0 0.5 0.002\n")
        f.write("#dx_dy_dz: 0.002 0.002 0.002\n")
        
        # Enough time for two-way travel: depth ~0.25m -> d_total ~0.5m. v=1.22e8 m/s -> t ~ 4.1 ns. Use 8 ns
        f.write("#time_window: 8e-9\n")
        f.write("#material: 6 0 1 0 ice\n\n")
        
        f.write("#waveform: ricker 1 1.5e9 my_ricker\n")
        f.write("#box: 0 0 0 1.0 0.5 0.002 ice\n\n")
        
        # Two air cavities (cylinders across the z slice to represent 2D circles)
        # 6cm spacing between centers, each 1cm radius
        f.write("#cylinder: 0.41 0.25 0 0.41 0.25 0.002 0.01 free_space\n")
        f.write("#cylinder: 0.47 0.25 0 0.47 0.25 0.002 0.01 free_space\n")
        f.write("#cylinder: 0.53 0.25 0 0.53 0.25 0.002 0.01 free_space\n\n")
        f.write("#cylinder: 0.59 0.25 0 0.59 0.25 0.002 0.01 free_space\n")
        
        # 20 Transmitters stepping across
        tx_start = 0.1
        rx_x = np.linspace(0.1, 0.9, 20)
        step_dx = rx_x[1] - rx_x[0]
        f.write(f"#hertzian_dipole: z {tx_start:.4f} 0.45 0.001 my_ricker\n")
        f.write(f"#src_steps: {step_dx:.6f} 0 0\n\n")
        
        # 20 Receivers spread evenly across top surface
        f.write("#rx_steps: 0 0 0\n") # Single shot SIMO array (rx doesn't move during steps)
        # Using a loop to place multiple static receivers
        for i, x in enumerate(rx_x):
            f.write(f"#rx: {x:.4f} 0.45 0.001\n")
            
    print(f"Generated {in_file}")
    print("Running GprMax (may take a minute or two)....")
    
    # Run the model for 20 steps (FMC)
    os.system(f"python -m gprMax {in_file} -n 20")

if __name__ == '__main__':
    generate_in_file()