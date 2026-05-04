import os
import numpy as np

def generate_bscan_array_in_file():
    in_file = "RTM/subwavelength_bscan_array.in"
    
    with open(in_file, "w") as f:
        f.write("#title: 2D Sub-wavelength scatterers - Moving Array B-Scan\n")
        f.write("#domain: 1.0 0.5 0.002\n")
        f.write("#dx_dy_dz: 0.002 0.002 0.002\n")
        
        f.write("#time_window: 8e-9\n")
        f.write("#material: 6 0 1 0 ice\n\n")
        
        f.write("#waveform: ricker 1 1.5e9 my_ricker\n")
        f.write("#box: 0 0 0 1.0 0.5 0.002 ice\n\n")
        
        # Scatterers
        f.write("#cylinder: 0.47 0.25 0 0.47 0.25 0.002 0.01 free_space\n")
        f.write("#cylinder: 0.53 0.25 0 0.53 0.25 0.002 0.01 free_space\n\n")
        
        # Moving array parameters
        num_steps = 40
        tx_start = 0.2
        tx_end = 0.8
        step_dx = (tx_end - tx_start) / (num_steps - 1)
        
        f.write(f"#hertzian_dipole: z {tx_start:.4f} 0.45 0.001 my_ricker\n")
        f.write(f"#src_steps: {step_dx:.6f} 0 0\n\n")
        
        # 11 Receivers moving with Tx (-10cm to +10cm offsets)
        offsets = np.linspace(-0.1, 0.1, 11)
        for rx_off in offsets:
            f.write(f"#rx: {tx_start + rx_off:.4f} 0.45 0.001\n")
        f.write(f"#rx_steps: {step_dx:.6f} 0 0\n")
            
    print(f"Generated {in_file}")
    print("Running GprMax (may take a minute or two)....")
    os.system(f"python -m gprMax {in_file} -n {num_steps}")

if __name__ == '__main__':
    generate_bscan_array_in_file()
