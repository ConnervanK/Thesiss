import numpy as np
import os
import sys
import matplotlib
matplotlib.use('Agg')


# Ensure project root on path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from functions.FWI_functions import FWIInversion

def main():
    print("--- Running Generalized FWI ---")
    nx, nz = 30, 30
    dx, dz = 0.2, 0.2      # 20 cm grid spacing
    nt = 300               # time steps
    dt = 2e-10             # 0.2 nanoseconds
    
    fwi = FWIInversion(nx, nz, dx, dz, nt, dt)
    
    # Receivers and sources
    # src_coords = [(1, 1), (8, 8)]
    src_coords = [(1, 1), (8, 8), (1, 8), (8, 1)]
    rec_coords = []

    for i in range(0, 10):
        rec_coords.extend([(i, 0), (i, 28), (0, i), (28, i)])
    rec_coords = list(set([r for r in rec_coords if r not in src_coords]))
    
    fwi.setup_geometry(src_coords, rec_coords)
    
    # Permittivity models for GPR (Electromagnetic formulation)
    # Background: Air/Vacuum or Dry sand (e.g., eps_r = 1.0 or 4.0)
    eps_true = np.full((nx, nz), 4.0)
    eps_true[4:17, 21:27] = 9.0  # Anomaly (e.g. wet sand or target)
    
    eps_background = np.full((nx, nz), 4.0)
    mu_bg = 1.25663706e-6  # Vacuum permeability in H/m
    eps_0 = 8.854187817e-12 # Vacuum permittivity in F/m
    
    # mask_inner = np.zeros((nx, nz), dtype=bool)
    # mask_inner[1:9, 1:9] = True
    mask_inner = None  

    fwi.setup_models(eps_true, eps_background, mu_bg, eps_0, mask_inner = mask_inner)
    
    # User Configuration
    misfit_choice = 'geometric' # Choose: 'L2' or 'geometric'
    optimizer_choice = 'LBFGS'  # Choose: 'LBFGS' or 'GD'
    
    print(f"Running Inversion with Misfit: {misfit_choice}, Optimizer: {optimizer_choice}")
    fwi.run_inversion(
        freq_schedule=[5e7, 1e8, 1.5e8, 2e8, 2.5e8],  # 50 MHz to 250 MHz
        max_iters=[25, 25, 20, 20, 15], 
        misfit=misfit_choice, 
        optimizer=optimizer_choice,
        output_folder=f"figures/fwi_{misfit_choice}_{optimizer_choice}"
    )

if __name__ == "__main__":
    main()
