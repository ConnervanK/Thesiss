import numpy as np
import os
import sys
import matplotlib
matplotlib.use('Agg')
import time
from joblib import Parallel, delayed
import imageio.v2 as imageio

# ensure project root on path
# Output folders
output_folder = "figures/fwi_L2_frechet_simple"
os.makedirs(output_folder, exist_ok=True)
output_dir = os.path.normpath(output_folder)
os.makedirs(output_dir, exist_ok=True)

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from functions.assembly import A_assembly, Binv_sqrt_staggered, get_dof_index, reduced_space_time_system_Ghat
from functions.visualization.plot_FWI import (
    plot_fwi_results,
    plot_convergence,
    plot_misfit_decomposition,
    plot_velocity_model_comparison,
)
from functions.simulation_utils import ricker_wavelet
from functions.testing.forward_model_unitary import run_forward_simulation_unitary
from functions.gradients import (
    construct_b_vec,
    compute_gradient_frechet_L2,
)

def main():
    print("Initializing Simple FWI with Geometric Misfit tracking...")
    
    # --- Parameters ---
    nx, nz = 10, 10
    dx, dz = 100.0, 100.0
    nt = 80 
    dt = 0.02

    Np = nx * nz
    N = [nx, nz]
    d = [dx, dz]

    print('Building geometry...')
    A_geom = A_assembly(N, d, ['Neumann', 'Dirichlet', 'Dirichlet', 'Neumann']).tocsr()

    n_rho = (A_geom.shape[0] - Np) // 2
    rho_bg = 1.0
    rho_list = [np.ones(n_rho) * rho_bg for _ in range(2)]

    # --- True Model ---
    v_true = np.full((nx, nz), 1500.0)
    v_true[5:9, 5:9] = 1800.0  # The Inclusion
    kappa_true = (v_true.flatten() ** 2) * rho_bg

    Binv_sqrt_true = Binv_sqrt_staggered(kappa_true, rho_list, sqrt=True)
    b_diag_true = Binv_sqrt_true.diagonal()

    # --- Sources & Receivers ---
    x_srcs_map = np.zeros(A_geom.shape[0], dtype=bool)
    # src_coords = [(2, 5), (6, 9)]
    src_coords = [(2, 5), (6, 9), (1,1), (8,5)]
    src_idxs = [get_dof_index('p', c, N, d) for c in src_coords]
    x_srcs_map[src_idxs] = True

    x_rec_map = np.zeros(A_geom.shape[0], dtype=bool)
    # rec_coords = [(2, 10), (2, 11), (1,2), (3,10), (3,11), (4,10), (4,11), (10, 10), (10, 2)]
    rec_coords = [(2, 1), (2, 8), (1,4), (9,9), (3,9), (4,4), (9,4), (7, 2), (6, 5)]
    rec_idxs = [get_dof_index('p', c, N, d) for c in rec_coords]
    x_rec_map[rec_idxs] = True

    # --- Source Signal Setup ---
    t_axis = np.arange(nt) * dt
    
    # Pre-compute True Ghat (Depends only on True Model which is constant if grid/dt/material is constant)
    # Note: Reduced space time system Ghat depends on A, B, and dt.
    # It does NOT depend on the source frequency directly.
    print('Assembling True Ghat...')
    Ghat_true = reduced_space_time_system_Ghat(A_geom * dt, Binv_sqrt_true, nt, x_rec_map, x_srcs_map, use_real=True, verbose=False)
    
    # Visualization Setup of a single forward simulation (for sanity check and to capture frames for GIF)
    f = 5.0  # Frequency for the source wavelet in the forward simulation check
    run_forward_simulation_unitary(
        nx, nz, dx, dz, nt, dt,
        A_geom, Binv_sqrt_true, src_idxs, ricker_wavelet(t_axis, f, 1/f),
        output_path=os.path.join(output_dir, 'forward_simulation_check.gif'),
        show_progress=True
    )
    
    # exit()  # Remove this after confirming forward simulation works and captures frames correctly  

    # --- Inversion Setup ---
    # Frequency Schedule: List of (frequency_hz, num_iterations)
    schedule = [(1.5, 20), (2, 15), (2.25, 10)]
    
    # Initial Guess: Homogeneous
    v_background = np.full((nx, nz), 1500.0)
    v_current = np.copy(v_background)
    
    learning_rate = 2.5e1
    delta_v = 15.0 
    
    inv_x_start, inv_x_end = 1, 9
    inv_z_start, inv_z_end = 1, 9
    
    history = {
        'total_loss': [],
        'velocity_models': []
    }
    
    print(f"Starting Multi-Stage FWI with schedule: {schedule}")
    start_time = time.time()
    total_iter_count = 0
    all_frames = []
    misfit_history = []
    # Stage boundaries for convergence shading (cumulative iteration indices)
    stage_boundaries = [0]


    for freq_hz, n_iter_stage in schedule:
        print(f"\n>>> STAGE: {freq_hz} Hz, {n_iter_stage} iterations <<<")
        start_frame_idx = len(all_frames)
        
        # 1. Update Wavelet for this stage
        src_time = ricker_wavelet(t_axis, freq_hz, 1.0/freq_hz)
        
        # 2. Update True Data for this frequency
        # Note: We assume observed data is generated from the true model at the current frequency
        b_vec_true = construct_b_vec(b_diag_true, src_time, src_idxs, nt, dt)
        d_obs = (Ghat_true @ b_vec_true).real
        
        for iter_sub in range(n_iter_stage):
            iter_start = time.time()
            total_iter_count += 1
            
            # --- Current State ---
            kappa_curr = (v_current.flatten('F') ** 2) * rho_bg
            Binv_sqrt_curr = Binv_sqrt_staggered(kappa_curr, rho_list, sqrt=True)
            # Use the current frequency source wavelet
            b_vec_curr = construct_b_vec(Binv_sqrt_curr.diagonal(), src_time, src_idxs, nt, dt)
            
            Ghat_curr = reduced_space_time_system_Ghat(A_geom * dt, Binv_sqrt_curr, nt, x_rec_map, x_srcs_map, use_real=True, verbose=False)
            d_pred = (Ghat_curr @ b_vec_curr).real
            total_sq = np.linalg.norm(d_obs - d_pred) ** 2
            
            history['total_loss'].append(total_sq)
            history['velocity_models'].append(v_current.copy())
            misfit_history.append(total_sq)
            
            print(f"Freq {freq_hz}Hz | Iter {total_iter_count:02d} | Loss (L2): {total_sq:.4e}")
            
            
            # --- Compute Gradient (Parallel Finite Difference) ---
            grad_v = np.zeros_like(v_current)

            
            tasks = []
            for ix in range(inv_x_start, inv_x_end):
                for iz in range(inv_z_start, inv_z_end):
                    tasks.append((ix, iz))
            
            results = Parallel(n_jobs=-1)(
                delayed(compute_gradient_frechet_L2)(
                    ix, iz, v_current.copy(), rho_bg, rho_list, A_geom, dt, nt, 
                    x_rec_map, x_srcs_map, d_obs, delta_v, src_time, src_idxs, total_sq
                ) for ix, iz in tasks
            )
            
            for ix, iz, val in results:
                grad_v[ix, iz] = val
            
            # --- Update Model ---
            grad_norm = np.linalg.norm(grad_v)
            
            update_mag = 0.0
            if grad_norm > 1e-20:
                update = learning_rate * grad_v
                v_current = v_current - update
                update_mag = np.linalg.norm(update)
            
            v_current = np.clip(v_current, 1000, 3000)
            
            # --- Visualization ---
            # Using plot_fwi_results to generate frame
            # Note: grad_v is likely small at first iteration or zero if not updated yet, but we are inside loop
            # We treat grad_v as the gradient to visualize
            
            # We also need 'grad_2d' which plot_fwi_results expects. 
            # In LBFGS code it was gaussian filtered. Here grad_v is pixel-wise. 
            # We can use it directly or filter it for visualization.
            
            img_data = plot_fwi_results(
                v_true, v_current, v_background, grad_v, src_coords, rec_coords, nx, nz, dx, dz,
                total_iter_count, save_path=None, source_frequency=freq_hz
            )
            all_frames.append(img_data)
            
            print(f"    Grad norm: {grad_norm:.4e}. Update magnitude: {update_mag:.2f}. Time: {time.time() - iter_start:.2f}s")
        
        # End of Stage Logic
        stage_boundaries.append(len(misfit_history))
        
        # Save Stage GIF
        stage_frames = all_frames[start_frame_idx:]
        if stage_frames:
             try:
                 imageio.mimsave(os.path.join(output_dir, f'fwi_stage_{freq_hz}Hz.gif'), stage_frames, fps=5)
             except Exception as e:
                 print(f"Error saving stage gif: {e}")

    total_time = time.time() - start_time
    print(f"FWI completed in {total_time:.2f}s = {total_time/60:.2f} minutes")

    # --- Final Visualization ---
    print("Generating final figures...")
    
    if all_frames:
        imageio.mimsave(os.path.join(output_dir, 'geometric_fwi_evolution.gif'), all_frames, fps=5)
        
    # 1. Convergence (Using dedicated function)
    plot_convergence(misfit_history, schedule=schedule, stage_boundaries=stage_boundaries,
                     save_path=os.path.join(output_dir, 'fwi_convergence_history.png'), show=False)

    # 1. Detailed Convergence (Standardized shared format)
    plot_misfit_decomposition(history, save_path=os.path.join(output_dir, 'convergence.png'))

    # 2. Model Comparison (Standardized shared format)
    plot_velocity_model_comparison(
        v_true,
        history['velocity_models'][0],
        v_current,
        nx,
        nz,
        dx,
        dz,
        total_iter_count,
        save_path=os.path.join(output_dir, 'velocity_result.png'),
    )

    print(f"Results saved to {output_dir}")

if __name__ == "__main__":
    main()
