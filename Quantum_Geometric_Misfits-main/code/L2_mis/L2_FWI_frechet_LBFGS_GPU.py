import numpy as np
import os
import sys
import matplotlib
matplotlib.use('Agg')
from scipy.optimize import minimize
import time
import imageio.v2 as imageio
import cupy as cp
import cupyx.scipy.sparse as cps
HAS_GPU = True
print("CuPy imported successfully. GPU acceleration enabled.")

# ensure project root on path
# Output folders
output_folder = "figures/fwi_L2_frechet_LBFGS_GPU"
os.makedirs(output_folder, exist_ok=True)
output_dir = os.path.normpath(output_folder)
os.makedirs(output_dir, exist_ok=True)

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from functions.assembly import A_assembly, Binv_sqrt_staggered, get_dof_index, reduced_space_time_system_Ghat_gpu as reduced_space_time_system_Ghat
from functions.visualization.plot_FWI import plot_fwi_results, plot_convergence, plot_misfit_decomposition, plot_final_velocity_comparison
from functions.simulation_utils import ricker_wavelet
from functions.testing.forward_model_unitary import run_forward_simulation_unitary
from functions.testing.FWI_utils import gpu_simulation_loss_l2 as gpu_simulation_loss, give_misfits
from scipy.optimize import minimize

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
    unique_src_indices = np.where(x_srcs_map)[0]
    num_srcs = len(unique_src_indices)
    
    def construct_b_vec(b_diag_local, wavelet_array):
        b_vec_out = np.zeros(nt * num_srcs)
        for it in range(nt):
            for i, s_idx in enumerate(unique_src_indices):
                b_vec_out[it * num_srcs + i] = wavelet_array[it] * b_diag_local[s_idx] * dt
        return b_vec_out
    
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
    schedule = [(1.5, 20), (2, 15), (2.25, 10), (3, 5)]
    
    # User controls
    # loss_change_tolerance = 2e-3 # Stop if change in loss is below this value
    loss_change_tolerance = None # Disable automatic stopping based on loss change for now

    # Initial Guess: Homogeneous
    v_background = np.full((nx, nz), 1500.0)
    v_current = np.copy(v_background)
    
    delta_v = 15.0 
    
    inv_x_start, inv_x_end = 1, 9
    inv_z_start, inv_z_end = 1, 9
    
    history = {
        'total_loss': [],
        'geom_misfit': [],
        'param_misfit': [],
        'velocity_models': []
    }
    
    print(f"Starting Multi-Stage FWI with schedule: {schedule}")
    start_time = time.time()
    total_iter_count = 0
    all_frames = []
    misfit_history = []
    # Stage boundaries for convergence shading (cumulative iteration indices)
    stage_boundaries = [0]
    
    # Objective function wrapper for minimize
    def fwi_objective_function(v_flat_inner, v_full_shape, mask_inner, freq_hz, n_iter_stage_max):
        nonlocal total_iter_count # Access outer scope counter
        
        iter_start = time.time()
        
        # Reconstruct full velocity grid
        v_temp = v_full_shape.copy()
        v_temp[mask_inner] = v_flat_inner
        
        # --- Current State ---
        kappa_curr = (v_temp.flatten() ** 2) * rho_bg
        Binv_sqrt_curr = Binv_sqrt_staggered(kappa_curr, rho_list, sqrt=True)
        # Use the current frequency source wavelet
        b_vec_curr = construct_b_vec(Binv_sqrt_curr.diagonal(), src_time)
        
        Ghat_curr = reduced_space_time_system_Ghat(A_geom * dt, Binv_sqrt_curr, nt, x_rec_map, x_srcs_map, use_real=True, verbose=False)
        
        # Misfit Calculation
        comp_A, comp_B, total_sq = give_misfits(Ghat_curr, d_obs, b_vec_curr)
        
        # Tracking
        # Note: minimize() might call this multiple times per iteration (line search)
        # We only want to log "major" iterations or every call?
        # Usually logging every call is fine for debugging, but total_iter_count increases monotonic.
        
        # --- Compute Gradient ---
        grad_v = np.zeros_like(v_temp)
        
        # Transfer static data to GPU
        A_geom_gpu = cps.csr_matrix(A_geom)
        d_obs_gpu = cp.array(d_obs)
        src_time_gpu = cp.array(src_time)
        x_rec_idx_gpu = cp.array(np.where(x_rec_map)[0])
        x_src_idx_gpu = cp.array(np.where(x_srcs_map)[0])

        # Current Loss (GPU) for consistency
        current_loss_gpu = gpu_simulation_loss(
                v_temp, rho_bg, rho_list, A_geom_gpu, dt, nt, 
                x_rec_idx_gpu, x_src_idx_gpu, d_obs_gpu, None, src_time_gpu, unique_src_indices
        )
        # print(f"    (GPU Loss: {current_loss_gpu:.4e} vs CPU: {total_sq:.4e})")

        for ix in range(inv_x_start, inv_x_end):
            for iz in range(inv_z_start, inv_z_end):
                v_temp[ix, iz] += delta_v
                loss_p = gpu_simulation_loss(
                        v_temp, rho_bg, rho_list, A_geom_gpu, dt, nt, 
                        x_rec_idx_gpu, x_src_idx_gpu, d_obs_gpu, None, src_time_gpu, unique_src_indices
                )
                grad_v[ix, iz] = (loss_p - current_loss_gpu) / delta_v
                v_temp[ix, iz] -= delta_v
            
        # Extract gradient for the inner region only
        grad_flat_inner = grad_v[mask_inner]
        
        # --- Visualization & Logging ---
        # Only log/visualize if this is a "valid" step (misfit decreased or accepted)
        # But L-BFGS calls this function. We can just log everything.
        
        history['total_loss'].append(total_sq)
        history['geom_misfit'].append(comp_A)
        history['param_misfit'].append(comp_B)
        history['velocity_models'].append(v_temp.copy())
        misfit_history.append(total_sq)
        
        total_iter_count += 1
        iter_time = time.time() - iter_start
        print(f"Freq {freq_hz}Hz | Call {total_iter_count} | Loss: {total_sq:.4e} | Geom: {comp_A:.4e} | Param: {comp_B:.4e} | Time: {iter_time:.2f}s")
            
        # Check tolerance (Change in loss)
        if loss_change_tolerance is not None and len(history['geom_misfit']) > 1:
            loss_change = abs(history['total_loss'][-2] - total_sq)
            if loss_change < loss_change_tolerance:
                print(f"Converged: Loss change {loss_change:.4e} < tolerance {loss_change_tolerance}")
                raise StopIteration

        img_data = plot_fwi_results(
            v_true, v_temp, v_background, grad_v, src_coords, rec_coords, nx, nz, dx, dz,
            total_iter_count, save_path=None, source_frequency=freq_hz
        )
        all_frames.append(img_data)
        
        return total_sq, grad_flat_inner


    for freq_hz, n_iter_stage in schedule:
        print(f"\n>>> STAGE: {freq_hz} Hz, {n_iter_stage} iterations (max) <<<")
        start_frame_idx = len(all_frames)
        
        # 1. Update Wavelet for this stage
        src_time = ricker_wavelet(t_axis, freq_hz, 1.0/freq_hz)
        
        # 2. Update True Data for this frequency
        b_vec_true = construct_b_vec(b_diag_true, src_time)
        d_obs = (Ghat_true @ b_vec_true).real
        
        # 3. Setup Mask for Optimization Region (Inner pixels)
        mask_inner = np.zeros_like(v_current, dtype=bool)
        mask_inner[inv_x_start:inv_x_end, inv_z_start:inv_z_end] = True
        
        v_initial_flat = v_current[mask_inner]
        
        # 4. Run L-BFGS-B
        # We need to wrap the objective to match scipy signature: func(x, *args) -> (f, g)
        
        # Bounds: 1000 to 3000
        bounds = [(1000.0, 3000.0) for _ in range(len(v_initial_flat))]
        
        try:
            res = minimize(
                fun=fwi_objective_function,
                x0=v_initial_flat,
                args=(v_current, mask_inner, freq_hz, n_iter_stage),
                method='L-BFGS-B',
                jac=True, # We return gradient
                bounds=bounds,
                options={'maxiter': n_iter_stage, 'disp': True}
            )
            
            # 5. Update v_current with result
            v_current[mask_inner] = res.x
        except StopIteration:
            print("Optimization stopped early due to tolerance.")
            if history['velocity_models']:
                 v_current[:] = history['velocity_models'][-1][:]
        
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

    # 1. Detailed Convergence (Decomposition)
    plot_misfit_decomposition(history, schedule=schedule, stage_boundaries=stage_boundaries,
                              save_path=os.path.join(output_dir, 'convergence.png'))
    
    # 2. Model Comparison
    plot_final_velocity_comparison(v_true, history['velocity_models'][0], v_current, nx, nz, dx, dz,
                                   save_path=os.path.join(output_dir, 'velocity_result.png'))

    print(f"Results saved to {output_dir}")

if __name__ == "__main__":
    main()
