import numpy as np
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import imageio.v2 as imageio

# Ensure project root on path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from functions.assembly import A_assembly, Binv_sqrt_staggered, get_dof_index, reduced_space_time_system_Ghat
from functions.visualization.plot_FWI import plot_fwi_results, plot_convergence
from functions.simulation_utils import ricker_wavelet
from functions.gradients import construct_b_vec, compute_gradient_adjoint_L2

def main():
    # Use a different output directory for Gradient Descent
    output_dir = "figures/fwi_L2_adjoint_simple"
    os.makedirs(output_dir, exist_ok=True)

    nx, nz = 14, 14
    dx, dz = 100.0, 100.0
    nt, dt = 80, 0.02
    Np = nx * nz
    N, d = [nx, nz], [dx, dz]

    A_geom = A_assembly(N, d, ['Neumann', 'Dirichlet', 'Dirichlet', 'Neumann']).tocsr()
    rho_bg = 1.0
    n_rho = (A_geom.shape[0] - Np) // 2
    rho_list = [np.ones(n_rho) * rho_bg for _ in range(2)]

    v_true = np.full((nx, nz), 1500.0)
    v_true[4:7, 3:7] = 1800.0
    
    v_flat_true = v_true.flatten('F')
    kappa_true = (v_flat_true ** 2) * rho_bg
    Binv_sqrt_true = Binv_sqrt_staggered(kappa_true, rho_list, sqrt=True)

    src_coords = [(1, 1), (8, 8), (1, 8), (8, 1)] 
    src_idxs = sorted([get_dof_index('p', c, N, d) for c in src_coords])
    x_srcs_map = np.zeros(A_geom.shape[0], dtype=bool); x_srcs_map[src_idxs] = True
    
    rec_coords = []
    for i in range(0, 14):
        rec_coords.extend([(i, 0), (i, 13), (0, i), (13, i)])
    rec_coords = list(set([r for r in rec_coords if r not in src_coords]))

    rec_idxs = sorted([get_dof_index('p', c, N, d) for c in rec_coords])
    x_rec_map = np.zeros(A_geom.shape[0], dtype=bool); x_rec_map[rec_idxs] = True

    Ghat_true = reduced_space_time_system_Ghat(A_geom * dt, Binv_sqrt_true, nt, x_rec_map, x_srcs_map, use_real=True, verbose=False)

    v_background = np.full((nx, nz), 1500.0)
    v_current = np.copy(v_background)
    mask_inner = np.zeros_like(v_current, dtype=bool)
    mask_inner[1:13, 1:13] = True
    
    history = {'total_loss': [], 'velocity_models': []}
    all_frames = []
    t_axis = np.arange(nt) * dt
    
    # --- Gradient Descent Settings ---
    learning_rate = 5.0 # Adjustable step size
    
    schedule = [(1.0, 35), (1.25, 25), (1.5, 20)] 
    
    total_iter = 0

    for freq, n_iter in schedule:
        print(f"\n>>> STAGE: {freq} Hz, {n_iter} iterations (max) - Gradient Descent <<<")
        start_frame_idx = len(all_frames)

        src_time = ricker_wavelet(t_axis, freq, 1.0/freq)
        b_vec_true_stage = construct_b_vec(Binv_sqrt_true.diagonal(), src_time, src_idxs, nt, dt)
        d_obs_stage = (Ghat_true @ b_vec_true_stage).real
        
        for k in range(n_iter):
            iter_start_time = time.time()
            total_iter += 1

            # 1. Update Material Properties for Current Model
            v_flat_curr = v_current.flatten('F')
            kappa_curr = (v_flat_curr ** 2) * rho_bg
            Binv_sqrt_curr = Binv_sqrt_staggered(kappa_curr, rho_list, sqrt=True)
            b_vec_curr = construct_b_vec(Binv_sqrt_curr.diagonal(), src_time, src_idxs, nt, dt)

            # 2. Compute Gradient (using shared function)
            loss, grad_full, _ = compute_gradient_adjoint_L2(
                v_current, rho_bg, A_geom, b_vec_curr, d_obs_stage, mask_inner, 
                src_coords, rec_coords, dx, dz, dt, nt, apply_smoothing=True
            )
            
            # 3. Descent Step
            grad_update = grad_full[mask_inner]
            
            # Normalize gradient to ensure physical step size
            # This makes 'learning_rate' equivalent to the maximum velocity change per iteration (m/s)
            grad_max = np.max(np.abs(grad_update))
            if grad_max > 0.0:
                 grad_update /= grad_max
                 
            v_current[mask_inner] -= learning_rate * grad_update
            
            # 4. Enforce Bounds
            v_current[mask_inner] = np.clip(v_current[mask_inner], 1200, 2500)
            
            # Logging
            iter_duration = time.time() - iter_start_time
            history['total_loss'].append(loss)
            history['velocity_models'].append(v_current.copy())
            
            print(f"Iter {total_iter} | Loss: {loss:.4e} | V_mean: {np.mean(v_current):.1f} | Time: {iter_duration:.2f}s")
            
            img_data = plot_fwi_results(
                v_true, v_current, v_background, grad_full, src_coords, rec_coords, nx, nz, dx, dz,
                total_iter, save_path=None, source_frequency=freq
            )
            all_frames.append(img_data)

        stage_frames = all_frames[start_frame_idx:]
        if stage_frames:
            imageio.mimsave(os.path.join(output_dir, f'gd_fwi_stage_{freq}Hz.gif'), stage_frames, fps=5)

    if all_frames:
        imageio.mimsave(os.path.join(output_dir, 'gd_fwi_evolution.gif'), all_frames, fps=5)
    
    plot_convergence(history['total_loss'], save_path=os.path.join(output_dir, 'convergence.png'))

if __name__ == "__main__":
    main()
