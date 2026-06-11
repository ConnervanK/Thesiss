import numpy as np
import os
import sys
import matplotlib
matplotlib.use('Agg')

from scipy.optimize import minimize
import time
import imageio.v2 as imageio

# Ensure project root on path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from functions.assembly import A_assembly, Binv_sqrt_staggered, get_dof_index, reduced_space_time_system_Ghat_no_block as reduced_space_time_system_Ghat
from functions.visualization.plot_FWI import plot_fwi_results, plot_convergence
from functions.simulation_utils import ricker_wavelet
from functions.gradients import construct_b_vec, compute_gradient_adjoint_L2

def main():
    output_dir = "figures/fwi_L2_adjoint_LBFGS"
    os.makedirs(output_dir, exist_ok=True)

    nx, nz = 10, 10
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

    src_coords = [(1, 1), (8, 8)] 
    src_idxs = sorted([get_dof_index('p', c, N, d) for c in src_coords])
    x_srcs_map = np.zeros(A_geom.shape[0], dtype=bool); x_srcs_map[src_idxs] = True
    
    rec_coords = []
    for i in range(0, 10):
        rec_coords.extend([(i, 0), (i, 9), (0, i), (9, i)])
    rec_coords = list(set([r for r in rec_coords if r not in src_coords]))

    rec_idxs = sorted([get_dof_index('p', c, N, d) for c in rec_coords])
    x_rec_map = np.zeros(A_geom.shape[0], dtype=bool); x_rec_map[rec_idxs] = True

    Ghat_true = reduced_space_time_system_Ghat(A_geom * dt, Binv_sqrt_true, nt, x_rec_map, x_srcs_map, use_real=True, verbose=False)

    v_background = np.full((nx, nz), 1500.0)
    v_current = np.copy(v_background)
    mask_inner = np.zeros_like(v_current, dtype=bool)
    mask_inner[1:9, 1:9] = True
    
    history = {'total_loss': [], 'velocity_models': []}
    all_frames = []
    total_iter_count = [0]
    t_axis = np.arange(nt) * dt
    last_iter_time = [time.time()]

    def fwi_objective(v_flat_inner, v_full_shape, freq_hz, src_time, d_obs_stage):
        current_time = time.time()
        iter_duration = current_time - last_iter_time[0]
        last_iter_time[0] = current_time
        
        v_temp = v_full_shape.copy()
        v_temp[mask_inner] = v_flat_inner
        
        v_flat_temp = v_temp.flatten('F')
        kappa_curr = (v_flat_temp ** 2) * rho_bg
        Binv_sqrt_curr = Binv_sqrt_staggered(kappa_curr, rho_list, sqrt=True)
        b_vec_curr = construct_b_vec(Binv_sqrt_curr.diagonal(), src_time, src_idxs, nt, dt)
        
        loss, grad_full, _ = compute_gradient_adjoint_L2(
            v_temp, rho_bg, A_geom, b_vec_curr, d_obs_stage, mask_inner, 
            src_coords, rec_coords, dx, dz, dt, nt, apply_smoothing=True
        )
        
        grad_flat = grad_full[mask_inner]
        
        total_iter_count[0] += 1
        history['total_loss'].append(loss)
        history['velocity_models'].append(v_temp.copy())
        
        # 'total_iter_count' counts function evaluations, which can be > maxiter due to line search
        print(f"FuncCall {total_iter_count[0]} | Loss: {loss:.4e} | V_mean: {np.mean(v_temp):.1f} | Time: {iter_duration:.2f}s")
        
        img_data = plot_fwi_results(
            v_true, v_temp, v_background, grad_full, src_coords, rec_coords, nx, nz, dx, dz,
            total_iter_count[0], save_path=None, source_frequency=freq_hz
        )
        all_frames.append(img_data)
        
        return loss, grad_flat
        

    schedule = [(1.0, 25), (1.25, 25), (1.5, 20), (1.75, 20), (2.0, 15), (2.5, 15), (3.0, 10)] 
    
    for freq, n_iter in schedule:
        print(f"\n>>> STAGE: {freq} Hz, {n_iter} iterations (max) <<<")
        start_frame_idx = len(all_frames)

        src_time = ricker_wavelet(t_axis, freq, 1.0/freq)
        b_vec_true_stage = construct_b_vec(Binv_sqrt_true.diagonal(), src_time, src_idxs, nt, dt)
        d_obs_stage = (Ghat_true @ b_vec_true_stage).real
        
        res = minimize(
            fun=fwi_objective,
            x0=v_current[mask_inner],
            args=(v_current, freq, src_time, d_obs_stage),
            method='L-BFGS-B',
            jac=True,
            # callback=optimization_callback,
            bounds=[(1200, 2500)] * np.sum(mask_inner),
            options={'maxiter': n_iter, 'gtol': 1e-6}
        )
        v_current[mask_inner] = res.x

        stage_frames = all_frames[start_frame_idx:]
        if stage_frames:
            imageio.mimsave(os.path.join(output_dir, f'adjoint_fwi_stage_{freq}Hz.gif'), stage_frames, fps=5)

    if all_frames:
        imageio.mimsave(os.path.join(output_dir, 'adjoint_fwi_evolution.gif'), all_frames, fps=5)
    
    plot_convergence(history['total_loss'], save_path=os.path.join(output_dir, 'convergence.png'))

if __name__ == "__main__":
    main()