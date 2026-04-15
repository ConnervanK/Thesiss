import numpy as np
import os
import sys
import matplotlib
matplotlib.use('Agg')

from scipy.optimize import minimize
import time
import imageio.v2 as imageio

# Ensure project root on path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from functions.assembly import A_assembly, Binv_sqrt_staggered, get_dof_index, reduced_space_time_system_Ghat_no_block as reduced_space_time_system_Ghat
from functions.visualization.plot_FWI import plot_fwi_results, plot_convergence
from functions.simulation_utils import ricker_wavelet
from functions.gradients import compute_gradient_adjoint_L2, compute_gradient_adjoint_geometric
from functions.geometric_mis import give_misfits

class FWIInversion:
    def __init__(self, nx, nz, dx, dz, nt, dt):
        self.nx = nx
        self.nz = nz
        self.dx = dx
        self.dz = dz
        self.nt = nt
        self.dt = dt
        self.N = [nx, nz]
        self.d = [dx, dz]
        self.Np = nx * nz
        
        self.A_geom = None
        self.src_coords = None
        self.rec_coords = None
        
        self.v_true = None
        self.v_background = None
        self.rho_bg = None
        self.mask_inner = None
        
        self.x_srcs_map = None
        self.x_rec_map = None
        self.src_idxs = None
        self.rec_idxs = None
        
        self.history = {'total_loss': [], 'geometric_misfit': [], 'param_misfit': [], 'velocity_models': []}
        self.all_frames = []
        self.total_iter_count = 0
        
    def setup_geometry(self, src_coords, rec_coords, boundaries=['Neumann', 'Dirichlet', 'Dirichlet', 'Neumann']):
        print('Building geometry...')
        self.A_geom = A_assembly(self.N, self.d, boundaries).tocsr()
        self.src_coords = src_coords
        self.rec_coords = rec_coords
        
        self.src_idxs = sorted([get_dof_index('p', c, self.N, self.d) for c in src_coords])
        self.x_srcs_map = np.zeros(self.A_geom.shape[0], dtype=bool)
        self.x_srcs_map[self.src_idxs] = True
        
        self.rec_idxs = sorted([get_dof_index('p', c, self.N, self.d) for c in rec_coords])
        self.x_rec_map = np.zeros(self.A_geom.shape[0], dtype=bool)
        self.x_rec_map[self.rec_idxs] = True
        
    def setup_models(self, eps_true, eps_background, mu_bg, eps_0, mask_inner=None):
        self.eps_true = eps_true
        self.eps_background = eps_background
        self.mu_bg = mu_bg
        self.eps_0 = eps_0 # Vacuum permittivity
        
        if mask_inner is None:
            self.mask_inner = np.zeros_like(self.eps_background, dtype=bool)
            self.mask_inner[1:-1, 1:-1] = True # Default skip 1 point boundary
        else:
            self.mask_inner = mask_inner
            
        n_mu = (self.A_geom.shape[0] - self.Np) // 2
        # Magnetic permeability components across the staggered grid
        self.mu_list = [np.ones(n_mu) * self.mu_bg for _ in range(2)]
        
    def generate_true_data(self, freq_hz):
        eps_flat_true = self.eps_true.flatten('F')
        # Here "kappa" corresponds to 1 / (epsilon_r * epsilon_0)
        inv_eps_true = 1.0 / (eps_flat_true * self.eps_0)
        Binv_sqrt_true = Binv_sqrt_staggered(inv_eps_true, self.mu_list, sqrt=True)
        
        t_axis = np.arange(self.nt) * self.dt
        wavelet = ricker_wavelet(t_axis, freq_hz, 1.5/freq_hz)
        
        unique_src_indices = np.where(self.x_srcs_map)[0]
        num_srcs = len(unique_src_indices)
        
        # Build source vector
        b_vec = np.zeros(self.nt * num_srcs)
        b_diag_local = Binv_sqrt_true.diagonal()
        for it in range(self.nt):
            for i, s_idx in enumerate(unique_src_indices):
                b_vec[it * num_srcs + i] = wavelet[it] * b_diag_local[s_idx] * self.dt
                
        Ghat_true = reduced_space_time_system_Ghat(
            self.A_geom * self.dt, Binv_sqrt_true, self.nt, 
            self.x_rec_map, self.x_srcs_map, use_real=True, verbose=False, auto_substeps=False, rk4_substeps=25
        )
        
        d_obs = (Ghat_true @ b_vec).real
        return d_obs, wavelet, b_vec
        
    def run_inversion(self, freq_schedule=[15.0, 25.0], max_iters=[10, 10], misfit='L2', optimizer='LBFGS', output_folder='figures/fwi_general', gd_learning_rate=500.0, loss_change_tolerance=1e-6):
        os.makedirs(output_folder, exist_ok=True)
        eps_current = np.copy(self.eps_background)
        self.all_frames = []
        self.history = {'total_loss': [], 'geometric_misfit': [], 'param_misfit': [], 'eps_models': []}
        self.total_iter_count = 0
        
        for stage_idx, (freq_hz, max_iter) in enumerate(zip(freq_schedule, max_iters)):
            print(f"--- Starting Stage {stage_idx + 1}/{len(freq_schedule)}: Frequency = {freq_hz} Hz ---")
            
            d_obs, wavelet, b_vec_true = self.generate_true_data(freq_hz)
            
            last_iter_time = time.time()
            
            def objective(eps_flat_inner):
                nonlocal last_iter_time
                current_time = time.time()
                iter_duration = current_time - last_iter_time
                last_iter_time = current_time
                
                eps_temp = np.copy(self.eps_background)
                eps_temp[self.mask_inner] = eps_flat_inner
                
                eps_flat_temp = eps_temp.flatten('F')
                # EM substitution for inverse permittivity
                inv_eps_curr = 1.0 / (eps_flat_temp * self.eps_0)
                Binv_sqrt_curr = Binv_sqrt_staggered(inv_eps_curr, self.mu_list, sqrt=True)
                
                # Source vector for current model
                unique_src_indices = np.where(self.x_srcs_map)[0]
                num_srcs = len(unique_src_indices)
                b_vec_curr = np.zeros(self.nt * num_srcs)
                b_diag_local = Binv_sqrt_curr.diagonal()
                for it in range(self.nt):
                    for i, s_idx in enumerate(unique_src_indices):
                        b_vec_curr[it * num_srcs + i] = wavelet[it] * b_diag_local[s_idx] * self.dt
                
                # Setup objective and gradient variables
                if misfit == 'L2':
                    loss, grad_full, _ = compute_gradient_adjoint_L2(
                        eps_temp, self.mu_bg, self.A_geom, b_vec_curr, d_obs, self.mask_inner, 
                        self.src_coords, self.rec_coords, self.dx, self.dz, self.dt, self.nt, apply_smoothing=True
                    )
                    
                    comp_A, comp_B, total_sq = loss, 0.0, loss
                    
                elif misfit == 'geometric':
                    Ghat_curr = reduced_space_time_system_Ghat(
                        self.A_geom * self.dt, Binv_sqrt_curr, self.nt, 
                        self.x_rec_map, self.x_srcs_map, use_real=True, verbose=False, auto_substeps=False, rk4_substeps=25
                    )
                    comp_A, comp_B, total_sq = give_misfits(Ghat_curr, d_obs, b_vec_curr)
                    loss = comp_A
                    
                    grad_full = compute_gradient_adjoint_geometric(
                        eps_temp, self.mu_bg, self.mu_list, self.A_geom, self.dt, self.nt, 
                        self.x_rec_map, self.x_srcs_map, d_obs
                    )
                    
                    # Applying taper / smoothing explicitly as it's missing in pure geometric standard call
                    from functions.gradients import apply_source_taper
                    grad_full = apply_source_taper(grad_full, self.src_coords, self.nx, self.nz, radius=3.0)
                    from scipy.ndimage import gaussian_filter
                    grad_full = gaussian_filter(grad_full, sigma=0.7)
                    
                else:
                    raise ValueError("misfit parameter must be 'L2' or 'geometric'")
                    
                grad_flat_inner = grad_full[self.mask_inner]
                
                # Logging
                self.history['total_loss'].append(total_sq)
                self.history['geometric_misfit'].append(comp_A)
                self.history['param_misfit'].append(comp_B)
                self.history['eps_models'].append(eps_temp.copy())
                self.total_iter_count += 1
                
                print(f"Freq {freq_hz}Hz | Call {self.total_iter_count} | Misfit Type: {misfit} | Loss: {loss:.4e} | |grad|: {np.linalg.norm(grad_flat_inner):.2e} | Time: {iter_duration:.2f}s")
                
                # Visuals
                img_data = plot_fwi_results(
                    self.eps_true, eps_temp, self.eps_background, grad_full, self.src_coords, self.rec_coords, 
                    self.nx, self.nz, self.dx, self.dz, self.total_iter_count, save_path=None, source_frequency=freq_hz
                )
                self.all_frames.append(img_data)
                
                return loss, grad_flat_inner

            # Inversion Iteration logic
            eps_flat_inner = eps_current[self.mask_inner]
            
            if optimizer == 'LBFGS':
                res = minimize(
                    objective,
                    eps_flat_inner,
                    method='L-BFGS-B',
                    jac=True,
                    options={'maxiter': max_iter, 'disp': True}
                )
                eps_current[self.mask_inner] = res.x
            
            elif optimizer == 'GD':
                for i in range(max_iter):
                    loss, grad = objective(eps_flat_inner)
                    eps_flat_inner = eps_flat_inner - gd_learning_rate * grad
                    eps_current[self.mask_inner] = eps_flat_inner
            else:
                raise ValueError("Optimizer should be 'LBFGS' or 'GD'")
                
        # Finalization
        # Plot converence
        if len(self.history['eps_models']) > 0:
            plot_convergence(self.history['geometric_misfit'], output_folder)
            
        gif_path = os.path.join(output_folder, f"fwi_inversion_{misfit}_{optimizer}.gif")
        imageio.mimsave(gif_path, self.all_frames, fps=2)
        print(f"Animation saved to {gif_path}")
        return eps_current
