'''
Here you find all the needed functions to run the full waveform inversion test. In partucular the gradient calculation
'''

import os
import sys

import numpy as np
import scipy.sparse as sp
import cupy as cp
import cupyx.scipy.sparse as cps
import cupyx.scipy.linalg as cpl
HAS_GPU = True

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from assembly import Binv_sqrt_staggered

def give_misfits(Ghat_guess, x0hat_vec, b): 
    """
    Computes the geometric and parametric misfit components.
    """       
    if sp.issparse(Ghat_guess):
        G_dense = Ghat_guess.toarray()
    else:
        G_dense = np.asarray(Ghat_guess)

    # R = G * pinv(G) is the projection matrix onto the 'guess' subspace
    # Note: Moore-Penrose pseudo-inverse
    R_guess = (G_dense @ np.linalg.pinv(G_dense)).real

    # 1. Component A: Irreducible/Geometric Misfit vector
    comp_A = (np.eye(R_guess.shape[0]) - R_guess) @ x0hat_vec
    comp_A_norm_sq = np.linalg.norm(comp_A)**2


    # 2. Component B: Parametric/Source Misfit vector
    comp_B = R_guess @ x0hat_vec - Ghat_guess @ b
    comp_B_norm_sq = np.linalg.norm(comp_B)**2

    # 3. Total Error in the unitary space
    total_error_vec = x0hat_vec - (Ghat_guess @ b)
    total_error_norm_sq = np.linalg.norm(total_error_vec)**2
    
    return comp_A_norm_sq, comp_B_norm_sq, total_error_norm_sq

def _compute_G_dense_gpu(v_test, rho_bg, rho_list, A_geom_gpu, dt, nt, x_rec_idx_gpu, x_src_idx_gpu):
    """
    Helper to compute the dense block Toeplitz matrix G on GPU.
    """
    # 1. Setup Material Properties (B matrix)
    kappa_test = (v_test.flatten() ** 2) * rho_bg
    B_cpu = Binv_sqrt_staggered(kappa_test, rho_list, sqrt=True)
    B_diag = B_cpu.diagonal() # numpy array
    
    # Move B to GPU
    B_gpu_diag = cp.array(B_diag)
    B_gpu = cps.diags(B_gpu_diag) # Sparse diagonal on GPU
    
    # 2. Construct Operator H = B @ A @ B * dt
    H_gpu = (B_gpu @ A_geom_gpu @ B_gpu) * dt
    
    # 3. Matrix Exponential on GPU
    H_dense = H_gpu.toarray()
    U_gpu = cpl.expm(H_dense)
    
    # 4. Simulation / Convolution
    N_sys = U_gpu.shape[0]
    dtype_U = U_gpu.dtype
    
    # Preallocate U_all [U^0, ..., U^nt-1]
    U_all = cp.zeros((nt, N_sys, N_sys), dtype=dtype_U)
    U_all[0] = cp.eye(N_sys, dtype=dtype_U)
    
    current_idx = 1
    P = U_gpu 
    
    while current_idx < nt:
        count = min(current_idx, nt - current_idx)
        # Vectorized multiplication: U[idx:idx+count] = P @ U[0:count]
        chunk = U_all[:count]
        U_all[current_idx : current_idx + count] = cp.matmul(P, chunk)
        
        current_idx += count
        if current_idx < nt:
            P = P @ P

    # Pre-extract reduced operators U^k[rec, src]
    temp = U_all[:, x_rec_idx_gpu, :]
    reduced_ops = temp[:, :, x_src_idx_gpu]
            
    # Construct full Block-Toeplitz Operator G_semidense
    n_rec = len(x_rec_idx_gpu)
    n_src = len(x_src_idx_gpu)
    
    G_dense = cp.zeros((nt * n_rec, nt * n_src), dtype=dtype_U)
    
    for i in range(nt):
        for j in range(i, nt):
             lag = j - i
             G_dense[j*n_rec:(j+1)*n_rec, i*n_src:(i+1)*n_src] = reduced_ops[lag]
             
    return G_dense, B_gpu_diag

def gpu_simulation_loss_l2(v_test, rho_bg, rho_list, A_geom_gpu, dt, nt, x_rec_idx_gpu, x_src_idx_gpu, d_obs_gpu, b_vec_skeleton_gpu, src_time_gpu, unique_src_indices):
    """
    Computes L2 loss on GPU using vectorized operations.
    Minimizes || G * b - d_obs ||^2
    """
    G_dense, B_gpu_diag = _compute_G_dense_gpu(v_test, rho_bg, rho_list, A_geom_gpu, dt, nt, x_rec_idx_gpu, x_src_idx_gpu)
    
    # Construct full source vector b
    B_src_gpu = B_gpu_diag[unique_src_indices] 
    scale_factors = B_src_gpu * dt 
    
    # Outer product: time_profile (nt,1) * scale_factors (1, num_srcs) -> (nt, num_srcs)
    source_grid = cp.outer(src_time_gpu, scale_factors)
    b_vec = source_grid.flatten()
    
    # Matrix-Vector Multiply (One big kernel)
    d_pred = G_dense @ b_vec
        
    diff = d_pred - d_obs_gpu
    loss = (diff ** 2).sum()
    
    return float(loss)

def gpu_simulation_loss_geometric(v_test, rho_bg, rho_list, A_geom_gpu, dt, nt, x_rec_idx_gpu, x_src_idx_gpu, d_obs_gpu, b_vec_skeleton_gpu, src_time_gpu, unique_src_indices):
    """
    Computes Geometric Misfit loss on GPU.
    Minimizes || (I - G * pinv(G)) * d_obs ||^2
    """
    G_dense, _ = _compute_G_dense_gpu(v_test, rho_bg, rho_list, A_geom_gpu, dt, nt, x_rec_idx_gpu, x_src_idx_gpu)
    
    # cupy.linalg.lstsq solves G * x = d_obs in least-squares sense
    sol, residuals, rank, s = cp.linalg.lstsq(G_dense, d_obs_gpu, rcond=None)
    
    if residuals.size > 0:
        loss = residuals.sum()
    else:
        proj = G_dense @ sol
        diff = d_obs_gpu - proj
        loss = (diff.flatten() ** 2).sum()

    return float(loss)