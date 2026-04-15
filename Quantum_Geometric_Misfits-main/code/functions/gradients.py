"""
Shared adjoint gradient computation utilities for FWI.
Contains common functions used across L2 and Geometric misfit implementations.
"""

import numpy as np

import scipy.sparse as sp
from scipy.ndimage import gaussian_filter
from functions.assembly import A_assembly, Binv_sqrt_staggered, get_dof_index, reduced_space_time_system_Ghat_no_block as reduced_space_time_system_Ghat
from functions.geometric_mis import give_misfits



def apply_source_taper(grad_v, src_coords, nx, nz, radius=3.0):
    """
    Applies a smooth tapered mute around source locations.
    The artifact often occurs if the radius is too small relative to the grid spacing.
    
    Parameters
    ----------
    grad_v : ndarray of shape (nx, nz)
        Gradient field to taper
    src_coords : list of tuples
        Source coordinates [(x1, z1), (x2, z2), ...]
    nx, nz : int
        Grid dimensions
    radius : float
        Tapering radius in grid points
        
    Returns
    -------
    ndarray
        Tapered gradient
    """
    taper = np.ones((nx, nz))
    x = np.arange(nx)
    z = np.arange(nz)
    X, Z = np.meshgrid(x, z, indexing='ij')
    
    for (sx, sz) in src_coords:
        dist_sq = (X - sx)**2 + (Z - sz)**2
        sigma = radius
        mask = 1.0 - np.exp(-dist_sq / (2 * sigma**2))
        taper *= mask
        
    return grad_v * taper


def construct_b_vec(b_diag_local, wavelet_array, src_idxs, nt, dt):
    """
    Constructs the source vector for a given wavelet and source locations.
    
    Parameters
    ----------
    b_diag_local : ndarray
        Diagonal of the Binv_sqrt matrix
    wavelet_array : ndarray of shape (nt,)
        Source time series (Ricker wavelet or equivalent)
    src_idxs : list or array
        Sorted indices of source locations on the full grid
    nt : int
        Number of time steps
    dt : float
        Time step size
        
    Returns
    -------
    ndarray of shape (nt * num_srcs,)
        Flattened source vector
    """
    num_srcs = len(src_idxs)
    b_vec = np.zeros(nt * num_srcs)
    for it in range(nt):
        for i in range(num_srcs):
            b_vec[it * num_srcs + i] = wavelet_array[it] * b_diag_local[src_idxs[i]] * dt
    return b_vec


def compute_src_indices_and_maps(src_coords, rec_coords, A_geom, N, d):
    """
    Compute source/receiver indices and boolean maps for grid operations.
    
    Parameters
    ----------
    src_coords : list of tuples
        Source coordinates
    rec_coords : list of tuples
        Receiver coordinates
    A_geom : sparse matrix
        Geometry matrix
    N : list of int
        Grid dimensions [nx, nz]
    d : list of float
        Grid spacing [dx, dz]
        
    Returns
    -------
    dict
        Dictionary with keys: 'src_idxs', 'x_srcs_map', 'rec_idxs', 'x_rec_map', 'p_all_map'
    """
    Np = N[0] * N[1]
    
    src_idxs = sorted([get_dof_index('p', c, N, d) for c in src_coords])
    x_srcs_map = np.zeros(A_geom.shape[0], dtype=bool)
    x_srcs_map[src_idxs] = True
    
    rec_idxs = sorted([get_dof_index('p', c, N, d) for c in rec_coords])
    x_rec_map = np.zeros(A_geom.shape[0], dtype=bool)
    x_rec_map[rec_idxs] = True
    
    p_all_map = np.zeros(A_geom.shape[0], dtype=bool)
    p_all_map[:Np] = True
    
    return {
        'src_idxs': src_idxs,
        'x_srcs_map': x_srcs_map,
        'rec_idxs': rec_idxs,
        'x_rec_map': x_rec_map,
        'p_all_map': p_all_map,
    }


def compute_gradient_frechet_L2(
    ix,
    iz,
    v_current,
    rho_bg,
    rho_list,
    A_geom,
    dt,
    nt,
    x_rec_map,
    x_srcs_map,
    d_obs,
    delta_v,
    src_time,
    src_idxs,
    current_loss,
):
    """
    Forward-difference Frechet gradient for a single pixel.
    """
    v_orig = v_current[ix, iz]

    v_current[ix, iz] = v_orig + delta_v
    k_p = (v_current.flatten("F") ** 2) * rho_bg
    B_p = Binv_sqrt_staggered(k_p, rho_list, sqrt=True)

    G_p = reduced_space_time_system_Ghat(
        A_geom * dt,
        B_p,
        nt,
        x_rec_map,
        x_srcs_map,
        use_real=True,
        verbose=False,
    )
    b_p_vec = construct_b_vec(B_p.diagonal(), src_time, src_idxs, nt, dt)
    d_p = (G_p @ b_p_vec).real
    loss_p = np.linalg.norm(d_obs - d_p) ** 2

    v_current[ix, iz] = v_orig
    grad_val = (loss_p - current_loss) / delta_v
    return ix, iz, grad_val



def compute_gradient_adjoint_L2(eps_current, mu_bg, A_geom, b_vec_curr, d_obs, mask_inner, 
                                  src_coords, rec_coords, dx, dz, dt, nt, apply_smoothing=True):
    """
    Computes L2 misfit loss and adjoint gradient for Electromagnetics.
    
    Parameters
    ----------
    eps_current : ndarray of shape (nx, nz)
        Current relative permittivity (eps_r)
    mu_bg : float
        Background magnetic permeability (mu_0)
    ...
    """
    eps_0 = 8.854187817e-12
    nx, nz = eps_current.shape
    Np = nx * nz
    N, d = [nx, nz], [dx, dz]

    # Material properties (Fortran order). kappa maps to 1/(eps*eps_0)
    eps_flat = eps_current.flatten('F')
    kappa_curr = 1.0 / (eps_flat * eps_0)
    n_mu = (A_geom.shape[0] - Np) // 2
    mu_list = [np.ones(n_mu) * mu_bg for _ in range(2)]
    Binv_sqrt_curr = Binv_sqrt_staggered(kappa_curr, mu_list, sqrt=True)
    
    # Maps
    maps = compute_src_indices_and_maps(src_coords, rec_coords, A_geom, N, d)
    src_idxs = maps['src_idxs']
    x_srcs_map = maps['x_srcs_map']
    rec_idxs = maps['rec_idxs']
    x_rec_map = maps['x_rec_map']
    p_all_map = maps['p_all_map']
    
    # Forward pass
    G_u = reduced_space_time_system_Ghat(
        A_geom * dt, Binv_sqrt_curr, nt, 
        p_x=p_all_map, p_b=x_srcs_map, 
        use_real=True, verbose=False
    )
    u_full_flat = (G_u @ b_vec_curr).real
    u_history = u_full_flat.reshape((nt, Np))
    
    # Residual
    local_rec_idxs = [idx for idx in rec_idxs if idx < Np] # Filter to keep only 'p' dofs
    d_pred_full = u_history[:, local_rec_idxs]
    d_pred_flat = d_pred_full.flatten()
    
    residual = d_pred_flat - d_obs
    loss = np.sum(residual**2)
    
    # Adjoint pass
    ################################################
    # Version 1
    # G_r = reduced_space_time_system_Ghat(
    #     A_geom * dt, Binv_sqrt_curr, nt, 
    #     p_x=x_rec_map, p_b=p_all_map, 
    #     use_real=True, verbose=False
    # )
    # lambda_full_flat = (G_r.T @ residual).real
    # lambda_history = lambda_full_flat.reshape((nt, Np))
    ################################################
    # Version 2
    num_recs = np.sum(x_rec_map)
    rec_res_reshaped = residual.reshape((nt, num_recs))
    rec_res_rev = rec_res_reshaped[::-1, :].flatten()
    
    G_adj_forward = reduced_space_time_system_Ghat(
        A_geom * dt, Binv_sqrt_curr, nt, 
        p_x=p_all_map, p_b=x_rec_map, use_real=True, verbose=False
    )
    lambda_fast = (G_adj_forward @ rec_res_rev).real
    
    lambda_fast_reshaped = lambda_fast.reshape((nt, Np))
    lambda_history = lambda_fast_reshaped[::-1, :]
    ################################################
    
    eps_0 = 8.854187817e-12
    # Gradient computation (first-order system)
    # Time derivative
    dudt = np.zeros_like(u_history)
    dudt[1:-1, :] = (u_history[2:, :] - u_history[:-2, :]) / (2 * dt)
    dudt[0, :] = (u_history[1, :] - u_history[0, :]) / dt
    dudt[-1, :] = (u_history[-1, :] - u_history[-2, :]) / dt
    
    scattering_term = np.sum(lambda_history * dudt, axis=0) * 2
    
    # Source term
    src_term_grid = np.zeros(Np)
    if b_vec_curr is not None:
        num_srcs = len(src_idxs)
        if num_srcs > 0:
            b_reshaped = b_vec_curr.reshape((nt, num_srcs))
            lambda_at_src = lambda_history[:, src_idxs]
            src_dot_time = np.sum(lambda_at_src * b_reshaped, axis=0)
            src_term_grid[src_idxs] = src_dot_time

    # derivative w.r.t permittivity eps_r using chain rule from acoustic formula D_j/D_v to D_j/D_eps_r
    # Acoustic used: dJ/dv = (scat + src) / v 
    # That means dJ/d(v^2) was proportional to (scat + src)/ (2 v^2)
    # But in EM: kappa = 1 / (eps_r * eps_0)
    # So dJ/deps_r = dJ/dkappa * dkappa/deps_r
    # We will use the direct transformation here:
    v_equiv = 1.0 / np.sqrt(eps_flat * eps_0 * mu_bg)
    raw_grad_v = (scattering_term + src_term_grid) / (v_equiv + 1e-6) * dt
    # d_eps_r / d_v = - 2 / (v^3 * eps_0 * mu_bg)
    raw_grad_eps = raw_grad_v * (-2.0 / (v_equiv**3 * eps_0 * mu_bg + 1e-12))

    grad_eps = raw_grad_eps.reshape((nx, nz), order='F')
    
    # Post-processing
    grad_eps = apply_source_taper(grad_eps, src_coords, nx, nz, radius=3.0)
    
    if apply_smoothing:
        grad_eps = gaussian_filter(grad_eps, sigma=0.7)

    return loss, grad_eps, d_pred_flat


def compute_gradient_frechet_geometric(ix, iz, eps_current, mu_bg, mu_list, eps_0, A_geom, dt, nt, x_rec_map, x_srcs_map, d_obs, delta_eps, src_time, unique_src_indices, num_srcs, current_loss):
    """
    Helper function to compute gradient for a single pixel using Forward Difference for permittivity.
    """
    eps_orig = eps_current[ix, iz]

    # --- Perturb ---
    eps_current[ix, iz] = eps_orig + delta_eps
    k_p = 1.0 / (eps_current.flatten('F') * eps_0)
    B_p = Binv_sqrt_staggered(k_p, mu_list, sqrt=True)
    G_p = reduced_space_time_system_Ghat(A_geom * dt, B_p, nt, p_x=x_rec_map, p_b=x_srcs_map, use_real=True, verbose=False)
    
    # --- Calculate Geometric Misfit ---
    b_p_vec = construct_b_vec(B_p.diagonal(), src_time, unique_src_indices, nt, dt)
    loss_p, _, _ = give_misfits(G_p, d_obs, b_p_vec)

    # Restore value
    eps_current[ix, iz] = eps_orig
    
    # Forward Difference: (f(x+h) - f(x)) / h
    grad_val = (loss_p - current_loss) / delta_eps
    return ix, iz, grad_val


def compute_gradient_adjoint_geometric(eps_current, mu_bg, mu_list, A_geom, dt, nt, x_rec_map, x_srcs_map, d_obs):
    """
    Computes the Adjoint Gradient for the Geometric Misfit based on the .tex derivation.
    """
    eps_0 = 8.854187817e-12
    nx, nz = eps_current.shape
    Np = nx * nz
    
    eps_flat = eps_current.flatten('F')
    kappa_curr = 1.0 / (eps_flat * eps_0)
    Binv_sqrt_curr = Binv_sqrt_staggered(kappa_curr, mu_list, sqrt=True)

    # 1. Base Operator G
    G_raw = reduced_space_time_system_Ghat(A_geom * dt, Binv_sqrt_curr, nt, p_x=x_rec_map, p_b=x_srcs_map, use_real=True, verbose=False)
    if sp.issparse(G_raw):
        G_dense = G_raw.toarray()
    else:
        G_dense = np.asarray(G_raw)
    
    # 2. Geometric Misfit Components
    # Regularize inversion to prevent huge w_vec and unstable gradients
    # rcond_val = 1e-4
    
    # Alternative to np.linalg.pinv 
    #####################################
    eps = 1e-4 * np.trace(G_dense.T @ G_dense) / G_dense.shape[1]
    eps = max(eps, 1e-12)
    G_pinv = np.linalg.inv(G_dense.T @ G_dense + eps * np.eye(G_dense.shape[1])) @ G_dense.T
    #####################################

    R_guess = (G_dense @ G_pinv).real
    comp_A_vec = (np.eye(R_guess.shape[0]) - R_guess) @ d_obs
    w_vec = G_pinv @ d_obs

    # 3. Forward Pass (inject reconstructed source w_vec and compute wavefield u in all the domain)
    p_all_map = np.zeros(A_geom.shape[0], dtype=bool)
    p_all_map[:Np] = True
    G_u = reduced_space_time_system_Ghat(A_geom * dt, Binv_sqrt_curr, nt, p_x=p_all_map, p_b=x_srcs_map, use_real=True, verbose=False)
    u_full_flat = (G_u @ w_vec).real
    u_history = u_full_flat.reshape((nt, Np))

    # 4. Adjoint Pass (inject Geometric Misfit Residual)
    # The gradient formula is dJ/dm = -2 Re( u^* (I - P) dG G^+ u ).
    # Note from Golub and Pereyra the variation in sign leads to:
    # dJ/dm = -2 u^* (I-P) * dG(m) * w
    # 
    # To match FD direction identically and standard notation conventions,
    # the adjoint source applied backwards must be: -2.0 * comp_A_vec (comp_A_vec is r on the .tex so the residual) 
    adj_src = -2.0 * comp_A_vec
    
    
    # We inject the time-reversed adjoint source at receiver locations
    # and perform a single forward simulation, then time-reverse the result.
    num_recs = np.sum(x_rec_map)
    adj_src_reshaped = adj_src.reshape((nt, num_recs))
    adj_src_rev = adj_src_reshaped[::-1, :].flatten()
    
    # Reverse-time forward simulation with adjoint source to get lambda(t)
    # Inject at receiver and record at all grid points (obtain lambda_history)
    G_adj_forward = reduced_space_time_system_Ghat(
        A_geom * dt, Binv_sqrt_curr, nt, 
        p_x=p_all_map, p_b=x_rec_map, 
        use_real=True, verbose=False, auto_substeps=False, rk4_substeps=25)
    lambda_fast = (G_adj_forward @ adj_src_rev).real
    lambda_fast_reshaped = lambda_fast.reshape((nt, Np))
    lambda_history = lambda_fast_reshaped[::-1, :]

    # 5. Adjoint Gradient Computation (Scattering term)
    # Using the correct first-order system gradient logic similar to L2
    dudt = np.zeros_like(u_history)
    dudt[1:-1, :] = (u_history[2:, :] - u_history[:-2, :]) / (2 * dt)
    dudt[0, :] = (u_history[1, :] - u_history[0, :]) / dt
    dudt[-1, :] = (u_history[-1, :] - u_history[-2, :]) / dt
    
    # dudt is basically dC/dm w if no sources are present.
    scattering_term = np.sum(lambda_history * dudt, axis=0) * 2.0
    
    src_term_grid = np.zeros(Np)
    if w_vec is not None:
        # w_vec has shape (nt * num_srcs)
        unique_src_indices_loc = np.where(x_srcs_map)[0]
        num_srcs_loc = len(unique_src_indices_loc)
        if num_srcs_loc > 0:
            b_reshaped = w_vec.reshape((nt, num_srcs_loc))
            lambda_at_src = lambda_history[:, unique_src_indices_loc]
            src_dot_time = np.sum(lambda_at_src * b_reshaped, axis=0)
            src_term_grid[unique_src_indices_loc] = src_dot_time

    # Taking into account of the source term and applying EM derivative transformation
    v_equiv = 1.0 / np.sqrt(eps_flat * eps_0 * mu_bg)
    raw_grad_v = (scattering_term + src_term_grid) / (v_equiv + 1e-6) * dt
    
    # d_eps_r / d_v mapping:
    raw_grad_eps = raw_grad_v * (-2.0 / (v_equiv**3 * eps_0 * mu_bg + 1e-12))
    
    return raw_grad_eps.reshape((nx, nz), order='F')


