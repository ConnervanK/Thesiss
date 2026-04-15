import numpy as np
import os
import sys
import matplotlib
matplotlib.use('Agg')
import scipy.sparse as sp

# ensure project root on path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def geometric_misfit(R, x0hat_vec):
    """Return geometric misfit ||(I - R) x||^2 for a projection matrix R."""
    x0hat_vec = np.asarray(x0hat_vec)
    geometric_mis = (np.eye(R.shape[0], dtype=np.asarray(R).dtype) - R) @ x0hat_vec
    return np.linalg.norm(geometric_mis) ** 2

def give_misfits(Ghat_guess, x0hat_vec, b, check_pythagorean=True, pythagorean_tol=1e-6): 
    """
    Computes the geometric and parametric misfit components.
    """       
    if sp.issparse(Ghat_guess):
        G_dense = Ghat_guess.toarray()
    else:
        G_dense = np.asarray(Ghat_guess)

    # R = G * pinv(G) is the projection matrix onto the 'guess' subspace
    # Note: Moore-Penrose pseudo-inverse with rcond parameter to regularize
    rcond_val = 1e-4
    R_guess = (G_dense @ np.linalg.pinv(G_dense, rcond=rcond_val)).real

    # 1. Component A: Irreducible/Geometric Misfit vector
    comp_A = (np.eye(R_guess.shape[0]) - R_guess) @ x0hat_vec
    comp_A_norm_sq = np.linalg.norm(comp_A)**2


    # 2. Component B: Amplitude Misfit vector
    comp_B = R_guess @ x0hat_vec - Ghat_guess @ b 
    comp_B_norm_sq = np.linalg.norm(comp_B)**2

    # 3. Total Error in the unitary space
    total_error_vec = x0hat_vec - (Ghat_guess @ b)
    total_error_norm_sq = np.linalg.norm(total_error_vec)**2
    
    if check_pythagorean:
        check = (comp_A_norm_sq + comp_B_norm_sq) - total_error_norm_sq
        scale = max(1.0, total_error_norm_sq)
        if abs(check) > pythagorean_tol * scale:
            print("Warning: Pythagorean relation does not hold within numerical precision.")
            print(f"Check Value (Should be close to 0): {check:.6e}")

    return comp_A_norm_sq, comp_B_norm_sq, total_error_norm_sq
