import numpy as np
import scipy.sparse as sp
from joblib import Parallel, delayed
from tqdm import tqdm

# GPU-accelerated versions of the full space-time system G computation and Ghat assembly.
try:
    import cupy as cp
    import cupyx.scipy.sparse as cps
    import cupyx.scipy.linalg as cpl
except ImportError:
    cp = None
    cps = None
    cpl = None

def diff_staggered_1d(N_cell, dh, bc_start, bc_end):
    N_in = N_cell      
    N_out = N_cell + 1 
    
    # 1. Standard Interior Stencil (COO)
    data_main = np.ones(N_in) / dh
    row_main = np.arange(N_in)
    col_main = np.arange(N_in)
    
    data_sub = -np.ones(N_in) / dh
    row_sub = np.arange(1, N_out)
    col_sub = np.arange(N_in)

    data = np.concatenate([data_main, data_sub])
    row = np.concatenate([row_main, row_sub])
    col = np.concatenate([col_main, col_sub])

    D1d = sp.coo_matrix((data, (row, col)), shape=(N_out, N_in)).tolil()
    
    # 2. Corrected Boundary Conditions for Maxwell's Equations (e.g. Ez, Hx, Hy)
    
    # Start Boundary (x=0)
    if bc_start in ['PMC', 'PML']: 
        # PMC (Perfect Magnetic Conductor): H[0] = 0. 
        # We zero out the row so dH/dt = 0
        D1d[0, :] = 0.0
    elif bc_start in ['PEC', 'PEC_Surface']: 
        # PEC (Perfect Electric Conductor): E=0 at x=0.
        # On a staggered grid, E[0] is half a cell in. 
        # To get E=0 at the boundary, the "ghost" E[-1] must be -E[0].
        # The derivative at the boundary is (E[0] - E[-1])/dh = (E[0] - (-E[0]))/dh = 2*E[0]/dh.
        D1d[0, 0] = 2.0 / dh

    # End Boundary (x=L)
    if bc_end in ['PMC', 'PML']:
        # PMC: H[N] = 0
        D1d[N_cell, :] = 0.0
    elif bc_end in ['PEC', 'PEC_Surface']:
        # PEC: E=0 at x=L.
        # E[N-1] is half a cell before the boundary. Ghost E[N] = -E[N-1].
        # derivative = (E[N] - E[N-1])/dh = (-E[N-1] - E[N-1])/dh = -2*E[N-1]/dh.
        D1d[N_cell, N_cell-1] = -2.0 / dh
        
    return D1d.asformat('csr')


def operators_staggered_2d(Nx, Ny, dx, dy, bc_x_start, bc_x_end, bc_y_start, bc_y_end):
    """
    Assembles the 2D Staggered Curl Operators (CurlE and CurlH) for Maxwell's equations.
    
    Returns:
    - CurlE, CurlH: The 2D Curl sparse matrices representing TMz or TEz formulations.
    - NEz, NHx, NHy: The sizes of the vectors (for TMz mode).
    """
    
    # --- 1. Construct 1D Operators with BCs ---
    
    Dx_Gr = diff_staggered_1d(Nx, dx, bc_x_start, bc_x_end)
    Dy_Gr = diff_staggered_1d(Ny, dy, bc_y_start, bc_y_end)
    
    Ix = sp.eye(Nx, format='csr')
    Iy = sp.eye(Ny, format='csr')

    # # Sizes
    # NEz = Nx * Ny          # Ez grid points
    # NHx = (Nx + 1) * Ny    # Hx grid points
    # NHy = Nx * (Ny + 1)    # Hy grid points
    
    # --- 2. Assemble 2D Curl Operator on E field (CurlE) ---
    # In TMz: dEz/dy -> Hx, -dEz/dx -> Hy
    
    # Gx: I_y * D_x_Gr. Maps Ez (Nx*Ny) -> Hx or Hy depending on formulation
    Grx = sp.kron(Iy, Dx_Gr, format='csr') 
    
    # Gy: D_y_Gr * I_x. Maps Ez -> Hx or Hy
    Gry = sp.kron(Dy_Gr, Ix, format='csr') 

    # For TMz mode specifically: CurlE maps Ez to (Hx, Hy) components natively.
    CurlE = sp.vstack([Gry, -Grx], format='csr')
    
    # --- 3. Assemble 2D Curl Operator on H fields (CurlH) from Transpose ---
    # CurlH maps (Hx, Hy) -> Ez; it's proportional to the transpose of CurlE
    CurlH = -CurlE.transpose().tocsr()
    
    return CurlE, CurlH
# Note: The create_1d_staggered_derivative function defined previously 
# will be reused here for Gr (Gradient).


# Assuming create_1d_staggered_derivative is defined and available
# (The COO-based version from our previous exchange)

def operators_staggered_3d(Nx, Ny, Nz, dx, dy, dz, 
                                    bc_x_start, bc_x_end, 
                                    bc_y_start, bc_y_end,
                                    bc_z_start, bc_z_end):
    """
    Assembles the 3D Staggered Gradient (Gr) and Divergence (Div) operators.
    
    Returns:
    - Gr, Div: The 3D Gradient and Divergence sparse matrices.
    - Np, Nvx, Nvy, Nvz: The sizes of the P, V_x, V_y, and V_z vectors.
    """
    
    # --- 1. Construct 1D Operators ---
    
    # 1D Derivatives (Gradient components)
    Dx_Gr = diff_staggered_1d(Nx, dx, bc_x_start, bc_x_end)
    Dy_Gr = diff_staggered_1d(Ny, dy, bc_y_start, bc_y_end)
    Dz_Gr = diff_staggered_1d(Nz, dz, bc_z_start, bc_z_end)
    
    # 1D Identity matrices
    Ix = sp.eye(Nx, format='csr')
    Iy = sp.eye(Ny, format='csr')
    Iz = sp.eye(Nz, format='csr')

    # Sizes
    # Np = Nx * Ny * Nz                     # P grid points
    # Nvx = (Nx + 1) * Ny * Nz              # V_x grid points
    # Nvy = Nx * (Ny + 1) * Nz              # V_y grid points
    # Nvz = Nx * Ny * (Nz + 1)              # V_z grid points
    
    # --- 2. Assemble 3D Gradient Operator (G) using Kronecker Product ---
    
    # Gx (d/dx): I_z * I_y * D_x_Gr. Maps P (Np) -> V_x (Nvx)
    # The order for Kronecker product is often Z-Y-X (slowest to fastest index)
    Grx = sp.kron(Iz, sp.kron(Iy, Dx_Gr, format='csr'), format='csr') 
    
    # Gy (d/dy): I_z * D_y_Gr * I_x. Maps P (Np) -> V_y (Nvy)
    Gry = sp.kron(Iz, sp.kron(Dy_Gr, Ix, format='csr'), format='csr') 

    # Gz (d/dz): D_z_Gr * I_y * I_x. Maps P (Np) -> V_z (Nvz)
    Grz = sp.kron(Dz_Gr, sp.kron(Iy, Ix, format='csr'), format='csr') 

    # Gr stacks Grx, Gry, Grz vertically
    Gr = sp.vstack([Grx, Gry, Grz], format='csr')
    
    # --- 3. Assemble 3D Divergence Operator (Div) from Transpose ---

    # Div maps V (Nvx + Nvy + Nvz) -> P (Np)
    # Div is the NEGATIVE transpose of Gr for skew-symmetry
    Div = -Gr.transpose().tocsr()
    
    return Gr, Div

def A_staggered_1d(N, dx, bc_start, bc_end):
    """
    Assembles the full sparse spatial operator matrix A_spatial for the 
    1D acoustic wave equation on a STAGGERED grid, EXCLUDING material properties.
    
    The resulting matrix is A_spatial = 
        | 0 | -D |
        | -G | 0 |
    
    Parameters:
    - N: Number of pressure grid points (N_p).
    - dx: Grid spacing.
    - bc_start: Boundary condition at the start.
    - bc_end: Boundary condition at the end.

    Returns:
    - A_spatial: The full (N_p + N_v x N_p + N_v) sparse system matrix.
    """
    
    # --- 1. Construct 1D Operators ---
    
    # G (Gradient): Maps P (N_p) -> V (N_v = N_p + 1). 
    G = diff_staggered_1d(N, dx, bc_start, bc_end)
    
    # D (Divergence): Maps V (N_v) -> P (N_p). 
    # Must satisfy D = -G^T for skew-symmetry and energy conservation.
    D = -G.transpose().tocsr()
    
    # Sizes
    Np = N          # Pressure points
    Nv = N + 1      # Velocity points
    
    # --- 2. Assemble the final block matrix ---

    # Zero matrix blocks
    Zero_P = sp.csc_matrix((Np, Np))
    Zero_V = sp.csc_matrix((Nv, Nv))
    
    # Assembly:
    # Row 1 (for dp/dt): [ 0 | -D ]
    # Row 2 (for dv/dt): [ -G | 0 ]
    
    # A = sp.bmat([
    #     [Zero_P, -D],     # Row 1: dp/dt
    #     [-G, Zero_V]      # Row 2: dv/dt
    # ], format='csc')

    return sp.bmat([[Zero_P, -D],[ -G, Zero_V]], format='csc')

def A_staggered_2d(Nx, Ny, dx, dy, bc_x_start, bc_x_end, bc_y_start, bc_y_end):
    """
    Assembles the full sparse spatial operator matrix A_spatial for 
    2D Maxwell equations (TMz) on a STAGGERED grid, EXCLUDING material properties.
    
    Returns:
    - A_spatial: The full sparse system matrix.
    """
    # Sizes
    NEz = Nx * Ny          # Ez grid points
    NHx = (Nx + 1) * Ny    # Hx grid points
    NHy = Nx * (Ny + 1)    # Hy grid points

    # 1. Assemble the Curl operators
    CurlE, CurlH = operators_staggered_2d(
        Nx, Ny, dx, dy, bc_x_start, bc_x_end, bc_y_start, bc_y_end)
    
    # 2. Identify the sub-matrices from CurlE and CurlH
    
    # CurlH: (NEz x N_h) where N_h = NHx + NHy
    # In TMz, dEz/dt = 1/e * (dHy/dx - dHx/dy)
    # The CurlH matrix inherently represents this combination mapping (Hx, Hy) -> Ez
    CurlHx = CurlH[:, :NHx] # Hx component -> (NEz x NHx)
    CurlHy = CurlH[:, NHx:] # Hy component -> (NEz x NHy)
    
    # CurlE: (N_h x NEz) where N_h = NHx + NHy
    CurlEx = CurlE[:NHx, :] # (NHx x NEz)
    CurlEy = CurlE[NHx:, :] # (NHy x NEz)

    # 3. Zero matrix blocks
    Zero_Ez = sp.csc_matrix((NEz, NEz))
    Zero_Hx = sp.csc_matrix((NHx, NHx))
    Zero_Hy = sp.csc_matrix((NHy, NHy))
    
    # 4. Assemble the final block matrix
    # A_spatial = 
    # [ 0,     CurlHx, CurlHy ]
    # [ CurlEx,    0,      0   ]
    # [ CurlEy,    0,      0   ]
    
    return sp.bmat([[Zero_Ez, CurlHx, CurlHy],[ CurlEx, Zero_Hx, None],[ CurlEy, None, Zero_Hy]], format='csc')

def A_staggered_3d(Nx, Ny, Nz, dx, dy, dz, 
                   bc_x_start, bc_x_end, 
                   bc_y_start, bc_y_end,
                   bc_z_start, bc_z_end):
    """
    Assembles the full sparse spatial operator matrix A_spatial for the 
    3D acoustic wave equation on a STAGGERED grid (material-free).
    
    The resulting matrix is A_spatial = 
        | 0 | -Dx | -Dy | -Dz |
        | -Gx | 0 | 0 | 0 |
        | -Gy | 0 | 0 | 0 |
        | -Gz | 0 | 0 | 0 |
    
    Returns:
    - A_spatial: The full (Np+Nvx+Nvy+Nvz x Np+Nvx+Nvy+Nvz) sparse system matrix.
    """
    Np = Nx * Ny * Nz                     # P grid points
    Nvx = (Nx + 1) * Ny * Nz              # V_x grid points
    Nvy = Nx * (Ny + 1) * Nz              # V_y grid points
    Nvz = Nx * Ny * (Nz + 1)              # V_z grid points
    
    # 1. Assemble the G and D operators
    Gr, Div = operators_staggered_3d(
        Nx, Ny, Nz, dx, dy, dz, 
        bc_x_start, bc_x_end, bc_y_start, bc_y_end, bc_z_start, bc_z_end
    )
    
    # 2. Identify the sub-matrices from D and G
    
    # D: (Np x N_v) where N_v = Nvx + Nvy + Nvz
    Divx = Div[:, :Nvx]          # d/dx component of Divergence (Np x Nvx)
    Divy = Div[:, Nvx:Nvx+Nvy]   # d/dy component of Divergence (Np x Nvy)
    Divz = Div[:, Nvx+Nvy:]      # d/dz component of Divergence (Np x Nvz)
    
    # G: (N_v x Np) where N_v = Nvx + Nvy + Nvz
    Grx = Gr[:Nvx, :]          # d/dx component of Gradient (Nvx x Np)
    Gry = Gr[Nvx:Nvx+Nvy, :]   # d/dy component of Gradient (Nvy x Np)
    Grz = Gr[Nvx+Nvy:, :]      # d/dz component of Gradient (Nvz x Np)

    # 3. Zero matrix blocks
    # We only need the diagonal zero blocks for the velocity partition
    Zero_P = sp.csc_matrix((Np, Np))
    Zero_Vx = sp.csc_matrix((Nvx, Nvx))
    Zero_Vy = sp.csc_matrix((Nvy, Nvy))
    Zero_Vz = sp.csc_matrix((Nvz, Nvz))
    
    # 4. Assemble the final block matrix using sp.bmat
    
    # A_spatial = sp.bmat([
    #     [Zero_P, -Dx, -Dy, -Dz],      # Row 1: dp/dt
    #     [-Grx, Zero_Vx, None, None],   # Row 2: dvx/dt
    #     [-Gry, None, Zero_Vy, None],   # Row 3: dvy/dt
    #     [-Grz, None, None, Zero_Vz]    # Row 4: dvz/dt
    # ], format='csc')

    return sp.bmat([[Zero_P,-Divx,-Divy,-Divz],[ -Grx, Zero_Vx, None, None],[ -Gry, None, Zero_Vy, None],[ -Grz, None, None, Zero_Vz]], format='csc')

def A_assembly(N,d,boundaries):
    """
    General N-dimensional assembly of the spatial operator A for the staggered-grid acoustic system.
    Parameters:
    - N: List of grid dimensions [Nx] or [Nx, Ny] or [Nx, Ny, Nz].
    - d: List of grid spacings [dx] or [dx, dy] or [dx, dy, dz].
    - boundaries: List of boundary condition types for each dimension, e.g. ['Dirichlet', 'Neumann', ...].
    """
    N_dim = len(N)
    
    if len(d) != N_dim:
        raise ValueError(f"Grid spacings list 'd' must have length {N_dim} for {N_dim}D simulation. Got {len(d)}.")
    if len(boundaries) != 2 * N_dim:
        raise ValueError(f"Boundaries list must have length {2 * N_dim} for {N_dim}D simulation. Got {len(boundaries)}.")
    
    if N_dim == 1:
        return A_staggered_1d(N[0], d[0], boundaries[0], boundaries[1])
    elif N_dim == 2:
        return A_staggered_2d(N[0], N[1], d[0], d[1], boundaries[0], boundaries[1], boundaries[2], boundaries[3])
    elif N_dim == 3:
        return A_staggered_3d(N[0], N[1], N[2], d[0], d[1], d[2],
                              boundaries[0], boundaries[1],
                              boundaries[2], boundaries[3],
                              boundaries[4], boundaries[5])
    else:
        raise ValueError("N-dimensional assembly only supports 1D, 2D, and 3D cases.")

def Binv_sqrt_staggered(kappa_p, rho_v_components, sqrt=True):
    """
    Assembles the general N-dimensional diagonal material properties matrix B 
    for the staggered-grid acoustic system.

    Parameters:
    - kappa_p: 1D NumPy array containing the bulk modulus values for the pressure grid points.
    - rho_v_components: Additional 1D NumPy arrays containing the density values for the velocity components, in order:
                       [rho_vx, rho_vy, rho_vz, ...]
                       
                       Note: The arrays must already be correctly ordered (lexicographically 
                       flattened) and correctly sized for the system.

    Returns:
    - B_spatial: The diagonal sparse matrix (in CSC format) of size 
                 (Total_DOF x Total_DOF), where the diagonal contains 
                 [1/kappa_p, 1/rho_vx, 1/rho_vy, 1/rho_vz, ...]
    """
    
    # 1. Calculate the inverse of each material property array
    # The first array is kappa (bulk modulus) for P, the rest are rho (density) for V.
    if sqrt:
        # Inverse Square Root
        inv_kappa_p = np.sqrt(kappa_p)
        inv_rho_v_components = [ 1.0 / np.sqrt(rho_v) for rho_v in rho_v_components]
    else:
        # Inverse
        inv_kappa_p = kappa_p
        inv_rho_v_components = [ 1.0 / rho_v for rho_v in rho_v_components]
    
    # Start with the P components
    data = [inv_kappa_p]
    
    # Add the V components (Vx, Vy, Vz, ...)
    data.extend(inv_rho_v_components)
    
    # Concatenate all parts into a single 1D vector
    combined_data = np.concatenate(data)
    
    # 3. Create the diagonal sparse matrix
    # Binv_sqrt_nd = sp.diags(combined_data, offsets=0, format='csc')
    
    return sp.diags(combined_data, offsets=0, format='csc')

def create_nd_inclusion_arrays(N_dims = None, kappa_bg=None, rho_bg=None, inclusion_start=None, inclusion_end=None, 
                               inclusion_kappa_val=None, inclusion_rho_val=None):
    """
    Generates flattened material property arrays (kappa_p and rho_v_components) 
    for an N-dimensional domain with a specified material inclusion.

    The inclusion properties (kappa and rho) are now defined directly by the user.

    The background properties (kappa_bg and rho_bg) can be either scalars (for a homogeneous 
    start) or existing arrays (to add multiple inclusions sequentially).

    Parameters:
    - N_dims: List of grid dimensions [Nx] or [Nx, Ny] or [Nx, Ny, Nz].
    - kappa_bg: Background bulk modulus. Can be a scalar (float) or an existing flattened 1D array.
    - rho_bg: Background density. Can be a scalar (float) or a list of flattened 1D arrays 
              [rho_vx, rho_vy, ...].
    - inclusion_start: List of start indices [x_start, y_start, z_start].
    - inclusion_end: List of end indices [x_end, y_end, z_end].
    - inclusion_kappa_val: The specific bulk modulus value (scalar) for the inclusion region.
    - inclusion_rho_val: The specific density value(s) for the inclusion region. 
                          Expected to be a scalar, or a list/array of length N_dim 
                          [rho_vx, rho_vy, rho_vz] for component-specific density assignment.

    Returns:
    - tuple: (kappa_p, rho_v_components) where kappa_p is the P-point bulk modulus array, 
             and rho_v_components is a list of V-point density arrays (one for each dimension).
    """
    
    N_dim = len(N_dims)
    
    # Pad N_dims and inclusion indices for easy slicing in the loop
    # For 1D/2D, we pad with 1/0 to make the structure uniform
    N = list(N_dims) + [1] * (3 - N_dim)
    
    Nx, Ny, Nz = N[0], N[1], N[2]
    
    # 1. Direct Assignment of Inclusion Properties
    
    # Handle kappa: If a list of length 1 is passed, extract the scalar.
    if isinstance(inclusion_kappa_val, (list, np.ndarray)) and len(inclusion_kappa_val) == 1:
        kappa_inclusion = inclusion_kappa_val[0]
    else:
        kappa_inclusion = inclusion_kappa_val

    # rho_vals handling is done inside the V-component loop (Step 3)

    
    # 2. Create Bulk Modulus (P-points) Array
    
    # P-points grid shape: (Nz, Ny, Nx) 
    P_shape = (Nz, Ny, Nx)
    
    # Initialize grid from scalar or existing array
    if np.ndim(kappa_bg) == 0:
        # Scalar background
        kappa_p_grid = np.full(P_shape, kappa_bg, dtype=float)
    else:
        # Array background
        if len(kappa_bg) != np.prod(P_shape):
             raise ValueError(f"Size of kappa_bg array ({len(kappa_bg)}) does not match grid size ({np.prod(P_shape)})")
        kappa_p_grid = np.array(kappa_bg).reshape(P_shape)

    # Calculate P-point slices (uses un-padded dimensions for indexing)
    # The slices are for [X, Y, Z] order in the inclusion_start/end lists
    p_slices_xyz = [slice(inclusion_start[i], inclusion_end[i]) for i in range(N_dim)]
    
    # Apply inclusion for P points (uses Z, Y, X index order)
    if N_dim == 1:
        # P_shape is (1, 1, Nx) -> use the X slice
        kappa_p_grid[..., p_slices_xyz[0]] = kappa_inclusion 
    elif N_dim == 2:
        # P_shape is (1, Ny, Nx) -> use the Y and X slices
        # Indexing: [Z=0, Y, X]
        kappa_p_grid[0, p_slices_xyz[1], p_slices_xyz[0]] = kappa_inclusion
    elif N_dim == 3:
        # P_shape is (Nz, Ny, Nx) -> use Z, Y, X slices
        # Indexing: [Z, Y, X]
        kappa_p_grid[p_slices_xyz[2], p_slices_xyz[1], p_slices_xyz[0]] = kappa_inclusion
    
    kappa_p = kappa_p_grid[:Nz, :Ny, :Nx].flatten() # Ensure we only flatten the active part
    
    
    # 3. Create Density (V-points) Arrays
    
    rho_v_components = []
    
    # V-components loop: i=0 (Vx), i=1 (Vy), i=2 (Vz)
    for i in range(N_dim):
        # Determine the inclusion rho value for this specific V-component
        if isinstance(inclusion_rho_val, (list, np.ndarray)):
            if i < len(inclusion_rho_val):
                 rho_inclusion = inclusion_rho_val[i]
            else:
                 rho_inclusion = rho_bg if np.ndim(rho_bg)==0 else rho_bg[i] # Fallback
        else:
            rho_inclusion = inclusion_rho_val
        # Determine the shape of the current V component grid
        v_shape_list = [Nx, Ny, Nz]
        v_shape_list[i] += 1
        v_shape_list = v_shape_list[:N_dim] 

        rho_v_flat_size = np.prod(v_shape_list)
        
        # --- Handle Scalar vs List/Array rho_bg ---
        current_bg_rho = None
        if isinstance(rho_bg, list) or (isinstance(rho_bg, np.ndarray) and rho_bg.ndim > 0 and len(rho_bg) == N_dim and hasattr(rho_bg[0], '__len__')):
             # It is a list of arrays (one per component)
             current_bg_rho = rho_bg[i]
        elif np.ndim(rho_bg) == 0:
             # Scalar
             current_bg_rho = rho_bg
        elif isinstance(rho_bg, (list, np.ndarray)) and len(rho_bg) == N_dim:
             # List of scalars (e.g. [1000, 1000])
             current_bg_rho = rho_bg[i]
        else:
             # Fallback for simple scalar wrapped in array
              if np.size(rho_bg) == 1:
                  current_bg_rho = np.array(rho_bg).item()
              else:
                   pass

        # Initialize background array
        if np.ndim(current_bg_rho) == 0:
             rho_v = np.full(rho_v_flat_size, current_bg_rho, dtype=float)
        else:
              if len(current_bg_rho) != rho_v_flat_size:
                   raise ValueError(f"Size of rho_bg component {i} ({len(current_bg_rho)}) does not match grid size ({rho_v_flat_size})")
              rho_v = np.array(current_bg_rho)

        
        # Reshape for 3D indexing (Z, Y, X order)
        if N_dim == 1:
            rho_v_grid = rho_v # Already 1D
            v_grid_shape = (v_shape_list[0],)
        elif N_dim == 2:
            v_grid_shape = (v_shape_list[1], v_shape_list[0]) # (Y, X)
            rho_v_grid = rho_v.reshape(v_grid_shape)
        elif N_dim == 3:
            v_grid_shape_3d = list(v_shape_list)[::-1]
            rho_v_grid = rho_v.reshape(v_grid_shape_3d)

        # Build slices for the inclusion: N_dim slices in X, Y, Z order
        v_slices_xyz = []
        for j in range(N_dim):
            start = inclusion_start[j]
            end = inclusion_end[j]
            if j == i: 
                v_slices_xyz.append(slice(start, end + 1))
            else: 
                v_slices_xyz.append(slice(start, end))

        # Apply inclusion to the reshaped grid
        if N_dim == 1:
             rho_v_grid[v_slices_xyz[0]] = rho_inclusion
        elif N_dim == 2:
            rho_v_grid[v_slices_xyz[1], v_slices_xyz[0]] = rho_inclusion
        elif N_dim == 3:
            rho_v_grid[v_slices_xyz[2], v_slices_xyz[1], v_slices_xyz[0]] = rho_inclusion
        
        rho_v_components.append(rho_v_grid.flatten())
    
    # We only return the components needed for the given N_dim
    return kappa_p, rho_v_components[:N_dim]

def U_assembly(H, dt):
    """
    Docstring for U_assembly
    
    :param A: Description
    :param Binv2: Description
    :param dt: Description
    """
    # H = -1j * Binv2 @ A @ Binv2
    return sp.linalg.expm(dt*H)

def get_dof_index(field, coords, N, d):
    """
    Get the global index in the state vector for a specific field and coordinate.
    
    Parameters:
    - field: str, one of 'p', 'vx', 'vy', 'vz'
    - coords: list or tuple of coordinates (x, y, z). 
              If int, treated as grid index. 
              If float, treated as physical position and converted to nearest index.
    - N: list of grid sizes [Nx, Ny, Nz] (or fewer for 1D/2D)
    - d: list of grid spacings [dx, dy, dz] (or fewer for 1D/2D)
    
    Returns:
    - int: The global index in the state vector.
    """
    field = field.lower()
    ndim = len(N)
    
    # 1. Parse Coordinates to Indices
    indices = []
    for i, c in enumerate(coords):
        if isinstance(c, int):
            idx = c
        else:
            idx = int(round(c / d[i]))
        indices.append(idx)
        
    # Pad indices for 3D logic if lower dim
    ix = indices[0] if ndim >= 1 else 0
    iy = indices[1] if ndim >= 2 else 0
    iz = indices[2] if ndim >= 3 else 0
    
    Nx = N[0]
    Ny = N[1] if ndim >= 2 else 1
    Nz = N[2] if ndim >= 3 else 1
    
    # 2. Calculate Offsets and Local Dimensions
    # State Vector Order: [P, Vx, Vy, Vz]
    
    # Sizes of each block
    size_p = Nx * Ny * Nz
    size_vx = (Nx + 1) * Ny * Nz
    size_vy = Nx * (Ny + 1) * Nz
    # size_vz = Nx * Ny * (Nz + 1) # Only needed if we go to Vz
    
    offset = 0
    
    # Determine which block we are in and the dimensions of that block
    if field == 'p':
        # P block is first
        offset = 0
        nx_local, ny_local, nz_local = Nx, Ny, Nz
        
    elif field == 'vx':
        offset = size_p
        nx_local, ny_local, nz_local = Nx + 1, Ny, Nz
        
    elif field == 'vy':
        if ndim < 2: raise ValueError("Vy not available in 1D")
        offset = size_p + size_vx
        nx_local, ny_local, nz_local = Nx, Ny + 1, Nz
        
    elif field == 'vz':
        if ndim < 3: raise ValueError("Vz not available in < 3D")
        offset = size_p + size_vx + size_vy
        nx_local, ny_local, nz_local = Nx, Ny, Nz + 1
        
    else:
        raise ValueError(f"Unknown field '{field}'")
        
    # 3. Check Bounds
    if not (0 <= ix < nx_local): raise ValueError(f"x index {ix} out of bounds for {field} (max {nx_local})")
    if ndim >= 2 and not (0 <= iy < ny_local): raise ValueError(f"y index {iy} out of bounds for {field} (max {ny_local})")
    if ndim >= 3 and not (0 <= iz < nz_local): raise ValueError(f"z index {iz} out of bounds for {field} (max {nz_local})")

    # 4. Calculate Flattened Index (Z-Y-X order)
    # index = iz * (ny_local * nx_local) + iy * nx_local + ix
    
    local_idx = 0
    if ndim == 1:
        local_idx = ix
    elif ndim == 2:
        local_idx = iy * nx_local + ix
    elif ndim == 3:
        local_idx = iz * (ny_local * nx_local) + iy * nx_local + ix
        
    return offset + local_idx


# Progress-tracked and optimized version of the full space-time system G computation and Ghat assembly.
def full_space_time_system_G(U, nt, n_jobs=-1, verbose=True):
    """
    Computes [U^0, U^1, ..., U^(T-1)] using parallel doubling with progress tracking.
    """
    # Initialize with Identity
    G = [sp.eye(U.shape[0], dtype=U.dtype, format='csc')]
    P = U.copy()
    block = 1

    # Initialize Progress Bar
    if verbose:
        pbar = tqdm(total=nt, desc="Parallel Power Doubling", unit="matrix")
        pbar.update(1)

    def multiply(P_op, M_op):
        return P_op @ M_op

    while len(G) < nt:
        chunk = G[-block:]
        need = min(block, nt - len(G))
        
        # Parallel computation of the next 'stripe'
        new = Parallel(n_jobs=n_jobs, prefer="threads")(
            delayed(multiply)(P, M) for M in chunk[:need]
        )
        
        G.extend(new)
        if verbose:
            pbar.update(len(new))
        
        # Update P for the next doubling step (P = U^block)
        if len(G) < nt:
            P = P @ P
            block <<= 1
    if verbose:
        pbar.close()
    return G

def reduced_space_time_system_Ghat(A, Binv_sqrt, nt, p_x, p_b, n_jobs=-1, use_real=True, verbose=True):
    """
    Constructs the Linear System (Ghat) with Progress Bars and Real-Valued optimization.
    """
    # 1. Construct the operator
    # OPTIMIZATION: If use_real is True, we simplify exp(-1j * (1j * AH)) to exp(AH)
    if use_real:
        # H is real, so expm(H) remains real. Memory usage is halved.
        H = Binv_sqrt @ A @ Binv_sqrt
    else:
        H = -1j * Binv_sqrt @ A @ Binv_sqrt
    
    if verbose:
        print(f"--- Operator Construction (Real Mode: {use_real}) ---")
    U = sp.linalg.expm(H) 
    
    # 2. Compute the full time evolution sequence
    G_list = full_space_time_system_G(U, nt, n_jobs=n_jobs, verbose=verbose)
    
    # 3. Pre-project the G matrices (Slicing)
    # This is much faster than slicing inside the nested block loop
    idx_x = p_x.astype(bool)
    idx_b = p_b.astype(bool)
    
    reduced_G = []
    if verbose:
        for g in tqdm(G_list, desc="Slicing Sub-blocks", unit="block"):
            reduced_G.append(g[idx_x, :][:, idx_b])
    else:
        for g in G_list:
            reduced_G.append(g[idx_x, :][:, idx_b])
    # 4. Construct the Block-Toeplitz System
    blocks = []
    # Outer progress bar for the block-row assembly
    if verbose:
        for j in tqdm(range(nt), desc="Assembling Block-Toeplitz", unit="row"):
            row = []
            for i in range(nt):
                if i <= j:
                    # Causal: Effect at time j from source at time i
                    row.append(reduced_G[j - i])
                else:
                    # Lower triangular: Future cannot affect the past
                    row.append(None)
            blocks.append(row)
    else:        
        for j in range(nt):
            row = []
            for i in range(nt):
                if i <= j:
                    row.append(reduced_G[j - i])
                else:
                    row.append(None)
            blocks.append(row)
    
    # 5. Final assembly into a large sparse matrix
    if verbose:
        print("--- Finalizing Sparse Assembly ---")
    Ghat = sp.block_array(blocks, format='csc')
    return Ghat


def full_space_time_system_G_gpu(U_gpu, nt, verbose=True):
    """
    GPU-accelerated doubling. U_gpu should be a cupyx.scipy.sparse matrix (or compatible).
    """
    G = [cps.eye(U_gpu.shape[0], dtype=U_gpu.dtype, format='csr')]
    P = U_gpu.copy()
    
    block = 1

    # Simple loop for doubling
    # (Parallelization is internal to CuPy kernels)
    
    if verbose:
        # We can implement a simple progress indicator if desired, 
        # but gpu is usually fast enough for small N.
        # For consistency with CPU version, we print.
        print(f"GPU Doubling: Target nt={nt}...", end="", flush=True)

    while len(G) < nt:
        need = min(block, nt - len(G))
        
        # Parallel computation of the next 'chunk'
        # On GPU, we execute these sequentially or batched. 
        # Simple loop is fine as operations are large.
        
        new_chunk = []
        for i in range(need):
             # G[len(G) - block + i] ? No, doubling logic is simpler:
             # multiply P (which is U^block) by the previous block of G
             # The existing logic in CPU version:
             # chunk = G[-block:] -> list of last 'block' matrices
             # M in chunk -> multiply P @ M
             
             # Actually, CPU version logic:
             # chunk = G[-block:]
             # need = min(block, T - len(G))
             # new = P @ M for M in chunk[:need]
             
             # Replicating that logic on GPU:
             M = G[-(block) + i] # The matrix from 'block' steps ago in the list? No.
             # CPU: chunk = G[-block:]. First element is G[len(G)-block].
             # so M corresponds to G[len(G)-block+i]
             
             M = G[len(G) - block + i]
             new_chunk.append(P @ M)
        
        G.extend(new_chunk)
        
        if len(G) < nt:
            P = P @ P
            block <<= 1
            
    if verbose:
        print(" Done.")
    return G

def reduced_space_time_system_Ghat_gpu(A, Binv_sqrt, nt, p_x, p_b, use_real=True, verbose=True):
    """
    GPU version of reduced system assembly.
    Returns a SciPy sparse matrix (CPU) for compatibility with downstream solvers.
    """
    # 1. Transfer to GPU if not already
    if sp.issparse(A):
        A_gpu = cps.csr_matrix(A)
    elif isinstance(A, np.ndarray):
        A_gpu = cps.csr_matrix(A) # Convert to sparse on GPU
    else: 
         # Assuming it might be already GPU array?
        A_gpu = A 
        if not cps.issparse(A_gpu):
             A_gpu = cps.csr_matrix(A_gpu)

    if sp.issparse(Binv_sqrt):
        B_gpu = cps.csr_matrix(Binv_sqrt)
    elif isinstance(Binv_sqrt, np.ndarray):
        B_gpu = cps.csr_matrix(Binv_sqrt)
    else:
        B_gpu = Binv_sqrt
        if not cps.issparse(B_gpu):
             B_gpu = cps.csr_matrix(B_gpu)

    
    H = B_gpu @ A_gpu @ B_gpu
    if not use_real:
        H = -1j * H
    
    if verbose:
        print(f"--- GPU Operator Construction (Real Mode: {use_real}) ---")

    # 2. Matrix Exponential on GPU
    # cupyx.scipy.linalg does not have expm for sparse matrices usually (maps to dense).
    # We convert to dense for the exponential if reasonable size.
    if verbose: print(f"Computing matrix exponential on GPU (Size: {H.shape})...")
    H_dense = H.toarray()
    U_dense = cpl.expm(H_dense)
    U_gpu = cps.csr_matrix(U_dense)

    
    # 3. Doubling on GPU
    G_list_gpu = full_space_time_system_G_gpu(U_gpu, nt, verbose=verbose)
    
    # 4. Slicing on GPU
    # Indices also need to be on device
    idx_x_gpu = cp.array(p_x.astype(bool))
    idx_b_gpu = cp.array(p_b.astype(bool))
    
    reduced_G_cpu = []
    if verbose:
        print("Slicing and transferring to CPU...")
        
    for g in G_list_gpu:
        # GPU slicing is much faster than CPU slicing
        # We slice on GPU, then transfer the small result to CPU
        sliced = g[idx_x_gpu, :][:, idx_b_gpu]
        reduced_G_cpu.append(sliced.get())
        
    # 5. Final assembly on CPU
    if verbose:
        print("Assembling Block-Toeplitz on CPU...")
        
    blocks = []
    for j in range(nt):
        row = []
        for i in range(nt):
            if i <= j:
                row.append(reduced_G_cpu[j - i])
            else:
                row.append(None)
        blocks.append(row)
    
    Ghat = sp.block_array(blocks, format='csc')
    return Ghat

def reduced_space_time_system_Ghat_no_block(
    A,
    Binv_sqrt,
    nt,
    p_x,
    p_b,
    use_real=True,
    verbose=True,
    rk4_substeps=1,
    integration_scheme="rk4",
    auto_substeps=True,
    max_auto_substeps=256,
):
    """
        Builds Ghat directly from impulse responses without constructing full G(t).

        This computes reduced kernels K_k = P_x U^k P_b^T for k=0..nt-1
        via RK4 time marching of impulse states, then assembles the
        lower-triangular block-Toeplitz matrix Ghat.

        Notes:
        - A is assumed to already include the physical time step scaling
            (same convention used in reduced_space_time_system_Ghat).
        - One macro time step is advanced per output sample; RK4 can be refined
            using rk4_substeps > 1.
        - integration_scheme options: 'rk4' (default), 'midpoint', 'heun', 'euler'.
        - If auto_substeps=True, rk4_substeps is automatically increased for
            stability using a conservative operator-norm heuristic.

        Parameters:
        - A: Sparse system matrix (n_state x n_state).
        - Binv_sqrt: Sparse diagonal matrix (n_state x n_state) of inverse square roots of material properties.
        - nt: Number of time steps (samples) to compute.
        - p_x: Boolean mask (length n_state) for receiver DOFs.
        - p_b: Boolean mask (length n_state) for source DOFs.
        - use_real: If True, computes with real-valued states and operator (H = Binv_sqrt @ A @ Binv_sqrt). If False, uses complex-valued states
            and operator (H = -1j * Binv_sqrt @ A @ Binv_sqrt) to capture oscillatory behavior directly.
        - verbose: If True, prints progress and diagnostic information.
        - rk4_substeps: Number of RK4 sub-steps per macro time step (default 1). Higher values increase stability at the cost of more computations.
        - integration_scheme: Time integration method to use for stepping the impulse response. Options are 'rk4', 'midpoint', 'heun', and 'euler'. RK4 is the default and most accurate, while the others are simpler but less stable.
        - auto_substeps: If True, automatically increases rk4_substeps based on a heuristic stability criterion derived from the infinity norm of the operator H. This can help ensure stability for stiff systems without manual tuning.
        - max_auto_substeps: Maximum allowed value for rk4_substeps when auto_substeps is enabled, to prevent excessive computation. Default is 256, which should be sufficient for most cases.

        Returns:
        - Ghat: Sparse block-Toeplitz matrix of size (nt * n_rec) x (nt * n_src), where n_rec and n_src are the number of active receivers and sources defined by p_x and p_b, respectively.    
    """
    if nt < 1:
        raise ValueError(f"nt must be >= 1, got {nt}.")

    n_state = A.shape[0]
    if A.shape[0] != A.shape[1]:
        raise ValueError(f"A must be square. Got shape {A.shape}.")
    if Binv_sqrt.shape != A.shape:
        raise ValueError(
            f"Binv_sqrt shape {Binv_sqrt.shape} must match A shape {A.shape}."
        )

    idx_x = np.asarray(p_x).astype(bool)
    idx_b = np.asarray(p_b).astype(bool)
    if idx_x.ndim != 1 or idx_x.size != n_state:
        raise ValueError(
            f"p_x must be a 1D mask of length {n_state}. Got shape {idx_x.shape}."
        )
    if idx_b.ndim != 1 or idx_b.size != n_state:
        raise ValueError(
            f"p_b must be a 1D mask of length {n_state}. Got shape {idx_b.shape}."
        )

    src_idx = np.flatnonzero(idx_b)
    n_src = src_idx.size
    n_rec = int(idx_x.sum())

    if n_src == 0 or n_rec == 0:
        if verbose:
            print("No active sources or receivers in masks. Returning empty Ghat.")
        return sp.csc_matrix((nt * n_rec, nt * n_src))

    if rk4_substeps < 1 or int(rk4_substeps) != rk4_substeps:
        raise ValueError(f"rk4_substeps must be a positive integer. Got {rk4_substeps}.")
    rk4_substeps = int(rk4_substeps)
    if max_auto_substeps < 1 or int(max_auto_substeps) != max_auto_substeps:
        raise ValueError(f"max_auto_substeps must be a positive integer. Got {max_auto_substeps}.")
    max_auto_substeps = int(max_auto_substeps)

    scheme = str(integration_scheme).lower()
    valid_schemes = {"rk4", "midpoint", "heun", "euler"}
    if scheme not in valid_schemes:
        raise ValueError(
            f"Unknown integration_scheme '{integration_scheme}'. "
            f"Choose one of {sorted(valid_schemes)}."
        )

    if use_real:
        H = Binv_sqrt @ A @ Binv_sqrt
        state_dtype = np.result_type(H.dtype, float)
    else:
        H = -1j * (Binv_sqrt @ A @ Binv_sqrt)
        state_dtype = np.result_type(H.dtype, complex)

    # Conservative stability heuristic for explicit methods.
    if auto_substeps:
        abs_row_sum = np.asarray(np.abs(H).sum(axis=1)).ravel()
        H_inf = float(abs_row_sum.max()) if abs_row_sum.size > 0 else 0.0
        target_mu = {
            "rk4": 1.5,
            "midpoint": 0.6,
            "heun": 0.6,
            "euler": 0.3,
        }[scheme]
        required = int(np.ceil(H_inf / max(target_mu, 1e-12))) if H_inf > 0 else 1
        required = max(1, required)
        stabilized = min(max_auto_substeps, max(rk4_substeps, required))
        if stabilized > rk4_substeps and verbose:
            print(
                f"Auto-substeps: increased rk4_substeps {rk4_substeps} -> {stabilized} "
                f"(||H||_inf~{H_inf:.3e}, scheme={scheme})."
            )
        rk4_substeps = stabilized

    if verbose:
        print(f"--- Direct Impulse Response Assembly (Real Mode: {use_real}) ---")
        print(f"State size={n_state}, receivers={n_rec}, sources={n_src}, nt={nt}")
        print(f"RK4 substeps per sample={rk4_substeps}")
        print(f"Integration scheme={scheme}")

    # Sparse selector whose columns are canonical basis vectors at source DOFs.
    selector_data = np.ones(n_src, dtype=state_dtype)
    B_sel = sp.csc_matrix(
        (selector_data, (src_idx, np.arange(n_src))),
        shape=(n_state, n_src),
        dtype=state_dtype,
    )

    # Time-march all source impulses simultaneously: V(:, s) is state for source s.
    V = B_sel.toarray()
    h = 1.0 / rk4_substeps

    def _step(V_curr, dt_step):
        if scheme == "euler":
            return V_curr + dt_step * (H @ V_curr)

        if scheme == "midpoint":
            k1 = H @ V_curr
            k2 = H @ (V_curr + 0.5 * dt_step * k1)
            return V_curr + dt_step * k2

        if scheme == "heun":
            k1 = H @ V_curr
            pred = V_curr + dt_step * k1
            k2 = H @ pred
            return V_curr + 0.5 * dt_step * (k1 + k2)

        # default: rk4
        k1 = H @ V_curr
        k2 = H @ (V_curr + 0.5 * dt_step * k1)
        k3 = H @ (V_curr + 0.5 * dt_step * k2)
        k4 = H @ (V_curr + dt_step * k3)
        return V_curr + (dt_step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    reduced_G = []
    if verbose:
        iterator = tqdm(range(nt), desc="RK4 Impulse Marching", unit="step")
    else:
        iterator = range(nt)

    for _ in iterator:
        # K_k = P_x V_k (receiver rows extracted from current impulse states).
        block_k = V[idx_x, :]
        reduced_G.append(sp.csc_matrix(block_k))

        # Advance one macro-step unless we are at the last snapshot.
        if len(reduced_G) < nt:
            for _sub in range(rk4_substeps):
                V = _step(V, h)
            if not np.isfinite(V).all():
                raise FloatingPointError(
                    "Non-finite values encountered during time marching. "
                    "Try increasing rk4_substeps, keeping auto_substeps=True, "
                    "or switch back to reduced_space_time_system_Ghat (expm-based)."
                )

    blocks = []
    if verbose:
        row_iter = tqdm(range(nt), desc="Assembling Block-Toeplitz", unit="row")
    else:
        row_iter = range(nt)

    for j in row_iter:
        row = []
        for i in range(nt):
            row.append(reduced_G[j - i] if i <= j else None)
        blocks.append(row)

    return sp.block_array(blocks, format='csc')


def reduced_space_time_system_Ghat_no_block_gpu(
    A,
    Binv_sqrt,
    nt,
    p_x,
    p_b,
    n_jobs=-1,
    use_real=True,
    verbose=True,
    rk4_substeps=1,
    integration_scheme="rk4",
    auto_substeps=True,
    max_auto_substeps=256,
):
    """
    GPU version of reduced_space_time_system_Ghat_no_block.

    Computes reduced impulse-response kernels with explicit time marching on GPU,
    then assembles the lower-triangular block-Toeplitz Ghat on CPU.

    Note:
    - n_jobs is accepted for API compatibility with CPU versions but is not used.
        - If auto_substeps=True, rk4_substeps is automatically increased for
            stability using a conservative operator-norm heuristic.
    """
    if nt < 1:
        raise ValueError(f"nt must be >= 1, got {nt}.")

    n_state = A.shape[0]
    if A.shape[0] != A.shape[1]:
        raise ValueError(f"A must be square. Got shape {A.shape}.")
    if Binv_sqrt.shape != A.shape:
        raise ValueError(
            f"Binv_sqrt shape {Binv_sqrt.shape} must match A shape {A.shape}."
        )

    idx_x = np.asarray(p_x).astype(bool)
    idx_b = np.asarray(p_b).astype(bool)
    if idx_x.ndim != 1 or idx_x.size != n_state:
        raise ValueError(
            f"p_x must be a 1D mask of length {n_state}. Got shape {idx_x.shape}."
        )
    if idx_b.ndim != 1 or idx_b.size != n_state:
        raise ValueError(
            f"p_b must be a 1D mask of length {n_state}. Got shape {idx_b.shape}."
        )

    src_idx = np.flatnonzero(idx_b)
    rec_idx = np.flatnonzero(idx_x)
    n_src = src_idx.size
    n_rec = rec_idx.size

    if n_src == 0 or n_rec == 0:
        if verbose:
            print("No active sources or receivers in masks. Returning empty Ghat.")
        return sp.csc_matrix((nt * n_rec, nt * n_src))

    if rk4_substeps < 1 or int(rk4_substeps) != rk4_substeps:
        raise ValueError(f"rk4_substeps must be a positive integer. Got {rk4_substeps}.")
    rk4_substeps = int(rk4_substeps)
    if max_auto_substeps < 1 or int(max_auto_substeps) != max_auto_substeps:
        raise ValueError(f"max_auto_substeps must be a positive integer. Got {max_auto_substeps}.")
    max_auto_substeps = int(max_auto_substeps)

    scheme = str(integration_scheme).lower()
    valid_schemes = {"rk4", "midpoint", "heun", "euler"}
    if scheme not in valid_schemes:
        raise ValueError(
            f"Unknown integration_scheme '{integration_scheme}'. "
            f"Choose one of {sorted(valid_schemes)}."
        )

    # Move operators to GPU sparse matrices.
    if sp.issparse(A):
        A_gpu = cps.csr_matrix(A)
    elif isinstance(A, np.ndarray):
        A_gpu = cps.csr_matrix(A)
    else:
        A_gpu = A
        if not cps.issparse(A_gpu):
            A_gpu = cps.csr_matrix(A_gpu)

    if sp.issparse(Binv_sqrt):
        B_gpu = cps.csr_matrix(Binv_sqrt)
    elif isinstance(Binv_sqrt, np.ndarray):
        B_gpu = cps.csr_matrix(Binv_sqrt)
    else:
        B_gpu = Binv_sqrt
        if not cps.issparse(B_gpu):
            B_gpu = cps.csr_matrix(B_gpu)

    if use_real:
        H_gpu = B_gpu @ A_gpu @ B_gpu
        state_dtype = cp.result_type(H_gpu.dtype, cp.float64)
    else:
        H_gpu = -1j * (B_gpu @ A_gpu @ B_gpu)
        state_dtype = cp.result_type(H_gpu.dtype, cp.complex128)

    if auto_substeps:
        # Estimate ||H||_inf from CPU copy for robust, backend-agnostic behavior.
        H_cpu = H_gpu.get()
        abs_row_sum = np.asarray(np.abs(H_cpu).sum(axis=1)).ravel()
        H_inf = float(abs_row_sum.max()) if abs_row_sum.size > 0 else 0.0
        target_mu = {
            "rk4": 1.5,
            "midpoint": 0.6,
            "heun": 0.6,
            "euler": 0.3,
        }[scheme]
        required = int(np.ceil(H_inf / max(target_mu, 1e-12))) if H_inf > 0 else 1
        required = max(1, required)
        stabilized = min(max_auto_substeps, max(rk4_substeps, required))
        if stabilized > rk4_substeps and verbose:
            print(
                f"Auto-substeps: increased rk4_substeps {rk4_substeps} -> {stabilized} "
                f"(||H||_inf~{H_inf:.3e}, scheme={scheme})."
            )
        rk4_substeps = stabilized

    if verbose:
        print(f"--- GPU Direct Impulse Response Assembly (Real Mode: {use_real}) ---")
        print(f"State size={n_state}, receivers={n_rec}, sources={n_src}, nt={nt}")
        print(f"RK4 substeps per sample={rk4_substeps}")
        print(f"Integration scheme={scheme}")

    src_idx_gpu = cp.asarray(src_idx)
    rec_idx_gpu = cp.asarray(rec_idx)

    # Dense state matrix on GPU: each column is one source impulse state.
    V_gpu = cp.zeros((n_state, n_src), dtype=state_dtype)
    V_gpu[src_idx_gpu, cp.arange(n_src)] = 1
    h = 1.0 / rk4_substeps

    def _step_gpu(V_curr, dt_step):
        if scheme == "euler":
            return V_curr + dt_step * (H_gpu @ V_curr)

        if scheme == "midpoint":
            k1 = H_gpu @ V_curr
            k2 = H_gpu @ (V_curr + 0.5 * dt_step * k1)
            return V_curr + dt_step * k2

        if scheme == "heun":
            k1 = H_gpu @ V_curr
            pred = V_curr + dt_step * k1
            k2 = H_gpu @ pred
            return V_curr + 0.5 * dt_step * (k1 + k2)

        # default: rk4
        k1 = H_gpu @ V_curr
        k2 = H_gpu @ (V_curr + 0.5 * dt_step * k1)
        k3 = H_gpu @ (V_curr + 0.5 * dt_step * k2)
        k4 = H_gpu @ (V_curr + dt_step * k3)
        return V_curr + (dt_step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    reduced_G = []
    if verbose:
        iterator = tqdm(range(nt), desc="GPU Impulse Marching", unit="step")
    else:
        iterator = range(nt)

    for _ in iterator:
        block_k_gpu = V_gpu[rec_idx_gpu, :]
        block_k_cpu = cp.asnumpy(block_k_gpu)
        reduced_G.append(sp.csc_matrix(block_k_cpu))

        if len(reduced_G) < nt:
            for _sub in range(rk4_substeps):
                V_gpu = _step_gpu(V_gpu, h)
            if not bool(cp.isfinite(V_gpu).all()):
                raise FloatingPointError(
                    "Non-finite values encountered during GPU time marching. "
                    "Try increasing rk4_substeps, keeping auto_substeps=True, "
                    "or switch back to reduced_space_time_system_Ghat (expm-based)."
                )

    if verbose:
        row_iter = tqdm(range(nt), desc="Assembling Block-Toeplitz", unit="row")
    else:
        row_iter = range(nt)

    blocks = []
    for j in row_iter:
        row = []
        for i in range(nt):
            row.append(reduced_G[j - i] if i <= j else None)
        blocks.append(row)

    return sp.block_array(blocks, format='csc')
