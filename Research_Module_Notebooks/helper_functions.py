import numpy as np

def kirchhoff_migration(
    data,        # seismic data d(x,t) [nx, nt]
    x,           # surface positions [nx]
    t,           # time axis [nt]
    z,           # depth axis [nz]
    v             # constant velocity
):
    """
    2D Zero-offset Kirchhoff Time Migration (constant velocity)

    Parameters
    ----------
    data : ndarray (nx, nt)
        Seismic section
    x : ndarray (nx,)
        Surface positions
    t : ndarray (nt,)
        Time samples
    z : ndarray (nz,)
        Depth samples
    v : float
        Migration velocity

    Returns
    -------
    image : ndarray (nx, nz)
        Migrated image
    """

    nt, nx = data.shape
    nz = len(z)

    image = np.zeros((nz, nx))

    dt = t[1] - t[0]

    for ix0 in range(nx):          # image location x0
        for iz in range(nz):       # depth
            z0 = z[iz]

            for ix in range(nx):   # summation over data traces
                dx = x[ix] - x[ix0]

                # two-way travel time
                tau = 2.0 * np.sqrt(z0**2 + dx**2) / v

                it = int(tau / dt)

                if 0 <= it < nt:
                    image[iz, ix0] += data[it, ix]

                    

    return image



def rpca_ialm(D, lam=None, mu=None, tol=1e-7, max_iter=100):
    """
    Robust PCA via Inexact Augmented Lagrange Multiplier (IALM)

    Parameters
    ----------
    D : ndarray (nt x nx)
        B-scan data (time x trace)
    lam : float
        Sparsity weight (default: 1/sqrt(max(nt, nx)))
    mu : float
        Augmented Lagrangian parameter
    tol : float
        Convergence tolerance
    max_iter : int
        Maximum iterations

    Returns
    -------
    L : ndarray
        Low-rank component (direct wave + clutter)
    S : ndarray
        Sparse component (reflections)
    """

    D = D.astype(float)
    m, n = D.shape

    if lam is None:
        lam = 1.0 / np.sqrt(max(m, n))

    norm_D = np.linalg.norm(D, ord='fro')

    # Initialize
    L = np.zeros_like(D)
    S = np.zeros_like(D)
    Y = np.zeros_like(D)

    if mu is None:
        mu = (m * n) / (4.0 * np.sum(np.abs(D)) + 1e-16)

    mu_inv = 1.0 / mu

    for it in range(max_iter):
        # --- Update L via singular value thresholding ---
        U, sigma, VT = np.linalg.svd(D - S + mu_inv * Y, full_matrices=False)

        U = np.asarray(U)
        sigma = np.asarray(sigma)
        VT = np.asarray(VT)

        sigma_thresh = np.maximum(sigma - mu_inv, 0)
        L = (U * sigma_thresh) @ VT

        # --- Update S via soft thresholding ---
        S = np.sign(D - L + mu_inv * Y) * np.maximum(
            np.abs(D - L + mu_inv * Y) - lam * mu_inv, 0
        )

        # --- Update Lagrange multiplier ---
        Z = D - L - S
        Y = Y + mu * Z

        # --- Check convergence ---
        
        err = np.linalg.norm(Z, ord='fro') / norm_D
        if err < tol:
            print(f"RPCA converged in {it} iterations.")
            break

    return L, S

import h5py
import json

def save_gpr_data_to_h5(gpr_data_obj, output_path, field_comp='Ez'):
    """
    Save GPRData object to HDF5 file in gprMax-compatible format.
    
    This creates the structure expected by ImpDAR's load_gprMax function:
    - /rxs/rx1/Ez: the data array
    - attrs['dt']: time step as attribute
    - Additional metadata stored as JSON
    
    Parameters:
    -----------
    gpr_data_obj : gdp.GPRData
        GPRData object to save
    output_path : str
        Path to save the HDF5 file
    field_comp : str
        Field component to save (e.g., 'Ez')
    """
    with h5py.File(output_path, 'w') as f:
        # Create the group structure expected by gprMax/ImpDAR
        rx_group = f.create_group('rxs/rx1')
        
        # Save the data array at the expected location
        rx_group.create_dataset(field_comp, data=gpr_data_obj.data, compression='gzip', compression_opts=6)
        
        # Save dt as a direct attribute (required by load_gprMax)
        dt_value = gpr_data_obj.info.get('dt (ns)', 0.0)
        if isinstance(dt_value, str):
            dt_value = float(dt_value)
        f.attrs['dt'] = dt_value
        
        # Save all metadata as individual attributes for easy access
        for key, value in gpr_data_obj.info.items():
            try:
                if isinstance(value, str):
                    f.attrs[key] = value
                else:
                    f.attrs[key] = float(value)
            except (TypeError, ValueError):
                # If conversion fails, store as string
                f.attrs[key] = str(value)
        
        # Also store metadata as JSON for Python dict reconstruction
        metadata_str = json.dumps({str(k): float(v) if isinstance(v, (int, float)) else str(v) 
                                   for k, v in gpr_data_obj.info.items()})
        f.attrs['metadata_json'] = metadata_str
    
    print(f"✓ Successfully saved processed GPR data to: {output_path}")
    print(f"  Data shape: {gpr_data_obj.data.shape}")
    print(f"  Data dtype: {gpr_data_obj.data.dtype}")
    print(f"  dt (ns): {gpr_data_obj.info.get('dt (ns)', 'N/A')}")
    print(f"  Field component: {field_comp}")


def load_gpr_data_from_h5(file_path, field_comp='Ez'):
    """
    Load GPRData object from HDF5 file.
    Compatible with both gprMax-format and stored metadata.
    
    Parameters:
    -----------
    file_path : str
        Path to the HDF5 file
    field_comp : str
        Field component to load (e.g., 'Ez')
        
    Returns:
    --------
    gdp.GPRData
        Loaded GPRData object
    """
    with h5py.File(file_path, 'r') as f:
        # Load data from the standard location
        if f'/rxs/rx1/{field_comp}' in f:
            data = f[f'/rxs/rx1/{field_comp}'][:]
        else:
            # Fallback for older format
            data = f['data'][:]
        
        # Load metadata
        if 'metadata_json' in f.attrs:
            info = json.loads(f.attrs['metadata_json'])
        else:
            # Extract metadata from attributes
            info = {}
            for key in f.attrs.keys():
                if key not in ['metadata_json']:
                    try:
                        info[key] = float(f.attrs[key])
                    except (TypeError, ValueError):
                        info[key] = f.attrs[key]
    
    # Create GPRData object
    gpr_data = gdp.GPRData(data, info)
    
    print(f"✓ Successfully loaded GPRData from: {file_path}")
    print(f"  Data shape: {gpr_data.data.shape}")
    print(f"  Data dtype: {gpr_data.data.dtype}")
    print(f"  Metadata: {dict(gpr_data.info)}")
    
    return gpr_data
