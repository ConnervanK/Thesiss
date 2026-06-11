import h5py
import numpy as np
import os
import glob
from scipy.interpolate import interp1d

def prepare_fwi_data(data_path="RTM/subwavelength_scatters*.out", num_shots=20, num_receivers=20, dt_new=None):
    """
    Reads GprMax output files, returns data in shape (num_shots, num_receivers, num_time_steps).
    Resamples time to dt_new if dt_new is specified (in nanoseconds).
    """
    out_files = []
    # Collect files 1 to num_shots
    for i in range(1, num_shots + 1):
        file_name = f"RTM/subwavelength_scatters{i}.out"
        if os.path.exists(file_name):
            out_files.append(file_name)
    
    if len(out_files) != num_shots:
        raise ValueError(f"Found {len(out_files)} files, expected {num_shots}")

    # Read first file to get dimensions
    with h5py.File(out_files[0], 'r') as f:
        # standard GPRMax .out file format 
        # normally inside group 'rxs' -> 'rx1' -> 'Ez'
        rx1_ez = f['rxs']['rx1']['Ez'][:]
        nt_original = len(rx1_ez)
        dt_original = f.attrs['dt'] * 1e9 # Convert to nanoseconds
        
    t_original = np.arange(nt_original) * dt_original
    
    # Devito often automatically computes its own internal dt to be CFL compliant.
    # We can pass dt_new to resample to typical numerical values. 
    # Or just return raw and handle in the Devito script.
    if dt_new is None:
        t_new = t_original
        nt_new = nt_original
    else:
        t_max = np.max(t_original)
        nt_new = int(t_max / dt_new) + 1
        t_new = np.linspace(0, t_max, nt_new)
        
    data = np.zeros((num_shots, num_receivers, nt_new))
    
    for shot_idx, fname in enumerate(out_files):
        with h5py.File(fname, 'r') as f:
            for rx_idx in range(1, num_receivers + 1):
                raw_ez = f['rxs'][f'rx{rx_idx}']['Ez'][:]
                if dt_new is not None:
                    interp_func = interp1d(t_original, raw_ez, kind='cubic', fill_value='extrapolate')
                    data[shot_idx, rx_idx - 1, :] = interp_func(t_new)
                else:
                    data[shot_idx, rx_idx - 1, :] = raw_ez
                    
    # The GprMax outputs units are typically Volts/meter. 
    # Return dt in ns, and the numpy data array
    return data, (dt_original if dt_new is None else dt_new), t_new

if __name__ == "__main__":
    d, dt, t = prepare_fwi_data()
    print(f"Data shape: {d.shape}, dt: {dt:.5f} ns, total time: {t[-1]:.2f} ns")
