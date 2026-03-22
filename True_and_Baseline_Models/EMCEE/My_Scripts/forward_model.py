import numpy as np
import os
import h5py
import subprocess
from scipy.signal import hilbert

# Global constants needed for simulation
tw = 7.96e-9
pec_Y = 0.04

def model_2L(p1, c1, p2, c2, sd2):
    sd1 = 0.15
    gpr_Y = pec_Y + sd1 + sd2
    
    # Generate unique filenames based on the Process ID
    pid = os.getpid()
    unique_input = f"Simulation_{pid}.in"
    unique_output = f"Simulation_{pid}.out"
    unique_geom = f"Sim_{pid}"

    # Write to unique input file
    with open(unique_input, "w") as f:
        f.write("#title: 15 cm soil + PEC + top layer of depth d2")
        f.write("\n#domain: 0.30 0.50 0.001")
        f.write("\n#dx_dy_dz: 0.001 0.001 0.001")  
        f.write("\n#time_window: {}".format(tw))
        f.write("\n#material: {} {} 1 0 mat1 \n#material: {} {} 1 0 mat2".format(p1, c1, p2, c2))
        f.write("\n#waveform: gaussian 1 {} my_wave".format(1.5778e9))
        f.write("\n#hertzian_dipole: z 0.12 {} 0 my_wave\n#rx: 0.18 {} 0".format(gpr_Y, gpr_Y))
        f.write("\n#box: 0 {} 0 0.3 {} 0.001 mat1".format(pec_Y, pec_Y + sd1))
        f.write("\n#box: 0 {} 0 0.3 {} 0.001 mat2".format(pec_Y+sd1, pec_Y+sd1+sd2))
        f.write("\n#box: 0 {} 0 0.30 {} 0.001 pec".format(pec_Y-0.005, pec_Y))
        # Update specific lines to use unique names if necessary, e.g. geometry view
        # It is HIGHLY recommended to remove #geometry_view for MCMC to save speed/disk
        f.close()
    
    # Run gprMax with unique input
    env = os.environ.copy()
    env['CUDA_VISIBLE_DEVICES'] = ''
    try:
        subprocess.run([
            'python', '-c',
            f'from gprMax.gprMax import api; api("{unique_input}", geometry_only=False)'
        ], env=env, check=True)
        
        # Read from unique output file
        hf2 = h5py.File(unique_output, 'r', libver='latest', swmr=True)
        Ascan = hf2['rxs/rx1/Ez'][:]
        hf2.close()
    except Exception as e:
        print(f"Simulation failed: {e}")
        # Return a high error penalty or handle gracefully
        # For now, return zeros or re-raise. Since this is MCMC, failures are bad.
        # But we must clean up.
        Ascan = None
        raise e
    finally:
        # Cleanup: Delete the temporary files to save space
        if os.path.exists(unique_input):
            os.remove(unique_input)
        if os.path.exists(unique_output):
            os.remove(unique_output)
            
    return Ascan

def mae_func(p1, c1, p2, c2, sd2, E_new_normalized):
    vector1 = E_new_normalized                  # Experimental signal (normalized amplitude envelope)
    vector2 = model_2L(p1, c1, p2, c2, sd2)     # Run gprMax Simulation
    vector2 = vector2/(np.max(abs(vector2)))    # Normalizing simulation signal
    vector2 = np.abs(hilbert(vector2))          # Simulation signal (normalized amplitude envelope)
    
    error = vector1-vector2     # compute error
    sigma = np.std(error)       # compute the standard deviation of the residual
    result = np.sum(0.5*(error/sigma)**2 + np.log(np.sqrt(2*3.1416)*sigma)) #calculate the log likelihood
    
    return -result # We want to maximize the negative log likelihood

# Define the log-posterior distribution function
def log_posterior(params, E_new_normalized):
    p1, c1, p2, c2, sd2 = params
    
    # Check if parameters are within bounds
    if (p1 < 1 or p1 > 23 or c1 < 0 or c1 > 0.1 or p2 < 1 or p2 > 9 or c2 < 0 or c2 > 0.1 or sd2 <0.02 or sd2 > 0.18):
        return -np.inf # Return -inf if the parameters are outside the bounds
    else:
        # Evaluate the log-posterior distribution as the negative of the objective function
        log_likelihood = mae_func(p1, c1, p2, c2, sd2, E_new_normalized)
        return log_likelihood       # No prior added as unifrom prior is assumed
