import numpy as np
import subprocess
from pathlib import Path
from gprMax.reader import read_output
import helper


WORKDIR = Path(__file__).resolve().parent
TEMP_INPUT = WORKDIR / "temp_run.in"
TEMP_OUTPUT = WORKDIR / "temp_run.out"


def build_input_file(fracture_epsilon):
    """Build a thin-fracture gprMax input file for the proposed permittivity."""
    helper.write_thin_fracture_gprmax_input(
        filename=str(TEMP_INPUT),
        fracture_epsilon_r=float(fracture_epsilon),
        background_epsilon_r=9.63,
        domain=(0.8, 0.6, 0.001),
        dx_dy_dz=(0.001, 0.001, 0.001),
        time_window=10e-9,
        center_frequency=1.5e9,
        source=(0.4, 0.55, 0.0),
        receiver=(0.41, 0.55, 0.0),
        src_steps=(0.0375, 0.0, 0.0),
        rx_steps=(0.0375, 0.0, 0.0),
        host_box=((0.0, 0.0, 0.0), (0.8, 0.5, 0.001)),
        fracture_y=(0.29375, 0.3),
        extend_fracture_into_pml=True,
        title="Thin Fracture Model",
        fracture_material_name="fracture",
    )

def run_gprmax(epsilon_val):
    """Builds input, runs gprMax, and returns the A-scan."""
    build_input_file(epsilon_val)

    # Execute gprMax (using -n 1 for a single A-scan)
    result = subprocess.run(
        ["python", "-m", "gprMax", str(TEMP_INPUT), "-n", "1"],
        capture_output=True,
        text=True,
        cwd=WORKDIR,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "gprMax simulation failed.\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    # Parse output (.out file is HDF5) using gprMax's built-in reader.
    outputdata, dt = read_output(str(TEMP_OUTPUT), 1, "Ez")
    return np.array(outputdata)

def log_likelihood(simulated, observed, sigma=0.01):
    """Calculates how well the simulation matches observations."""
    # Assuming Gaussian noise
    ssr = np.sum((simulated - observed)**2)
    return -0.5 * ssr / (sigma**2)

# --- Simple Metropolis-Hastings Loop ---
current_eps = 4.0  # Initial guess
observed_data = np.load(WORKDIR / "real_scan.npy")
trace = []

for i in range(100):
    # Propose new value from a normal distribution
    proposal = current_eps + np.random.normal(0, 0.2)
    
    # Constraints: permittivity can't be < 1
    if proposal < 1: 
        continue

    sim_data = run_gprmax(proposal)
    
    log_p_current = log_likelihood(run_gprmax(current_eps), observed_data)
    log_p_proposal = log_likelihood(sim_data, observed_data)
    
    # Acceptance ratio
    if np.log(np.random.rand()) < (log_p_proposal - log_p_current):
        current_eps = proposal
    
    trace.append(current_eps)
    print(f"Iteration {i}: Epsilon = {current_eps}")