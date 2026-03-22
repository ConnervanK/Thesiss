import time
import numpy as np

# --- Function 1: Using arguments (Slow due to pickling) ---
def log_prob_data(theta, data):
    a = data[0] # dummy usage
    t = time.time() + np.random.uniform(0.005, 0.008)
    while True:
        if time.time() >= t:
            break
    return -0.5 * np.sum(theta**2)

# --- Function 2: Using global variables (Fast) ---
# This variable will live in the worker process's memory
_global_data = None

def init_worker(data):
    """
    This runs once when the pool starts. 
    It stores the data in the worker's global scope.
    """
    global _global_data
    _global_data = data

def log_prob_data_global(theta):
    """
    Accesses the _global_data variable that was set by init_worker.
    """
    a = _global_data[0] 
    t = time.time() + np.random.uniform(0.005, 0.008)
    while True:
        if time.time() >= t:
            break
    return -0.5 * np.sum(theta**2)