# Acoustic Wave Simulation - Thesis Code

This repository contains the Python implementation for the simulation of acoustic wave propagation using Hamiltonian formulation and Unitary operators, as described in the thesis.

## Project Structure

The code is organized as follows:

```
code/
├── acustic_wave/
│   ├── main.py                 # Main simulation script (CPU)
│   ├── freq_cpu.py             # Frequency domain solver (CPU)
│   ├── freq_gpu.py             # Frequency domain solver (GPU - requires CuPy)
│   ├── simulation_utils.py     # Helper functions for time-stepping and plotting
│   ├── assembly/
│   │   ├── assembly.py         # Matrix assembly logic (Hamiltonian, Material properties)
│   │   └── Px_Pb.py            # Projection operators (Source/Receiver)
│   └── visualization/
│       ├── generates_fig_report.py  # Script to generate thesis figures and performance tables
│       ├── test_material_visualization.py # Debug script for material properties
│       └── visual.py           # Visualization helpers
└── requirements.txt            # Python dependencies
```

## Installation

1.  **Environment Setup**: It is recommended to use a virtual environment.
    ```bash
    python3 -m venv venv
    source venv/bin/activate
    ```

2.  **Install Dependencies**:
    ```bash
    pip install -r requirements.txt
    ```

### Optional GPU Support
If you have an NVIDIA GPU and want to run the accelerated scripts (like `freq_gpu.py`), you must install `cupy`. The installation command depends on your CUDA version (e.g., for CUDA 12.x):
```bash
pip install cupy-cuda12x
```
Refer to the [CuPy installation guide](https://docs.cupy.dev/en/stable/install.html) for details.

## Usage

### Running the Main Simulation
To run the standard time-domain simulation (CPU):

1.  Navigate to the `acustic_wave` directory:
    ```bash
    cd acustic_wave
    ```
2.  Run the main script:
    ```bash
    python3 main.py
    ```

### Reproducing Thesis Figures
To generate the figures and performance tables used in the "Results" section of the thesis:

1.  Navigate to the visualization directory:
    ```bash
    cd acustic_wave/visualization
    ```
2.  Run the report generation script:
    ```bash
    python3 generates_fig_report.py
    ```
    *Note: The script will attempt to save figures to `latency_doc/figures` if it exists. Otherwise, it will create a local `figures` directory.*

## Modules Description

*   **`main.py`**: Performs time evolution of the wave equation on a 2D grid using staggered grids.
*   **`assembly`**: Contains the core logic for building the Hamiltonian matrices ($\mathbf{H}$), spatial derivatives ($\mathbf{A}$), and material matrices ($\mathbf{B}$).
*   **`visualization`**: Scripts dedicated to producing publication-quality plots and verifying material geometry.
