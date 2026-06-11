import numpy as np
import matplotlib.pyplot as plt

# --- Time Integration Schemes ---
def step_euler(state, M, dt, Np=None):
    """Forward Euler (1st order)"""
    return state + dt * (M @ state)

def step_rk2(state, M, dt, Np=None):
    """Runge-Kutta 2nd order (Heun's method)"""
    k1 = M @ state
    state_pred = state + dt * k1
    k2 = M @ state_pred
    return state + (dt/2.0) * (k1 + k2)

def step_rk4(state, M, dt, Np=None):
    """Runge-Kutta 4th order"""
    k1 = M @ state
    k2 = M @ (state + 0.5 * dt * k1)
    k3 = M @ (state + 0.5 * dt * k2)
    k4 = M @ (state + dt * k3)
    return state + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

def step_symplectic_euler(state, M, dt, Np):
    """
    Symplectic Euler (1st order, Energy Conserving).
    Updates V then P. Better for long-term energy stability than standard Euler.
    Requires 'Np' (number of pressure points).
    """
    if Np is None:
        raise ValueError("step_symplectic_euler requires 'Np' argument")
        
    # 1. Update Velocity (V) using current Pressure (P)
    # dV/dt is in the second half of the result vector (indices Np:)
    k1 = M @ state
    state[Np:] += dt * k1[Np:]
    
    # 2. Update Pressure (P) using new Velocity (V)
    # dP/dt is in the first half of the result vector (indices :Np)
    k2 = M @ state
    state[:Np] += dt * k2[:Np]
    
    return state

def step_velocity_verlet(state, M, dt, Np):
    """
    Velocity Verlet (2nd order, Symplectic).
    Excellent energy conservation and stability for wave equations.
    Requires 'Np' (number of pressure points).
    """
    if Np is None:
        raise ValueError("step_velocity_verlet requires 'Np' argument")

    # 1. Half-step Velocity update
    k1 = M @ state
    state[Np:] += 0.5 * dt * k1[Np:]
    
    # 2. Full-step Pressure update (using V_half)
    k2 = M @ state
    state[:Np] += dt * k2[:Np]
    
    # 3. Second Half-step Velocity update (using P_new)
    k3 = M @ state
    state[Np:] += 0.5 * dt * k3[Np:]
    
    return state

# --- SOURCE FUNCTIONS ---
def ricker_wavelet(t, f_peak, t_delay):
    """
    Generates a Ricker wavelet (second derivative of a Gaussian).
    
    Parameters:
    - t: Time array or scalar.
    - f_peak: Peak frequency (Hz). Controls the width of the pulse.
              Ensure f_peak <= c / (10 * dx) to avoid dispersion.
    - t_delay: Time shift to center the wavelet (s).
    """
    arg = (np.pi * f_peak * (t - t_delay)) ** 2
    return (1.0 - 2.0 * arg) * np.exp(-arg)

# --- PLOTTING ---

def plot_model(kappa_p, rho_v, N, d, source_pos=None, receiver_indices=None):
    """
    Visualizes the material properties (Bulk Modulus and Density) of the 2D domain,
    overlaid with source and receiver positions.
    """
    if len(N) != 2:
        print("plot_model currently only supports 2D grids.")
        return

    Nx, Ny = N[0], N[1]
    dx, dy = d[0], d[1]
    extent = [0, Nx*dx, 0, Ny*dy]
    
    # Prepare Kappa Map (Pressure Grid)
    # kappa_p is flat, reshape to (Ny, Nx)
    kappa_map = kappa_p.reshape((Ny, Nx))
    
    # Prepare Density Map (Velocity Grid)
    # rho_v might be a list of arrays [rho_vx, rho_vy] or a single flat array
    # We'll visualize one component (e.g., rho_vx) or an average if possible, 
    # but usually plotting rho_vx (staggered in X) is sufficient representation.
    if isinstance(rho_v, list):
        rho_map_flat = rho_v[0] # Take rho_vx
        # rho_vx shape is usually (Ny, Nx+1), let's crop to (Ny, Nx) for easy visualization matching P grid
        rho_map = rho_map_flat.reshape((Ny, Nx+1))[:, :-1] 
    else:
        # Assume it's somehow already compatible or we just take what we get
        rho_map = rho_v.reshape((Ny, Nx))

    fig, ax = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: Bulk Modulus
    im1 = ax[0].imshow(kappa_map, origin='lower', extent=extent, cmap='viridis')
    ax[0].set_title("Bulk Modulus ($\kappa$)")
    ax[0].set_xlabel("x (m)")
    ax[0].set_ylabel("y (m)")
    fig.colorbar(im1, ax=ax[0], label="Pa")
    
    # Plot 2: Density
    im2 = ax[1].imshow(rho_map, origin='lower', extent=extent, cmap='plasma')
    ax[1].set_title("Density ($\\rho$)")
    ax[1].set_xlabel("x (m)")
    ax[1].set_ylabel("y (m)")
    fig.colorbar(im2, ax=ax[1], label="kg/m$^3$")

    # Overlay Source and Receivers on both plots
    for axis in ax:
        # Receivers
        if receiver_indices:
            # Convert indices to physical coordinates
            rx_x = []
            rx_y = []
            Np = Nx * Ny
            for idx in receiver_indices:
                 if idx < Np: # Pressure index
                    iy = idx // Nx
                    ix = idx % Nx
                    rx_x.append(ix * dx)
                    rx_y.append(iy * dy)
            
            if rx_x:
                axis.scatter(rx_x, rx_y, c='white', marker='o', s=20, edgecolors='black', linewidth=0.5, label='Receivers')

        # Source
        if source_pos:
            axis.scatter(source_pos[0], source_pos[1], c='red', marker='*', s=150, edgecolors='white', linewidth=1.0, label='Source')
        
        axis.legend(loc='upper right')

    plt.tight_layout()
    plt.show()

def plot_results(wavefield_frames, seismogram, receiver_indices, dt, field_name="Pressure", N=None, d=None, snapshot_interval=10, source_pos=None):
    if not wavefield_frames:
        print("No frames to plot.")
        return

    # Find the time step where the receiver(s) record the maximum amplitude
    if seismogram.shape[1] > 0:
        t_max_step = np.argmax(np.max(np.abs(seismogram), axis=1))
        max_frame_idx = t_max_step // snapshot_interval
        max_frame_idx = min(max_frame_idx, len(wavefield_frames) - 1)
    else:
        # Fallback: use max amplitude of wavefield if no receivers
        max_amp = -1.0
        max_frame_idx = 0
        for i, frame in enumerate(wavefield_frames):
            current_max = np.max(np.abs(frame))
            if current_max > max_amp:
                max_amp = current_max
                max_frame_idx = i

    time_snapshot = max_frame_idx * snapshot_interval * dt
    
    # Calculate common scale for both plots
    vlim_snapshot = np.max(np.abs(wavefield_frames[max_frame_idx]))
    vlim_seis = np.max(np.abs(seismogram)) if seismogram.size > 0 else 0
    common_vlim = max(vlim_snapshot, vlim_seis)
    if common_vlim == 0: common_vlim = 1e-10
    
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # Plot 1: Wavefield Snapshot
    im = ax[0].imshow(wavefield_frames[max_frame_idx], cmap='RdBu', vmin=-common_vlim, vmax=common_vlim, origin='lower')
    ax[0].set_title(f"{field_name} Snapshot (t = {time_snapshot:.3f} s)")
    fig.colorbar(im, ax=ax[0], label=f"{field_name} Amplitude")
    
    # Add receiver dots
    if N is not None:
        Nx = N[0]
        Np = Nx * N[1]
        rx_xs = []
        rx_ys = []
        for idx in receiver_indices:
            # Only plot if it's a pressure index (since we plot pressure snapshot)
            if idx < Np: 
                iy = idx // Nx
                ix = idx % Nx
                rx_xs.append(ix)
                rx_ys.append(iy)
        
        if rx_xs:
            ax[0].scatter(rx_xs, rx_ys, c='green', marker='o', s=30, edgecolors='black', linewidth=0.5, label='Receiver')
            
        # Add Source dot if provided
        if source_pos is not None and d is not None:
             # Convert physical position (meters) to grid indices
             src_x_idx = source_pos[0] / d[0]
             src_y_idx = source_pos[1] / d[1]
             ax[0].scatter([src_x_idx], [src_y_idx], c='red', marker='*', s=100, edgecolors='black', linewidth=0.5, label='Source')

        if rx_xs or source_pos:
             ax[0].legend(loc='upper right')

    # Plot 2: Seismogram or Single Trace
    if len(receiver_indices) == 1:
        # 1D Plot: Quantity over time
        time_axis = np.arange(seismogram.shape[0]) * dt
        ax[1].plot(time_axis, seismogram[:, 0])
        
        title = f"{field_name} Receiver Signal"
        # Try to decode position if N and d are provided and it's a Pressure index
        if N is not None and d is not None:
            idx = receiver_indices[0]
            Nx = N[0]
            Np = Nx * N[1]
            if idx < Np: # Assuming Pressure index
                iy = idx // Nx
                ix = idx % Nx
                x_pos = ix * d[0]
                y_pos = iy * d[1]
                title += f"\n(x={x_pos:.1f}m, y={y_pos:.1f}m)"
        
        ax[1].set_title(title)
        ax[1].set_xlabel("Time (s)")
        ax[1].set_ylabel(f"{field_name} Amplitude")
        ax[1].set_ylim(-common_vlim, common_vlim)
        ax[1].grid(True)
    else:
        # 2D Plot: Seismogram
        # Y-axis is time (s), X-axis is receiver index.
        # extent = [x_min, x_max, y_min, y_max]
        extent = [0, len(receiver_indices), 0, seismogram.shape[0] * dt]
        
        im2 = ax[1].imshow(seismogram, aspect='auto', cmap='RdBu', origin='lower', extent=extent, vmin=-common_vlim, vmax=common_vlim)
        ax[1].set_title(f"{field_name} Seismogram")
        ax[1].set_xlabel("Receiver Index")
        ax[1].set_ylabel("Time (s)")
        fig.colorbar(im2, ax=ax[1], label=f"{field_name} Amplitude")

    plt.show()
