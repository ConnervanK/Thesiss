import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.colors import Normalize
from pathlib import Path

# Try importing helper functions, or define them inline if helper.py is missing
try:
    from helper import horizontal_resolution, central_wavelength, fracture_depth as calc_min_far_field_depth
except ImportError:
    try:
        from .helper import horizontal_resolution, central_wavelength, fracture_depth as calc_min_far_field_depth
    except ImportError:
        def central_wavelength(f_central, epsilon_r):
            c = 299792458.0
            c_medium = c / np.sqrt(epsilon_r)
            return c_medium / f_central

        def horizontal_resolution(f_central, epsilon_r, z):
            lam = central_wavelength(f_central, epsilon_r)
            I = (z + lam / 4)**2
            val = I - z**2
            return np.sqrt(val) if val > 0 else 0.0

        def calc_min_far_field_depth(f_central, epsilon_r, source_receiver_distance):
            lam = central_wavelength(f_central, epsilon_r)
            return 2 * source_receiver_distance**2 / lam

def analyze_and_plot_model_geometry(
    length=1.0,
    depth=0.5,
    f_central=1.5e9,
    epsilon_r=6.0,
    fracture_depth_lambda=None,     # Depth as fraction of wavelength (distance from receiver)
    fracture_thickness_lambda=None, # Thickness as fraction of wavelength
    fracture_depth_m=0.3,           # Absolute depth (default if lambda not specified)
    fracture_thickness_m=0.01,      # Absolute thickness (default if lambda not specified)
    x_receiver=0.5,
    y_receiver=0.1,
    source_offset=0.01,             # Distance between source and receiver for far-field calc
    show_plot=True
):
    """
    Analyzes and plots model geometry, calculating resolution and far-field depth.
    Allows specifying fracture properties as fractions of wavelength.
    """
    
    # Calculate Wavelength
    lam = central_wavelength(f_central, epsilon_r)

    # Determine Fracture Properties
    if fracture_depth_lambda is not None:
        dist_from_receiver = fracture_depth_lambda * lam
    else:
        dist_from_receiver = fracture_depth_m

    if fracture_thickness_lambda is not None:
        thickness = fracture_thickness_lambda * lam
    else:
        thickness = fracture_thickness_m

    fracture_y0 = y_receiver + dist_from_receiver
    fracture_y1 = fracture_y0 + thickness
    
    # Calculate Metrics
    h_res = horizontal_resolution(f_central, epsilon_r, dist_from_receiver)
    min_depth_far_field = calc_min_far_field_depth(f_central, epsilon_r, source_receiver_distance=source_offset) 

    results = {
        'wavelength_m': lam,
        'fracture_depth_m': dist_from_receiver,
        'fracture_depth_lambda': dist_from_receiver / lam,
        'fracture_thickness_m': thickness,
        'fracture_thickness_lambda': thickness / lam,
        'horizontal_resolution_m': h_res,
        'min_far_field_depth_m': min_depth_far_field,
        'first_fresnel_radius_m': h_res, # Same as horizontal resolution
        'is_in_far_field': dist_from_receiver >= min_depth_far_field
    }

    if not show_plot:
        return results

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Draw Domain
    outer = Rectangle((0, 0), length, depth, fill=False, edgecolor='black', linewidth=2)
    ax.add_patch(outer)
    
    # Draw Upper Layer (for visualization, assuming y_receiver implies an interface)
    upper_layer_thickness = y_receiver
    # Draw background first
    # Draw "Air" zone
    air_zone = Rectangle((0, 0), length, upper_layer_thickness, facecolor='#cfe8ff', edgecolor='none', alpha=0.5)
    ax.add_patch(air_zone)
    ax.text(length/2, upper_layer_thickness/2, "Medium 1 (e.g. Air)", ha='center', va='center', fontsize=9, color='blue')
    
    # Draw Fracture
    fracture = Rectangle((0, fracture_y0), length, thickness,
                         facecolor='#ffcccb', edgecolor='red',
                         linewidth=1.5, alpha=0.9, label='Fracture')
    ax.add_patch(fracture)

    # Draw First Fresnel Zone on Fracture
    fresnel_x_left = x_receiver - h_res
    fresnel_x_right = x_receiver + h_res
    
    # Clip to domain
    f_x_min = max(0, fresnel_x_left)
    f_x_max = min(length, fresnel_x_right)
    
    if f_x_max > f_x_min:
        fresnel_zone = Rectangle((f_x_min, fracture_y0), f_x_max - f_x_min, thickness,
                                 facecolor='gold', edgecolor='darkorange', linewidth=2, alpha=1.0, 
                                 zorder=5, label='1st Fresnel Zone')
        ax.add_patch(fresnel_zone)

    # Draw Receiver
    ax.scatter(x_receiver, y_receiver, marker='v', s=150, color='black', zorder=10, label='Receiver')
    
    # Annotations
    ax.annotate(f'h = {dist_from_receiver:.3f} m\n({dist_from_receiver/lam:.2f} $\\lambda$)', 
                xy=(x_receiver, (y_receiver + fracture_y0)/2), 
                xytext=(x_receiver + 0.1, (y_receiver + fracture_y0)/2),
                arrowprops=dict(arrowstyle='<->', color='black'), va='center')

    # Draw Semicircles (Fresnel zones)
    theta = np.linspace(0, np.pi, 200)
    
    # First Fresnel radius (dist_from_receiver)
    r1 = dist_from_receiver
    arc1_x = x_receiver + r1 * np.cos(theta)
    arc1_y = y_receiver + r1 * np.sin(theta)
    ax.plot(arc1_x, arc1_y, color='blue', linestyle='--', linewidth=1.5, alpha=0.7, label=f'R = h ({r1:.3f} m)')

    # Second Fresnel radius (dist + lambda/4)
    r2 = dist_from_receiver + lam / 4
    arc2_x = x_receiver + r2 * np.cos(theta)
    arc2_y = y_receiver + r2 * np.sin(theta)
    ax.plot(arc2_x, arc2_y, color='navy', linestyle='-', linewidth=1.5, alpha=0.7, label=f'R = h + $\\lambda$/4 ({r2:.3f} m)')

    ax.set_title("Model Geometry & Fresnel Zone", fontsize=14)
    ax.set_xlabel("Length (m)")
    ax.set_ylabel("Depth (m)")
    ax.set_ylim(depth, 0) # Invert Y for depth
    ax.set_xlim(0, length)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='lower left')
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.show()

    return results


def generate_alternating_block_model(
    filename,
    length=1.0,
    depth=0.5,
    f_central=1.5e9,
    epsilon_r=6.0,
    fracture_depth_lambda=None,
    fracture_thickness_lambda=None,
    fracture_depth_m=0.3,
    fracture_thickness_m=0.01,
    block_width_lambda=0.1,  # Width of alternating blocks as fraction of wavelength
    percent_delta=5.0,       # Percentage change in PERMITTIVITY
    baseline_conductivity=0.01,
    x_receiver=0.5,
    y_receiver=0.1,
    src_pos=None,
    rx_pos=None,
    show_plot=True,
    # gprMax specific params
    dx_dy_dz=(0.002, 0.002, 0.002),
    time_window=10e-9,
    background_material='ice', # Name of background material in gprMax
    fracture_material_base='fracture_alt'
):
    """
    Generates a gprMax model with alternating PERMITTIVITY blocks in the fracture.
    Plots the model and writes the input file.
    """
    # 1. Resolve Geometry (Reuse logic)
    lam = central_wavelength(f_central, epsilon_r)
    
    if fracture_depth_lambda is not None:
        dist_from_receiver = fracture_depth_lambda * lam
    else:
        dist_from_receiver = fracture_depth_m

    if fracture_thickness_lambda is not None:
        thickness = fracture_thickness_lambda * lam
    else:
        thickness = fracture_thickness_m
        
    fracture_y0 = y_receiver + dist_from_receiver
    fracture_y1 = fracture_y0 + thickness
    
    # 2. Block generation logic
    # Calculate approximate block width from lambda
    approx_block_width = block_width_lambda * lam
    
    # Calculate number of blocks
    # Round to nearest even integer to ensure equal number of + and - blocks
    # This guarantees the average permittivity change across the entire length is 0.
    num_blocks_temp = length / approx_block_width
    num_blocks = int(round(num_blocks_temp))
    
    if num_blocks % 2 != 0:
        num_blocks += 1 # specific adjustment to ensure even number
        
    if num_blocks < 2:
        num_blocks = 2 # Minimum 2 blocks for alternation
        
    # Recalculate block width to fit exactly into length
    block_width = length / num_blocks
    
    # Generate Alternating Permittivity Values
    # Pattern: +delta, -delta, +delta... centering around epsilon_r
    
    blocks = []
    current_x = 0.0
    
    # Multipliers for permutation
    # +5% -> 1.05, -5% -> 0.95
    perm_factors = [1 + percent_delta/100.0, 1 - percent_delta/100.0]
    
    for i in range(num_blocks):
        x_start = current_x
        x_end = current_x + block_width
        
        # Ensure last block aligns perfectly with length (floating point safety)
        if i == num_blocks - 1:
            x_end = length
            
        epsilon_val = epsilon_r * perm_factors[i % 2]
        
        blocks.append({
            'x_start': x_start,
            'x_end': x_end,
            'epsilon': epsilon_val,
            'index': i
        })
        current_x = x_end

    # 3. Plotting
    if show_plot:
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Domain
        outer = Rectangle((0, 0), length, depth, fill=False, edgecolor='black', linewidth=2)
        ax.add_patch(outer)
        
        # Normalize color map based on epsilon values
        all_eps = [b['epsilon'] for b in blocks]
        vmin, vmax = min(all_eps), max(all_eps)
        # Handle case where delta is 0
        if vmin == vmax:
             vmin -= 0.1
             vmax += 0.1
             
        norm = Normalize(vmin=vmin, vmax=vmax)
        cmap = plt.cm.coolwarm
        
        for b in blocks:
            color = cmap(norm(b['epsilon']))
            rect = Rectangle((b['x_start'], fracture_y0), b['x_end'] - b['x_start'], thickness,
                             facecolor=color, edgecolor='none', alpha=1.0)
            ax.add_patch(rect)
            
        # Colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, orientation='horizontal', pad=0.15, fraction=0.05)
        cbar.set_label(f'Relative Permittivity ($\\epsilon_r$) [Base: {epsilon_r}]')
        
        # Receiver
        ax.scatter(x_receiver, y_receiver, marker='v', s=150, color='black', zorder=10, label='Receiver')
        
        # Upper layer (visual)
        upper_zone = Rectangle((0, 0), length, y_receiver, facecolor='#cfe8ff', edgecolor='none', alpha=0.3)
        ax.add_patch(upper_zone)

        # Draw Semicircles (Fresnel zones)
        theta = np.linspace(0, np.pi, 200)
        
        # First Fresnel radius (dist_from_receiver)
        r1 = dist_from_receiver
        arc1_x = x_receiver + r1 * np.cos(theta)
        arc1_y = y_receiver + r1 * np.sin(theta)
        ax.plot(arc1_x, arc1_y, color='blue', linestyle='--', linewidth=1.5, alpha=0.7, label=f'R = h ({r1:.3f} m)')

        # Second Fresnel radius (dist + lambda/4)
        r2 = dist_from_receiver + lam / 4
        arc2_x = x_receiver + r2 * np.cos(theta)
        arc2_y = y_receiver + r2 * np.sin(theta)
        ax.plot(arc2_x, arc2_y, color='navy', linestyle='-', linewidth=1.5, alpha=0.7, label=f'R = h + $\\lambda$/4 ({r2:.3f} m)')

        ax.set_title(f"Alternating Permittivity Model (Block Width: {block_width_lambda:.2f}$\\lambda$)", fontsize=14)
        ax.set_xlabel("Length (m)")
        ax.set_ylabel("Depth (m)")
        ax.set_ylim(depth, 0)
        ax.set_xlim(0, length)
        ax.set_aspect('equal')
        ax.grid(True, linestyle=':', alpha=0.5)
        
        plt.tight_layout()
        plt.show()

    # 4. Export gprMax Input File
    if src_pos is None:
        src_pos = (x_receiver, y_receiver - 0.05, 0.0) # slightly above rx?
    if rx_pos is None:
        rx_pos = (x_receiver + 0.01, y_receiver, 0.0)
        
    lines = []
    lines.append(f"#title: Alternating Permittivity Model")
    lines.append(f"#domain: {length} {depth} {dx_dy_dz[2]}")
    lines.append(f"#dx_dy_dz: {dx_dy_dz[0]} {dx_dy_dz[1]} {dx_dy_dz[2]}")
    lines.append(f"#time_window: {time_window}")
    lines.append("")
    # Define Materials
    lines.append(f"#material: {epsilon_r} {baseline_conductivity} 1 0 {background_material}")
    lines.append(f"#material: 1 0 1 0 air")
    lines.append("")
    
    # We define specific materials for +delta and -delta to save definitions
    eps_plus = epsilon_r * (1 + percent_delta/100.0)
    eps_minus = epsilon_r * (1 - percent_delta/100.0)
    
    # Only create two materials for alternating blocks
    mat_plus = f"{fracture_material_base}_plus"
    mat_minus = f"{fracture_material_base}_minus"
    
    lines.append(f"#material: {eps_plus:.4f} {baseline_conductivity} 1 0 {mat_plus}")
    lines.append(f"#material: {eps_minus:.4f} {baseline_conductivity} 1 0 {mat_minus}")
    lines.append("")
    
    # Define Geometry
    # 1. Background fills up to surface (y_receiver assumed surface or interface depth)
    # Actually, simplistic approach: Fill everything below y_receiver with background material.
    # Fill everything above y_receiver with air.
    # The user didn't specify strict layering, but usually y_receiver is near surface.
    
    # Assuming standard model:
    # Air: 0 to y_receiver
    # Background: y_receiver to depth
    
    # However, standard gprMax coordinates: Origin (0,0,0). Y increases up in 2D plot? 
    # Wait, usually plotting uses (0, height) at top for implementation, but gprMax uses (0,0,0) at bottom left.
    # The plot functions above used `ax.set_ylim(depth, 0)` -> implying 0 is top, depth is bottom in visual.
    # BUT gprMax uses Cartesian coordinates where Y=0 is bottom.
    
    # We must be careful!
    # In `plot_model_geometry`, `fracture_y0 = 0.4` and `depth` was confusingly handled.
    # Let's assume standard gprMax coordinate system for the input file: (0,0,0) is bottom-left.
    # Length is X, Depth is Y.
    # If `depth` = 0.5m (total height of domain), then Y goes from 0 to 0.5.
    # If plotting uses `ylim(depth, 0)`, then 0 is Top visually.
    # So `y_receiver = 0.1` means 0.1m from TOP in Plot, which is Y = 0.4 in gprMax?
    # Or is y_receiver calculated from bottom?
    
    # Let's check `plot_model_geometry` in old file:
    # `upper_layer_thickness = 0.1`
    # `fracture_y0 = 0.4` (for depth=0.5 presumably)
    # If `fracture_y0` is 0.4 and `upper_layer` is 0.1, and receiver at 0.1?
    # Actually, in the plot it creates rectangles.
    # `Rectangle((0, fracture_y0), ...)` -> Standard matplotlib uses bottom-left origin for Rect.
    # If `ylim(depth, 0)` flips Y-axis, then (0,0) is top-left visual. 
    # But Matplotlib coords don't change, just the view.
    # So `Rectangle((0, 0), ...)` is at "Top" visual (0).
    # If `fracture_y0` = 0.4, it is 0.4 down from top (0).
    
    # HOWEVER, gprMax uses Y=0 at bottom.
    # We should convert "Visual Depth" to "gprMax Y".
    # gprMax Y = Domain Depth - Visual Depth.
    
    # Let's confirm coordinate system preference.
    # Most GPR users think in Depth (0 at surface, increasing downwards).
    # gprMax input requires Y (0 at bottom).
    
    # The inputs `fracture_depth_m` (distance from receiver) imply geometry relative to receiver.
    # If receiver is at Y_visual = 0.1 (0.1m depth), and fracture is at 0.3m distance -> Y_visual = 0.4.
    # Total Depth = 0.5.
    # In gprMax coords (Height=0.5):
    # Surface is at Y=0.4 (If 0.1m air layer).
    # Receiver is at Y=0.4.
    # Fracture is at Y = 0.4 - 0.3 = 0.1.
    
    # I will implement this conversion for the gprMax file output transparently.
    # Visualization uses "Depth" (0 at top).
    # gprMax file uses "Height" (0 at bottom).
    
    # Conversion:
    # gprMax_Y = depth - Visual_Y
    
    gpr_ysurf = depth - y_receiver # Interface height
    gpr_yfrac_top = depth - fracture_y0 # Fracture top height
    gpr_yfrac_bot = depth - fracture_y1 # Fracture bottom height
    
    # Ensure background filling
    # Box 1: Air (from gpr_ysurf to depth)
    lines.append(f"#box: 0 {gpr_ysurf:.6f} 0 {length} {depth} {dx_dy_dz[2]} air")
    
    # Box 2: Background (from 0 to gpr_ysurf)
    lines.append(f"#box: 0 0 0 {length} {gpr_ysurf:.6f} {dx_dy_dz[2]} {background_material}")
    
    lines.append("")
    lines.append(f"// Fracture Blocks ({num_blocks} blocks)")
    
    for b in blocks:
        mat_suffix = "plus" if b['index'] % 2 == 0 else "minus"
        mat_name = f"{fracture_material_base}_{mat_suffix}"
        # Note: gprMax needs min_y, max_y. 
        # Since gpr_yfrac_bot is lower than gpr_yfrac_top, we order them correctly.
        y_min = min(gpr_yfrac_top, gpr_yfrac_bot)
        y_max = max(gpr_yfrac_top, gpr_yfrac_bot)
        
        lines.append(f"#box: {b['x_start']:.6f} {y_min:.6f} 0 {b['x_end']:.6f} {y_max:.6f} {dx_dy_dz[2]} {mat_name}")

    lines.append("")
    # Source/Receiver conversion
    gpr_src_y = depth - src_pos[1] if src_pos else gpr_ysurf + 0.01 # Slightly above interface?
    gpr_rx_y = depth - rx_pos[1] if rx_pos else gpr_ysurf + 0.01
    
    # If implicit positions, use receiver location as reference
    if src_pos is None:
        gpr_src_y = depth - y_receiver # At interface
    if rx_pos is None:
        gpr_rx_y = depth - y_receiver # At interface

    # Default logic: Receiver at specified (x, y) visual.
    # gprMax Y = depth - y_receiver
    
    final_src_y = depth - y_receiver # Collocated for simplicity if not separated
    final_rx_y = depth - y_receiver
    
    # But usually source/rx are slightly offset or at specific coords.
    # I'll use the passed x/y_receiver directly converted for gprMax Y.
    
    lines.append("#waveform: ricker 1 1.5e9 my_ricker")
    lines.append(f"#hertzian_dipole: z {x_receiver} {final_src_y:.4f} 0 my_ricker")
    lines.append(f"#rx: {x_receiver + 0.02} {final_rx_y:.4f} 0") 
    lines.append("#messages: y")
    lines.append(f"#geometry_view: 0 0 0 {length} {depth} {dx_dy_dz[2]} {dx_dy_dz[0]} {dx_dy_dz[1]} {dx_dy_dz[2]} Model_Geometry n")

    # Write file
    out_path = Path(filename)
    if not out_path.parent.exists():
        out_path.parent.mkdir(parents=True, exist_ok=True)
        
    with open(out_path, 'w') as f:
        f.write("\n".join(lines))
        
    return {
        'filename': str(out_path),
        'block_width_m': block_width,
        'epsilon_plus': eps_plus,
        'epsilon_minus': eps_minus,
        'fracture_depth_m': dist_from_receiver,
        'fracture_thickness_m': thickness,
        'gprMax_y_fracture_center': (gpr_yfrac_top + gpr_yfrac_bot)/2
    }
