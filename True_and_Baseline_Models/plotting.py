import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

def plot_model_geometry(length, depth, f_central, epsilon_r, fracture_thickness, x_receiver=1.0, y_receiver=0.2):
    upper_layer_thickness = 0.1  # m
    fracture_y0 = 0.4  # m
    fracture_y1 = fracture_y0 + fracture_thickness  # m
    pml_thickness = 0.02  # m

    receiver_x = x_receiver  # m
    receiver_y = y_receiver  # m

    f0 = f_central
    c = 299792458.0  # m/s
    lam = (c / np.sqrt(epsilon_r)) / f0
    radius_1 = 0.4
    radius_2 = radius_1 + lam / 4.0

    theta = np.linspace(0, np.pi, 400)
    arc1_x = receiver_x + radius_1 * np.cos(theta)
    arc1_y = receiver_y + radius_1 * np.sin(theta)
    arc2_x = receiver_x + radius_2 * np.cos(theta)
    arc2_y = receiver_y + radius_2 * np.sin(theta)

    dy_frac = fracture_y0 - receiver_y
    delta_x = np.sqrt(max(radius_2**2 - dy_frac**2, 0.0))
    fresnel_x_min = max(0.0, receiver_x - delta_x)
    fresnel_x_max = min(length, receiver_x + delta_x)

    fig, ax = plt.subplots(figsize=(10, 5))

    outer = Rectangle((0, 0), length, depth, fill=False, edgecolor='black', linewidth=2)
    ax.add_patch(outer)

    upper = Rectangle((0, 0), length, upper_layer_thickness, facecolor='#cfe8ff', edgecolor='none', alpha=0.8)
    ax.add_patch(upper)

    fracture_color = '#ffcccb' if fracture_thickness < 0.05 else '#ff9999'
    fracture_edge_color = 'red' if fracture_thickness < 0.05 else 'darkred'
    
    fracture = Rectangle((0, fracture_y0), length, fracture_thickness,
                         facecolor=fracture_color, edgecolor=fracture_edge_color,
                         linewidth=1.5, alpha=0.9)
    ax.add_patch(fracture)

    fresnel_segment = Rectangle((fresnel_x_min, fracture_y0), fresnel_x_max - fresnel_x_min, fracture_thickness,
                                facecolor='gold', edgecolor='darkorange', linewidth=1.5, alpha=0.9,
                                zorder=5, label='Fracture within 1st Fresnel zone')
    ax.add_patch(fresnel_segment)
    pml_style = dict(facecolor='#dddddd', edgecolor='gray', alpha=0.7, hatch='///')
    ax.add_patch(Rectangle((0, 0), pml_thickness, depth, **pml_style))
    ax.add_patch(Rectangle((length - pml_thickness, 0), pml_thickness, depth, **pml_style))
    ax.add_patch(Rectangle((0, 0), length, pml_thickness, **pml_style))
    ax.add_patch(Rectangle((0, depth - pml_thickness), length, pml_thickness, **pml_style))
    ax.text(x_receiver, upper_layer_thickness / 2, 'Air',
            ha='center', va='center', fontsize=10, color='black')
    ax.text(x_receiver, fracture_y0 + 0.1, f'Fracture (thickness={fracture_thickness:.3f} m)',
            ha='center', va='center', fontsize=10, color=fracture_edge_color)
    ax.text(length - pml_thickness / 2, depth / 2, 'PML', ha='center', va='center', fontsize=16, color='gray', rotation=90)
    ax.text(pml_thickness / 2, depth / 2, 'PML', ha='center', va='center', fontsize=16, color='gray', rotation=90)
    ax.scatter(receiver_x, receiver_y, marker='v', s=140, color='black', zorder=7, label='Receiver')
    ax.plot(arc1_x, arc1_y, color='blue', linewidth=2.0, label='r = 1.0 m')
    ax.plot(arc2_x, arc2_y, color='navy', linestyle='--', linewidth=2.0,
            label=f'r = 1 + $\\lambda/4$')
    arrow1 = FancyArrowPatch((receiver_x, receiver_y), (receiver_x, fracture_y0),
                            arrowstyle='->', mutation_scale=20, color='black', linewidth=1.5)   
    ax.add_patch(arrow1)
    ax.text(receiver_x - 0.12, (receiver_y + fracture_y0) / 2, 'h', fontsize=14, fontweight='bold', color='black')
    dy_frac_from_receiver = fracture_y0 - receiver_y
    delta_x_2 = np.sqrt(max(radius_2**2 - dy_frac_from_receiver**2, 0.0))
    x_left = receiver_x - delta_x_2
    x_right = receiver_x + delta_x_2
    arrow2_left = FancyArrowPatch((receiver_x, receiver_y), (x_left, fracture_y0),
                                 arrowstyle='->', mutation_scale=20, color='navy', linewidth=1.5, linestyle='--')
    ax.add_patch(arrow2_left)
    arrow2_right = FancyArrowPatch((receiver_x, receiver_y), (x_right, fracture_y0),
                                   arrowstyle='->', mutation_scale=20, color='navy', linewidth=1.5, linestyle='--')
    ax.add_patch(arrow2_right)
    ax.text((receiver_x + x_left) / 2 - 0.18, (receiver_y + fracture_y0) / 2 - 0.08, 'h + $\\lambda$/4', fontsize=11, color='navy', fontweight='bold')
    ax.text((receiver_x + x_right) / 2 + 0.08, (receiver_y + fracture_y0) / 2 - 0.08, 'h + $\\lambda$/4', fontsize=11, color='navy', fontweight='bold')
    ax.set_title('Model Geometry and 1st Fresnel Zone on Fracture', fontsize=16)
    ax.set_xlabel('Length (m)', fontsize=16)
    ax.set_ylabel('Depth (m)', fontsize=16)
    ax.set_xlim(0, length)
    ax.set_ylim(depth, 0)  # invert to show depth downward
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), borderaxespad=0.0, fontsize=11)
    fig.subplots_adjust(right=0.78)
    plt.show()
