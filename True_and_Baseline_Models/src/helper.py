import re
import numpy as np

def central_wavelength(f_central, epsilon_r):
    c = 299792458.0  # m/s
    c_medium = c / np.sqrt(epsilon_r)  # m/s
    return c_medium / f_central  # m

def max_fracture_thickness(f_central, epsilon_r, sub_wavelength_factor=4):
    lam = central_wavelength(f_central, epsilon_r)
    thickness = lam / sub_wavelength_factor  # m
    return thickness

def fracture_depth(f_central, epsilon_r, source_receiver_distance):
    lam = central_wavelength(f_central, epsilon_r)
    return 2 * source_receiver_distance**2 / lam # Minimum distance to get to the far field in metre

def horizontal_resolution(f_central, epsilon_r, z):
    lam = central_wavelength(f_central, epsilon_r)
    I = (z + lam / 4)**2
    return np.sqrt( I - z**2 ) # radius in metres
    



def write_thin_fracture_gprmax_input(
    filename,
    fracture_epsilon_r,
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
):
    """Write a gprMax input file for a thin-fracture model.

    The fracture can span the full domain in x so it extends into side PML regions.
    """
    if fracture_epsilon_r < 1.0:
        raise ValueError("fracture_epsilon_r must be >= 1.0")

    x_len, y_len, z_len = domain
    (host_min, host_max) = host_box
    frac_y0, frac_y1 = fracture_y

    if frac_y1 <= frac_y0:
        raise ValueError("fracture_y upper bound must be greater than lower bound")

    if extend_fracture_into_pml:
        frac_x0 = 0.0
        frac_x1 = x_len
    else:
        frac_x0 = host_min[0]
        frac_x1 = host_max[0]

    lines = [
        f"#title: {title}",
        "#num_threads: 8",
        f"#domain: {x_len:.3f} {y_len:.3f} {z_len:.3f}",
        f"#dx_dy_dz: {dx_dy_dz[0]:.3f} {dx_dy_dz[1]:.3f} {dx_dy_dz[2]:.3f}",
        f"#time_window: {time_window:g}",
        "",
        f"#material: {background_epsilon_r:.6g} 0 1 0 marble",
        f"#material: {fracture_epsilon_r:.6g} 0 1 0 {fracture_material_name}",
        "",
        f"#waveform: ricker 1 {center_frequency:.6g} my_ricker",
        f"#hertzian_dipole: z {source[0]:.3f} {source[1]:.3f} {source[2]:.3f} my_ricker",
        f"#rx: {receiver[0]:.3f} {receiver[1]:.3f} {receiver[2]:.3f}",
        "",
        f"#src_steps: {src_steps[0]:.4f} {src_steps[1]:.4f} {src_steps[2]:.4f}",
        f"#rx_steps: {rx_steps[0]:.4f} {rx_steps[1]:.4f} {rx_steps[2]:.4f}",
        "",
        f"#box: {host_min[0]:.3f} {host_min[1]:.3f} {host_min[2]:.3f} "
        f"{host_max[0]:.3f} {host_max[1]:.3f} {host_max[2]:.3f} marble",
        f"#box: {frac_x0:.3f} {frac_y0:.5f} 0 {frac_x1:.3f} {frac_y1:.5f} {z_len:.3f} {fracture_material_name}",
        "",
        "#geometry_view: 0 0 0 "
        f"{x_len:.3f} {y_len:.3f} {z_len:.3f} "
        f"{dx_dy_dz[0]:.3f} {dx_dy_dz[1]:.3f} {dx_dy_dz[2]:.3f} Thin_fracture_model n",
        "#messages: y",
        "",
    ]

    with open(filename, "w", encoding="ascii") as f:
        f.write("\n".join(lines))


def update_fracture_permittivity_in_input(
    filename,
    fracture_epsilon_r,
    fracture_material_name="fracture",
):
    """Update only the fracture material permittivity in an existing gprMax file."""
    if fracture_epsilon_r < 1.0:
        raise ValueError("fracture_epsilon_r must be >= 1.0")

    with open(filename, "r", encoding="ascii") as f:
        content = f.read()

    pattern = (
        r"(^\s*#material:\s*)"
        r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
        r"(\s+0\s+1\s+0\s+"
        + re.escape(fracture_material_name)
        + r"\s*$)"
    )

    replacement = r"\g<1>" + f"{fracture_epsilon_r:.6g}" + r"\g<3>"
    new_content, n_subs = re.subn(pattern, replacement, content, flags=re.MULTILINE)

    if n_subs == 0:
        raise ValueError(
            "Could not find fracture material line. "
            f"Expected '#material: <eps> 0 1 0 {fracture_material_name}'."
        )

    with open(filename, "w", encoding="ascii") as f:
        f.write(new_content)


def adjust_fracture_aperture(filename, new_aperture, fixed_reference='top', material_name='water'):
    """
    Adjusts the aperture of the thin fracture (box with specific material) in the gprMax input file.
    
    Args:
        filename (str): Path to the .in file.
        new_aperture (float): The new aperture (thickness) in meters.
        fixed_reference (str): 'top', 'bottom', or 'center'. 
                               'top' (default) keeps y_max constant.
                               'bottom' keeps y_min constant.
                               'center' keeps the center constant.
        material_name (str): The material name to look for (default: 'water').
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    new_lines = []
    fracture_found = False
    
    # Regex to match the box line with specific material
    # Expecting: #box: x_min y_min z_min x_max y_max z_max material
    pattern_str = r'^#box:\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+' + re.escape(material_name)
    box_pattern = re.compile(pattern_str)
    
    for line in lines:
        match = box_pattern.match(line.strip())
        if match:
            x_min, y_min_str, z_min, x_max, y_max_str, z_max = match.groups()
            y_min = float(y_min_str)
            y_max = float(y_max_str)
            
            # Decide new coordinates
            if fixed_reference == 'top':
                new_y_max = y_max
                new_y_min = y_max - new_aperture
            elif fixed_reference == 'bottom':
                new_y_min = y_min
                new_y_max = y_min + new_aperture
            elif fixed_reference == 'center':
                center = (y_min + y_max) / 2.0
                new_y_min = center - new_aperture / 2.0
                new_y_max = center + new_aperture / 2.0
            else:
                 raise ValueError(f"Unknown fixed_reference: {fixed_reference}")

            # Format the new line. Using 6 decimal places for precision.
            new_line = f"#box: {x_min} {new_y_min:.6f} {z_min} {x_max} {new_y_max:.6f} {z_max} {material_name}\n"
            new_lines.append(new_line)
            fracture_found = True
            print(f"Updated fracture in {filename}:")
            print(f"  Old Y: {y_min} to {y_max} (Aperture: {y_max - y_min:.6f})")
            print(f"  New Y: {new_y_min:.6f} to {new_y_max:.6f} (Aperture: {new_aperture:.6f})")
        else:
            new_lines.append(line)
            
    if not fracture_found:
        print(f"Warning: No box with material '{material_name}' found in {filename}.")
        return

    with open(filename, 'w') as f:
        f.writelines(new_lines)

