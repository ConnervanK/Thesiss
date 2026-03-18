import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle, FancyArrowPatch
from pathlib import Path

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
    radius_1 = 0.3
    radius_2 = radius_1 + lam / 4.0

    theta = np.linspace(0, np.pi, 400)
    arc1_x = receiver_x + radius_1 * np.cos(theta)
    arc1_y = receiver_y + radius_1 * np.sin(theta)
    arc2_x = receiver_x + radius_2 * np.cos(theta)
    arc2_y = receiver_y + radius_2 * np.sin(theta)

    # Calculate first Fresnel zone intersection with fracture
    dy_frac_top = fracture_y0 - receiver_y
    delta_x_top = np.sqrt(max(radius_1**2 - dy_frac_top**2, 0.0))
    fresnel_x_left_top = receiver_x - delta_x_top
    fresnel_x_right_top = receiver_x + delta_x_top

    dy_frac_bottom = fracture_y1 - receiver_y
    delta_x_bottom = np.sqrt(max(radius_1**2 - dy_frac_bottom**2, 0.0))
    fresnel_x_left_bottom = receiver_x - delta_x_bottom
    fresnel_x_right_bottom = receiver_x + delta_x_bottom

    # Calculate second Fresnel zone intersection with fracture
    delta_x_top_2 = np.sqrt(max(radius_2**2 - dy_frac_top**2, 0.0))
    fresnel_x_left_top_2 = receiver_x - delta_x_top_2
    fresnel_x_right_top_2 = receiver_x + delta_x_top_2

    delta_x_bottom_2 = np.sqrt(max(radius_2**2 - dy_frac_bottom**2, 0.0))
    fresnel_x_left_bottom_2 = receiver_x - delta_x_bottom_2
    fresnel_x_right_bottom_2 = receiver_x + delta_x_bottom_2

    # For visualization: second Fresnel zone at fracture top
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
    ax.plot(arc1_x, arc1_y, color='blue', linewidth=2.0, label=f'r = {radius_1:.1f} m')
    ax.plot(arc2_x, arc2_y, color='navy', linestyle='--', linewidth=2.0,
            label=f'r = {radius_2:.1f} + 1/4 $\\lambda$')
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

    # Return first and second Fresnel zone x-coordinates where they intersect the fracture
    return {
        'fresnel_x_left_top': fresnel_x_left_top,
        'fresnel_x_right_top': fresnel_x_right_top,
        'fresnel_x_left_bottom': fresnel_x_left_bottom,
        'fresnel_x_right_bottom': fresnel_x_right_bottom,
        'fresnel_x_left_top_2': fresnel_x_left_top_2,
        'fresnel_x_right_top_2': fresnel_x_right_top_2,
        'fresnel_x_left_bottom_2': fresnel_x_left_bottom_2,
        'fresnel_x_right_bottom_2': fresnel_x_right_bottom_2,
    }


def plot_alternating_fracture_conductivity(
        length,
        depth,
        fracture_thickness,
        baseline_conductivity,
        fracture_y0=0.4,
        upper_layer_thickness=0.1,
        x_receiver=0.5,
        y_receiver=0.1,
        percent_delta=5.0,
):
        """
        Plot fracture blocks with alternating conductivity perturbations.

        Block width is set equal to the fracture aperture (fracture_thickness).
        Conductivity alternates as +percent_delta, -percent_delta, +percent_delta, ...
        to keep the average perturbation close to zero.
        """
        if fracture_thickness <= 0:
                raise ValueError('fracture_thickness must be positive.')
        if baseline_conductivity <= 0:
                raise ValueError('baseline_conductivity must be positive.')

        fracture_y1 = fracture_y0 + fracture_thickness
        block_width = fracture_thickness

        n_blocks = int(np.ceil(length / block_width))
        pct_values = np.zeros(n_blocks, dtype=float)
        n_alternating = n_blocks if n_blocks % 2 == 0 else n_blocks - 1
        for i in range(n_alternating):
                pct_values[i] = percent_delta if i % 2 == 0 else -percent_delta
        conductivity_values = baseline_conductivity * (1.0 + pct_values / 100.0)

        fig, ax = plt.subplots(figsize=(11, 5))

        outer = Rectangle((0, 0), length, depth, fill=False, edgecolor='black', linewidth=2)
        ax.add_patch(outer)

        upper = Rectangle((0, 0), length, upper_layer_thickness, facecolor='#cfe8ff', edgecolor='none', alpha=0.75)
        ax.add_patch(upper)

        cmap = plt.get_cmap('coolwarm')
        norm = Normalize(vmin=-percent_delta, vmax=percent_delta)

        for i in range(n_blocks):
                x0 = i * block_width
                width_i = min(block_width, length - x0)
                if width_i <= 0:
                        break

                pct_i = pct_values[i]
                color_i = cmap(norm(pct_i))
                block = Rectangle(
                        (x0, fracture_y0),
                        width_i,
                        fracture_thickness,
                        facecolor=color_i,
                        edgecolor='black',
                        linewidth=0.6,
                        alpha=0.95,
                        zorder=4,
                )
                ax.add_patch(block)

        ax.scatter(x_receiver, y_receiver, marker='v', s=120, color='black', zorder=6, label='Receiver')
        ax.text(x_receiver, upper_layer_thickness / 2, 'Air', ha='center', va='center', fontsize=10)
        ax.text(
                min(length * 0.5, length - 0.02),
                fracture_y0 + fracture_thickness + 0.02,
                f'Alternating conductivity blocks (block width = aperture = {fracture_thickness:.4f} m)',
                ha='center',
                va='bottom',
                fontsize=10,
                color='black',
        )

        sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.02)
        cbar.set_label('Conductivity perturbation [%]')

        ax.set_title('Fracture with Alternating Conductivity Perturbations', fontsize=15)
        ax.set_xlabel('Length (m)', fontsize=12)
        ax.set_ylabel('Depth (m)', fontsize=12)
        ax.set_xlim(0, length)
        ax.set_ylim(depth, 0)
        ax.set_aspect('equal', adjustable='box')
        ax.grid(True, linestyle='--', alpha=0.35)
        ax.legend(loc='upper right')
        plt.tight_layout()
        plt.show()

        return {
                'block_width': block_width,
                'n_blocks': n_blocks,
                'percent_delta_sequence': pct_values,
                'conductivity_sequence': conductivity_values,
                'mean_percent_delta': float(np.mean(pct_values)),
                'mean_conductivity': float(np.mean(conductivity_values)),
        }


def export_alternating_fracture_conductivity_gprmax(
        length,
        fracture_thickness,
        baseline_conductivity,
        relative_permittivity,
        fracture_y0=0.4,
        z0=0.0,
        z1=0.0002,
        percent_delta=5.0,
        base_material_name='fracture_alt',
        write_path=None,
):
        """
        Export alternating fracture conductivity blocks as gprMax commands.

        Block width is set to fracture_thickness and conductivity alternates
        +percent_delta/-percent_delta. If the block count is odd, the final
        block is set to 0% so the overall mean perturbation is exactly zero.
        """
        if fracture_thickness <= 0:
                raise ValueError('fracture_thickness must be positive.')
        if baseline_conductivity <= 0:
                raise ValueError('baseline_conductivity must be positive.')
        if relative_permittivity <= 0:
                raise ValueError('relative_permittivity must be positive.')
        if length <= 0:
                raise ValueError('length must be positive.')

        fracture_y1 = fracture_y0 + fracture_thickness
        block_width = fracture_thickness
        n_blocks = int(np.ceil(length / block_width))

        pct_values = np.zeros(n_blocks, dtype=float)
        n_alternating = n_blocks if n_blocks % 2 == 0 else n_blocks - 1
        for i in range(n_alternating):
                pct_values[i] = percent_delta if i % 2 == 0 else -percent_delta

        conductivity_values = baseline_conductivity * (1.0 + pct_values / 100.0)

        material_lines = []
        box_lines = []
        rows = []

        for i in range(n_blocks):
                x0 = i * block_width
                x1 = min(length, x0 + block_width)
                if x1 <= x0:
                        continue

                material_name = f'{base_material_name}_{i + 1:03d}'
                sigma_i = float(conductivity_values[i])
                pct_i = float(pct_values[i])

                material_line = (
                        f'#material: {relative_permittivity:.6g} {sigma_i:.6g} 1 0 {material_name}'
                )
                box_line = (
                        f'#box: {x0:.6f} {fracture_y0:.6f} {z0:.6f} '
                        f'{x1:.6f} {fracture_y1:.6f} {z1:.6f} {material_name}'
                )

                material_lines.append(material_line)
                box_lines.append(box_line)
                rows.append({
                        'block_index': i,
                        'x0': x0,
                        'x1': x1,
                        'percent_delta': pct_i,
                        'conductivity': sigma_i,
                        'material_name': material_name,
                })

        header = [
                '// Alternating fracture conductivity blocks generated by plotting.py',
                f'// block_width = {block_width:.6f} m (equals fracture aperture)',
                f'// mean percent perturbation = {float(np.mean(pct_values)):.6f} %',
                f'// baseline conductivity = {baseline_conductivity:.6g} S/m',
                f'// number of blocks = {len(rows)}',
                '',
                '// Materials',
        ]
        command_lines = header + material_lines + ['', '// Geometry boxes'] + box_lines
        command_text = '\n'.join(command_lines) + '\n'

        if write_path is not None:
                out_path = Path(write_path)
                out_path.parent.mkdir(parents=True, exist_ok=True)
                out_path.write_text(command_text, encoding='utf-8')

        return {
                'block_width': block_width,
                'n_blocks': len(rows),
                'percent_delta_sequence': pct_values,
                'conductivity_sequence': conductivity_values,
                'mean_percent_delta': float(np.mean(pct_values)),
                'mean_conductivity': float(np.mean(conductivity_values)),
                'material_lines': material_lines,
                'box_lines': box_lines,
                'command_lines': command_lines,
                'command_text': command_text,
                'rows': rows,
                'write_path': str(write_path) if write_path is not None else None,
        }


def create_base_gpr_max_input_template(
        title='Baseline Model',
        domain=(1.0, 0.5, 0.0002),
        dx_dy_dz=(0.0002, 0.0002, 0.0002),
        time_window=7e-9,
        num_threads=8,
        upper_air_depth=0.5,
        background_depth=0.4,
        fracture_y0=0.0918,
        fracture_y1=0.1,
        baseline_conductivity=0.01,
        relative_permittivity=6.0,
        percent_delta=5.0,
        base_material_name='fracture_alt',
        background_material='ice',
        air_material='air',
        waveform_type='ricker',
        waveform_amp=1.0,
        f_central=1.5e9,
        waveform_name='my_ricker',
        src_axis='z',
        src_pos=(0.5, 0.4, 0.0),
        rx_pos=(0.51, 0.4, 0.0),
        src_steps=(0.0375, 0.0, 0.0),
        rx_steps=(0.0375, 0.0, 0.0),
        snapshot_count=9,
        snapshot_dt=6e-10,
        geometry_view_name='Baseline_Model',
        messages='y',
        write_path=None,
):
        """
        Create a base gprMax .in template using block-wise fracture conductivity.

        The fracture interval [fracture_y0, fracture_y1] is discretized into
        x-blocks whose width equals fracture aperture (fracture_y1 - fracture_y0).
        Each block gets an alternating +percent_delta/-percent_delta conductivity
        perturbation via generated #material and #box commands.
        """
        lx, ly, lz = domain
        dx, dy, dz = dx_dy_dz
        sx, sy, sz = src_pos
        rx, ry, rz = rx_pos
        ssx, ssy, ssz = src_steps
        rsx, rsy, rsz = rx_steps

        if lx <= 0 or ly <= 0 or lz <= 0:
                raise ValueError('domain values must be positive.')
        if dx <= 0 or dy <= 0 or dz <= 0:
                raise ValueError('dx_dy_dz values must be positive.')
        if time_window <= 0:
                raise ValueError('time_window must be positive.')
        if fracture_y1 <= fracture_y0:
                raise ValueError('fracture_y1 must be greater than fracture_y0.')
        if snapshot_count < 1:
                raise ValueError('snapshot_count must be >= 1.')
        if snapshot_dt <= 0:
                raise ValueError('snapshot_dt must be positive.')

        fracture_thickness = fracture_y1 - fracture_y0
        fracture_export = export_alternating_fracture_conductivity_gprmax(
                length=lx,
                fracture_thickness=fracture_thickness,
                baseline_conductivity=baseline_conductivity,
                relative_permittivity=relative_permittivity,
                fracture_y0=fracture_y0,
                z0=0.0,
                z1=lz,
                percent_delta=percent_delta,
                base_material_name=base_material_name,
                write_path=None,
        )

        lines = [
                f'#title: {title}',
                f'#num_threads: {int(num_threads)}',
                f'#domain: {lx:.6f} {ly:.6f} {lz:.6f}',
                f'#dx_dy_dz: {dx:.6f} {dy:.6f} {dz:.6f}',
                f'#time_window: {time_window:.6g}',
                '',
                '#material: 6 1e-6 1 0 ice',
                '#material: 1 0 1 0 air',
                '#material: 80 0.01 1 0 water',
                '// Block-wise fracture materials',
                *fracture_export['material_lines'],
                '',
                '// Above surface',
                f'#box: 0 0 0 {lx:.6f} {upper_air_depth:.6f} {lz:.6f} {air_material}',
                '',
                '// Background Model',
                f'#box: 0 0 0 {lx:.6f} {background_depth:.6f} {lz:.6f} {background_material}',
                '',
                '// Thin fracture (block-wise alternating conductivity)',
                *fracture_export['box_lines'],
                '',
                f'#waveform: {waveform_type} {waveform_amp:.6g} {f_central:.6g} {waveform_name}',
                f'#hertzian_dipole: {src_axis} {sx:.6f} {sy:.6f} {sz:.6f} {waveform_name}',
                f'#rx: {rx:.6f} {ry:.6f} {rz:.6f}',
                '',
                f'#src_steps: {ssx:.6f} {ssy:.6f} {ssz:.6f}',
                f'#rx_steps: {rsx:.6f} {rsy:.6f} {rsz:.6f}',
                '',
                '#python:',
                f'for i in range(1, {int(snapshot_count) + 1}):',
                (
                        "    print('#snapshot: 0 0 0 "
                        f"{lx:.6f} {ly:.6f} {lz:.6f} "
                        f"{dx:.6f} {dy:.6f} {dz:.6f} "
                        f"{{}} snapshot{{}}'.format(i*{snapshot_dt:.6g}, i))"
                ),
                '#end_python:',
                '',
                (
                        f'#geometry_view: 0 0 0 {lx:.6f} {ly:.6f} {lz:.6f} '
                        f'{dx:.6f} {dy:.6f} {dz:.6f} {geometry_view_name} n'
                ),
                f'#messages: {messages}',
                '',
        ]

        template_text = '\n'.join(lines)

        if write_path is not None:
                out_path = Path(write_path)
                out_path.parent.mkdir(parents=True, exist_ok=True)
                out_path.write_text(template_text, encoding='utf-8')

        return {
                'template_lines': lines,
                'template_text': template_text,
                'fracture_export': fracture_export,
                'write_path': str(write_path) if write_path is not None else None,
        }


def create_base_gprmax_input_template(**kwargs):
        """Backward-compatible alias for create_base_gpr_max_input_template."""
        return create_base_gpr_max_input_template(**kwargs)
