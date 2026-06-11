
import json, sys
sys.stdout.reconfigure(encoding='utf-8')

nb_path = r'c:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\CWT_playground.ipynb'
with open(nb_path, encoding='utf-8') as f:
    nb = json.load(f)

new_code = """\
# =============================================================================
# 2D cross-spectrum estimator: error vs shift magnitude
# Row 0: absolute error |est - true| [mm]
# Row 1: relative error |est - true| / ||true shift||  [%]
#         normalised by the total applied shift magnitude so that small and
#         large shifts are on a comparable footing.
# =============================================================================
from scipy.ndimage import shift as nd_shift

shift_fracs = np.linspace(0.02, 0.60, 40)

case_defs = [
    ('Vertical',
     lambda s: (s * lam,              0.0              ),
     1.0 / 3.0,
     'steelblue', 'cornflowerblue'),
    ('Horizontal',
     lambda s: (0.0,                  s * lam          ),
     1.0 / 3.0,
     'tomato', 'lightsalmon'),
    ('Diagonal',
     lambda s: (s * lam / np.sqrt(2), s * lam / np.sqrt(2)),
     1.0 / (3.0 * np.sqrt(2)),
     'seagreen', 'mediumaquamarine'),
]

fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharey='row')
fig.suptitle(
    f'2D cross-spectrum: error vs shift magnitude  '
    f'(noiseless, band-limited to 1.5*kz_c)\\n'
    f'lambda = {lam*1e3:.1f} mm  |  '
    f'dz_sample = {dz_mig*1e3:.2f} mm  |  '
    f'dx_sample = {dx_mig*1e3:.2f} mm',
    fontsize=10
)

for col, (name, shift_fn, wrap_lim, col_z, col_x) in enumerate(case_defs):
    ax_abs = axes[0, col]
    ax_rel = axes[1, col]

    err_z_abs, err_x_abs = [], []
    err_z_rel, err_x_rel = [], []

    for sf in shift_fracs:
        true_dz, true_dx = shift_fn(sf)
        norm_mm = sf * lam * 1e3          # total shift magnitude [mm]

        mon = nd_shift(mig_v_base,
                       (true_dz / dz_mig, true_dx / dx_mig), order=3)
        dz_est, dx_est, *_ = estimate_shift_2d(
            mig_v_base, mon, dz_mig, dx_mig, kz_c
        )

        ez = abs(dz_est - true_dz) * 1e3   # mm
        ex = abs(dx_est - true_dx) * 1e3

        err_z_abs.append(max(ez, 1e-5))
        err_x_abs.append(max(ex, 1e-5))
        err_z_rel.append(max(ez / norm_mm * 100, 1e-4))   # percent
        err_x_rel.append(max(ex / norm_mm * 100, 1e-4))

    err_z_abs = np.array(err_z_abs)
    err_x_abs = np.array(err_x_abs)
    err_z_rel = np.array(err_z_rel)
    err_x_rel = np.array(err_x_rel)

    for ax, y_z, y_x, ylabel, ref_lines in [
        (ax_abs, err_z_abs, err_x_abs,
         'Absolute error [mm]',
         [(dz_mig * 1e3,   'dimgray',    f'Sample spacing ({dz_mig*1e3:.2f} mm)'),
          (lam * 1e3 / 10, 'darkorange', f'lambda/10 ({lam*1e3/10:.1f} mm)'),
          (lam * 1e3 / 100,'goldenrod',  f'lambda/100 ({lam*1e3/100:.2f} mm)')]),
        (ax_rel, err_z_rel, err_x_rel,
         'Relative error [%]',
         [(100,   'dimgray',    '100 % (error = shift)'),
          (10,    'darkorange', '10 %'),
          (1,     'goldenrod',  '1 %'),
          (0.1,   'lightgreen', '0.1 %')]),
    ]:
        ax.semilogy(shift_fracs, y_z, '-o',  ms=4, lw=1.5,
                    color=col_z, label='|err dz|')
        ax.semilogy(shift_fracs, y_x, '--s', ms=4, lw=1.5,
                    color=col_x, label='|err dx|')

        for level, color, label in ref_lines:
            ax.axhline(level, color=color, ls=':', lw=1.1, alpha=0.8,
                       label=label)

        ax.axvline(wrap_lim, color='red', ls='--', lw=1.3,
                   label=f'Wrap limit ({wrap_lim:.2f} lambda)')

        ax.set_xlabel('Shift magnitude [lambda]')
        if col == 0:
            ax.set_ylabel(ylabel)
        ax.set_xlim(shift_fracs[0] - 0.01, 0.62)
        ax.legend(fontsize=7, loc='upper left')
        ax.grid(True, alpha=0.2, which='both')

    ax_abs.set_title(f'{name} shift')

plt.tight_layout()
plt.show()
"""

# Replace cell 46
nb['cells'][46]['source']          = new_code
nb['cells'][46]['outputs']         = []
nb['cells'][46]['execution_count'] = None

with open(nb_path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f'Updated cell 46. Total cells: {len(nb["cells"])}')
