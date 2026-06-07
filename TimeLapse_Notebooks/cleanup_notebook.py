import json, sys
sys.stdout.reconfigure(encoding='utf-8')

nb_path = r'c:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\R_PhasePlaneFit.ipynb'
with open(nb_path, encoding='utf-8') as f:
    nb = json.load(f)

def cell_by_id(nb, cid):
    for c in nb['cells']:
        if c.get('id') == cid:
            return c
    raise KeyError(cid)

cell_by_id(nb, 'ef61c5a0')['source'] = """\
# Dataset and shift index to analyse  (<-- tune)
LABEL     = '⅛λ'
SC1_IDX   = 5          # index into x_scatterer_1 for the chosen monitor dataset

dl = 0.002   # gprMax FDTD cell size [m]; positions snap to this grid

# Nominal scatterer positions rounded to mm, then snapped to FDTD grid.
# For 1/16-lambda: nominal 1.669 m = 834.5 cells -> gprMax snaps to 834*0.002 = 1.668 m.
# Using the snapped positions gives the true simulated displacement.
x_s1_base_g = np.round(x_scatterer_1[0]       / dl) * dl
x_s1_mon_g  = np.round(x_scatterer_1[SC1_IDX] / dl) * dl
true_dz = 0.0
true_dx = float(x_s1_mon_g - x_s1_base_g)

dz_mig = z_img[1]    - z_img[0]
dx_mig = x_traces[1] - x_traces[0]
kz_c   = 2 * np.pi / wavelength

# ── Crop to the moving scatterer only ─────────────────────────────────────────
# Two scatterers are present: one moves (+dx), one is fixed at x_centre + 3λ.
# Without cropping the cross-spectrum averages both contributions and the phase
# plane slope converges to dx/2 instead of dx  (one fixed + one moving → average).
x_crop_cen = 0.5 * (x_s1_base_g + x_s1_mon_g)
x_crop_hw  = 2.5 * wavelength          # half-width: fixed scatterer is >4λ away

ix_lo = np.searchsorted(x_traces, x_crop_cen - x_crop_hw)
ix_hi = np.searchsorted(x_traces, x_crop_cen + x_crop_hw)
x_crop = x_traces[ix_lo:ix_hi]

base_crop = migrated_gz['Baseline'][:, ix_lo:ix_hi]
mon_crop  = migrated_gz[LABEL][:,     ix_lo:ix_hi]

# ── Estimate ──────────────────────────────────────────────────────────────────
dz_est, dx_est, XS, kz_ax, kx_ax = estimate_shift_2d(
    base_crop, mon_crop, dz_mig, dx_mig, kz_c
)

XS_s     = np.fft.fftshift(XS)
kz_s     = np.fft.fftshift(kz_ax)
kx_s     = np.fft.fftshift(kx_ax)
energy   = np.abs(XS_s)
phi_show = np.where(energy > 0.01 * energy.max(),
                    np.degrees(np.angle(XS_s)), np.nan)

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 6))
fig.suptitle(
    f'2D cross-spectrum — Horizontal shift  |  dataset: {LABEL}\\n'
    f'Cropped to moving scatterer  x = [{x_crop[0]:.3f}, {x_crop[-1]:.3f}] m  '
    f'(FDTD grid {dl*1e3:.0f} mm)',
    fontsize=11
)

# Col 0: difference image (cropped)
ax   = axes[0]
diff = mon_crop - base_crop
vmax = np.max(np.abs(diff))
ax.pcolormesh(x_crop * 100, z_img * 100, diff,
              cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto')
ax.axvline(x_s1_base_g * 100, color='g', ls='--', lw=1, label='x_s1 base (snapped)')
ax.axvline(x_s1_mon_g  * 100, color='r', ls='--', lw=1, label=f'x_s1 {LABEL} (snapped)')
ax.set_title(f'Difference (mon − base)\\ndx_true = {true_dx*1e3:.2f} mm  (FDTD-snapped)', fontsize=9)
ax.set_xlabel('x [cm]');  ax.set_ylabel('z [cm]')
ax.invert_yaxis();  ax.legend(fontsize=8)

# Col 1: cross-spectrum phase (cropped patch)
ax2 = axes[1]
im2 = ax2.pcolormesh(kx_s, kz_s, phi_show,
                     cmap='RdBu_r', vmin=-180, vmax=180, shading='auto')
klim = 1.4 * kz_c
for sgn in [-1, 1]:
    ax2.axhline(sgn * klim, color='k', lw=0.8, ls='--', alpha=0.5)
    ax2.axvline(sgn * klim, color='k', lw=0.8, ls='--', alpha=0.5)
ax2.set_title('Cross-spectrum phase [deg]\\ndashed = fit band', fontsize=9)
ax2.set_xlabel('kx [rad/m]');  ax2.set_ylabel('kz [rad/m]')
ax2.set_xlim(kx_s.min() / 2, kx_s.max() / 2)
ax2.set_ylim(kz_s.min() / 2, kz_s.max() / 2)
plt.colorbar(im2, ax=ax2, fraction=0.046)

# Col 2: results
ax3 = axes[2]
ax3.axis('off')
txt = (
    f'True:  dz = {true_dz * 1e3:+7.3f} mm\\n'
    f'       dx = {true_dx * 1e3:+7.3f} mm\\n\\n'
    f'Est:   dz = {dz_est  * 1e3:+7.3f} mm\\n'
    f'       dx = {dx_est  * 1e3:+7.3f} mm\\n\\n'
    f'Err:   dz = {(dz_est - true_dz) * 1e3:+.4f} mm\\n'
    f'       dx = {(dx_est - true_dx) * 1e3:+.4f} mm'
)
ax3.text(0.05, 0.92, txt, transform=ax3.transAxes,
         fontsize=10, family='monospace', verticalalignment='top')
ax3.set_title('Results', fontsize=9)

plt.tight_layout()
plt.show()
"""

with open(nb_path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print('Updated cell ef61c5a0 with FDTD-snapped true_dx.')
