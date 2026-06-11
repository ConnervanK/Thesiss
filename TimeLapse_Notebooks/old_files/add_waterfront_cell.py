import json, sys
sys.stdout.reconfigure(encoding='utf-8')

nb_path = r'c:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\R_PhasePlaneFit.ipynb'
with open(nb_path, encoding='utf-8') as f:
    nb = json.load(f)

# ---------- Markdown header cell ----------
md_cell = {
    "cell_type": "markdown",
    "id": "c3e7a1f2",
    "metadata": {},
    "source": "# Moving Water Front Test"
}

# ---------- Code cell ----------
code_src = '''\
from scipy.signal import fftconvolve
from scipy.signal import hilbert as sp_hilbert
from matplotlib.patches import Rectangle

# ── Physics ──────────────────────────────────────────────────────────────────────
c_light_wf  = 0.299792458
eps_ice_wf  = 3.15;  n_ice_wf = np.sqrt(eps_ice_wf)
n_air_wf    = 1.0;   n_water_wf = 9.0
v_ice_wf    = c_light_wf / n_ice_wf   # ≈ 0.1689 m/ns
f_c_wf      = 1.5                      # GHz
lam_wf      = v_ice_wf / f_c_wf       # ≈ 0.1126 m
kz_c_wf     = 2 * np.pi / lam_wf

# ── Grid ─────────────────────────────────────────────────────────────────────────
Nx_wf     = 64;   dx_wf = lam_wf / 4
x_wf      = np.arange(Nx_wf) * dx_wf  # [m]

# dt chosen so one time-sample == one depth-sample after migration
dz_wf     = lam_wf / 20
v_mig_wf  = v_ice_wf / 2
dt_wf     = dz_wf / v_mig_wf          # ≈ 0.0667 ns
Nt_wf     = 512
t_wf      = np.arange(Nt_wf) * dt_wf  # [ns]

Nz_wf     = 256
z_wf      = np.arange(Nz_wf) * dz_wf  # [m]

# ── Fracture geometry ─────────────────────────────────────────────────────────────
iz0_wf       = Nz_wf // 3              # fracture top index in depth axis
z0_wf        = z_wf[iz0_wf]           # ≈ 0.48 m
d_frac_wf    = lam_wf / 20            # fracture thickness ≈ 5.6 mm

x_front_base = x_wf[Nx_wf // 2]      # water front in baseline (domain midpoint)
delta_x_true = lam_wf / 8             # subwavelength shift
x_front_mon  = x_front_base + delta_x_true

# ── Reflection coefficients ───────────────────────────────────────────────────────
def _refl(n1, n2): return (n1 - n2) / (n1 + n2)
R_top_air_wf   = _refl(n_ice_wf, n_air_wf)
R_bot_air_wf   = _refl(n_air_wf,   n_ice_wf)
R_top_water_wf = _refl(n_ice_wf, n_water_wf)
R_bot_water_wf = _refl(n_water_wf, n_ice_wf)
dR_wf          = R_top_water_wf - R_top_air_wf  # edge diffraction amplitude

# Integer sample delays through fracture
# Since dt_wf = dz_wf/v_mig_wf and the thin-layer delay = d*n_fill/n_ice / dz_wf,
# time delay in samples == depth delay in samples (1:1 mapping).
iz_del_air_wf   = max(1, int(round(d_frac_wf * n_air_wf   / n_ice_wf / dz_wf)))
iz_del_water_wf = max(1, int(round(d_frac_wf * n_water_wf / n_ice_wf / dz_wf)))

# ── Ricker wavelet (time domain) ──────────────────────────────────────────────────
Nw_t_wf   = 2 * int(round(3.5 / f_c_wf / dt_wf)) + 1
t_wl_wf   = (np.arange(Nw_t_wf) - Nw_t_wf // 2) * dt_wf
w_rk_t_wf = (1 - 2*(np.pi*f_c_wf*t_wl_wf)**2) * np.exp(-(np.pi*f_c_wf*t_wl_wf)**2)

# ── B-scan synthesis: flat reflector + edge diffraction from fill boundary ────────
def _make_bscan_wf(x_front):
    """
    Per-trace reflectivity determined by fill at that x (water left, air right).
    Edge diffraction from the fill step at x_front is added via point_bscan,
    amplitude-normalised so the hyperbola apex == dR_wf (the fill impedance step).
    """
    bscan = np.zeros((Nt_wf, Nx_wf))
    it0 = int(round(2.0 * z0_wf / v_ice_wf / dt_wf))  # == iz0_wf
    for ix in range(Nx_wf):
        if x_wf[ix] < x_front:
            Rt, Rb, it_d = R_top_water_wf, R_bot_water_wf, iz_del_water_wf
        else:
            Rt, Rb, it_d = R_top_air_wf,   R_bot_air_wf,   iz_del_air_wf
        r = np.zeros(Nt_wf)
        r[it0]        += Rt
        r[it0 + it_d] += Rb
        bscan[:, ix] = fftconvolve(r, w_rk_t_wf, mode='same')
    # Edge diffraction: point_bscan gives 1/sqrt(r) spreading;
    # multiply by sqrt(z0) so amplitude at apex (r=z0) equals dR_wf.
    edge = dR_wf * np.sqrt(z0_wf) * point_bscan(
        x_front, z0_wf, x_wf, t_wf, v_ice_wf, f_c_wf
    )
    return bscan + edge

bscan_base_wf = _make_bscan_wf(x_front_base)
bscan_mon_wf  = _make_bscan_wf(x_front_mon)

# ── Gazdag migration ──────────────────────────────────────────────────────────────
print('Migrating baseline ...')
mig_base_wf = gazdag_migration(bscan_base_wf, x_wf, t_wf, z_wf, v_ice_wf)
print('Migrating monitor  ...')
mig_mon_wf  = gazdag_migration(bscan_mon_wf,  x_wf, t_wf, z_wf, v_ice_wf)
diff_mig_wf = mig_mon_wf - mig_base_wf

# ── Cross-spectrum ROI ────────────────────────────────────────────────────────────
# z: wavelet half-width on each side of the fracture
Nw_z_wf  = 2 * int(round(3.5 * lam_wf / dz_wf)) + 1
n_mz_wf  = Nw_z_wf // 2 + 5
iz_lo_wf = max(0,     iz0_wf - n_mz_wf)
iz_hi_wf = min(Nz_wf, iz0_wf + iz_del_water_wf + n_mz_wf + 1)

# x: ±4λ centred on baseline water front (edge diffraction dominates there)
n_mx_wf  = int(round(4.0 * lam_wf / dx_wf)) + 1
ix_c_wf  = int(np.argmin(np.abs(x_wf - x_front_base)))
ix_lo_wf = max(0,     ix_c_wf - n_mx_wf)
ix_hi_wf = min(Nx_wf, ix_c_wf + n_mx_wf + 1)

roi_base_wf = mig_base_wf[iz_lo_wf:iz_hi_wf, ix_lo_wf:ix_hi_wf]
roi_mon_wf  = mig_mon_wf [iz_lo_wf:iz_hi_wf, ix_lo_wf:ix_hi_wf]

# Analytic signal along z (breaks Hermitian symmetry → enables phi_0 estimation)
a_base_wf = sp_hilbert(roi_base_wf, axis=0)
a_mon_wf  = sp_hilbert(roi_mon_wf,  axis=0)

dz_est_wf, dx_est_wf, XS_wf, kz_ax_wf, kx_ax_wf = estimate_shift_2d(
    a_base_wf, a_mon_wf, dz_wf, dx_wf, kz_c_wf
)

# Re-run WLS for phase intercept phi_0
KZ_wf, KX_wf_g = np.meshgrid(kz_ax_wf, kx_ax_wf, indexing='ij')
w_wf   = np.abs(XS_wf);  phi_wf = np.angle(XS_wf)
band_wf = (np.abs(KZ_wf) < 1.4*kz_c_wf) & (np.abs(KX_wf_g) < 1.4*kz_c_wf)
mask_wf = (w_wf > 0.1*w_wf.max()) & band_wf & ((np.abs(KZ_wf)+np.abs(KX_wf_g)) > 0)
A_wf    = np.column_stack([KZ_wf[mask_wf], KX_wf_g[mask_wf], np.ones(mask_wf.sum())])
W_wf    = w_wf[mask_wf]
c_wf    = np.linalg.lstsq(A_wf*W_wf[:,None], phi_wf[mask_wf]*W_wf, rcond=None)[0]
phi_0_wf = c_wf[2]

# ── Figure: 3 × 2 layout ─────────────────────────────────────────────────────────
t_ns = t_wf.copy()                 # already in ns
x_cm = x_wf * 100
z_cm = z_wf * 100

z0_cm    = z0_wf * 100
xfb_cm   = x_front_base * 100
xfm_cm   = x_front_mon  * 100
roi_zlo  = z_wf[iz_lo_wf]         * 100
roi_zhi  = z_wf[iz_hi_wf - 1]     * 100
roi_xlo  = x_wf[ix_lo_wf]         * 100
roi_xhi  = x_wf[ix_hi_wf - 1]     * 100
t0_ns_v  = 2.0 * z0_wf / v_ice_wf

fig, axes = plt.subplots(3, 2, figsize=(12, 14))
fig.suptitle(
    f'Moving Water Front Test — fracture full lateral extent,  d = \\u03bb/20 = {d_frac_wf*1e3:.1f} mm\\n'
    f'Baseline: water front at x = {xfb_cm:.1f} cm   \\u2192   '
    f'Monitor: front shifted by \\u03b4x = \\u03bb/8 = {delta_x_true*1e3:.1f} mm',
    fontsize=11
)

# Row 0: raw time-domain B-scans
vmax_raw = max(np.max(np.abs(bscan_base_wf)), np.max(np.abs(bscan_mon_wf)))
for col, (bscan, title, xf_cm) in enumerate([
    (bscan_base_wf, 'Raw baseline B-scan (time domain)\\nHyperbola centred at water front', xfb_cm),
    (bscan_mon_wf,  'Raw monitor B-scan (time domain)\\nHyperbola shifted by \\u03b4x = \\u03bb/8', xfm_cm),
]):
    ax = axes[0, col]
    ax.pcolormesh(x_cm, t_ns, bscan, cmap='RdBu_r',
                  vmin=-vmax_raw, vmax=vmax_raw, shading='auto')
    ax.axhline(t0_ns_v, color='gold', lw=1.2, ls='--',
               label=f'TWT(z\\u2080) = {t0_ns_v:.2f} ns')
    ax.axvline(xf_cm, color='cyan', lw=1.2, ls=':',
               label=f'x_front = {xf_cm:.1f} cm')
    ax.set_xlabel('x [cm]');  ax.set_ylabel('Two-way time [ns]')
    ax.set_title(title, fontsize=9);  ax.legend(fontsize=7)

# Row 1: migrated images (common colour scale)
vmax_mig = max(np.max(np.abs(mig_base_wf)), np.max(np.abs(mig_mon_wf)))
for col, (img, title, xf_cm) in enumerate([
    (mig_base_wf, 'Migrated baseline\\nHyperbola collapsed to focused edge', xfb_cm),
    (mig_mon_wf,  'Migrated monitor\\nFocused edge shifted to new front',    xfm_cm),
]):
    ax = axes[1, col]
    ax.pcolormesh(x_cm, z_cm, img, cmap='RdBu_r',
                  vmin=-vmax_mig, vmax=vmax_mig, shading='auto')
    ax.axhline(z0_cm, color='gold', lw=1.2, ls='--', label='z\\u2080 (fracture)')
    ax.axvline(xf_cm, color='cyan', lw=1.2, ls=':',
               label=f'x_front = {xf_cm:.1f} cm')
    ax.set_xlabel('x [cm]');  ax.set_ylabel('Depth [cm]')
    ax.invert_yaxis();  ax.set_title(title, fontsize=9);  ax.legend(fontsize=7)

# Row 2, left: migrated difference with ROI box
ax = axes[2, 0]
vmax_d = np.max(np.abs(diff_mig_wf)) or 1.0
ax.pcolormesh(x_cm, z_cm, diff_mig_wf, cmap='RdBu_r',
              vmin=-vmax_d, vmax=vmax_d, shading='auto')
ax.axhline(z0_cm,  color='gold',  lw=1.2, ls='--')
ax.axvline(xfb_cm, color='cyan',  lw=1.2, ls=':', label=f'x_front base = {xfb_cm:.1f} cm')
ax.axvline(xfm_cm, color='lime',  lw=1.2, ls=':', label=f'x_front mon  = {xfm_cm:.1f} cm')
roi_rect = Rectangle(
    (roi_xlo, roi_zlo), roi_xhi - roi_xlo, roi_zhi - roi_zlo,
    edgecolor='yellow', facecolor='none', lw=1.8
)
ax.add_patch(roi_rect)
ax.set_xlabel('x [cm]');  ax.set_ylabel('Depth [cm]')
ax.invert_yaxis()
ax.set_title(
    'Migrated difference  (Monitor \\u2212 Baseline)\\n'
    '[yellow box = ROI for phase-plane fit]', fontsize=9
)
ax.legend(fontsize=7)

# Row 2, right: results text
ax3 = axes[2, 1];  ax3.axis('off')
err_dx = (dx_est_wf - delta_x_true) * 1e3
txt = (
    'Model:\\n'
    f'  Fill: water (left)  |  air (right)\\n'
    f'  R_top water = {R_top_water_wf:+.3f}   R_top air = {R_top_air_wf:+.3f}\\n'
    f'  \\u0394R (edge) = {dR_wf:+.3f}\\n'
    f'  d = \\u03bb/20 = {d_frac_wf*1e3:.1f} mm\\n\\n'
    'True displacement:\\n'
    f'  \\u03b4x = \\u03bb/8 = {delta_x_true*1e3:.2f} mm\\n'
    f'  \\u03b4z = 0 mm\\n\\n'
    'Estimated (phase-plane fit on ROI):\\n'
    f'  \\u0394z = {dz_est_wf*1e3:+.1f} mm   [true: 0]\\n'
    f'  \\u0394x = {dx_est_wf*1e3:+.1f} mm   [true: {delta_x_true*1e3:.1f} mm]\\n'
    f'  Error = {err_dx:+.1f} mm\\n\\n'
    f'  \\u03c6\\u2080 = {np.degrees(phi_0_wf):+.1f}\\u00b0\\n'
    f'  [partial polarity change in \\u03b4x = {delta_x_true/lam_wf:.3f}\\u03bb strip]'
)
ax3.text(0.04, 0.97, txt, transform=ax3.transAxes,
         fontsize=10, family='monospace', verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.6))
ax3.set_title('Results', fontsize=9)

plt.tight_layout()
plt.show()

print(f'True \\u03b4x = {delta_x_true*1e3:.2f} mm  |  '
      f'Estimated \\u0394x = {dx_est_wf*1e3:.2f} mm  |  '
      f'Error = {err_dx:+.2f} mm')
'''

code_cell = {
    "cell_type": "code",
    "execution_count": None,
    "id": "d9f5c2a8",
    "metadata": {},
    "outputs": [],
    "source": code_src
}

# ---------- Insert after cell 505c92c1 ----------
idx = next(i for i, c in enumerate(nb['cells']) if c.get('id') == '505c92c1')
nb['cells'].insert(idx + 1, md_cell)
nb['cells'].insert(idx + 2, code_cell)

with open(nb_path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f'Inserted markdown cell c3e7a1f2 and code cell d9f5c2a8 after 505c92c1.')
print(f'Notebook now has {len(nb["cells"])} cells.')
