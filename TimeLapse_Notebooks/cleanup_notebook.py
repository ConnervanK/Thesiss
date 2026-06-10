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

cell_by_id(nb, '505c92c1')['source'] = """\
from scipy.signal import fftconvolve
from scipy.signal import hilbert as sp_hilbert
from scipy.signal.windows import tukey

# ── Physics ────────────────────────────────────────────────────────────────────
c_light   = 0.299792458          # m / ns
eps_ice   = 3.15
n_ice     = np.sqrt(eps_ice)
n_air     = 1.0                  # baseline fill
n_water   = 9.0                  # sqrt(81), monitor fill

v_ice     = c_light / n_ice      # ≈ 0.1689 m/ns
f_c       = 1.5                  # GHz
lam       = v_ice / f_c          # ≈ 0.1126 m
v_mig     = v_ice / 2
d         = lam / 20             # fracture thickness ≈ 5.6 mm

dz_mig = lam / 20
dx_mig = lam / 4
kz_c   = 2 * np.pi / lam

Nz, Nx = 256, 64
z_mc   = np.arange(Nz) * dz_mig
x_mc   = np.arange(Nx) * dx_mig

# ── Fresnel reflection coefficients ───────────────────────────────────────────
def refl(n1, n2):
    return (n1 - n2) / (n1 + n2)

R_top_air   = refl(n_ice, n_air)    # ice → air   ≈ +0.279
R_bot_air   = refl(n_air, n_ice)    # air → ice   ≈ -0.279
R_top_water = refl(n_ice, n_water)  # ice → water ≈ -0.671  (polarity reversal)
R_bot_water = refl(n_water, n_ice)  # water → ice ≈ +0.671

dz_del_air   = d * n_air   / n_ice
dz_del_water = d * n_water / n_ice

iz_del_air   = max(1, int(round(dz_del_air   / dz_mig)))
iz_del_water = max(1, int(round(dz_del_water / dz_mig)))

# ── Ricker wavelet (migration-depth domain) ────────────────────────────────────
Nw   = 2 * int(round(3.5 * lam / dz_mig)) + 1
z_wl = (np.arange(Nw) - Nw // 2) * dz_mig
u_wl = (np.pi * f_c * z_wl / v_mig) ** 2
w_rk = (1 - 2 * u_wl) * np.exp(-u_wl)

# ── Reflectivity series & 1D traces ───────────────────────────────────────────
iz0 = Nz // 3

r_base = np.zeros(Nz)
r_base[iz0]              += R_top_air
r_base[iz0 + iz_del_air] += R_bot_air

r_mon = np.zeros(Nz)
r_mon[iz0]                += R_top_water
r_mon[iz0 + iz_del_water] += R_bot_water

trace_base = fftconvolve(r_base, w_rk, mode='same')
trace_mon  = fftconvolve(r_mon,  w_rk, mode='same')

# ── 2D migrated B-scans ────────────────────────────────────────────────────────
# Outer product of 1D depth trace with a lateral Gaussian (sigma_x = lam/2).
# This models the focused PSF of a Gazdag-migrated horizontal reflector.
# A Ricker wavelet in depth  x  Gaussian in x  =  the migrated image.
sigma_x = lam / 2
gauss_x = np.exp(-0.5 * ((x_mc - x_mc.mean()) / sigma_x) ** 2)

base_bscan = np.outer(trace_base, gauss_x)   # shape (Nz, Nx)
mon_bscan  = np.outer(trace_mon,  gauss_x)

# ── ROI: include full Ricker wavelet on both sides of fracture ─────────────────
n_margin = Nw // 2 + 5
iz_lo    = max(0,  iz0 - n_margin)
iz_hi    = min(Nz, iz0 + iz_del_water + n_margin + 1)
z_roi    = z_mc[iz_lo:iz_hi]

base_roi = base_bscan[iz_lo:iz_hi, :]
mon_roi  = mon_bscan[iz_lo:iz_hi, :]
diff_roi = mon_roi - base_roi

# ── Analytic signal along z ────────────────────────────────────────────────────
# Cross-spectra of real images are Hermitian: phi(-k)=-phi(k),
# forcing phi_0=0 in the WLS fit regardless of material change.
# Hilbert along z removes negative-kz components, enabling phi_0 != 0.
a_base = sp_hilbert(base_roi, axis=0)
a_mon  = sp_hilbert(mon_roi,  axis=0)

# ── Shift estimate ─────────────────────────────────────────────────────────────
dz_est, dx_est, XS, kz_ax, kx_ax = estimate_shift_2d(
    a_base, a_mon, dz_mig, dx_mig, kz_c
)

# Re-run WLS for phase intercept c[2]
KZ_f, KX_f = np.meshgrid(kz_ax, kx_ax, indexing='ij')
w_f    = np.abs(XS);  phi_f = np.angle(XS)
band_f = (np.abs(KZ_f) < 1.4 * kz_c) & (np.abs(KX_f) < 1.4 * kz_c)
mask_f = (w_f > 0.1 * w_f.max()) & band_f & ((np.abs(KZ_f) + np.abs(KX_f)) > 0)
A_f    = np.column_stack([KZ_f[mask_f], KX_f[mask_f], np.ones(mask_f.sum())])
W_f    = w_f[mask_f]
c_all  = np.linalg.lstsq(A_f * W_f[:, None], phi_f[mask_f] * W_f, rcond=None)[0]
phi_0  = c_all[2]

# ── Figure: 2x3 layout ────────────────────────────────────────────────────────
XS_s    = np.fft.fftshift(XS)
kz_s    = np.fft.fftshift(kz_ax)
kx_s    = np.fft.fftshift(kx_ax)
energy  = np.abs(XS_s)
phi_vis = np.where(energy > 0.01 * energy.max(),
                   np.degrees(np.angle(XS_s)), np.nan)

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle(
    f'Material Change Test — horizontal fracture  d = λ/20 = {d*1e3:.1f} mm\n'
    f'Baseline: Air (ε = {n_air**2:.0f})   →   Monitor: Water (ε = {n_water**2:.0f})'
    '   |   no physical movement',
    fontsize=11
)

x_cm = x_mc * 100;  z_cm = z_roi * 100
z_frac_top = z_mc[iz0] * 100
z_frac_bot = z_mc[iz0 + iz_del_water] * 100

def _bscan_ax(ax, img, title, vmax=None):
    vmax = vmax or np.max(np.abs(img))
    ax.pcolormesh(x_cm, z_cm, img,
                  cmap='RdBu_r', vmin=-vmax, vmax=vmax, shading='auto')
    ax.axhline(z_frac_top, color='gold',  lw=1.2, ls='--', label='fracture top')
    ax.axhline(z_frac_bot, color='olive', lw=1.2, ls=':',  label='fracture bot (H₂O)')
    ax.set_xlabel('x [cm]');  ax.set_ylabel('Depth [cm]')
    ax.invert_yaxis()
    ax.set_title(title, fontsize=9)

# (0,0) Baseline migrated B-scan
vmax_b = np.max(np.abs(base_roi))
_bscan_ax(axes[0, 0], base_roi,
          f'Baseline migrated B-scan\nAir fill  R_top = {R_top_air:+.3f}', vmax_b)
axes[0, 0].legend(fontsize=7, loc='lower right')

# (0,1) Monitor migrated B-scan (same colour scale for direct comparison)
_bscan_ax(axes[0, 1], mon_roi,
          f'Monitor migrated B-scan\nWater fill  R_top = {R_top_water:+.3f}', vmax_b)

# (0,2) Difference: monitor – baseline
vmax_d = np.max(np.abs(diff_roi))
_bscan_ax(axes[0, 2], diff_roi,
          'Difference  (Monitor − Baseline)\nReveals polarity & depth change', vmax_d)

# (1,0) Central trace overlay
ax = axes[1, 0]
norm_b = np.max(np.abs(trace_base));  norm_m = np.max(np.abs(trace_mon))
ax.plot(trace_base[iz_lo:iz_hi] / norm_b, z_cm,
        color='C0', lw=1.8, label=f'Baseline Air   (R_top={R_top_air:+.3f})')
ax.plot(trace_mon[iz_lo:iz_hi]  / norm_m, z_cm,
        color='C1', lw=1.8, ls='--',
        label=f'Monitor Water  (R_top={R_top_water:+.3f})')
ax.axhline(z_frac_top, color='gold',  lw=1.2, ls='--')
ax.axhline(z_frac_bot, color='olive', lw=1.2, ls=':')
ax.set_xlabel('Normalised amplitude');  ax.set_ylabel('Depth [cm]')
ax.set_title('Central trace (normalised)\nNote polarity reversal & pulse-shape change', fontsize=9)
ax.invert_yaxis();  ax.legend(fontsize=8);  ax.grid(alpha=0.2)

# (1,1) Cross-spectrum phase
ax2 = axes[1, 1]
im  = ax2.pcolormesh(kx_s, kz_s, phi_vis,
                     cmap='RdBu_r', vmin=-180, vmax=180, shading='auto')
klim = 1.4 * kz_c
for sgn in [-1, 1]:
    ax2.axhline(sgn * klim, color='k', lw=0.8, ls='--', alpha=0.6)
    ax2.axvline(sgn * klim, color='k', lw=0.8, ls='--', alpha=0.6)
ax2.set_title(
    'Analytic cross-spectrum phase [deg]\n'
    'Approx. flat ≈ constant polarity offset\n'
    'Small kz slope → apparent shift from Δv', fontsize=9
)
ax2.set_xlabel('kx [rad/m]');  ax2.set_ylabel('kz [rad/m]')
ax2.set_xlim(-2 * kz_c, 2 * kz_c);  ax2.set_ylim(-2 * kz_c, 2 * kz_c)
plt.colorbar(im, ax=ax2, fraction=0.046)

# (1,2) Results
ax3 = axes[1, 2]
ax3.axis('off')
txt = (
    'Fracture model:\n'
    f'  d          = λ/20 = {d*1e3:.2f} mm\n'
    f'  Baseline: Air    R_top = {R_top_air:+.3f}\n'
    f'  Monitor:  Water  R_top = {R_top_water:+.3f}\n\n'
    'Depth delay in migrated image:\n'
    f'  Air:   Δz_bot = {dz_del_air*1e3:.1f} mm  (< 1 sample)\n'
    f'  Water: Δz_bot = {dz_del_water*1e3:.1f} mm  (≈ {iz_del_water} samples)\n\n'
    'Estimated (true displacement = 0):\n'
    f'  Δz = {dz_est*1e3:+.1f} mm   [apparent shift from Δv]\n'
    f'  Δx = {dx_est*1e3:+.1f} mm   [correct: 0]\n\n'
    'Material change signature:\n'
    f'  c = {np.degrees(phi_0):+.1f}°   [polarity reversal + thin-layer phase]'
)
ax3.text(0.04, 0.97, txt, transform=ax3.transAxes,
         fontsize=10, family='monospace', verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.6))
ax3.set_title('Results', fontsize=9)

plt.tight_layout()
plt.show()

print(f'dz={dz_est*1e3:+.1f} mm  dx={dx_est*1e3:+.1f} mm  phi_0={np.degrees(phi_0):+.1f} deg')
"""

with open(nb_path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print('Updated cell 505c92c1 with 2x3 layout showing B-scans.')
