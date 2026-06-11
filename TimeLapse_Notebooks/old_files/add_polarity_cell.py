import json, sys
sys.stdout.reconfigure(encoding='utf-8')

nb_path = r'c:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\R_PhasePlaneFit.ipynb'
with open(nb_path, encoding='utf-8') as f:
    nb = json.load(f)

md_cell = {
    "cell_type": "markdown",
    "id": "e2f4b8a1",
    "metadata": {},
    "source": "# Two-Step Polarity Retrieval"
}

code_src = '''\
from scipy.signal import hilbert as sp_hilbert

# Reuses from Material Change Test cell:
#   lam, dz_mig, dx_mig, kz_c, x_mc, sigma_x, iz_lo, iz_hi
#   trace_base (air fill), trace_mon (water fill)
#   R_top_air, R_top_water, estimate_shift_2d

# ── Scenario: fill change (air->water) + lateral shift delta_x = lam/4 ──────────
delta_x_p = lam / 4                 # shift = 1 spatial sample (dx_mig)
gauss_ctr = x_mc.mean()

gauss_b_p = np.exp(-0.5 * ((x_mc - gauss_ctr            ) / sigma_x)**2)
gauss_m_p = np.exp(-0.5 * ((x_mc - gauss_ctr - delta_x_p) / sigma_x)**2)

base_p = np.outer(trace_base, gauss_b_p)   # air fill, centred
mon_p  = np.outer(trace_mon,  gauss_m_p)   # water fill, shifted right

roi_b_p = base_p[iz_lo:iz_hi, :]
roi_m_p = mon_p [iz_lo:iz_hi, :]

a_b_p = sp_hilbert(roi_b_p, axis=0)
a_m_p = sp_hilbert(roi_m_p, axis=0)

# ── Step 1: joint phase-plane fit  ->  Delta_z, Delta_x, phi_0 (naive) ───────────
dz_p, dx_p, XS_p, kz_p, kx_p = estimate_shift_2d(a_b_p, a_m_p, dz_mig, dx_mig, kz_c)

KZ_p, KX_p = np.meshgrid(kz_p, kx_p, indexing='ij')
w_p    = np.abs(XS_p);  phi_p = np.angle(XS_p)
band_p = (np.abs(KZ_p) < 1.4*kz_c) & (np.abs(KX_p) < 1.4*kz_c)
mask_p = (w_p > 0.1*w_p.max()) & band_p & ((np.abs(KZ_p)+np.abs(KX_p)) > 0)
A_p    = np.column_stack([KZ_p[mask_p], KX_p[mask_p], np.ones(mask_p.sum())])
W_p    = w_p[mask_p]
c_p    = np.linalg.lstsq(A_p*W_p[:,None], phi_p[mask_p]*W_p, rcond=None)[0]
phi_0_s1 = c_p[2]

# ── Step 2: divide out the estimated shift ramp -> isolate phi_0 ──────────────────
# Multiply XS by exp(-i*(kz*Dz + kx*Dx)) to remove the displacement phase.
# What remains is the constant polarity offset phi_0.
XS_corr  = XS_p * np.exp(-1j * (KZ_p * dz_p + KX_p * dx_p))
phi_corr = np.angle(XS_corr)
phi_0_s2 = np.average(phi_corr[mask_p], weights=W_p)

# ── Figure: 1x3 ───────────────────────────────────────────────────────────────────
kz_s = np.fft.fftshift(kz_p)
kx_s = np.fft.fftshift(kx_p)
energy_s = np.fft.fftshift(np.abs(XS_p))
thresh   = 0.01 * energy_s.max()

vis_raw  = np.where(energy_s > thresh,
                    np.degrees(np.angle(np.fft.fftshift(XS_p))),    np.nan)
vis_corr = np.where(energy_s > thresh,
                    np.degrees(np.angle(np.fft.fftshift(XS_corr))), np.nan)

fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=120)
fig.suptitle(
    'Two-step polarity retrieval  \\u2014  fill change (air\\u2192water, \\u03c6\\u2080 = \\u2212180\\u00b0) '
    '+ lateral shift \\u03b4x = \\u03bb/4\\n'
    'Step 1: joint phase-plane fit extracts shift  |  '
    'Step 2: remove shift ramp \\u2192 flat phase \\u2248 \\u2212180\\u00b0',
    fontsize=10
)

kw   = dict(cmap='RdBu_r', vmin=-180, vmax=180, shading='auto')
klim = 1.4 * kz_c

for ax, vis, ttl in [
    (axes[0], vis_raw,
     'Step 1: raw cross-spectrum phase\\n'
     'Tilted pattern: kx slope from shift,  offset from fill change'),
    (axes[1], vis_corr,
     f'Step 2: shift ramp removed  (\\u0394x = {dx_p*1e3:.1f} mm)\\n'
     'Flat pattern \\u2248 \\u2212180\\u00b0  \\u2192  polarity reversal confirmed'),
]:
    im = ax.pcolormesh(kx_s, kz_s, vis, **kw)
    for sgn in [-1, 1]:
        ax.axhline(sgn * klim, color='k', lw=0.8, ls='--', alpha=0.5)
        ax.axvline(sgn * klim, color='k', lw=0.8, ls='--', alpha=0.5)
    ax.set_xlim(-2*kz_c, 2*kz_c);  ax.set_ylim(-2*kz_c, 2*kz_c)
    ax.set_xlabel('kx [rad/m]');    ax.set_ylabel('kz [rad/m]')
    ax.set_title(ttl, fontsize=9)
    plt.colorbar(im, ax=ax, fraction=0.046, label='Phase [deg]')

ax3 = axes[2];  ax3.axis('off')
txt = (
    'Scenario:\\n'
    f'  Baseline: air fill    R_top = {R_top_air:+.3f}\\n'
    f'  Monitor:  water fill  R_top = {R_top_water:+.3f}\\n'
    f'            + shift \\u03b4x = \\u03bb/4 = {delta_x_p*1e3:.1f} mm\\n\\n'
    'Step 1 \\u2014 joint phase-plane fit:\\n'
    f'  \\u0394z = {dz_p*1e3:+.1f} mm   [true: 0]\\n'
    f'  \\u0394x = {dx_p*1e3:+.1f} mm   [true: {delta_x_p*1e3:.1f}]\\n'
    f'  \\u03c6\\u2080 = {np.degrees(phi_0_s1):+.1f}\\u00b0   (simultaneous fit)\\n\\n'
    'Step 2 \\u2014 remove shift ramp from XS:\\n'
    f'  \\u03c6\\u2080 = {np.degrees(phi_0_s2):+.1f}\\u00b0\\n\\n'
    'Discriminator:\\n'
    '  |\\u03c6\\u2080| \\u2248 180\\u00b0  \\u2192  fill change\\n'
    '  |\\u03c6\\u2080| \\u2248   0\\u00b0  \\u2192  pure displacement\\n'
    '  (cf. water-front test: \\u03c6\\u2080 \\u2248 +0.6\\u00b0)'
)
ax3.text(0.04, 0.97, txt, transform=ax3.transAxes,
         fontsize=10, family='monospace', verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.6))
ax3.set_title('Results', fontsize=9)

plt.tight_layout()
plt.show()

print(f'Step 1:  \\u0394x = {dx_p*1e3:.1f} mm   \\u03c6\\u2080 = {np.degrees(phi_0_s1):.1f}\\u00b0')
print(f'Step 2:  \\u03c6\\u2080 = {np.degrees(phi_0_s2):.1f}\\u00b0   (after shift correction)')
'''

code_cell = {
    "cell_type": "code",
    "execution_count": None,
    "id": "f3a9d1c7",
    "metadata": {},
    "outputs": [],
    "source": code_src
}

# Insert after the Moving Water Front Test code cell (d9f5c2a8)
idx = next(i for i, c in enumerate(nb['cells']) if c.get('id') == 'd9f5c2a8')
nb['cells'].insert(idx + 1, md_cell)
nb['cells'].insert(idx + 2, code_cell)

with open(nb_path, 'w', encoding='utf-8') as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f'Inserted cells e2f4b8a1 (md) and f3a9d1c7 (code) after d9f5c2a8.')
print(f'Notebook now has {len(nb["cells"])} cells.')
