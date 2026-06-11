import numpy as np, h5py, time as _time
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from pathlib import Path

DATA = Path(r'C:\Users\Administrator\OneDrive\Thesis\TimeLapse_Notebooks\B_Scan_Migration_data')
eps_r=3.15; v_ice=0.299792458/np.sqrt(eps_r); v_mig=v_ice/2; t0_ns=0.953
dx_grid=7.454e-4; x_survey_start=0.5; trace_spacing=15*dx_grid; n_traces=269
scatterer_spacing_m=0.20; x_s1=2.0; x_s2=2.20; z_scat=0.5
wavelength=v_ice/1.5; f_c_GHz=1.5; rayleigh_m=wavelength/2
x_traces=x_survey_start+np.arange(n_traces)*trace_spacing
x_img=np.linspace(0.5,3.5,300); z_img=np.linspace(0.0,0.8,120)

def load(path, comp='Ey'):
    with h5py.File(path,'r') as f:
        dt_s=float(f.attrs['dt'])
        rxs=sorted(f['rxs'].keys(), key=lambda s: int(s.replace('rx','')))
        bs=np.stack([f['rxs'][r][comp][:] for r in rxs], axis=0)
    return bs, dt_s*1e9

bscan_full, dt_ns = load(DATA/'psf_bscan_merged.out')
bscan_bg, _       = load(DATA/'bg_bscan_merged.out')
bscan_sc = bscan_full - bscan_bg
print('Data: {}, dt={:.5f} ns'.format(bscan_sc.shape, dt_ns))


def kirchhoff(bscan, dt_ns, x_tr, v_mig, x_img, z_img, t0_ns):
    n_tr, n_t = bscan.shape
    tr_idx = np.arange(n_tr)
    img = np.zeros((len(z_img), len(x_img)))
    for iz, z in enumerate(z_img):
        if z == 0:
            continue
        dx = x_img[:, None] - x_tr[None, :]
        r  = np.sqrt(dx**2 + z**2)
        t  = r / v_mig + t0_ns
        i_frac = t / dt_ns
        i0 = np.clip(i_frac.astype(int), 0, n_t - 2)
        frac = i_frac - i0
        samp = ((1 - frac) * bscan[tr_idx[None, :], i0]
               + frac * bscan[tr_idx[None, :], np.minimum(i0 + 1, n_t - 1)])
        img[iz] = samp.sum(axis=1)
    return img / (np.max(np.abs(img)) + 1e-30)


def gazdag(bscan, dt_ns, x_tr, v_mig, x_img, z_img, t0_ns, f_min=0.3, f_max=3.0):
    n_tr, n_t = bscan.shape
    dx_tr = float(x_tr[1] - x_tr[0])
    n_pad = 2 * n_tr
    tap = np.ones(n_tr)
    nb = max(1, int(0.05 * n_tr))
    tap[:nb]  = 0.5 * (1 - np.cos(np.pi * np.arange(nb) / nb))
    tap[-nb:] = tap[:nb][::-1]
    D = np.fft.fft2(bscan * tap[:, None], s=(n_pad, n_t))
    f_arr  = np.fft.fftfreq(n_t, dt_ns)
    omega  = 2 * np.pi * f_arr
    kx     = 2 * np.pi * np.fft.fftfreq(n_pad, dx_tr)
    band   = (np.abs(f_arr) >= f_min) & (np.abs(f_arr) <= f_max)
    D[:, ~band] = 0.0
    print('  Bandpass: {}/{} bins ({}-{} GHz)'.format(band.sum(), n_t, f_min, f_max))
    D *= np.exp(1j * omega[None, :] * t0_ns)
    kz2 = (omega[None, :] / v_mig)**2 - kx[:, None]**2
    kz  = np.where(kz2 > 0, np.sqrt(np.maximum(kz2, 0.0)), 0.0)
    D[:, ~band] = 0.0
    D[kz2 <= 0] = 0.0
    img = np.zeros((len(z_img), n_pad))
    z_prev = 0.0
    for iz, z in enumerate(z_img):
        dz = z - z_prev
        if dz > 0:
            D *= np.exp(1j * kz * dz)
            D[kz2 <= 0] = 0.0
        D_x = np.fft.ifft(D, axis=0)
        img[iz] = np.real(D_x.sum(axis=1))
        z_prev = z
    x_fft   = x_tr[0] + np.arange(n_pad) * dx_tr
    img_out = np.zeros((len(z_img), len(x_img)))
    for iz in range(len(z_img)):
        img_out[iz] = np.interp(x_img, x_fft, img[iz])
    return img_out / (np.max(np.abs(img_out)) + 1e-30)


def backprop(bscan, dt_ns, x_tr, v, x_img, z_img, t0_ns):
    n_tr, n_t = bscan.shape
    t0_samp = int(round(t0_ns / dt_ns))
    tr_idx = np.arange(n_tr)
    bscan_sh = np.roll(bscan, -t0_samp, axis=1)
    bscan_sh[:, -t0_samp:] = 0.0
    D_c = hilbert(bscan_sh, axis=1)
    img = np.zeros((len(z_img), len(x_img)))
    for iz, z in enumerate(z_img):
        if z == 0:
            continue
        dx = x_img[:, None] - x_tr[None, :]
        r  = np.sqrt(dx**2 + z**2)
        t2 = 2.0 * r / v
        i_frac = t2 / dt_ns
        i0   = np.clip(i_frac.astype(int), 0, n_t - 2)
        frac = i_frac - i0
        samp = ((1 - frac) * D_c[tr_idx[None, :], i0]
               + frac * D_c[tr_idx[None, :], np.minimum(i0 + 1, n_t - 1)])
        img[iz] = np.abs(samp.sum(axis=1))
    return img / (np.max(img) + 1e-30)


t0 = _time.perf_counter()
img_kh = kirchhoff(bscan_sc, dt_ns, x_traces, v_mig, x_img, z_img, t0_ns)
print('Kirchhoff: {:.1f}s'.format(_time.perf_counter() - t0))

t0 = _time.perf_counter()
img_ps = gazdag(bscan_sc, dt_ns, x_traces, v_mig, x_img, z_img, t0_ns)
print('Gazdag: {:.1f}s'.format(_time.perf_counter() - t0))

t0 = _time.perf_counter()
img_bp = backprop(bscan_sc, dt_ns, x_traces, v_ice, x_img, z_img, t0_ns)
print('Backprop: {:.1f}s'.format(_time.perf_counter() - t0))


def hilbert_env(img):
    return np.abs(hilbert(img / (np.max(np.abs(img)) + 1e-30), axis=1))

def fwhm(x, p):
    half  = np.max(p) / 2
    above = p >= half
    c = np.where(np.diff(above.astype(int)))[0]
    return (x[c[-1]] - x[c[0]]) if len(c) >= 2 else float('nan')

def peak_depth(env, z_img, z0, hw=0.10):
    win = (z_img >= z0 - hw) & (z_img <= z0 + hw)
    return np.where(win)[0][np.argmax(np.max(env[win], axis=1))]


# Hilbert envelopes
env_kh = hilbert_env(img_kh)
env_ps = hilbert_env(img_ps)
# Backprop is already magnitude-based — normalise directly
env_bp = img_bp / (np.max(img_bp) + 1e-30)

# Per-method peak depth in both raw and envelope (use envelope for depth finding)
iz_kh = peak_depth(env_kh, z_img, z_scat)
iz_ps = peak_depth(env_ps, z_img, z_scat)
iz_bp = peak_depth(env_bp, z_img, z_scat)

# 1D slices — envelope
psf_env_kh = env_kh[iz_kh] / (np.max(env_kh[iz_kh]) + 1e-30)
psf_env_ps = env_ps[iz_ps] / (np.max(env_ps[iz_ps]) + 1e-30)
psf_env_bp = env_bp[iz_bp] / (np.max(env_bp[iz_bp]) + 1e-30)

# 1D slices — raw (signed, normalised to peak abs)
raw_kh = img_kh[iz_kh] / (np.max(np.abs(img_kh[iz_kh])) + 1e-30)
raw_ps = img_ps[iz_ps] / (np.max(np.abs(img_ps[iz_ps])) + 1e-30)
raw_bp = img_bp[iz_bp] / (np.max(img_bp[iz_bp]) + 1e-30)  # already >= 0

fw_env_kh = fwhm(x_img, psf_env_kh)
fw_env_ps = fwhm(x_img, psf_env_ps)
fw_env_bp = fwhm(x_img, psf_env_bp)
fw_raw_kh = fwhm(x_img, np.abs(raw_kh))
fw_raw_ps = fwhm(x_img, np.abs(raw_ps))
fw_raw_bp = fwhm(x_img, raw_bp)

print('Peak depths: KH={:.4f}m  PS={:.4f}m  BP={:.4f}m'.format(
    z_img[iz_kh], z_img[iz_ps], z_img[iz_bp]))
print('Envelope FWHM  -> KH={:.1f}mm  PS={:.1f}mm  BP={:.1f}mm'.format(
    fw_env_kh*1e3, fw_env_ps*1e3, fw_env_bp*1e3))
print('Raw |PSF| FWHM -> KH={:.1f}mm  PS={:.1f}mm  BP={:.1f}mm'.format(
    fw_raw_kh*1e3, fw_raw_ps*1e3, fw_raw_bp*1e3))

# ── Figure: 3 rows ─────────────────────────────────────────────────────────────
# Row 1: raw 2D migrated images (seismic, signed)
# Row 2: Hilbert envelope 2D images (hot_r)
# Row 3 (split): raw 1D PSF  |  envelope 1D PSF
# ──────────────────────────────────────────────────────────────────────────────
titles  = ['Kirchhoff', 'Gazdag (Phase-Shift)', 'Backpropagation']
raws    = [img_kh,  img_ps,  img_bp]
envs    = [env_kh,  env_ps,  env_bp]
iz_list = [iz_kh,   iz_ps,   iz_bp]
colors  = ['tab:blue', 'tab:red', 'tab:green']
ls_list = ['-', '--', '-.']

fig = plt.figure(figsize=(18, 14))
gs  = fig.add_gridspec(3, 3, hspace=0.38, wspace=0.28,
                        height_ratios=[1, 1, 0.85])

ax_raw1d = fig.add_subplot(gs[2, :2])   # raw PSF (left 2/3)
ax_env1d = fig.add_subplot(gs[2, 2])    # envelope PSF (right 1/3)

for col, (title, raw, env, iz_pk, color, ls) in enumerate(
        zip(titles, raws, envs, iz_list, colors, ls_list)):

    # --- Row 1: raw signed image ---
    ax_r = fig.add_subplot(gs[0, col])
    vmax_r = np.percentile(np.abs(raw), 98)
    ax_r.imshow(raw, aspect='auto', cmap='seismic',
                vmin=-vmax_r, vmax=vmax_r,
                extent=[x_img[0], x_img[-1], z_img[-1], z_img[0]], origin='upper')
    ax_r.axhline(z_img[iz_pk], color='lime', ls='--', lw=0.8)
    ax_r.plot([x_s1, x_s2], [z_scat, z_scat], 'y*', ms=8, zorder=5)
    ax_r.set_title('Raw: ' + title, fontsize=10)
    ax_r.set_xlabel('x [m]'); ax_r.set_ylabel('Depth [m]')

    # --- Row 2: envelope image ---
    ax_e = fig.add_subplot(gs[1, col])
    ax_e.imshow(env, aspect='auto', cmap='hot_r', vmin=0, vmax=1,
                extent=[x_img[0], x_img[-1], z_img[-1], z_img[0]], origin='upper')
    ax_e.axhline(z_img[iz_pk], color='cyan', ls='--', lw=0.8)
    ax_e.plot([x_s1, x_s2], [z_scat, z_scat], 'c*', ms=8, zorder=5)
    ax_e.set_title('Envelope: ' + title, fontsize=10)
    ax_e.set_xlabel('x [m]'); ax_e.set_ylabel('Depth [m]')

    # --- Row 3 left: raw 1D PSF (absolute value, to compare peaks) ---
    raw_1d   = raws[col][iz_pk]
    raw_norm = raw_1d / (np.max(np.abs(raw_1d)) + 1e-30)
    fw_r = fwhm(x_img, np.abs(raw_norm))
    fw_r_str = '{:.1f} mm'.format(fw_r * 1e3) if not np.isnan(fw_r) else 'N/A'
    ax_raw1d.plot(x_img, raw_norm, color=color, ls=ls, lw=1.5, alpha=0.85,
                  label='{} (FWHM={})'.format(title, fw_r_str))

    # --- Row 3 right: envelope 1D PSF ---
    env_1d   = envs[col][iz_pk]
    env_norm = env_1d / (np.max(env_1d) + 1e-30)
    fw_e = fwhm(x_img, env_norm)
    fw_e_str = '{:.1f} mm'.format(fw_e * 1e3) if not np.isnan(fw_e) else 'N/A'
    ax_env1d.plot(x_img, env_norm, color=color, ls=ls, lw=1.5,
                  label='{} {}'.format(title, fw_e_str))

# Raw 1D panel formatting
ax_raw1d.axvline(x_s1, color='k', ls=':', lw=0.9, label='True scatterer x')
ax_raw1d.axvline(x_s2, color='k', ls=':', lw=0.9)
ax_raw1d.axhline(0.0,  color='k',    ls='-',  lw=0.4, alpha=0.3)
ax_raw1d.axhline(0.5,  color='gray', ls='--', lw=0.7, label='+0.5 level')
ax_raw1d.axhline(-0.5, color='gray', ls=':',  lw=0.7, label='-0.5 level')
ax_raw1d.set_xlim(x_s1 - 0.35, x_s2 + 0.35)
ax_raw1d.set_ylim(-1.05, 1.05)
ax_raw1d.set_xlabel('x [m]'); ax_raw1d.set_ylabel('Norm. amplitude (signed)')
ax_raw1d.set_title(
    'Raw migrated PSF at focus depth  |  Dx={:.0f}mm={:.2f}lam'.format(
        scatterer_spacing_m * 1e3, scatterer_spacing_m / wavelength),
    fontsize=10)
ax_raw1d.legend(fontsize=8, loc='upper right')

# Envelope 1D panel formatting
ax_env1d.axvline(x_s1, color='k', ls=':', lw=0.9)
ax_env1d.axvline(x_s2, color='k', ls=':', lw=0.9)
ax_env1d.axhline(0.5,  color='gray', ls='--', lw=0.7, label='-6 dB level')
ax_env1d.set_xlim(x_s1 - 0.35, x_s2 + 0.35)
ax_env1d.set_ylim(-0.05, 1.1)
ax_env1d.set_xlabel('x [m]'); ax_env1d.set_ylabel('Norm. envelope')
ax_env1d.set_title('Hilbert envelope PSF\nRayleigh={:.1f}mm'.format(rayleigh_m * 1e3), fontsize=10)
ax_env1d.legend(fontsize=7, loc='upper right')

fig.suptitle(
    'Zero-Offset B-Scan Migration  |  v_ice={:.4f} m/ns  |  f_c={:.1f} GHz  |'
    '  Dx={:.0f}mm={:.2f}lam'.format(v_ice, f_c_GHz,
                                      scatterer_spacing_m * 1e3,
                                      scatterer_spacing_m / wavelength),
    fontsize=12
)

out_png = DATA / 'psf_comparison_200mm.png'
plt.savefig(out_png, dpi=150, bbox_inches='tight')
print('Saved {}'.format(out_png))
