"""
Migration functions shared by Resolution_Playground and TimeLapse_Playground notebooks.
"""

import numpy as np
import pylops


def PylopsKirchoffMigration(data, t, x, vel_model, z, wav='Ricker', wavcenter='Center',
                             f0=0.08, recs=None, srcs=None, mode='analytic', dynamic=False,
                             aperture=None, angleaperture=None, return_op=False):
    """
    Perform Kirchhoff migration using PyLops zero offset operator. pass ns, gHz and meter
    """
    import importlib
    import helper_functions.KirchhoffPylopsZeroOffset as KirchhoffPylopsZeroOffset
    importlib.reload(KirchhoffPylopsZeroOffset)
    from pylops.utils.wavelets import ricker

    # For zero-offset, srcs = recs
    rx = x
    rz = np.zeros_like(rx)
    recs = np.vstack((rx, rz))
    srcs = recs

    # Handle wavelet selection
    if wav == 'Ricker':
        n_cycles = 6
        period_ns = 1.0 / f0
        wav_duration = n_cycles * period_ns
        n_wav_samples = int(np.ceil(wav_duration / (t[1] - t[0])))
        if n_wav_samples % 2 == 0:
            n_wav_samples += 1
        wav_arr, _, wcenter = ricker(t[0:n_wav_samples], f0=f0)
    else:
        wav_arr = wav

    if wavcenter == 'Center':
        wavcenter = wcenter

    K = KirchhoffPylopsZeroOffset.Kirchhoff(
        z=z, x=x, t=t,
        srcs=srcs, recs=recs,
        vel=vel_model,
        wav=wav_arr, wavcenter=wavcenter,
        mode=mode, dynamic=dynamic,
        aperture=aperture, angleaperture=angleaperture
    )

    d = data.flatten()
    print("d.shape:", d.shape, "K dimension:", K.shape)

    m = K.H @ d
    m = m.reshape(len(x), len(z)).T

    if return_op:
        return m, K
    return m


def gazdag_migration(data, x, t, z, vel):
    """
    Gazdag (1978) phase-shift migration for zero-offset post-stack data.

    data : (n_t, n_x)  B-scan, time axis first, t0-shifted
    x, t : 1-D arrays [m, ns]
    z    : 1-D depth axis [m], uniformly spaced from 0
    vel  : full medium velocity [m/ns]; v_mig = vel/2 used internally

    Returns image (n_z, n_x), peak-normalised to ±1.
    """
    n_t, n_x = data.shape
    dt    = float(t[1] - t[0])
    dx    = float(x[1] - x[0])
    dz    = float(z[1] - z[0])
    v_mig = vel / 2.0

    # Spatial cosine taper (5%) + 100% zero-pad each side → nx_pad = 3*n_x
    n_xtap  = max(3, int(0.05 * n_x))
    pad     = n_x
    sp_tap  = pylops.utils.tapers.taper2d(n_t, n_x, n_xtap)   # (n_x, n_t)
    data_sp = np.pad(
        data.T * sp_tap,   # (n_x, n_t)
        ((pad, pad), (0, 0)),
        mode='constant'
    )                      # (nx_pad, n_t)
    nx_pad = n_x + 2 * pad

    # f-kx dip filter: zero evanescent bins before migration to prevent smile artefacts
    D_fk       = np.fft.fft(np.fft.rfft(data_sp, axis=1), axis=0)   # (nx_pad, n_f)
    freq_arr   = np.fft.rfftfreq(n_t, dt)          # [GHz]
    kx_arr     = np.fft.fftfreq(nx_pad, dx)        # [1/m]
    evanescent = np.abs(kx_arr[:, None]) > np.abs(freq_arr[None, :]) / v_mig
    D_fk[evanescent] = 0.0
    data_sp = np.real(np.fft.irfft(np.fft.ifft(D_fk, axis=0), n=n_t, axis=1))

    print(f'  nx={n_x} → nx_pad={nx_pad},  evanescent bins zeroed: '
          f'{evanescent.sum()}/{evanescent.size} ({100*evanescent.mean():.1f}%)')

    # PhaseShift depth loop: downward continuation + t=0 imaging condition
    freq  = np.fft.rfftfreq(n_t, dt)
    kx    = np.fft.fftshift(np.fft.fftfreq(nx_pad, dx))
    Pop   = pylops.waveeqprocessing.PhaseShift(v_mig, dz, n_t, freq, kx)
    field = data_sp.T.ravel()   # (n_t * nx_pad,) row-major

    image = np.zeros((len(z), n_x))
    for iz in range(len(z)):
        image[iz] = np.real(field.reshape(n_t, nx_pad)[0, pad:-pad])
        if iz < len(z) - 1:
            field = Pop.H * field   # adjoint PhaseShift = downward continuation

    peak = np.max(np.abs(image))
    return image / (peak + 1e-30)


def write_backprop_files(study_root, label, slug, tapered_ntr_nt, dt_ns, x_midpoints,
                         t0_ns, eps_r, v_ice, stride=1, n_snap=30, snap_win=1.0,
                         sign_bit=False):
    """
    Write gprMax excitation file and .in file for back-propagation migration.
    Files go in study_root/backprop/<slug>/.

    Sources are placed at the Tx-Rx midpoints (x_midpoints), consistent with the
    exploding-reflector half-velocity (v_mig = v/2) assumption.

    sign_bit=True — inject sign(u) instead of the peak-normalised time-reversed
    wavefield. Spatial focusing during back-propagation is governed by phase
    (zero-crossings), not amplitude; sign-bit reversal keeps every phase trend
    of the wavelet but squashes impulsive noise spikes down to the same +-1 as
    the coherent signal, stripping them of the outsized amplitude they'd
    otherwise inject as competing point sources. See TimeLapse_Processing.ipynb.

    Returns (in_path, n_src, n_snaps, t_focus_ns).
    """
    import pathlib
    study_root = pathlib.Path(study_root)
    out_dir = study_root / 'backprop' / slug
    out_dir.mkdir(parents=True, exist_ok=True)

    _, n_t = tapered_ntr_nt.shape
    dt_s  = dt_ns * 1e-9
    T_ns  = n_t * dt_ns

    # Stride, time-reverse, then either sign-bit or peak-normalise per trace
    data_s   = tapered_ntr_nt[::stride]
    x_src    = x_midpoints[::stride]
    n_src    = len(x_src)
    data_rev = data_s[:, ::-1].copy()
    if sign_bit:
        data_rev = np.sign(data_rev)
    else:
        peak     = np.max(np.abs(data_rev), axis=1, keepdims=True)
        peak[peak == 0] = 1.0
        data_rev /= peak

    # Snapshot timing: all depths focus simultaneously at t_focus = T - t0
    t_focus_ns = T_ns - t0_ns
    t_start_ns = max(0.0, t_focus_ns - snap_win)
    t_start_s  = t_start_ns * 1e-9
    snap_step  = max(1, int((T_ns * 1e-9 - t_start_s) / (max(1, n_snap - 1) * dt_s)))
    n_snaps    = max(1, int((n_t * dt_s - t_start_s) / (snap_step * dt_s)))

    # Excitation file: rows = time steps, cols = [time, src_0, ..., src_n]
    exc_path = out_dir / 'excitation.txt'
    time_s   = np.arange(n_t) * dt_s
    N_PAD    = 10
    pad_t    = np.arange(n_t, n_t + N_PAD) * dt_s
    exc_arr  = np.vstack([
        np.column_stack([time_s, data_rev.T]),
        np.column_stack([pad_t,  np.zeros((N_PAD, n_src))]),
    ])
    header = 'time ' + ' '.join(f'bp_{i}' for i in range(n_src))
    np.savetxt(exc_path, exc_arr, fmt='%.6e', header=header, comments='')

    # Source x-positions (loaded at gprMax runtime)
    npy_path = out_dir / 'src_positions.npy'
    np.save(npy_path, x_src)

    # Write .in file
    eps_r_half = 4.0 * eps_r
    v_half     = v_ice / 2
    in_path    = out_dir / f'backprop_{slug}.in'

    excitation_mode = 'sign-bit (sign(u), amplitude stripped)' if sign_bit else 'peak-normalised'
    in_lines = [
        f'#title: Back-Propagation -- {label}',
        f'// Excitation mode: {excitation_mode}',
        '#domain: 4.000 1.000 0.001',
        '#dx_dy_dz: 0.001 0.001 0.001',
        f'#time_window: {n_t * dt_s:.6e}',
        '#pml_cells: 10 10 0 10 10 0',
        '',
        f'// Half-velocity: eps_r={eps_r_half:.2f} (=4x{eps_r}) -> v={v_half:.5f} m/ns',
        f'#material: {eps_r_half:.2f} 1e-6 1.0 0 ice',
        '',
        f'#box: 0 0 0 4.000 0.900 0.001 ice',
        '',
        f'#excitation_file: {exc_path.name}',
        '',
        f'// {n_src} time-reversed sources at Tx-Rx midpoints (zero-offset convention, stride={stride})',
        '#python:',
        'from gprMax.input_cmd_funcs import *',
        'import numpy as np',
        f"x_sources = np.load(r'{npy_path}')",
        'for i, x in enumerate(x_sources):',
        "    hertzian_dipole('z', float(x), 0.900, 0.0, 'bp_{}'.format(i))",
        '#end_python:',
        '',
        f'// {n_snaps} snapshots  window=[{t_start_ns:.1f}, {T_ns:.1f}] ns  focus at {t_focus_ns:.2f} ns',
        '#python:',
        'from gprMax.input_cmd_funcs import *',
        f't_start = {t_start_s:.8e}',
        f'step    = {snap_step}',
        f'dt_gpr  = {dt_s:.8e}',
        f'n_snap  = {n_snaps}',
        'for k in range(n_snap):',
        '    t_snap = t_start + k * step * dt_gpr',
        "    print('#snapshot: 0 0 0 4.000 1.000 0.001 0.001 0.001 0.001 %.8e bp_snap%04d' % (t_snap, k+1))",
        '#end_python:',
        '',
        '#messages: y',
    ]
    in_path.write_text('\n'.join(in_lines) + '\n', encoding='utf-8')
    return in_path, n_src, n_snaps, t_focus_ns


def dispersion_limited_cutoff(eps_r, dx, min_cells_per_wavelength=3, safety_factor=0.7):
    """
    Highest frequency (Hz) that a gprMax grid can resolve in a material,
    with a safety margin below the solver's own numerical-dispersion limit.

    gprMax refuses to run (GeneralError: "Non-physical wave propagation")
    whenever a source's spectral content samples the medium's wavelength
    with fewer than `min_cells_per_wavelength` cells (gprMax's own default
    is 3 -- see gprMax.grid.FDTDGrid.mingridsampling and dispersion_analysis()
    in gprMax/grid.py). For #excitation_file sources this is evaluated per
    source column, so a single noisy/sharp trace among hundreds can trip the
    check even though most of the model is well resolved.

    safety_factor < 1 leaves headroom for two reasons: (1) gprMax's check
    measures where the spectrum drops 40 dB below its peak, not a hard
    band-edge, and (2) the Butterworth filter applied in
    lowpass_filter_excitation() has a gradual roll-off, so its actual -40 dB
    point sits noticeably above the nominal cutoff passed to it. The default
    of 0.7 was tuned empirically against gprMax's own dispersion_analysis()
    for this study's eps_r=12.6, dx=1mm grid (nominal cutoff ~20 GHz landed
    the real post-filter significant frequency around ~22-23 GHz, safely
    under the ~28 GHz hard limit). Re-validate with a real gprMax run if you
    change eps_r, dx, or the filter order substantially.

    Args:
        eps_r (float): Relative permittivity of the fastest-attenuating /
            highest-permittivity material in the model (sets the minimum
            wave velocity, which is the worst case for dispersion).
        dx (float): Spatial step (m). Assumes a cubic/uniform grid as used
            by write_backprop_files (#dx_dy_dz all equal).
        min_cells_per_wavelength (int): gprMax's mingridsampling threshold.
        safety_factor (float): Fraction of the theoretical limit to target.

    Returns:
        cutoff_hz (float): Suggested low-pass cutoff frequency in Hz.
    """
    v = 299792458.0 / np.sqrt(eps_r)
    return safety_factor * v / (min_cells_per_wavelength * dx)


def lowpass_filter_excitation(exc_path, cutoff_hz, order=8, edge_exclude=0, keep_backup=True):
    """
    Zero-phase low-pass filter every source column of a gprMax excitation
    file in place, to bring its spectral content under gprMax's
    numerical-dispersion limit (see dispersion_limited_cutoff).

    Needed when write_backprop_files() produces virtual sources with
    sharper/sub-wavelength focusing kernels than the simulation grid can
    resolve -- gprMax then raises "Non-physical wave propagation" for that
    .in file instead of running. Filtering trims only the unresolvable
    high-frequency tail; the time column and overall pulse shape/timing are
    left intact.

    Filtering alone is not always enough: the outermost virtual sources sit
    at the edge of the migration aperture and can carry near-Nyquist content
    baked into the raw time-reversed data (observed: a single-sample swing of
    ~0.95 out of a +-1 range, i.e. close to the simulation's own sample
    rate). No low-pass filter can safely remove that -- too gentle and it
    survives, too aggressive and the filter's own ringing on such a sharp
    signal makes the measured significant frequency *worse*, not better
    (observed maxfreq climbing past 100GHz at very low cutoffs). Dropping a
    few edge sources entirely via `edge_exclude` (standard aperture-limiting
    practice in migration) is the reliable fix for that failure mode; the
    frequency-domain filter alone cannot resolve it.

    Args:
        exc_path (str or Path): Path to the excitation.txt file (format:
            header row 'time bp_0 bp_1 ...', then one row per time step,
            as written by write_backprop_files).
        cutoff_hz (float): Low-pass cutoff frequency in Hz. Use
            dispersion_limited_cutoff() to derive this from the model's
            material and grid spacing.
        order (int): Butterworth filter order (higher = sharper roll-off,
            closer to the nominal cutoff at the cost of more ringing).
        edge_exclude (int): Zero out this many source columns at each end
            of the array (outermost sources of the migration aperture)
            before filtering. 0 disables this (default).
        keep_backup (bool): If True, save the pre-filter file as
            '<exc_path>' with 'excitation' replaced by 'excitation_orig'
            (skipped if that backup already exists).

    Returns:
        exc_path (Path): The (now filtered) excitation file path.
    """
    import pathlib
    import shutil
    from scipy.signal import butter, sosfiltfilt

    exc_path = pathlib.Path(exc_path)
    with open(exc_path) as f:
        header = f.readline()

    data = np.loadtxt(exc_path, skiprows=1)
    dt = data[1, 0] - data[0, 0]
    nyquist = 0.5 / dt
    sos = butter(order, cutoff_hz / nyquist, btype='low', output='sos')

    filtered = data.copy()
    for col in range(1, data.shape[1]):
        filtered[:, col] = sosfiltfilt(sos, data[:, col])

    if edge_exclude:
        filtered[:, 1:1 + edge_exclude] = 0.0
        filtered[:, -edge_exclude:] = 0.0

    if keep_backup:
        backup_path = exc_path.with_name(exc_path.name.replace('excitation', 'excitation_orig'))
        if not backup_path.exists():
            shutil.copy(exc_path, backup_path)

    with open(exc_path, 'w') as f:
        f.write(header)
        np.savetxt(f, filtered, fmt='%.6e')

    return exc_path
