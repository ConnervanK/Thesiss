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
    import TimeLapse_Notebooks.helper_functions.KirchhoffPylopsZeroOffset as KirchhoffPylopsZeroOffset
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
                         t0_ns, eps_r, v_ice, stride=1, n_snap=30, snap_win=1.0):
    """
    Write gprMax excitation file and .in file for back-propagation migration.
    Files go in study_root/backprop/<slug>/.

    Sources are placed at the Tx-Rx midpoints (x_midpoints), consistent with the
    exploding-reflector half-velocity (v_mig = v/2) assumption.

    Returns (in_path, n_src, n_snaps, t_focus_ns).
    """
    import pathlib
    study_root = pathlib.Path(study_root)
    out_dir = study_root / 'backprop' / slug
    out_dir.mkdir(parents=True, exist_ok=True)

    _, n_t = tapered_ntr_nt.shape
    dt_s  = dt_ns * 1e-9
    T_ns  = n_t * dt_ns

    # Stride, time-reverse, normalise per trace
    data_s   = tapered_ntr_nt[::stride]
    x_src    = x_midpoints[::stride]
    n_src    = len(x_src)
    data_rev = data_s[:, ::-1].copy()
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

    in_lines = [
        f'#title: Back-Propagation -- {label}',
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
    in_path.write_text('\n'.join(in_lines) + '\n')
    return in_path, n_src, n_snaps, t_focus_ns
