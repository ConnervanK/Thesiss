"""
Migration functions shared by Resolution_Playground and TimeLapse_Playground notebooks.
"""

from dataclasses import dataclass, field

import numpy as np
import pylops
from scipy.fft import fft, ifft, fftfreq


def PylopsKirchoffMigration(data, t, x, vel_model, z, wav='Ricker', wavcenter='Center',
                             f0=0.08, recs=None, srcs=None, mode='analytic', dynamic=False,
                             aperture=None, angleaperture=None, return_op=False):
    """
    Kirchhoff migration for zero-offset GPR data using the PyLops adjoint operator.

    Applies the adjoint (H^T) of a zero-offset Kirchhoff demigration operator to
    map the B-scan from (time, x) to a reflectivity image in (depth, x). Sources
    and receivers are co-located (exploding-reflector assumption). All spatial
    units must be metres and all time units nanoseconds; frequency in GHz.

    Args:
        data (ndarray): B-scan of shape (n_t, n_x), time axis first.
        t (ndarray): 1-D time axis [ns], uniformly spaced.
        x (ndarray): 1-D along-profile axis [m], uniformly spaced.
        vel_model (ndarray or float): Velocity model [m/ns]. Can be a scalar
            (homogeneous) or a 2-D array of shape (n_z, n_x).
        z (ndarray): 1-D depth axis [m], uniformly spaced from 0.
        wav (str or ndarray): Wavelet to use. 'Ricker' generates a Ricker wavelet
            with centre frequency f0; otherwise pass the wavelet array directly.
        wavcenter (str or int): Index of the wavelet centre sample. 'Center' uses
            the centre returned by pylops.utils.wavelets.ricker.
        f0 (float): Centre frequency of the Ricker wavelet [GHz]. Ignored when
            wav is an array.
        recs (ndarray): Receiver positions — overridden internally (zero-offset
            forces recs = srcs = x). Retained for API compatibility.
        srcs (ndarray): Source positions — overridden internally (see recs).
        mode (str): Travel-time computation mode passed to the Kirchhoff operator
            ('analytic', 'eikonal', etc.).
        dynamic (bool): If True, apply amplitude (geometric-spreading) corrections.
        aperture (float or None): Maximum lateral migration aperture [m]. None
            uses the full aperture.
        angleaperture (float or None): Maximum dip angle included in migration
            [degrees]. None uses the full angle range.
        return_op (bool): If True, also return the Kirchhoff operator K.

    Returns:
        m (ndarray): Migrated reflectivity image of shape (n_z, n_x).
        K (LinearOperator): The Kirchhoff operator (only when return_op=True).
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
    Gazdag (1978) phase-shift migration for zero-offset post-stack GPR data.

    Downward-continues the wavefield depth-by-depth in the f-kx domain and
    extracts the t=0 imaging condition at each depth level. The half-velocity
    convention (v_mig = vel/2) implements the exploding-reflector assumption.

    Pre-processing applied internally before the depth loop:
      - 5 % spatial cosine taper + 100 % zero-pad on each side to suppress
        spatial wrap-around and Gibbs ringing.
      - f-kx evanescent filter: bins where |kx| > f/v_mig are zeroed to
        prevent migration smiles caused by energy outside the propagation cone.

    Args:
        data (ndarray): B-scan of shape (n_t, n_x), time axis first. Should be
            t0-shifted so that t=0 corresponds to the surface (z=0).
        x (ndarray): 1-D along-profile axis [m], uniformly spaced.
        t (ndarray): 1-D two-way travel-time axis [ns], uniformly spaced.
        z (ndarray): 1-D depth axis [m], uniformly spaced from 0. Sets the
            number of downward-continuation steps and the step size dz.
        vel (float): Full (round-trip) medium velocity [m/ns]. Halved internally
            to v_mig = vel/2 for the one-way phase-shift operator.

    Returns:
        image (ndarray): Migrated reflectivity image of shape (n_z, n_x),
            peak-normalised to the range [-1, 1].
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
                         sign_bit=False, dx=0.001, domain_y=1.0, src_y=0.9,
                         pml_cells=10):
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

    Returns (in_path, n_src, n_snaps, t_focus_ns). If min(x_midpoints[::stride])
    is closer to 0 than pml_cells*dx, every source (and the domain) is silently
    shifted by x_offset = pml_cells*dx - min(x_midpoints[::stride]) so sources
    clear the x-min PML band -- gprMax x=0 then corresponds to true x_midpoints
    value -x_offset. x_offset is not returned (kept out of the return tuple for
    backward compatibility with existing 4-way unpacking call sites); recompute
    it from the same inputs if you need to map gprMax x back to physical
    coordinates, e.g. `x_offset = max(0.0, pml_cells*dx - min(x_midpoints))`.
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
    # Extra polarity flip: gprMax's hertzian_dipole is a CURRENT source, not a direct
    # field source, so re-injecting a time-reversed E-field recording as-is does not
    # correctly implement EM time-reversal (Maxwell's equations are first-order in
    # time; naive field time-reversal needs a compensating sign flip on the source
    # term -- the standard result in time-reversal-mirror EM literature, unlike the
    # self-dual acoustic case). Verified empirically: a synthetic single-reflector
    # gprMax test showed the back-propagated focus correlates at -0.92 with the
    # correct polarity and +0.92 with the sign-flipped one, ruling out a phase/
    # quadrature error (Hilbert-transform correlation ~0) in favour of a clean
    # 180-degree inversion.
    data_rev = -data_s[:, ::-1].copy()
    if sign_bit:
        data_rev = np.sign(data_rev)
    else:
        peak     = np.max(np.abs(data_rev), axis=1, keepdims=True)
        peak[peak == 0] = 1.0
        data_rev /= peak

    # Snapshot timing: all depths focus simultaneously at t_focus = T - t0
    # t_start must be > 0 (gprMax rejects snapshot time = 0)
    t_focus_ns = T_ns - t0_ns
    t_start_ns = max(dt_ns, t_focus_ns - snap_win)
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

    # Domain geometry derived from source positions and grid parameters.
    # x_offset shifts every source away from x=0 only when the minimum source
    # position would otherwise land inside the x-min PML band (e.g. a profile
    # that starts at or near 0 m -- #pml_cells applies a pml_pad-thick layer on
    # *both* x edges, not just x_max, so sources near x=0 need the same margin
    # sources near x_max already get from the padded domain_x below). For the
    # normal case (min source position already well clear of pml_pad, as with
    # borehole depths of tens of metres) x_offset is 0 and behaviour is
    # unchanged from before.
    pml_pad  = pml_cells * dx
    x_min    = float(np.min(x_src))
    x_max    = float(np.max(x_src))
    x_offset = max(0.0, pml_pad - x_min)
    x_src    = x_src + x_offset
    domain_x = np.ceil((x_max + x_offset + pml_pad) / dx) * dx + pml_pad
    dz       = dx

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
        f'#domain: {domain_x:.6f} {domain_y:.6f} {dz:.6f}',
        f'#dx_dy_dz: {dx:.6f} {dx:.6f} {dz:.6f}',
        f'#time_window: {n_t * dt_s:.6e}',
        f'#pml_cells: {pml_cells} {pml_cells} 0 {pml_cells} {pml_cells} 0',
        '',
        f'// Half-velocity: eps_r={eps_r_half:.2f} (=4x{eps_r}) -> v={v_half:.5f} m/ns',
        f'#material: {eps_r_half:.2f} 1e-6 1.0 0 ice',
        '',
        f'#box: 0 0 0 {domain_x:.6f} {domain_y:.6f} {dz:.6f} ice',
        '',
        f'#excitation_file: {exc_path.name}',
        '',
        f'// {n_src} time-reversed sources at Tx-Rx midpoints (zero-offset convention, stride={stride})',
        '#python:',
        'from gprMax.input_cmd_funcs import *',
        'import numpy as np',
        f"x_sources = np.load(r'{npy_path}')",
        'for i, x in enumerate(x_sources):',
        f"    hertzian_dipole('z', float(x), {src_y:.6f}, 0.0, 'bp_{{}}'.format(i))",
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
        f"    print('#snapshot: 0 0 0 {domain_x:.6f} {domain_y:.6f} {dz:.6f} {dx:.6f} {dx:.6f} {dz:.6f} %.8e bp_snap%04d' % (t_snap, k+1))",
        '#end_python:',
        '',
        '#messages: y',
    ]
    in_path.write_text('\n'.join(in_lines) + '\n', encoding='utf-8')
    return in_path, n_src, n_snaps, t_focus_ns


def write_borehole_backprop_files(study_root, label, slug, tapered_ntr_nt, dt_ns, x_midpoints,
                                   t0_ns, eps_r, v_ice, eps_r_water=81.0, sigma_water=0.01,
                                   borehole_width=0.10, left_buffer=1.0, imaging_range=12.0,
                                   src_offset=3.2, stride=1, n_snap=30, snap_win=1.0,
                                   sign_bit=False, dx=0.001, pml_cells=10, scale_water_eps=True,
                                   normalize_mode='peak', norm_scale=None):
    """
    Write gprMax excitation file and .in file for back-propagation migration through an
    explicit single-borehole geometry, instead of assuming a homogeneous background
    medium starting right at the source (as write_backprop_files does). The borehole is
    modelled as a vertical, borehole_width-wide rectangle of water spanning the full
    depth extent of the domain; source/receiver sit on its lateral centreline, matching
    a real single-hole VRP survey where both antennas are inside the fluid-filled hole.
    Files go in study_root/backprop/<slug>/.

    Radial (y) geometry, left to right:
        [PML] [left_buffer]  [borehole_width, water]  [>= imaging_range, background]  [PML]
    Source/receiver sit at the borehole's lateral centreline (the middle of the
    borehole_width rectangle).

    Depth (x) geometry: spans min(x_midpoints) to max(x_midpoints) (the recorded
    interval), padded by src_offset/2 on each end so the domain comfortably contains the
    physical Tx-Rx antenna pair even at the shallowest/deepest recorded position, plus
    PML beyond that. The whole depth axis is then shifted so gprMax x=0 sits exactly
    pml_pad + src_offset/2 before min(x_midpoints) -- unlike write_backprop_files (which
    leaves real depths like 60-85 m unshifted, wasting a 0-60 m stretch of empty
    domain), this always shifts, so the modelled domain is only as tall as it needs to
    be. True depth = gprMax x - x_shift (both returned in geom, see below).

    ASSUMPTIONS WORTH CONFIRMING (flagged because they're judgment calls, not given
    directly by the problem spec):
      1. Sources are still placed one per trace at x_midpoints (the Tx-Rx midpoint),
         with the background medium halved in velocity -- the same exploding-reflector
         convention used everywhere else in this codebase (Kirchhoff, Gazdag,
         write_backprop_files). src_offset is used only to size the depth buffer above,
         NOT to place two separate offset Tx/Rx sources. If a real two-antenna forward
         model (true offset, true velocity, no exploding-reflector trick) was intended
         instead, this function needs a different source scheme.
      2. scale_water_eps=True (default) ALSO quadruples the borehole water's
         permittivity (eps_r_water_half = 4*eps_r_water), exactly like the background
         medium. This is required for internal consistency of the exploding-reflector
         trick: one-way travel time at half-velocity only equals the true two-way
         travel time if EVERY material in the model is scaled the same way -- leaving
         the borehole at its true permittivity while halving the background would
         distort the relative delay through the borehole vs. the surrounding medium.
         sigma_water is always left unscaled regardless (conductivity sets
         attenuation, not the travel-time equivalence the halving trick relies on).

         scale_water_eps=False uses eps_r_water as-is (unscaled). This breaks that
         travel-time consistency for the (short) borehole segment specifically, but
         also weakens the borehole's dielectric-waveguide effect: at the scaled
         permittivity the guided wavelength (~17 cm at eps=324, f0=100 MHz) is
         comparable to a 10 cm borehole, i.e. strong modal confinement; at the true
         permittivity (eps=81) the wavelength is ~33 cm, giving markedly weaker
         confinement. Worth testing as a source of excess clutter in the back-
         propagated image before assuming it's real physics.

    sign_bit=True -- see write_backprop_files' docstring; identical behaviour here.

    normalize_mode : {'peak', 'minmax'}
        How the reversed traces are scaled before injection (ignored if sign_bit=True).
        'peak' (default) divides by a single scalar max(|x|), sign-preserving -- the
        convention used everywhere else in this codebase. 'minmax' instead applies
        Eq 7 of Santos & Teixeira (2017): (x - x_min) / (x_max - x_min), per trace,
        range [0, 1]. Per that paper's own text and Fig. 2, Eq 7 most likely normalises
        the *output* TR wavefield for the std-based Mode 1/2/X12 statistics, not the
        pre-injection excitation -- so 'minmax' here is an explicit experiment, not a
        reproduction of the paper's actual injection step. It is NOT sign-preserving:
        every trace's background sits at a non-zero, non-physical DC level instead of
        zero, which can inject spurious low-frequency energy from a current source.
        Confirm this is what you want before trusting the result.

    norm_scale : float or None
        Only used when normalize_mode='peak'. If None (default), the scale is
        max(|data_rev|) over the WHOLE array for this one call -- i.e. a single global
        scalar per profile, not per trace. Divide every trace by a single number, not
        by its own individual peak: per-trace normalisation was tried first and found
        to destroy exactly the amplitude information a time-lapse study depends on --
        both the real reflectivity change between profiles (an observed ~1.8x RMS
        difference between two profiles was completely erased) and the real amplitude
        structure between traces within one profile (per-trace peaks varied by
        80-160x). Pass an explicit norm_scale (e.g. the peak found across an entire
        multi-profile time-lapse set) to keep multiple write_borehole_backprop_files()
        calls on the same amplitude scale, which is required for their back-propagated
        outputs to be meaningfully differenced against each other.

    Returns (in_path, n_src, n_snaps, t_focus_ns, geom). geom is a dict of the computed
    domain layout (domain_x, domain_y, dz, dx, pml_pad, x_shift, x_min_true, x_max_true,
    depth_buffer, y_bh_start, y_bh_end, y_img_end, src_y, borehole_width, left_buffer,
    imaging_range) -- also consumed by plot_borehole_domain() for the confirmation
    figure.
    """
    import pathlib
    study_root = pathlib.Path(study_root)
    out_dir = study_root / 'backprop' / slug
    out_dir.mkdir(parents=True, exist_ok=True)

    _, n_t = tapered_ntr_nt.shape
    dt_s  = dt_ns * 1e-9
    T_ns  = n_t * dt_ns
    dz    = dx

    # --- Time-reversal + polarity flip (identical to write_backprop_files -- see its
    # docstring/comments for why the extra 180-degree flip is needed for a current-
    # source hertzian_dipole) ---
    data_s     = tapered_ntr_nt[::stride]
    x_src_true = np.asarray(x_midpoints)[::stride]
    n_src      = len(x_src_true)
    data_rev   = -data_s[:, ::-1].copy()
    if sign_bit:
        data_rev = np.sign(data_rev)
    elif normalize_mode == 'minmax':
        x_min = np.min(data_rev, axis=1, keepdims=True)
        x_max = np.max(data_rev, axis=1, keepdims=True)
        span  = x_max - x_min
        span[span == 0] = 1.0
        data_rev = (data_rev - x_min) / span
    elif normalize_mode == 'peak':
        scale = norm_scale if norm_scale is not None else np.max(np.abs(data_rev))
        if scale == 0:
            scale = 1.0
        data_rev /= scale
    else:
        raise ValueError(f"normalize_mode must be 'peak' or 'minmax', got {normalize_mode!r}")

    # Snapshot timing: all depths focus simultaneously at t_focus = T - t0
    t_focus_ns = T_ns - t0_ns
    t_start_ns = max(dt_ns, t_focus_ns - snap_win)
    t_start_s  = t_start_ns * 1e-9
    snap_step  = max(1, int((T_ns * 1e-9 - t_start_s) / (max(1, n_snap - 1) * dt_s)))
    n_snaps    = max(1, int((n_t * dt_s - t_start_s) / (snap_step * dt_s)))

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

    # --- Depth (x) axis geometry: shift so gprMax x=0 sits pml_pad+depth_buffer before
    # the shallowest recorded position, and size the domain to just cover the recorded
    # interval plus buffers on both ends. ---
    pml_pad      = pml_cells * dx
    depth_buffer = src_offset / 2.0
    x_min_true   = float(np.min(x_src_true))
    x_max_true   = float(np.max(x_src_true))
    x_shift      = pml_pad + depth_buffer - x_min_true
    x_src        = x_src_true + x_shift
    domain_x     = np.ceil((x_max_true - x_min_true + 2 * depth_buffer) / dx) * dx + 2 * pml_pad

    # --- Radial (y) axis geometry: [PML][left_buffer][borehole][imaging_range][PML] ---
    y_bh_start = pml_pad + left_buffer
    y_bh_end   = y_bh_start + borehole_width
    y_img_end  = y_bh_end + imaging_range
    domain_y   = y_img_end + pml_pad
    src_y      = (y_bh_start + y_bh_end) / 2.0

    npy_path = out_dir / 'src_positions.npy'
    np.save(npy_path, x_src)

    eps_r_half       = 4.0 * eps_r
    eps_r_water_half = 4.0 * eps_r_water if scale_water_eps else eps_r_water
    v_half           = v_ice / 2
    in_path          = out_dir / f'backprop_{slug}.in'

    if sign_bit:
        excitation_mode = 'sign-bit (sign(u), amplitude stripped)'
    elif normalize_mode == 'minmax':
        excitation_mode = 'min-max normalised (Eq 7 style, range [0,1], NOT sign-preserving)'
    else:
        excitation_mode = 'peak-normalised (range [-1,1])'
    in_lines = [
        f'#title: Borehole Back-Propagation -- {label}',
        f'// Excitation mode: {excitation_mode}',
        f'#domain: {domain_x:.6f} {domain_y:.6f} {dz:.6f}',
        f'#dx_dy_dz: {dx:.6f} {dx:.6f} {dz:.6f}',
        f'#time_window: {n_t * dt_s:.6e}',
        f'#pml_cells: {pml_cells} {pml_cells} 0 {pml_cells} {pml_cells} 0',
        '',
        f'// Half-velocity background: eps_r={eps_r_half:.2f} (=4x{eps_r:.4f}) -> v={v_half:.5f} m/ns',
        f'#material: {eps_r_half:.4f} 1e-6 1.0 0 ice',
        f'// {"Half-velocity" if scale_water_eps else "TRUE, unscaled"} borehole fluid: '
        f'eps_r={eps_r_water_half:.2f} '
        f'({"=4x" + str(eps_r_water) if scale_water_eps else "unscaled, true permittivity"}), '
        f'sigma={sigma_water} S/m (always unscaled -- see docstring point 2)',
        f'#material: {eps_r_water_half:.4f} {sigma_water:.6f} 1.0 0 water',
        '',
        f'#box: 0 0 0 {domain_x:.6f} {domain_y:.6f} {dz:.6f} ice',
        f'// Borehole: {borehole_width * 100:.0f} cm wide water-filled rectangle, full depth extent',
        f'#box: 0 {y_bh_start:.6f} 0 {domain_x:.6f} {y_bh_end:.6f} {dz:.6f} water',
        '',
        f'#excitation_file: {exc_path.name}',
        '',
        f'// {n_src} time-reversed sources at Tx-Rx midpoints, on the borehole centreline (y={src_y:.4f} m)',
        '#python:',
        'from gprMax.input_cmd_funcs import *',
        'import numpy as np',
        f"x_sources = np.load(r'{npy_path}')",
        'for i, x in enumerate(x_sources):',
        f"    hertzian_dipole('z', float(x), {src_y:.6f}, 0.0, 'bp_{{}}'.format(i))",
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
        f"    print('#snapshot: 0 0 0 {domain_x:.6f} {domain_y:.6f} {dz:.6f} {dx:.6f} {dx:.6f} {dz:.6f} %.8e bp_snap%04d' % (t_snap, k+1))",
        '#end_python:',
        '',
        '#messages: y',
    ]
    in_path.write_text('\n'.join(in_lines) + '\n', encoding='utf-8')

    geom = dict(
        domain_x=domain_x, domain_y=domain_y, dz=dz, dx=dx, pml_pad=pml_pad,
        x_shift=x_shift, x_min_true=x_min_true, x_max_true=x_max_true,
        depth_buffer=depth_buffer, y_bh_start=y_bh_start, y_bh_end=y_bh_end,
        y_img_end=y_img_end, src_y=src_y, borehole_width=borehole_width,
        left_buffer=left_buffer, imaging_range=imaging_range,
    )
    return in_path, n_src, n_snaps, t_focus_ns, geom


def plot_borehole_domain(geom, x_src_true=None, save_path=None):
    """
    Draw a schematic of the borehole back-propagation domain built by
    write_borehole_backprop_files() -- PML bands, background medium, the water-filled
    borehole rectangle, and source/receiver positions -- so the geometry can be sanity-
    checked before spending any gprMax runtime on it. Pure matplotlib, no gprMax needed.

    Args:
        geom (dict): The geom dict returned by write_borehole_backprop_files().
        x_src_true (array-like or None): True (unshifted) source depths to mark on the
            depth axis. If None, only the modelled depth range is shown.
        save_path (str, Path, or None): If given, the figure is also saved there.

    Returns:
        fig (matplotlib.figure.Figure)
    """
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches

    domain_x, domain_y = geom['domain_x'], geom['domain_y']
    pml   = geom['pml_pad']
    y0    = geom['y_bh_start']
    y1    = geom['y_bh_end']
    y2    = geom['y_img_end']
    src_y = geom['src_y']
    x_shift = geom['x_shift']

    zoom_x_max = y1 + 1.0     # a bit past the borehole's right edge (radial)
    zoom_y_max = 4.0          # only show the top few metres of depth -- full 25+ m
                               # depth in a <3 m-wide panel would force an extreme
                               # aspect ratio that leaves no room for dimension labels

    fig, (ax_full, ax_zoom) = plt.subplots(1, 2, figsize=(13, 6.5),
                                            gridspec_kw={'width_ratios': [1.3, 1]})

    for ax in (ax_full, ax_zoom):
        ax.add_patch(patches.Rectangle((0, 0), domain_y, domain_x,
                                        facecolor='none', edgecolor='black', lw=1.2))
        pml_kw = dict(facecolor='lightgray', edgecolor='gray', hatch='//', lw=0.4)
        ax.add_patch(patches.Rectangle((0, 0), pml, domain_x, **pml_kw))
        ax.add_patch(patches.Rectangle((domain_y - pml, 0), pml, domain_x, **pml_kw))
        ax.add_patch(patches.Rectangle((0, 0), domain_y, pml, **pml_kw))
        ax.add_patch(patches.Rectangle((0, domain_x - pml), domain_y, pml, **pml_kw))
        ax.add_patch(patches.Rectangle((pml, pml), domain_y - 2 * pml, domain_x - 2 * pml,
                                        facecolor='#cfe8f3', edgecolor='none', zorder=0))
        ax.add_patch(patches.Rectangle((y0, 0), y1 - y0, domain_x,
                                        facecolor='#1f77b4', edgecolor='none', zorder=1))
        ax.axvline(src_y, color='red', lw=1, ls='--', zorder=2)
        if x_src_true is not None:
            xs = np.atleast_1d(x_src_true) + x_shift
            ax.plot(np.full_like(xs, src_y, dtype=float), xs, 'r.', ms=3, zorder=3)
        else:
            ax.plot([src_y], [domain_x / 2], 'r.', ms=6, zorder=3)
        ax.set_xlabel('Radial distance, y (m)')
        ax.set_aspect('equal', adjustable='box')

    ax_full.set_xlim(-0.3, domain_y)
    ax_full.set_ylim(domain_x + 1.0, -1.0)   # depth increases downward
    ax_full.set_ylabel(f'gprMax depth axis, x (m)  [true depth = x - {x_shift:.2f} m]')
    ax_full.set_title('Full domain')

    ax_zoom.set_xlim(-0.2, zoom_x_max)
    ax_zoom.set_ylim(zoom_y_max, -1.6)       # cropped depth window near the top PML
    ax_zoom.set_ylabel('gprMax depth axis, x (m)')
    ax_zoom.set_title(f'Near-borehole detail\n(radial 0-{zoom_x_max:.1f} m, depth 0-{zoom_y_max:.0f} m of {domain_x:.0f} m shown)')

    def _dim(ax, x0, x1, y, text, open_ended=False):
        style = '-|>' if open_ended else '<->'
        ax.annotate('', xy=(x1, y), xytext=(x0, y),
                    arrowprops=dict(arrowstyle=style, color='black', lw=0.9))
        ax.text((x0 + x1) / 2, y - 0.15, text, ha='center', va='bottom', fontsize=8)

    _dim(ax_zoom, pml, y0, -0.6, f'{geom["left_buffer"]:.1f} m buffer')
    _dim(ax_zoom, y0, y1, -1.3, f'{geom["borehole_width"] * 100:.0f} cm borehole')
    _dim(ax_zoom, y1, zoom_x_max, -0.6, f'>= {geom["imaging_range"]:.0f} m imaging', open_ended=True)
    ax_zoom.axvline(pml, color='gray', lw=0.6, ls=':')

    fig.suptitle('Borehole back-propagation domain (schematic)')
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    return fig


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

    if cutoff_hz >= nyquist:
        # The excitation's own sample rate (fixed by the trace dt, independent of the
        # grid dx cutoff_hz was derived from) already band-limits it below the
        # requested cutoff -- e.g. a finer grid raises the dispersion-safe cutoff past
        # the data's Nyquist frequency. Nothing to filter; butter() would otherwise
        # raise ValueError for a normalised frequency >= 1.
        print(f'  cutoff {cutoff_hz/1e9:.3f} GHz >= data Nyquist {nyquist/1e9:.3f} GHz -- '
              f'already band-limited, skipping filter.')
        filtered = data
    else:
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



def apply_3d_to_2d_correction(b_scan, dt, velocity, time_zero_idx=0,
                               phase_sign=-1.0, water_level=0.0):
    """
    Converts 3D recorded GPR data to a 2D equivalent format for 2D FDTD back-propagation.
    Applies the mathematical corrections outlined in the G_2D / G_3D Green's function ratio.

    Parameters:
    -----------
    b_scan : numpy.ndarray
        The 2D array of GPR traces (shape: [num_traces, num_time_samples]).
    dt : float
        The time step between samples in seconds.
    velocity : float
        The electromagnetic velocity of the background medium (m/s).
    time_zero_idx : int
        The index of the time zero (t=0) in the trace. Used to properly scale sqrt(t).
    phase_sign : float
        Sign of the pi/4 phase rotation applied to positive frequencies, i.e. the
        filter uses exp(phase_sign * 1j * pi/4). The Hankel-asymptotic derivation
        (see the docstring's Green's-function ratio) gives +1; the default of -1
        was chosen empirically by comparing full-pipeline back-propagation output
        against Kirchhoff/Gazdag images, which is a confounded test (SVD, tapers,
        the dip filter, and the separate hertzian-dipole polarity flip in
        write_backprop_files() all sit between this filter and that comparison).
        Exposed as a parameter so it can be validated in isolation against an
        analytic single-reflector synthetic instead.
    water_level : float
        Regularises the 1/sqrt(omega) pole near DC, expressed as a fraction of the
        Nyquist angular frequency: omega_eff = max(omega, water_level * omega_nyquist).
        0 (default) reproduces the original behaviour of exactly zeroing only the
        DC bin, leaving the next few bins with very large (unregularised) gain.
        A small positive value (e.g. 0.01-0.05) floors the gain across the whole
        low-frequency range instead, which should suppress the disproportionate
        low-frequency noise amplification that a hard-zero-only DC bin permits.

    Returns:
    --------
    corrected_b_scan : numpy.ndarray
        The pre-conditioned B-scan ready to be time-reversed and injected into gprMax.
    """
    num_traces, num_samples = b_scan.shape
    corrected_b_scan = np.zeros_like(b_scan)

    # 1. Spatial Amplitude Correction: sqrt(r)
    # Since r = (v * t) / 2 for a reflection, sqrt(r) is proportional to sqrt(t).
    # We apply a sqrt(t) gain to correct 1/r 3D spreading to 1/sqrt(r) 2D spreading.
    t_array = np.arange(num_samples) * dt
    # Shift time array so time-zero is actually t=0
    t_array = t_array - (time_zero_idx * dt)
    # Prevent negative times or zero (to avoid divide-by-zero or complex numbers)
    t_array[t_array <= 0] = 1e-12

    # The amplitude scalar proportional to sqrt(r)
    spatial_gain = np.sqrt(velocity * t_array / 2.0)

    # 2. Phase and Frequency Correction: (1 / sqrt(w)) * e^(phase_sign * i * pi / 4)
    # Prepare the frequency axis
    freqs = fftfreq(num_samples, d=dt)
    omega = 2.0 * np.pi * np.abs(freqs)

    if water_level > 0:
        # Water-level regularisation: floor omega everywhere (including DC) instead
        # of only zeroing bin 0, to avoid amplifying the low-frequency bins next to it.
        omega_nyquist = np.pi / dt
        omega_eff = np.maximum(omega, water_level * omega_nyquist)
    else:
        # Original behaviour: only guard against the exact DC divide-by-zero here;
        # the DC bin itself is hard-zeroed below.
        omega_eff = omega.copy()
        omega_eff[0] = 1e-12

    H_filter = (1.0 / np.sqrt(omega_eff)) * np.exp(phase_sign * 1j * np.pi / 4.0)

    # Ensure Hermitian symmetry for a real-valued time signal
    # Negative frequencies must be the complex conjugate of positive frequencies
    H_filter[freqs < 0] = np.conj(H_filter[freqs < 0])

    if water_level <= 0:
        # Zero out the DC component to prevent massive baseline drift from the
        # 1/sqrt(w) integration (skipped when water-level regularisation is active,
        # since omega_eff already gives DC a finite, bounded gain).
        H_filter[0] = 0.0 + 0.0j

    # Apply corrections trace by trace
    for i in range(num_traces):
        trace = b_scan[i, :]

        # Step A: Apply spatial gain sqrt(r)
        trace_gained = trace * spatial_gain

        # Step B: Apply the half-integration phase shift filter in frequency domain
        trace_fft = fft(trace_gained)
        trace_filtered = np.real(ifft(trace_fft * H_filter))

        corrected_b_scan[i, :] = trace_filtered

    return corrected_b_scan


@dataclass
class MigratedImage:
    """Unified container returned by load_migrated_image() for a single migrated /
    back-propagated image. image/depth_axis/radial_axis always follow ONE convention
    -- row index increases with physical depth (row 0 = shallowest), column index
    increases with radial distance from the borehole (col 0 = nearest) -- regardless
    of which native convention the underlying migration technique used. See
    load_migrated_image's docstring for why this matters.
    """
    image: np.ndarray
    depth_axis: np.ndarray
    radial_axis: np.ndarray
    dz: float
    dx: float
    kz_cent: float
    method: str
    run: int
    meta: dict = field(default_factory=dict)


def load_migrated_image(method, run, *, migrated_dir=None, borehole_dir=None,
                         depth_gk=None, x_img_gk=None, dL_gk=None,
                         f0_mig=None, v=None):
    """
    Single, technique-agnostic loader for the three migration techniques used in
    FieldData_Playground.ipynb's ROI / displacement-estimation workflow (Gazdag,
    Kirchhoff-BP, back-propagation through the borehole-geometry gprMax pipeline).

    Resolves, architecturally, a sign-convention bug that was previously patched
    one-off per call site (e.g. a manual dz_g_bh negation in the notebook's
    borehole-displacement cell): Gazdag/Kirchhoff-BP's cached 'migrated/
    {method}_{run}.npy' arrays and the shared `depth` axis DECREASE with row index
    (row 0 = 85 m, deepest), while back-propagation's per-run
    'borehole_prof_{run}_dx02_final_processed.npz' INCREASES with row index (row 0 =
    shallowest). Feeding a WLS phase-plane fit two images on opposite row-direction
    conventions makes its raw Delta_z sign mean opposite physical directions
    depending on which technique produced the image, even though the fit code itself
    never changes. This function is the single place that normalises every
    technique's image + axes to the same convention before anything downstream (ROI
    selection, cropping, the WLS fit) ever sees it, so no per-technique sign
    correction is needed anywhere else.

    Parameters
    ----------
    method : {'gazdag', 'kirchhoff_bp', 'backprop'}
    run : int
        Profile number.
    migrated_dir : Path
        Directory holding f'{method}_{run}.npy' (gazdag / kirchhoff_bp only).
    borehole_dir : Path
        Directory holding f'borehole_prof_{run}_dx02_final_processed.npz'
        (method='backprop' only -- this is `bh_out` in the notebook).
    depth_gk, x_img_gk, dL_gk : ndarray, ndarray, float
        The notebook's shared `depth` (decreasing, m), `x_img` (increasing, m) and
        `dL` (trace spacing, m) globals -- gazdag/kirchhoff_bp only; the cached .npy
        arrays don't carry their own axes.
    f0_mig, v : float
        Centre frequency [GHz] and full (round-trip) velocity [m/ns].
        gazdag/kirchhoff_bp: kz_cent = 2*pi*f0_mig/v (full velocity).
        backprop: kz_cent = 2*pi*f0_mig/(v/2) (half-velocity exploding-reflector
        convention -- gprMax's back-propagation snapshots live in that domain,
        unlike the already-migrated Kirchhoff/Gazdag .npy images).

    Returns
    -------
    MigratedImage, or None if the underlying file does not exist (mirrors the
    notebook's previous per-cell `_load_img` contract).
    """
    import pathlib

    if method in ('gazdag', 'kirchhoff_bp'):
        p = pathlib.Path(migrated_dir) / f'{method}_{run}.npy'
        if not p.exists():
            return None
        image_dec = np.load(p)
        n = image_dec.shape[0]
        radial_axis = np.asarray(x_img_gk)
        return MigratedImage(
            image=image_dec[::-1, :].copy(),
            depth_axis=np.asarray(depth_gk)[:n][::-1].copy(),
            radial_axis=radial_axis,
            dz=float(dL_gk),
            dx=float(radial_axis[1] - radial_axis[0]),
            kz_cent=2.0 * np.pi * f0_mig / v,
            method=method, run=run, meta={},
        )
    elif method == 'backprop':
        p = pathlib.Path(borehole_dir) / f'borehole_prof_{run}_dx02_final_processed.npz'
        if not p.exists():
            return None
        d = np.load(p)
        depth_axis = d['depth_axis']
        radial_axis = d['radial_axis']
        meta = {k: (d[k].item() if d[k].ndim == 0 else d[k])
                for k in d.files if k not in ('image', 'depth_axis', 'radial_axis')}
        return MigratedImage(
            image=d['image'], depth_axis=depth_axis, radial_axis=radial_axis,
            dz=float(depth_axis[1] - depth_axis[0]),
            dx=float(radial_axis[1] - radial_axis[0]),
            kz_cent=2.0 * np.pi * f0_mig / (v / 2.0),
            method=method, run=run, meta=meta,
        )
    else:
        raise ValueError(
            f"method must be 'gazdag', 'kirchhoff_bp', or 'backprop', got {method!r}")


def wls_phase_plane_fit(base, mon, dz, dx, kz_cent, *,
                         roi_px=None, data_pad=8, taper='edge', tukey_alpha=0.15,
                         kz_band_fac=0.5, kx_band_fac=2.0, amp_thr=0.20, wls_pow=1,
                         pad_fac=10, force_dz_zero=False, mask=None,
                         pos_kz_only=True, return_diagnostics=False):
    """
    Weighted least-squares (WLS) cross-spectrum phase-plane fit -- the displacement
    estimator behind every ROI-workflow strategy in FieldData_Playground.ipynb
    (rectangular window, sliding window, k-space / depth-radial manual picking).
    Fits phi(kz, kx) = kz*dz_est + kx*dx_est + phi_0 to the phase of the cross-
    spectrum XS = FFT(base) * conj(FFT(mon)), weighted by |XS|**wls_pow and
    restricted to either an automatic band + amplitude-threshold mask or an explicit
    caller-supplied k-space mask (manual picking).

    Consolidates what were previously ~4 separately-maintained, quietly-diverged
    copies of this same fit (two rectangular-window cells, one sliding-window
    helper, one bare mask-only helper) into one function, with what divergence
    turned out to be genuine (rather than accidental drift) exposed as explicit,
    per-technique-configurable parameters -- see WLS_FIT_DEFAULTS_BY_METHOD in the
    notebook.

    Parameters
    ----------
    base, mon : ndarray (n_depth, n_radial)
        base = earlier profile, mon = later profile, on a shared grid already in the
        row-index-increases-with-depth convention (e.g. from load_migrated_image).
        If roi_px is None, base/mon are used exactly as passed -- crop to the region
        of interest yourself first (sliding-window, or a Gaussian-softened painted
        mask already baked in for picking). If roi_px is given, pass the full image
        and this function does the cropping.
    dz, dx : float
        Grid spacing [m] (positive).
    kz_cent : float
        Dominant wavenumber [rad/m] (2*pi*f0/v, or the half-velocity equivalent for
        back-propagation); sets the default kz/kx band via kz_band_fac/kx_band_fac.
    roi_px : (z0, z1, x0, x1) or None
        Pixel-index ROI within base/mon. If given, crops with `data_pad` extra
        pixels of real data on each side before tapering, so the taper rolls off
        through genuine surrounding data instead of attenuating the ROI signal
        itself (the rectangular-window recipe). If None, no cropping is done here.
    data_pad : int
        Real-data margin pixels added around roi_px before tapering. Only used when
        roi_px is given and taper='edge'.
    taper : {'edge', 'tukey', 'none'}
        'edge': sin^2 ramp over exactly the data_pad margin -- requires roi_px.
        'tukey': tukey(alpha=tukey_alpha) window over the whole array as passed.
        'none': no taper applied here -- use when the caller already tapered/masked
        base/mon before calling (e.g. picking cells bake a Gaussian-softened painted
        mask into base/mon directly).
    kz_band_fac, kx_band_fac : float
        Automatic mask band: |KZ| < kz_band_fac*kz_cent, |KX| < kx_band_fac*kz_cent.
        Ignored if `mask` is given.
    amp_thr : float
        Automatic mask amplitude gate: |XS| > amp_thr * max(|XS|). Ignored if `mask`
        is given.
    wls_pow : float
        Weight exponent -- W = |XS|**wls_pow. A higher power (e.g. 3) suppresses
        low-energy peripheral k-cells that would otherwise bias the fit; whether a
        given technique needs this is empirical, not universal -- see
        WLS_FIT_DEFAULTS_BY_METHOD.
    pad_fac : int
        Zero-padding factor before the FFT2 (denser kz/kx sampling).
    force_dz_zero : bool
        If True, fit only (dx_est, phi_0) -- for a known lateral-only displacement.
    mask : ndarray (bool) or None
        Explicit k-space mask over the (padded) KZ/KX grid. If given, the fit is
        restricted to exactly these cells and kz_band_fac/kx_band_fac/amp_thr/
        pos_kz_only are ignored entirely (manual picking).
    pos_kz_only : bool
        If True (default), the automatic mask additionally requires KZ > 0.
        base/mon are real-valued, so XS = FFT(base)*conj(FFT(mon)) is exactly
        Hermitian: XS(-kz,-kx) = conj(XS(kz,kx)), i.e. phi(-kz,-kx) = -phi(kz,kx).
        The symmetric band |KZ| < kz_band_fac*kz_cent therefore always admits a
        point cloud together with its exact mirror through the origin. Fitting
        phi = kz*dz + kx*dx + phi_0 through both at once forces a single shared
        phi_0 onto data that actually needs phi_0 and -phi_0 on the two halves;
        the least-squares fit resolves that conflict by driving phi_0 towards
        zero and biasing the fitted slope (dz_est/dx_est) to compensate -- i.e.
        it averages the tangent between the two mirrored point clouds instead of
        fitting either one correctly. Restricting to KZ > 0 keeps only one
        physical point cloud (the positive-frequency one) and removes its
        Hermitian mirror, eliminating that bias. Only kz_cent's carrier axis
        (KZ) needs this restriction -- KX is not itself a Hermitian-paired
        frequency axis once KZ's sign is fixed, so kx_band_fac's range is left
        symmetric. Ignored if `mask` is given. Set False to recover the old
        symmetric-band behaviour.
    return_diagnostics : bool
        If True, also return a dict with KZ, KX, w (amplitude), phi (phase), the
        selected mask, and the fitted-plane array -- everything a 5-panel diagnostic
        figure needs; 1-D kz/kx slice extraction for that figure is left to the
        caller (a plotting concern, not a fit concern).

    Returns
    -------
    (dz_est, dx_est, phi_0, n_mask) or, if return_diagnostics, with a 5th dict
    element `diag` containing KZ, KX, w, phi, mask, fitted, kz_ax, kx_ax.
    """
    from scipy.signal.windows import tukey as _tukey

    if roi_px is not None:
        z0, z1, x0, x1 = roi_px
        Nz_img, Nx_img = base.shape
        z0p, z1p = max(0, z0 - data_pad), min(Nz_img, z1 + data_pad)
        x0p, x1p = max(0, x0 - data_pad), min(Nx_img, x1 + data_pad)
        base_crop = base[z0p:z1p, x0p:x1p]
        mon_crop = mon[z0p:z1p, x0p:x1p]
        pad_lo_hi = (z0 - z0p, z1p - z1, x0 - x0p, x1p - x1)
    else:
        base_crop, mon_crop = base, mon
        pad_lo_hi = None

    Nz, Nx = base_crop.shape
    Nz_pad, Nx_pad = Nz * pad_fac, Nx * pad_fac

    if taper == 'edge':
        if pad_lo_hi is None:
            raise ValueError(
                "taper='edge' requires roi_px (needs a real-data margin to ramp through)")

        def _edge_taper_1d(N, n_lo, n_hi):
            win = np.ones(N)
            if n_lo > 0:
                win[:n_lo] = np.sin(np.linspace(0, np.pi / 2, n_lo)) ** 2
            if n_hi > 0:
                win[-n_hi:] = np.sin(np.linspace(np.pi / 2, 0, n_hi)) ** 2
            return win

        nz_lo, nz_hi, nx_lo, nx_hi = pad_lo_hi
        win2d = np.outer(_edge_taper_1d(Nz, nz_lo, nz_hi), _edge_taper_1d(Nx, nx_lo, nx_hi))
    elif taper == 'tukey':
        win2d = np.outer(_tukey(Nz, alpha=tukey_alpha), _tukey(Nx, alpha=tukey_alpha))
    elif taper == 'none':
        win2d = 1.0
    else:
        raise ValueError(f"taper must be 'edge', 'tukey', or 'none', got {taper!r}")

    kz_ax = np.fft.fftfreq(Nz_pad, d=dz) * 2 * np.pi
    kx_ax = np.fft.fftfreq(Nx_pad, d=dx) * 2 * np.pi
    KZ, KX = np.meshgrid(kz_ax, kx_ax, indexing='ij')

    XS = (np.fft.fft2(base_crop * win2d, s=(Nz_pad, Nx_pad)) *
          np.conj(np.fft.fft2(mon_crop * win2d, s=(Nz_pad, Nx_pad))))
    w_amp, phi = np.abs(XS), np.angle(XS)

    if mask is not None:
        kmask = mask
    else:
        band = (np.abs(KZ) < kz_band_fac * kz_cent) & (np.abs(KX) < kx_band_fac * kz_cent)
        if pos_kz_only:
            band = band & (KZ > 0)
        kmask = (w_amp > amp_thr * w_amp.max()) & band & ((np.abs(KZ) + np.abs(KX)) > 0)

    n_mask = int(kmask.sum())
    if n_mask < 3:
        dz_est = dx_est = phi_0 = 0.0
    else:
        W = w_amp[kmask] ** wls_pow
        if force_dz_zero:
            A = np.column_stack([KX[kmask], np.ones(n_mask)])
            c = np.linalg.lstsq(A * W[:, None], phi[kmask] * W, rcond=None)[0]
            dz_est, dx_est, phi_0 = 0.0, float(c[0]), float(c[1])
        else:
            A = np.column_stack([KZ[kmask], KX[kmask], np.ones(n_mask)])
            c = np.linalg.lstsq(A * W[:, None], phi[kmask] * W, rcond=None)[0]
            dz_est, dx_est, phi_0 = float(c[0]), float(c[1]), float(c[2])

    if not return_diagnostics:
        return dz_est, dx_est, phi_0, n_mask

    fitted = KX * dx_est + KZ * dz_est + phi_0
    diag = dict(KZ=KZ, KX=KX, w=w_amp, phi=phi, mask=kmask, fitted=fitted,
                kz_ax=kz_ax, kx_ax=kx_ax)
    return dz_est, dx_est, phi_0, n_mask, diag


def _wrap_to_pi(angle):
    """Wrap an array of angles [rad] to (-pi, pi]."""
    return (angle + np.pi) % (2 * np.pi) - np.pi


def ransac_phase_plane_fit(base, mon, dz, dx, kz_cent, *,
                            roi_px=None, data_pad=8, taper='edge', tukey_alpha=0.15,
                            kz_band_fac=0.5, kx_band_fac=2.0, amp_thr=0.20,
                            pad_fac=10, force_dz_zero=False, mask=None,
                            pos_kz_only=True,
                            n_iter=1000, residual_thr=0.35, weight_by_amp=True,
                            final_wls_pow=1, random_state=0,
                            return_diagnostics=False):
    """
    RANSAC cross-spectrum phase-plane fit -- a robust alternative to
    wls_phase_plane_fit for the same model phi(kz, kx) = kz*dz_est + kx*dx_est + phi_0.

    wls_phase_plane_fit minimises squared phase residual over every k-cell it is
    given, on the raw np.angle(XS) values (no unwrapping). Any masked cell --
    cross-spectrum noise, side-lobe energy from a reflector other than the one
    being tracked, or a cell where the true phase ramp exceeds +-pi across the
    fitting window and so wraps -- pulls the fitted plane, and a high-amplitude
    wrapped cell can pull it hard even under |XS|**wls_pow weighting. This is the
    suspected mechanism behind the larger, harder-to-trust Delta_z on the
    wider-band pairs (Chase/Wait/Pull; c.f. Push's "sub-wavelength increments,
    zero wrapping" in WLS_FIT_DEFAULTS_BY_METHOD's notebook comment).

    This function fits the same plane with RANSAC instead: repeatedly (a) fits the
    exact plane through a minimal random sample of 2 (force_dz_zero) or 3 k-cells,
    (b) scores it by the amplitude-weighted count of cells whose CIRCULAR residual
    (phase wrapped to (-pi, pi] via _wrap_to_pi, so a cell that is a whole 2*pi off
    the candidate plane is correctly scored as a perfect inlier, not an outlier)
    falls under `residual_thr`, and (c) keeps the best-scoring consensus set across
    `n_iter` trials. The final estimate is an amplitude-weighted least-squares
    refit (`final_wls_pow`, same |XS|**power weighting as wls_phase_plane_fit)
    restricted to that best inlier set only -- cells that don't lie on any
    consistent plane (the actual outliers) never enter the estimate, rather than
    being merely down-weighted alongside genuine signal.

    Shares its crop/taper/mask/cross-spectrum construction with
    wls_phase_plane_fit -- same base/mon/dz/dx/kz_cent/roi_px/data_pad/taper/
    kz_band_fac/kx_band_fac/amp_thr/mask/pos_kz_only semantics (see that
    docstring -- in particular pos_kz_only=True's fix for the Hermitian-mirror
    tangent-averaging bias, since the same automatic-mask construction is shared
    here), so the two are directly comparable: call both on the same (base, mon,
    roi_px, ...), or both with the same explicit `mask` (e.g. from a napari
    manual-picking cell), to isolate the fitting method as the only difference
    between the two estimates.

    Parameters (in addition to wls_phase_plane_fit's shared ones above)
    ----------
    n_iter : int
        Number of random minimal-sample trials.
    residual_thr : float
        Circular phase-residual threshold [rad] for a k-cell to count as an inlier
        of a candidate plane.
    weight_by_amp : bool
        If True (default), a candidate plane's consensus score is the summed |XS|
        of its inliers rather than the raw inlier count -- a plane explained by a
        few strong reflections outscores one explained by many weak/noisy cells,
        consistent with wls_phase_plane_fit's own amplitude weighting.
    final_wls_pow : float
        Weight exponent (|XS|**final_wls_pow) for the final WLS refit on the best
        inlier set. Set to 0 for an unweighted refit.
    random_state : int or None
        Seed for the minimal-sample draws (np.random.default_rng). Fixed by
        default so repeated notebook runs reproduce the same fit; pass None for a
        fresh draw each call.
    return_diagnostics : bool
        If True, also return a dict like wls_phase_plane_fit's diag, plus
        `inlier_mask` / `outlier_mask` (2D boolean, full KZ/KX grid, split from
        `mask`) and `best_score`.

    Returns
    -------
    (dz_est, dx_est, phi_0, n_inliers) or, if return_diagnostics, with a 5th dict
    element `diag`.
    """
    from scipy.signal.windows import tukey as _tukey

    if roi_px is not None:
        z0, z1, x0, x1 = roi_px
        Nz_img, Nx_img = base.shape
        z0p, z1p = max(0, z0 - data_pad), min(Nz_img, z1 + data_pad)
        x0p, x1p = max(0, x0 - data_pad), min(Nx_img, x1 + data_pad)
        base_crop = base[z0p:z1p, x0p:x1p]
        mon_crop = mon[z0p:z1p, x0p:x1p]
        pad_lo_hi = (z0 - z0p, z1p - z1, x0 - x0p, x1p - x1)
    else:
        base_crop, mon_crop = base, mon
        pad_lo_hi = None

    Nz, Nx = base_crop.shape
    Nz_pad, Nx_pad = Nz * pad_fac, Nx * pad_fac

    if taper == 'edge':
        if pad_lo_hi is None:
            raise ValueError(
                "taper='edge' requires roi_px (needs a real-data margin to ramp through)")

        def _edge_taper_1d(N, n_lo, n_hi):
            win = np.ones(N)
            if n_lo > 0:
                win[:n_lo] = np.sin(np.linspace(0, np.pi / 2, n_lo)) ** 2
            if n_hi > 0:
                win[-n_hi:] = np.sin(np.linspace(np.pi / 2, 0, n_hi)) ** 2
            return win

        nz_lo, nz_hi, nx_lo, nx_hi = pad_lo_hi
        win2d = np.outer(_edge_taper_1d(Nz, nz_lo, nz_hi), _edge_taper_1d(Nx, nx_lo, nx_hi))
    elif taper == 'tukey':
        win2d = np.outer(_tukey(Nz, alpha=tukey_alpha), _tukey(Nx, alpha=tukey_alpha))
    elif taper == 'none':
        win2d = 1.0
    else:
        raise ValueError(f"taper must be 'edge', 'tukey', or 'none', got {taper!r}")

    kz_ax = np.fft.fftfreq(Nz_pad, d=dz) * 2 * np.pi
    kx_ax = np.fft.fftfreq(Nx_pad, d=dx) * 2 * np.pi
    KZ, KX = np.meshgrid(kz_ax, kx_ax, indexing='ij')

    XS = (np.fft.fft2(base_crop * win2d, s=(Nz_pad, Nx_pad)) *
          np.conj(np.fft.fft2(mon_crop * win2d, s=(Nz_pad, Nx_pad))))
    w_amp, phi = np.abs(XS), np.angle(XS)

    if mask is not None:
        kmask = mask
    else:
        band = (np.abs(KZ) < kz_band_fac * kz_cent) & (np.abs(KX) < kx_band_fac * kz_cent)
        if pos_kz_only:
            band = band & (KZ > 0)
        kmask = (w_amp > amp_thr * w_amp.max()) & band & ((np.abs(KZ) + np.abs(KX)) > 0)

    n_mask = int(kmask.sum())
    min_pts = 2 if force_dz_zero else 3

    def _empty_result():
        if not return_diagnostics:
            return 0.0, 0.0, 0.0, 0
        diag = dict(KZ=KZ, KX=KX, w=w_amp, phi=phi, mask=kmask,
                    inlier_mask=np.zeros_like(kmask), outlier_mask=kmask.copy(),
                    fitted=np.zeros_like(KZ), kz_ax=kz_ax, kx_ax=kx_ax, best_score=0.0)
        return 0.0, 0.0, 0.0, 0, diag

    if n_mask < min_pts:
        return _empty_result()

    kz_pts, kx_pts, phi_pts, w_pts = KZ[kmask], KX[kmask], phi[kmask], w_amp[kmask]
    n_pts = phi_pts.size
    A_full = (np.column_stack([kx_pts, np.ones(n_pts)]) if force_dz_zero else
              np.column_stack([kz_pts, kx_pts, np.ones(n_pts)]))

    rng = np.random.default_rng(random_state)
    best_inliers = np.zeros(n_pts, dtype=bool)
    best_score = -1.0
    for _ in range(n_iter):
        sample = rng.choice(n_pts, size=min_pts, replace=False)
        c = np.linalg.lstsq(A_full[sample], phi_pts[sample], rcond=None)[0]
        resid = _wrap_to_pi(phi_pts - A_full @ c)
        inliers = np.abs(resid) < residual_thr
        score = float(w_pts[inliers].sum()) if weight_by_amp else float(inliers.sum())
        if score > best_score:
            best_score, best_inliers = score, inliers

    n_inliers = int(best_inliers.sum())
    if n_inliers < min_pts:
        return _empty_result()

    W = w_pts[best_inliers] ** final_wls_pow
    c = np.linalg.lstsq(A_full[best_inliers] * W[:, None], phi_pts[best_inliers] * W,
                         rcond=None)[0]
    if force_dz_zero:
        dz_est, dx_est, phi_0 = 0.0, float(c[0]), float(c[1])
    else:
        dz_est, dx_est, phi_0 = float(c[0]), float(c[1]), float(c[2])

    if not return_diagnostics:
        return dz_est, dx_est, phi_0, n_inliers

    inlier_full = np.zeros_like(kmask)
    outlier_full = np.zeros_like(kmask)
    idx2d = np.where(kmask)
    inlier_full[idx2d[0][best_inliers], idx2d[1][best_inliers]] = True
    outlier_full[idx2d[0][~best_inliers], idx2d[1][~best_inliers]] = True
    fitted = KX * dx_est + KZ * dz_est + phi_0
    diag = dict(KZ=KZ, KX=KX, w=w_amp, phi=phi, mask=kmask,
                inlier_mask=inlier_full, outlier_mask=outlier_full,
                fitted=fitted, kz_ax=kz_ax, kx_ax=kx_ax, best_score=best_score)
    return dz_est, dx_est, phi_0, n_inliers, diag


def plot_phase_slice(ax, k_pts, other_pts, phi_pts, w_pts, slope_est, other_slope_est, phi_0, *,
                      xlabel='k [rad/m]', color='tab:blue', outlier=None, line_label=None):
    """1-D phase-plane-fit diagnostic: scatter the fit's own input points against
    one k-axis (with the OTHER axis's estimated contribution subtracted out, so the
    fitted plane collapses to a single line), and overlay that line -- lets a viewer
    see at a glance whether a fitted plane (wls_phase_plane_fit or
    ransac_phase_plane_fit) actually tracks its input points, the same question
    Figure 3's 1-D kz/kx slice panels answer, but built from the fit's own
    (kz_pts, kx_pts, phi_pts, w_pts) directly rather than a row/column-averaged band
    -- so this shows exactly the points a given fit was computed from (all of them
    for a WLS call, inliers-vs-rejected for a RANSAC call), not a broader summary.

    No phase unwrapping is applied (matches Figure 3 / _wls_diag's convention) --
    valid as long as the fitted plane's total phase excursion across the plotted
    k-range stays within about a cycle, which holds for this notebook's
    sub-wavelength-to-moderate displacement regime.

    Parameters
    ----------
    ax : matplotlib Axes
    k_pts, other_pts, phi_pts, w_pts : ndarray, 1-D, same shape
        This axis's k-coordinate, the other axis's k-coordinate, phase, and
        amplitude weight, for every point being plotted (e.g. KZ[mask], KX[mask],
        phi[mask], w[mask] for a kz-slice of the fit's full input; or the same
        indexed by an inlier/outlier split for a RANSAC inlier-only slice).
    slope_est, other_slope_est, phi_0 : float
        The fitted (dz_est or dx_est), the OTHER axis's fitted slope (used only to
        subtract its contribution from phi_pts), and phi_0 -- from the same fit
        call that produced k_pts's population.
    outlier : ndarray (bool) or None
        Same shape as k_pts. True marks a point as rejected (plotted as a red '×'
        instead of the amplitude-coloured dot) -- pass RANSAC's outlier split to
        show what got excluded and where it sits relative to the fitted line.
    line_label : str or None
        Legend label for the fitted line (e.g. 'WLS' / 'RANSAC'); None omits it.

    Returns
    -------
    The scatter PathCollection (for an external colorbar), or None if k_pts is empty.
    """
    if outlier is None:
        outlier = np.zeros(k_pts.shape, dtype=bool)
    inlier = ~outlier
    y = phi_pts - other_pts * other_slope_est

    sc = None
    if inlier.any():
        sc = ax.scatter(k_pts[inlier], y[inlier], c=w_pts[inlier], cmap='viridis', s=16, zorder=3)
    if outlier.any():
        ax.scatter(k_pts[outlier], y[outlier], marker='x', c='red', s=30, linewidths=1.3,
                   zorder=4, label='rejected')
    if k_pts.size:
        k_line = np.array([k_pts.min(), k_pts.max()])
        ax.plot(k_line, k_line * slope_est + phi_0, color=color, lw=1.6, zorder=5, label=line_label)
    ax.axhline(0, color='k', lw=0.5, ls='--', zorder=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('phase (other-axis contribution removed) [rad]')
    return sc