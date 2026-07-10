"""
Standardized plotting functions shared by the five TimeLapse study notebooks
(Diagonal_TimeLapse_Playground, FluidFlow_Playground, Resolution_Playground,
TimeLapse_Playground, Vertical_TimeLapse_Playground).

This module is presentation-only: every public function builds a Matplotlib
figure and returns it (``(fig, ax)`` or ``(fig, axes)``, or a list of such
pairs for the ``zooms=`` grid variants) but never calls ``plt.show()`` or
``fig.savefig()``. Saving is handled entirely by the existing autosave
pipeline in :mod:`helper_functions.figures` (``setup_autosave`` monkey-patches
``plt.show()`` so every displayed figure is written to the central
``TimeLapse_Figures/`` store, named from its title). Because that pipeline
derives the saved filename from the figure's suptitle/title, functions here
still set the title text from the caller-supplied string — the calling
notebook is only responsible for assembling the data, choosing the title
content, and calling ``plt.show()`` once per figure.

Style constants below are the single source of truth for colormaps, marker
styles, and amplitude-scaling conventions. Several of these standardize
choices that had drifted across notebooks during copy-paste (four different
vmax-scaling conventions, hand-toggled/sometimes-missing empty-panel hiding,
and at least two places where a computed vmax was silently never passed to
``imshow``) — see ``_compute_symmetric_vmax`` and ``_hide_unused_axes``.
"""

import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert

# ---- Colormaps (one name per data semantics, never a bare string in a call site) ----
CMAP_SIGNED = 'seismic'      # signed field / migrated / diff data (Ez, images, diffs)
CMAP_MAGNITUDE = 'inferno'   # |E| wavefield magnitude snapshots
CMAP_ENVELOPE = 'hot'        # Hilbert-envelope display mode

# ---- Line/marker colors and styles ----
COLOR_AMPLITUDE = 'steelblue'   # PSF/trace amplitude line
COLOR_ENVELOPE = 'tomato'       # PSF/trace Hilbert-envelope line
MARKER_BASELINE = 'g*'          # baseline/reference position marker (ax.plot fmt string)
MARKER_CURRENT = 'g^'           # current/shifted/timelapsed position marker
MARKER_SIZE = 10
TX_COLOR, RX_COLOR = 'navy', 'darkorange'
TX_MARKER, RX_MARKER = '^', 'v'

# ---- Placeholder panel (missing-data cells in comparison/PSF grids) ----
PLACEHOLDER_FACECOLOR = '#f5f5f5'
PLACEHOLDER_TEXT_COLOR = '#999'
PLACEHOLDER_TEXT = 'N/A'

# ---- vmax scaling: single source of truth replacing 4 legacy conventions ----
# (np.max(np.abs(x))*0.8, np.max(np.abs(x))*0.5, np.percentile(np.abs(x),100),
# np.percentile(np.abs(x),99)) are all expressible as percentile*headroom;
# defaults below reproduce the dominant "*0.8" convention exactly, since
# percentile(x, 100) == max(x).
DEFAULT_VMAX_PERCENTILE = 100
DEFAULT_VMAX_HEADROOM = 0.8
SIGNBIT_VMAX = 1.0   # sign-bit data is legitimately fixed at +-1

# ---- Typography / sizing ----
SUPTITLE_FONTSIZE = 12
TITLE_FONTSIZE = 11
LABEL_FONTSIZE = 11
TICK_FONTSIZE = 7
LEGEND_FONTSIZE = 8
PANEL_WIDTH_IN = 2.8
PANEL_HEIGHT_IN = 5.0

# Standardize titles on the em dash (cosmetic only: figures.py._safe_stem
# strips non-word characters, so this has no effect on saved filenames).
DASH = '—'


# =============================================================================
# Private helpers
# =============================================================================

def _grid_shape(n_items, ncols):
    """(n_rows, n_cols) for a grid holding n_items panels, ncols wide."""
    ncols = min(ncols, n_items) if n_items else ncols
    nrows = math.ceil(n_items / ncols)
    return nrows, ncols


def _hide_unused_axes(axes, n_used):
    """Hide any panel past n_used. Called unconditionally by every grid
    function so an unfilled grid can never leave a blank panel visible
    (this was previously hand-toggled per cell, and in at least one notebook
    the hide call had been commented out, leaving a visible empty panel in
    saved thesis figures)."""
    flat = np.atleast_1d(axes).ravel()
    for ax in flat[n_used:]:
        ax.set_visible(False)


def _compute_symmetric_vmax(data, percentile=DEFAULT_VMAX_PERCENTILE,
                             headroom=DEFAULT_VMAX_HEADROOM):
    """Symmetric vmax for a signed colormap: percentile(|data|, percentile) * headroom.

    Args:
        data: array or iterable of arrays; combined if iterable.
        percentile (float): Percentile of |data| to use as the base scale.
            100 (default) reproduces plain np.max(np.abs(data)).
        headroom (float): Multiplier applied after the percentile, e.g. 0.8
            reproduces the dominant "*0.8" convention used across the
            notebooks' migrated-image grids.

    Returns:
        float: vmax such that imshow(..., vmin=-vmax, vmax=vmax) is symmetric.
    """
    if isinstance(data, (list, tuple)):
        combined = np.concatenate([np.abs(np.asarray(d)).ravel() for d in data])
    else:
        combined = np.abs(np.asarray(data)).ravel()
    return float(np.percentile(combined, percentile) * headroom)


def _placeholder_panel(ax, text=PLACEHOLDER_TEXT, facecolor=PLACEHOLDER_FACECOLOR,
                        text_color=PLACEHOLDER_TEXT_COLOR, fontsize=12):
    """Render an empty/'N/A' panel for missing data in comparison-style grids."""
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha='center', va='center',
             fontsize=fontsize, color=text_color)
    ax.set_facecolor(facecolor)
    ax.set_xticks([])
    ax.set_yticks([])


def _default_figsize(n_cols, n_rows=1, panel_w=PANEL_WIDTH_IN, panel_h=PANEL_HEIGHT_IN):
    return (panel_w * n_cols, panel_h * n_rows)


def _resolve_vmax(data, vmax, vmax_percentile, vmax_headroom):
    if vmax is not None:
        return vmax
    return _compute_symmetric_vmax(data, percentile=vmax_percentile, headroom=vmax_headroom)


# =============================================================================
# 1. Domain geometry scaffold
# =============================================================================

def plot_domain_geometry(domain_x, domain_y, y_surface, pml_t, *, x_src=None,
                          x_rx=None, tx_stride=5, rx_offset=0.0, eps_r=None,
                          ax=None, figsize=(14, 4)):
    """Draw the ice/air/PML cross-section + Tx/Rx markers common to all five studies.

    Depth convention (matches every audited notebook): y = 0 is the ice
    surface, positive y is depth into the ice, negative y is height into the
    air. Domain top is at ``d_top = -(domain_y - y_surface)``.

    Draws, in order: ice region, air region, PML on the left/right/top edges
    only (no bottom PML, matching every gprMax input file in this project),
    a domain border, the y=0 surface line, and Tx/Rx scatter markers every
    ``tx_stride``-th trace. Does NOT set aspect ratio, legend, grid, title,
    axis limits, tight_layout, or call plt.show() — scenario-specific markers
    (plasma-colored circles, broken_barh rows, annotated cylinders, etc.)
    differ enough between notebooks that they stay in the calling notebook,
    layered on the returned axes, followed by the notebook's own
    ax.set_xlim/ylim, ax.legend(), ax.grid(), ax.set_title(),
    plt.tight_layout(), and plt.show().

    Args:
        domain_x (float): Domain width [m].
        domain_y (float): Domain height [m] (ice + air + PML).
        y_surface (float): Ice thickness below the surface shown in this
            cross-section [m] (i.e. how far down the ice region is drawn).
        pml_t (float): PML thickness [m].
        x_src (ndarray, optional): 1-D Tx x-positions [m]. Pass None to skip
            Tx/Rx markers entirely.
        x_rx (ndarray, optional): 1-D Rx x-positions [m]; defaults to
            ``x_src + rx_offset`` when x_src is given and x_rx is None.
        tx_stride (int): Draw every tx_stride-th Tx/Rx pair (avoids
            overplotting for profiles with many traces).
        rx_offset (float): Tx-Rx offset [m], used only when x_rx is None.
        eps_r (float, optional): Ice relative permittivity, included in the
            "Ice" region label when given.
        ax (Axes, optional): Draw into this axes instead of creating a new
            figure.
        figsize (tuple): Only used when ax is None.

    Returns:
        (fig, ax): fig is None when ax was supplied by the caller.
    """
    from matplotlib.patches import Rectangle

    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    air_h = domain_y - y_surface
    d_top = -air_h

    # Material regions
    ax.add_patch(Rectangle((0, 0), domain_x, y_surface,
                            facecolor='#cce5ff', edgecolor='none'))          # ice
    ax.add_patch(Rectangle((0, d_top), domain_x, air_h,
                            facecolor='#f2faff', edgecolor='none'))          # air

    # PML (left, right, top only)
    pml_kw = dict(facecolor='#d0d0d0', edgecolor='#888', hatch='///', alpha=0.80, lw=0.4)
    ax.add_patch(Rectangle((0, d_top), pml_t, domain_y, **pml_kw))
    ax.add_patch(Rectangle((domain_x - pml_t, d_top), pml_t, domain_y, **pml_kw))
    ax.add_patch(Rectangle((0, d_top), domain_x, pml_t, **pml_kw))

    ax.add_patch(Rectangle((0, d_top), domain_x, domain_y,
                            fill=False, edgecolor='black', lw=1.5))
    ax.axhline(0, color='steelblue', lw=1.0, ls='--', alpha=0.7)

    # Tx / Rx markers
    if x_src is not None:
        x_src = np.asarray(x_src)
        x_rx_arr = np.asarray(x_rx) if x_rx is not None else x_src + rx_offset
        ax.scatter(x_src[::tx_stride], np.zeros_like(x_src[::tx_stride]),
                   marker=TX_MARKER, s=10, color=TX_COLOR, zorder=6, label='Tx')
        ax.scatter(x_rx_arr[::tx_stride], np.zeros_like(x_rx_arr[::tx_stride]),
                   marker=RX_MARKER, s=10, color=RX_COLOR, zorder=6, label='Rx')

    # Region text
    ice_label = f'Ice  ({DASH} ε_r = {eps_r})' if eps_r is not None else 'Ice'
    ax.text(domain_x / 2, d_top + air_h / 2, 'Air  (ε_r = 1)',
            ha='center', va='center', fontsize=10, color='#444')
    ax.text(domain_x / 2, y_surface / 2, ice_label,
            ha='center', va='center', fontsize=10, color='#1a4a6e')
    ax.text(pml_t / 2, d_top + domain_y / 2, 'PML',
            ha='center', va='center', fontsize=8, color='#555', rotation=90)
    ax.text(domain_x - pml_t / 2, d_top + domain_y / 2, 'PML',
            ha='center', va='center', fontsize=8, color='#555', rotation=90)

    return fig, ax


# =============================================================================
# 2. B-scan panel grids
# =============================================================================

def plot_bscan_grid(data, x_traces, time_ns, *, ncols=None, title=None,
                     cmap=CMAP_SIGNED, vmax=None, vmax_percentile=DEFAULT_VMAX_PERCENTILE,
                     vmax_headroom=DEFAULT_VMAX_HEADROOM, shared_colorbar=True,
                     colorbar_label='Ez [V/m]', markers=None, axes=None, figsize=None):
    """Row (or grid) of imshow B-scan panels sharing one vmin/vmax.

    Used for raw / background-subtracted / noisy / sign-bit B-scan grids
    (the single most duplicated pattern across all five notebooks), and for
    the tapering-diagnostic imshow triptych via ``shared_colorbar=False``.

    Args:
        data (list[tuple[str, ndarray]]): (title, bscan) pairs, each bscan of
            shape (n_t, n_traces), time axis first.
        x_traces (ndarray): 1-D along-profile axis [m].
        time_ns (ndarray): 1-D time axis [ns].
        ncols (int, optional): Panels per row; defaults to len(data) (a
            single row, the dominant layout in every audited notebook).
        title (str, optional): Figure suptitle.
        cmap (str): Colormap for every panel.
        vmax (float, optional): Explicit vmin=-vmax/vmax=vmax override — set
            this to SIGNBIT_VMAX for sign-bit data (legitimately fixed +-1).
        vmax_percentile, vmax_headroom (float): Passed to
            _compute_symmetric_vmax when vmax is None; computed once across
            the combined data of all panels.
        shared_colorbar (bool): True draws one fig.colorbar spanning all axes
            (the dominant convention). False skips the colorbar and instead
            computes/display each panel's own vmax independently (matches
            the tapering diagnostic's B-scan triptych).
        colorbar_label (str): Label for the shared colorbar.
        markers (list[callable] or None): Optional per-panel callbacks
            ``marker_fn(ax, i)`` invoked after imshow for panel i, e.g. to
            draw axvline/axhline reference lines — this content genuinely
            differs by notebook (axis of motion under study) so it stays a
            caller-supplied callback.
        axes (ndarray of Axes, optional): Draw into this pre-existing 1-D
            axes array instead of creating a new figure (used to compose a
            B-scan row into a larger multi-row figure, e.g. paired with a
            spectrum row from plot_spectrum_grid).
        figsize (tuple, optional): Defaults to _default_figsize(len(data)).

    Returns:
        (fig, axes): fig is None when axes was supplied by the caller.
    """
    n = len(data)
    ncols = ncols or n
    extent = [x_traces[0], x_traces[-1], time_ns[-1], 0]

    fig = None
    if axes is None:
        figsize = figsize or _default_figsize(ncols)
        fig, axes = plt.subplots(1, ncols, figsize=figsize, facecolor='w', sharey=True)
    axes = np.atleast_1d(axes)

    if shared_colorbar:
        vmax_shared = _resolve_vmax([d for _, d in data], vmax, vmax_percentile, vmax_headroom)

    im = None
    for i, (ax, (panel_title, d)) in enumerate(zip(axes, data)):
        panel_vmax = vmax_shared if shared_colorbar else _resolve_vmax(
            d, vmax, vmax_percentile, vmax_headroom)
        im = ax.imshow(d, aspect='auto', cmap=cmap, extent=extent,
                        interpolation='nearest', vmin=-panel_vmax, vmax=panel_vmax)
        ax.set_title(panel_title, fontsize=TITLE_FONTSIZE)
        ax.set_xlabel('x [m]', fontsize=LABEL_FONTSIZE)
        if ax is axes[0]:
            ax.set_ylabel('Time [ns]', fontsize=LABEL_FONTSIZE)
        ax.grid(linestyle='-.', alpha=0.4)
        if markers and i < len(markers) and markers[i] is not None:
            markers[i](ax, i)

    _hide_unused_axes(axes, n)

    if title and fig is not None:
        fig.suptitle(title, fontsize=SUPTITLE_FONTSIZE, fontweight='bold', y=1.02)
    if fig is not None:
        fig.tight_layout()

    # Added after tight_layout (not before): fig.colorbar(ax=...) shrinks the
    # panel axes' gridspec to carve out its own space, but a later
    # tight_layout() call doesn't know about that reserved space and expands
    # the panels back into it, overlapping the last one. Adding the colorbar
    # last keeps its reserved space intact.
    if shared_colorbar and im is not None:
        fig_for_cbar = fig if fig is not None else axes[0].figure
        fig_for_cbar.colorbar(im, ax=list(axes[:n]), fraction=0.02, pad=0.02,
                               label=colorbar_label)

    return fig, axes


# =============================================================================
# 3 & 4. Migrated / time-lapse-difference grids (with optional zoom pairs)
# =============================================================================

def _plot_image_grid(images, extent, *, title, ncols, cmap, origin, vmin,
                      vmax, vmax_percentile, vmax_headroom, marker_x, marker_y,
                      marker_x_baseline, y_surface, xlim, ylim, figsize):
    """Shared body for plot_migrated_grid / plot_wavefield_grid (one figure).

    vmin=None means "symmetric about zero" (-vmax_shared) — the convention
    for signed field data. Pass vmin=0 explicitly for magnitude data (e.g.
    |E| wavefield snapshots), where a symmetric scale would waste half the
    colormap on negative values that never occur.
    """
    if isinstance(images, dict):
        items = list(images.items())
    else:
        items = list(images)
    n = len(items)
    nrows, ncols = _grid_shape(n, ncols)
    figsize = figsize or _default_figsize(ncols, nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)
    axes_flat = np.atleast_1d(axes).ravel()

    vmax_shared = _resolve_vmax([img for _, img in items], vmax, vmax_percentile, vmax_headroom)
    vmin_shared = vmin if vmin is not None else -vmax_shared

    def _marker_value(spec, i, label):
        if callable(spec):
            return spec(i, label)
        if isinstance(spec, dict):
            return spec.get(label)
        return spec

    for i, (ax, (label, img)) in enumerate(zip(axes_flat, items)):
        ax.imshow(img, aspect='auto', cmap=cmap, vmin=vmin_shared, vmax=vmax_shared,
                  extent=extent, origin=origin)

        mx = _marker_value(marker_x, i, label) if marker_x is not None else None
        my = _marker_value(marker_y, i, label) if marker_y is not None else None
        if mx is not None and my is not None:
            ax.plot(mx, my, MARKER_CURRENT, ms=MARKER_SIZE, zorder=5, label='current')
            if marker_x_baseline is not None:
                mxb = _marker_value(marker_x_baseline, i, label)
                ax.plot(mxb, my, MARKER_BASELINE, ms=MARKER_SIZE, zorder=5, label='baseline')

        if y_surface is not None:
            ax.axhline(y_surface, color='cyan', lw=0.8, ls='--')

        ax.set_title(label, fontsize=TITLE_FONTSIZE)
        ax.set_xlabel('x [m]', fontsize=LABEL_FONTSIZE)
        if i % ncols == 0:
            ax.set_ylabel('Depth z [m]' if origin == 'upper' else 'y [m]',
                          fontsize=LABEL_FONTSIZE)
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)

    _hide_unused_axes(axes_flat, n)
    if n and marker_x is not None:
        axes_flat[0].legend(fontsize=LEGEND_FONTSIZE, loc='upper right')

    if title:
        fig.suptitle(title, fontsize=SUPTITLE_FONTSIZE, fontweight='bold', y=1.01)
    fig.tight_layout()

    return fig, axes


def plot_migrated_grid(images, extent, *, title, ncols, marker_x=None, marker_z=None,
                        marker_x_baseline=None, cmap=CMAP_SIGNED,
                        vmax=None, vmax_percentile=DEFAULT_VMAX_PERCENTILE,
                        vmax_headroom=DEFAULT_VMAX_HEADROOM, origin='upper',
                        zooms=None, figsize=None):
    """Grid of migrated (or time-lapse-difference) image panels with baseline/
    current position markers, replacing both a notebook's duplicated
    "full extent" and "zoomed" cells with one call.

    Args:
        images (dict[str, ndarray] or list[tuple[str, ndarray]]): label ->
            image of shape (n_z, n_x). For time-lapse-difference grids, pass
            the already-differenced images (e.g.
            ``{label: migrated[label] - migrated['Baseline'] for label in labels}``
            — that subtraction is a two-line dict comprehension and stays in
            the notebook).
        extent (list[float]): [x0, x1, z1, z0] (depth axis reversed), i.e.
            ``[x_traces[0], x_traces[-1], z_img[-1], z_img[0]]``.
        title (str): Figure suptitle (caller pre-formats technique/f_c/
            aperture metadata into this string).
        ncols (int): Panels per row.
        marker_x, marker_z: Per-panel current-position marker coordinates.
            Each accepts a flat value, a dict keyed like images, or a
            callable ``(i, label) -> value`` (needed when marker depth
            differs by column, e.g. comparison grids where a back-prop
            column uses a different z reference than migration columns).
        marker_x_baseline: Same accepted types as marker_x; if given, plots a
            second MARKER_BASELINE star using marker_z as the shared depth.
        cmap (str): Colormap.
        vmax, vmax_percentile, vmax_headroom: see _compute_symmetric_vmax;
            vmax computed once across all panels when not given explicitly.
        origin (str): 'upper' for migrated/diff grids (z increases downward
            in the extent) — the only value used by this function; kept as
            an explicit parameter for symmetry with plot_wavefield_grid,
            whose gprMax-grid data uses origin='lower'.
        zooms (list[tuple[tuple, tuple, str]], optional): List of
            (xlim, ylim, suffix) specs. When given, produces one (fig, axes)
            per spec (title gets " " + suffix appended) instead of a single
            full-extent figure — generalizes the duplicated full+zoomed cell
            pairs into one call, mirroring the zoom-pair loop already used in
            Resolution_Playground.ipynb's LSM/ILSM cells.
        figsize (tuple, optional): Defaults from panel count.

    Returns:
        (fig, axes) if zooms is None, else list[(fig, axes)], one per spec.
    """
    if zooms is None:
        return _plot_image_grid(images, extent, title=title, ncols=ncols, cmap=cmap,
                                 origin=origin, vmin=None, vmax=vmax,
                                 vmax_percentile=vmax_percentile,
                                 vmax_headroom=vmax_headroom, marker_x=marker_x,
                                 marker_y=marker_z, marker_x_baseline=marker_x_baseline,
                                 y_surface=None, xlim=None, ylim=None, figsize=figsize)

    results = []
    for xlim, ylim, suffix in zooms:
        panel_title = f'{title} {suffix}'.strip() if suffix else title
        results.append(_plot_image_grid(
            images, extent, title=panel_title, ncols=ncols, cmap=cmap, origin=origin,
            vmin=None, vmax=vmax, vmax_percentile=vmax_percentile, vmax_headroom=vmax_headroom,
            marker_x=marker_x, marker_y=marker_z, marker_x_baseline=marker_x_baseline,
            y_surface=None, xlim=xlim, ylim=ylim, figsize=figsize))
    return results


# =============================================================================
# 5. Back-propagation wavefield snapshot grids
# =============================================================================

def plot_wavefield_grid(frames, extent, *, field, title, ncols, y_surface=None,
                         marker_x=None, marker_y=None, marker_x_baseline=None,
                         vmax=None, vmax_percentile=DEFAULT_VMAX_PERCENTILE,
                         vmax_headroom=DEFAULT_VMAX_HEADROOM, zooms=None, figsize=None):
    """Grid of back-propagation wavefield snapshots (|E| magnitude or Ez signed
    field), origin='lower' (gprMax's y-axis runs bottom-up — the opposite
    convention from plot_migrated_grid's origin='upper').

    This is the structural fix for two bugs found during the audit: a
    computed-but-never-applied vmax in one notebook's back-propagation cell,
    and a commented-out empty-panel-hide call in another notebook that left
    a visible blank panel in saved thesis figures. Both become impossible
    here since vmin/vmax always flows into the single imshow call site and
    _hide_unused_axes always runs.

    Args:
        frames (dict[str, ndarray]): label -> 2-D array, e.g.
            ``{label: frame['mag'] for label, frame in focus_frames.items()}``.
        extent (list[float]): [0, domain_x, 0, domain_y], the full gprMax
            domain extent (not reversed, unlike migrated-grid extents).
        field ({'magnitude', 'signed'}): 'magnitude' selects CMAP_MAGNITUDE
            (matches |E| panels, vmin=0); 'signed' selects CMAP_SIGNED with a
            symmetric vmin/vmax (matches Ez / Ez-diff panels). Selecting the
            colormap via this one semantic flag (rather than a raw cmap=
            string) prevents a copy-pasted cell from mixing the two up.
        title (str): Figure suptitle.
        ncols (int): Panels per row.
        y_surface (float, optional): If given, draws
            ``ax.axhline(y_surface, color='cyan', ls='--')`` on every panel
            (the ice-surface marker used on focus-frame grids; omitted on
            diff-frame grids in the audited notebooks, hence optional).
        marker_x, marker_y, marker_x_baseline: same semantics as
            plot_migrated_grid's marker_x/marker_z/marker_x_baseline.
        vmax, vmax_percentile, vmax_headroom: see _compute_symmetric_vmax.
        zooms: same (xlim, ylim, suffix) list contract as plot_migrated_grid.
        figsize (tuple, optional): Defaults from panel count.

    Returns:
        (fig, axes) or list[(fig, axes)] (same zoom-spec contract as
        plot_migrated_grid).
    """
    cmap = CMAP_MAGNITUDE if field == 'magnitude' else CMAP_SIGNED
    vmin = 0 if field == 'magnitude' else None

    if zooms is None:
        return _plot_image_grid(frames, extent, title=title, ncols=ncols, cmap=cmap,
                                 origin='lower', vmin=vmin, vmax=vmax,
                                 vmax_percentile=vmax_percentile,
                                 vmax_headroom=vmax_headroom, marker_x=marker_x,
                                 marker_y=marker_y, marker_x_baseline=marker_x_baseline,
                                 y_surface=y_surface, xlim=None, ylim=None, figsize=figsize)

    results = []
    for xlim, ylim, suffix in zooms:
        panel_title = f'{title} {suffix}'.strip() if suffix else title
        results.append(_plot_image_grid(
            frames, extent, title=panel_title, ncols=ncols, cmap=cmap, origin='lower',
            vmin=vmin, vmax=vmax, vmax_percentile=vmax_percentile, vmax_headroom=vmax_headroom,
            marker_x=marker_x, marker_y=marker_y, marker_x_baseline=marker_x_baseline,
            y_surface=y_surface, xlim=xlim, ylim=ylim, figsize=figsize))
    return results


# =============================================================================
# 6. Multi-method comparison grid
# =============================================================================

def plot_method_comparison_grid(images, extent, methods, row_labels, *, title,
                                 envelope=False, marker_x=None, marker_z=None,
                                 xlim=None, ylim=None, vmax_percentile=100, figsize=None):
    """n_scenarios x n_methods grid; each cell is an image, None (rendered as
    an N/A placeholder), signed-amplitude, or Hilbert-envelope display.

    Hilbert-envelope computation happens inside this function (gated by
    envelope=True) rather than in the caller: it is a pure, parameter-free
    display transform identical at every call site, and keeping it here is
    what fixes the "dead vmax" bug found in both audited comparison-grid
    cells, where vmin/vmax were computed but never actually passed to
    imshow — here there is exactly one imshow call site per panel and it
    always receives vmin/vmax.

    Args:
        images (list[list[ndarray or None]]): images[i][j], i = scenario row,
            j = method column. Assembling this structure (including any
            reprojection of back-propagation data onto the migration grid)
            is notebook-side data wrangling and stays in the notebook — this
            function only renders the assembled 2-D list.
        extent (list[float]): Shared [x0, x1, z1, z0] extent for every panel.
        methods (list[str]): Column titles.
        row_labels (list[str]): Row labels (e.g. scenario names).
        title (str): Figure suptitle.
        envelope (bool): False -> raw signed data, CMAP_SIGNED, vmin/vmax
            symmetric about 0. True -> np.abs(hilbert(img, axis=0)),
            CMAP_ENVELOPE, vmin=0.
        marker_x, marker_z: per-(row,col) marker position. Each accepts a
            flat value or a callable ``(i, j) -> value`` (marker depth
            differs by method column in the audited code, e.g. back-prop vs.
            migration columns use different z references).
        xlim, ylim (tuple, optional): Shared axis window for every panel.
        vmax_percentile (float): Percentile passed to _compute_symmetric_vmax
            (headroom=1.0), computed once per method column from all
            available (non-None) images in that column so every scenario for
            a given method shares one color scale. Colorscales are not
            shared across columns/methods, since different methods can
            produce genuinely different amplitude scales and sharing across
            columns lets the largest method dominate the others.
        figsize (tuple, optional): Defaults from grid shape.

    Returns:
        (fig, axes)
    """
    n_rows, n_cols = len(images), len(methods)
    figsize = figsize or _default_figsize(n_cols, n_rows, panel_w=4.5, panel_h=3.0)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, sharex=True, sharey=True)
    axes = np.atleast_2d(axes)

    def _marker_value(spec, i, j):
        if callable(spec):
            return spec(i, j)
        return spec

    if envelope:
        display_images = [[np.abs(hilbert(img, axis=0)) if img is not None else None
                            for img in row_imgs] for row_imgs in images]
        cmap, vmin = CMAP_ENVELOPE, 0
    else:
        display_images = images
        cmap, vmin = CMAP_SIGNED, None

    col_vmax = []
    for j in range(n_cols):
        available = [display_images[i][j] for i in range(n_rows)
                     if display_images[i][j] is not None]
        col_vmax.append(_compute_symmetric_vmax(
            available, percentile=vmax_percentile, headroom=1.0) if available else 1.0)

    for i, row_imgs in enumerate(display_images):
        for j, img in enumerate(row_imgs):
            ax = axes[i, j]
            if img is None:
                _placeholder_panel(ax)
                continue
            vmax_col = col_vmax[j]
            panel_vmin = -vmax_col if vmin is None else vmin
            ax.imshow(img, aspect='auto', cmap=cmap, vmin=panel_vmin, vmax=vmax_col,
                      extent=extent, origin='upper')

            mx = _marker_value(marker_x, i, j) if marker_x is not None else None
            mz = _marker_value(marker_z, i, j) if marker_z is not None else None
            if mx is not None and mz is not None:
                ax.plot(mx, mz, MARKER_CURRENT, ms=MARKER_SIZE, zorder=5)

            if xlim is not None:
                ax.set_xlim(xlim)
            if ylim is not None:
                ax.set_ylim(ylim)
            if i == 0:
                ax.set_title(methods[j], fontsize=TITLE_FONTSIZE)
            if j == 0:
                ax.set_ylabel(row_labels[i], fontsize=LABEL_FONTSIZE)

    fig.suptitle(
        f'{title}  {DASH}  ★ = baseline,  ▲ = timelapsed',
        fontsize=SUPTITLE_FONTSIZE, fontweight='bold', y=1.01)
    fig.tight_layout()

    return fig, axes


# =============================================================================
# 7. PSF line-profile grids
# =============================================================================

def plot_psf_grid(profiles, axis_values, methods, row_labels, *, orientation,
                   marker_positions=None, title, xlim=None, ylim=None, figsize=None):
    """n_scenarios x n_methods grid of amplitude + Hilbert-envelope line
    profiles, mirroring plot_method_comparison_grid's layout.

    Args:
        profiles (list[list[ndarray or None]]): profiles[i][j], already
            peak-normalised 1-D slices. Extraction (which needs
            method-dependent index lookups) stays in the notebook.
        axis_values (ndarray): Shared spatial axis paired with each profile
            (x_traces for orientation='lateral', z_img for 'vertical').
        methods (list[str]): Column titles.
        row_labels (list[str]): Row labels.
        orientation ({'lateral', 'vertical'}): 'lateral' plots
            ax.plot(axis_values, profile) with column titles on row 0 and
            ax.set_ylabel(row_label) on column 0. 'vertical' plots
            ax.plot(profile, axis_values) (axes swapped, since the spatial
            axis is now the y-axis) with the row label drawn via rotated
            ax.text instead of set_ylabel. Required and explicit — forcing
            callers to pre-swap arrays would silently break the row-label
            placement logic.
        marker_positions (list[list[list[float]]], optional): Per-(row,col)
            list of reference positions, drawn as axvline in 'lateral' mode
            or axhline in 'vertical' mode.
        title (str): Figure suptitle.
        xlim, ylim (tuple, optional): Shared axis window for every panel.
        figsize (tuple, optional): Defaults from grid shape (orientation
            determines whether panels are wider or taller).

    Returns:
        (fig, axes)
    """
    n_rows, n_cols = len(profiles), len(methods)
    if figsize is None:
        figsize = (_default_figsize(n_cols, n_rows, panel_w=4.0, panel_h=2.6) if
                   orientation == 'lateral' else
                   _default_figsize(n_cols, n_rows, panel_w=2.6, panel_h=4.0))

    sharey, sharex = (True, False) if orientation == 'lateral' else (False, True)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, sharex=sharex, sharey=sharey)
    axes = np.atleast_2d(axes)

    for i, row_profiles in enumerate(profiles):
        for j, profile in enumerate(row_profiles):
            ax = axes[i, j]
            if profile is None:
                _placeholder_panel(ax)
                continue

            envelope = np.abs(hilbert(profile))
            if orientation == 'lateral':
                ax.plot(axis_values, profile, color=COLOR_AMPLITUDE, lw=1.0, label='amplitude')
                ax.plot(axis_values, envelope, color=COLOR_ENVELOPE, lw=1.0, label='envelope')
            else:
                ax.plot(profile, axis_values, color=COLOR_AMPLITUDE, lw=1.0, label='amplitude')
                ax.plot(envelope, axis_values, color=COLOR_ENVELOPE, lw=1.0, label='envelope')

            if marker_positions is not None and marker_positions[i][j] is not None:
                for pos in marker_positions[i][j]:
                    if orientation == 'lateral':
                        ax.axvline(pos, color='green', lw=0.8, ls='--')
                    else:
                        ax.axhline(pos, color='green', lw=0.8, ls='--')

            if xlim is not None:
                ax.set_xlim(xlim)
            if ylim is not None:
                ax.set_ylim(ylim)
            if i == 0:
                ax.set_title(methods[j], fontsize=TITLE_FONTSIZE)
            if j == 0:
                if orientation == 'lateral':
                    ax.set_ylabel(row_labels[i], fontsize=LABEL_FONTSIZE, rotation=0,
                                  ha='right', va='center')
                else:
                    ax.text(-0.45, 0.5, row_labels[i], transform=ax.transAxes,
                            ha='center', va='center', fontsize=LABEL_FONTSIZE, rotation=90)

    fig.suptitle(title, fontsize=SUPTITLE_FONTSIZE, fontweight='bold', y=1.01)
    fig.tight_layout()

    return fig, axes


# =============================================================================
# 8. Frequency-spectrum comparison plots
# =============================================================================

def plot_spectrum_grid(spectra, freqs_ghz, *, title, row_colors=('C0', 'C1'),
                        xlim=None, axes=None, figsize=None):
    """Row(s) of |FFT| line-plot panels (e.g. clean vs. noisy spectra, or a
    single difference-spectrum row).

    FFT computation stays in the notebook (it needs dt_ns/n_t from the raw
    B-scans), but must use ``np.fft.rfft(d, axis=0, norm='forward')`` so
    amplitudes are directly comparable across scenarios and notebooks —
    this standardizes a real discrepancy found during the audit, where one
    notebook's spectrum cells omitted norm='forward' entirely.

    Args:
        spectra (list[list[tuple[str, ndarray]]]): One inner list per row,
            each (label, spectrum) pair pre-computed by the caller as
            ``np.abs(np.fft.rfft(d, axis=0, norm='forward')).mean(axis=1)``.
        freqs_ghz (ndarray): Shared x-axis, ``np.fft.rfftfreq(n_t, d=dt_ns)``.
        title (str): Figure suptitle.
        row_colors (tuple[str]): Line color per row; default reproduces the
            clean/noisy 'C0'/'C1' convention. Pass a single-color tuple for a
            1-row difference-spectrum plot.
        xlim (tuple, optional): Shared frequency window, e.g. (0, 7).
        axes (ndarray of Axes, optional): Draw into this pre-existing 2-D (or
            1-D for a single row) axes array instead of creating a new
            figure (used to compose a spectrum row underneath a B-scan row
            from plot_bscan_grid into one combined figure).
        figsize (tuple, optional): Defaults from panel count.

    Returns:
        (fig, axes): fig is None when axes was supplied by the caller.
    """
    n_rows = len(spectra)
    n_cols = len(spectra[0]) if n_rows else 0

    fig = None
    if axes is None:
        figsize = figsize or _default_figsize(n_cols, n_rows, panel_w=2.5, panel_h=3.0)
        fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, sharex=True, sharey='row')
    axes2d = np.atleast_2d(axes)

    for i, row in enumerate(spectra):
        color = row_colors[i] if i < len(row_colors) else row_colors[-1]
        for j, (label, spec) in enumerate(row):
            ax = axes2d[i, j]
            ax.plot(freqs_ghz, spec, lw=0.9, color=color)
            if i == 0:
                ax.set_title(label, fontsize=TITLE_FONTSIZE)
            if j == 0:
                ax.set_ylabel('|FFT| [-]', fontsize=LABEL_FONTSIZE)
            if i == n_rows - 1:
                ax.set_xlabel('Frequency [GHz]', fontsize=LABEL_FONTSIZE)
            if xlim is not None:
                ax.set_xlim(xlim)

    if title and fig is not None:
        fig.suptitle(title, fontsize=SUPTITLE_FONTSIZE, fontweight='bold', y=1.02)
    if fig is not None:
        fig.tight_layout()

    return fig, axes


# =============================================================================
# 9. Tapering / t0-shift diagnostic
# =============================================================================

def plot_trace_comparison(traces, time_ns, *, title, taper_curve=None, vline=None,
                           figsize=(16, 4)):
    """1x3 triptych: raw-vs-tapered, tapered-vs-shifted, and the taper weight
    curve. Kept as its own small function (not folded into the N-panel grid
    function) since the three panels show genuinely different content pairs.

    Args:
        traces (list[list[tuple[str, ndarray, str]]]): traces[panel_idx] =
            list of (label, y_values, color) line specs to plot in that
            panel. Only used for panels 0 and 1 (the raw/tapered/shifted
            trace comparisons); pass an empty list for panel 2 when using
            taper_curve instead.
        time_ns (ndarray): Shared x-axis for panels 0-1.
        title (str): Figure suptitle.
        taper_curve (ndarray, optional): If given, drawn in panel 2 as
            ``ax.plot(time_ns, taper_curve, 'k', lw=1.5)`` with
            ylim=(0, 1.05) and ylabel 'Weight [-]'.
        vline (float, optional): Shared reference x-position (e.g.
            taper_end_ns) drawn as axvline on every panel.
        figsize (tuple): Figure size.

    Returns:
        (fig, axes): axes is a length-3 array.
    """
    fig, axes = plt.subplots(1, 3, figsize=figsize)

    for panel_idx in range(2):
        ax = axes[panel_idx]
        for label, y, color in traces[panel_idx]:
            ax.plot(time_ns, y, color, lw=0.8, label=label)
        ax.set_xlabel('Time [ns]', fontsize=LABEL_FONTSIZE)
        ax.legend(fontsize=LEGEND_FONTSIZE)
        if vline is not None:
            ax.axvline(vline, color='grey', lw=0.8, ls=':')

    ax2 = axes[2]
    if taper_curve is not None:
        ax2.plot(time_ns, taper_curve, 'k', lw=1.5)
        ax2.set_ylim(0, 1.05)
        ax2.set_ylabel('Weight [-]', fontsize=LABEL_FONTSIZE)
    ax2.set_xlabel('Time [ns]', fontsize=LABEL_FONTSIZE)
    if vline is not None:
        ax2.axvline(vline, color='grey', lw=0.8, ls=':')

    fig.suptitle(title, fontsize=SUPTITLE_FONTSIZE, fontweight='bold')
    fig.tight_layout()

    return fig, axes
