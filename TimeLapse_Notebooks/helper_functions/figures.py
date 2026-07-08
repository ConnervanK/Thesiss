"""Standardized figure saving for all thesis notebooks.

Every notebook saves its figures into the central store
``TimeLapse_Figures/<study>/<category>/`` at the repository root, as
``<PREFIX><NNN>_<sanitized-title>.png`` at 300 dpi, where the category
(Kirchhoff, Gazdag, Back-Propagation, Noisy/<technique>, or General) is
derived from the filename by :func:`classify`. The Typst thesis
recomputes the same routing from the bare filename (see ``img()`` in
``TimeLapse_LaTeX/Thesis/30_June/template.typ``).

Full protocol: ``.wiki/FIGURES_PROTOCOL.md``.

Usage (top of a notebook)::

    from helper_functions.figures import setup_autosave
    setup_autosave(study="TimeLapse_Study", prefix="TL_")

After this call every ``plt.show()`` also writes the shown figure(s)
to the central store. For headless/batch notebooks (e.g. FieldData),
save individual thesis-bound figures explicitly::

    from helper_functions.figures import save_fig
    save_fig(fig, "strat3_stage_anchored_summary",
             study="FieldData_Study", prefix="FD_")
"""

import re
from pathlib import Path

import matplotlib.pyplot as plt

# Repository root = two levels above this file
# (<root>/TimeLapse_Notebooks/helper_functions/figures.py)
FIGURES_ROOT = Path(__file__).resolve().parents[2] / "TimeLapse_Figures"

# Project standard: print-quality raster; Typst cannot import PDF and
# SVG handles the imshow-heavy panels poorly, so high-dpi PNG it is.
DPI = 300

# ---- Technique taxonomy (mirrored by img() in template.typ) -----------------
# Every figure lands in TimeLapse_Figures/<study>/<category>/, where the
# category is derived from the figure's filename so that the Typst router
# can recompute it from the bare name. Keep these rules in sync with
# template.typ and .wiki/FIGURES_PROTOCOL.md.
_NOISY = re.compile(r"nois[ey]|laplace", re.I)
_TECH_RULES = [
    ("Kirchhoff",        re.compile(r"kirchh?off|\bi?lsm\b", re.I)),
    ("Gazdag",           re.compile(r"gazdag", re.I)),
    ("Back-Propagation", re.compile(r"back-?prop|sign-?bit|time-?reversal", re.I)),
]
_PREFIX_OVERRIDES = {
    "TLC_": "Noisy/Kirchhoff",  # cleaning notebook: exclusively noisy Kirchhoff data
    "FD_":  "Gazdag",           # field strategies run on Gazdag-migrated images
}


def classify(name):
    """Filename (or prefix+stem) -> category subpath within a study folder.

    Returns one of: 'Kirchhoff', 'Gazdag', 'Back-Propagation',
    'Noisy/<technique>', or 'General'.
    """
    for pre, cat in _PREFIX_OVERRIDES.items():
        if name.startswith(pre):
            return cat
    tech = next((t for t, rx in _TECH_RULES if rx.search(name)), None)
    if tech and _NOISY.search(name):
        return f"Noisy/{tech}"
    return tech or "General"

_orig_show = None          # original plt.show, kept for chaining/idempotency
_fig_counter = [0]


def _safe_stem(title_line):
    """First line of a figure title -> filesystem-safe stem (<=60 chars)."""
    safe = re.sub(r"[^\w\s-]", "", title_line)[:60].strip().replace(" ", "_")
    return safe


def _figure_title(fig):
    if fig._suptitle:
        return fig._suptitle.get_text()
    if fig.axes:
        return fig.axes[0].get_title()
    return ""


def save_fig(fig, name, study, prefix, dpi=DPI, numbered=False, category=None):
    """Save one figure to the central store and return its path.

    ``name`` may be a descriptive stem ("strat1_consecutive_summary") or,
    with ``numbered=True``, gets the shared run counter prepended the same
    way the autosave hook numbers figures. The technique subfolder is
    derived from the filename via :func:`classify`; pass ``category=`` to
    override (must match what the Typst router computes, or the figure
    won't resolve in the thesis).
    """
    stem = name
    if numbered:
        _fig_counter[0] += 1
        stem = f"{_fig_counter[0]:03d}_{name}" if name else f"{_fig_counter[0]:03d}"
    fname = f"{prefix}{stem}.png"
    out_dir = FIGURES_ROOT / study / (category or classify(fname))
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / fname
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    print(f"  [saved] {path.relative_to(FIGURES_ROOT)}")
    return path


def setup_autosave(study, prefix, dpi=DPI):
    """Wrap plt.show so every shown figure is saved to the central store.

    Naming is identical to the historical per-notebook hooks
    (``<PREFIX><NNN>_<safe-title>.png``, counter in cell-execution order)
    because the thesis references these exact filenames. Idempotent:
    calling it again just resets the counter instead of double-wrapping.
    """
    global _orig_show
    out_dir = FIGURES_ROOT / study
    out_dir.mkdir(parents=True, exist_ok=True)

    if _orig_show is None:
        _orig_show = plt.show
    _fig_counter[0] = 0

    def _autosave_show(*args, **kwargs):
        for fignum in plt.get_fignums():
            fig = plt.figure(fignum)
            _fig_counter[0] += 1
            title_1line = _figure_title(fig).split("\n")[0]
            safe = _safe_stem(title_1line)
            stem = f"{_fig_counter[0]:03d}_{safe}" if safe else f"{_fig_counter[0]:03d}"
            fname = prefix + stem + ".png"
            cat_dir = out_dir / classify(fname)
            cat_dir.mkdir(parents=True, exist_ok=True)
            path = cat_dir / fname
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
            print(f"  [saved] {path.relative_to(FIGURES_ROOT)}")
        _orig_show(*args, **kwargs)

    plt.show = _autosave_show
    print(f"Figure autosave -> {out_dir}  (prefix={prefix!r}, dpi={dpi})")
