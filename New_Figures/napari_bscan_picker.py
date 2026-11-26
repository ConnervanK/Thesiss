"""napari picker for Natural_Movement.ipynb: paint, for every consecutive profile pair, the part of
one migration method's section that the cross-spectrum fit should use.

The notebook runs this script in the gdp environment (which has napari): it writes the image stack,
this script opens it with the saved picks, and the picks are saved back when the window closes.

Usage: python napari_bscan_picker.py <stack.npz> <picks.npz> [--check]
    --check  build the viewer without showing it, run the key actions once and exit (a test)

In the viewer: the slider steps through the pairs; paint the 'pick' layer (paint brush, fill
bucket, eraser). Shift-C copies the previous pair's pick to this pair, Shift-S saves.
"""
import sys
from pathlib import Path

import numpy as np
import napari

stack_path, picks_path = Path(sys.argv[1]), Path(sys.argv[2])
check = '--check' in sys.argv[3:]
with np.load(stack_path) as f:
    stack = {k: f[k] for k in f.files}
with np.load(picks_path) as f:
    saved = {k: f[k] for k in f.files}

imgs, runs, method = stack['imgs'].astype(np.float32), stack['runs'], str(stack['method'])
scale, translate = tuple(stack['scale']), tuple(stack['translate'])   # (dz, dx), (z0, x0) [m]
wz0, wz1, wx0, wx1 = stack['window']                                  # analysis window [m]
pairs = list(zip(runs[:-1], runs[1:]))
picks = saved['picks'].astype(np.uint8)
assert picks.shape == (len(pairs), *imgs.shape[1:]), 'picks do not match the image stack'

# Display: every slice scaled to its peak inside the analysis window, so one colour scale fits all
z_ax = translate[0] + np.arange(imgs.shape[1]) * scale[0]
x_ax = translate[1] + np.arange(imgs.shape[2]) * scale[1]
in_win = np.ix_((z_ax >= wz0) & (z_ax <= wz1), (x_ax >= wx0) & (x_ax <= wx1))


def normalised(s):
    peak = np.array([np.abs(sl[in_win]).max() for sl in s])[:, None, None]
    return s / np.where(peak > 0, peak, 1)


viewer = napari.Viewer(title=f'B-scan picks: {method}', show=not check)
geo = dict(scale=(1, *scale), translate=(0, *translate))
img_kw = dict(colormap='RdBu_r', contrast_limits=(-1, 1), **geo)
viewer.add_image(normalised(imgs[:-1]), name='T1', visible=False, **img_kw)
viewer.add_image(normalised(imgs[1:]), name='T2', visible=False, **img_kw)
viewer.add_image(normalised(imgs[1:] - imgs[:-1]), name='T2 - T1', **img_kw)
viewer.add_shapes([[[wz0, wx0], [wz1, wx1]]], shape_type='rectangle', name='analysis window',
                  edge_color='lime', face_color='transparent', edge_width=0.05)
labels = viewer.add_labels(picks, name='pick', opacity=0.4, **geo)
labels.selected_label = 1
labels.brush_size = max(1, round(0.4 / min(scale)))       # ~0.4 m on every grid
labels.mode = 'paint'
viewer.layers.selection.active = labels
viewer.dims.axis_labels = ('pair', 'depth (m)', 'radial (m)')
viewer.text_overlay.visible = True


def show_pair(event=None):
    i = viewer.dims.current_step[0]
    a, b = pairs[i]
    viewer.text_overlay.text = (f'{method}: pair {a} -> {b} ({i + 1}/{len(pairs)}), '
                                f'{int((labels.data[i] > 0).sum())} px picked at slice change')


def save(_viewer=None):
    np.savez_compressed(picks_path, **{**saved, 'picks': labels.data > 0})
    print(f'saved {int((labels.data > 0).any(axis=(1, 2)).sum())}/{len(pairs)} picked pairs '
          f'to {picks_path}')


def copy_previous(_viewer=None):
    i = viewer.dims.current_step[0]
    if i > 0:
        labels.data[i] = labels.data[i - 1]
        labels.refresh()
        show_pair()


viewer.bind_key('Shift-S', save, overwrite=True)
viewer.bind_key('Shift-C', copy_previous, overwrite=True)
viewer.dims.events.current_step.connect(show_pair)
viewer.dims.set_current_step(0, 0)
show_pair()

if check:                           # exercise the slider and Shift-C without saving anything
    viewer.dims.set_current_step(0, 1)
    copy_previous()
    print('check:', viewer.text_overlay.text)
    viewer.close()
else:
    napari.run()                    # blocks until the window is closed
    save()
