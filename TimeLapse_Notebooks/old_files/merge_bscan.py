"""
Merge N individual gprMax B-scan output files into one HDF5 file.

Usage: python merge_bscan.py <base_path> <n_traces> [--remove-files]

Each per-trace file is:  <base_path>{i}.out   (i = 1 .. n_traces)
Merged output file:      <base_path>_merged.out

The merged file contains rxs/rx1 … rxs/rxN (one per trace), copying all
field components (Ey, Ex, Ez, Hx, Hy, Hz) from each source file.
"""

import sys
import h5py
import numpy as np
from pathlib import Path


def merge_bscan(base_path, n_traces, remove_files=False):
    base = Path(base_path)
    merged_path = Path(str(base) + '_merged.out')

    # --- Read header info from file 1 ---
    first = Path(f'{base}1.out')
    with h5py.File(first, 'r') as f0:
        dt         = float(f0.attrs['dt'])
        iterations = int(f0.attrs['Iterations'])
        title      = str(f0.attrs.get('Title', ''))
        dx_dy_dz   = f0.attrs['dx_dy_dz'][:]
        nx_ny_nz   = f0.attrs['nx_ny_nz'][:]
        gprmax_ver = str(f0.attrs.get('gprMax', ''))
        components = list(f0['rxs']['rx1'].keys())

    print(f'Merging {n_traces} traces -> {merged_path.name}')
    print(f'  dt = {dt:.4e} s,  iterations = {iterations},  '
          f'components = {components}')

    with h5py.File(merged_path, 'w') as fout:
        # Root attributes
        fout.attrs['dt']         = dt
        fout.attrs['Iterations'] = iterations
        fout.attrs['Title']      = title
        fout.attrs['dx_dy_dz']   = dx_dy_dz
        fout.attrs['nx_ny_nz']   = nx_ny_nz
        fout.attrs['gprMax']     = gprmax_ver
        fout.attrs['nrx']        = n_traces
        fout.attrs['nsrc']       = n_traces

        rxs_group = fout.create_group('rxs')

        for i in range(1, n_traces + 1):
            fpath = Path(f'{base}{i}.out')
            if not fpath.exists():
                raise FileNotFoundError(f'Missing: {fpath}')

            with h5py.File(fpath, 'r') as fi:
                rx_grp = rxs_group.create_group(f'rx{i}')
                for comp in components:
                    data = fi['rxs']['rx1'][comp][:]
                    rx_grp.create_dataset(comp, data=data)

            if i % 50 == 0 or i == n_traces:
                print(f'  {i}/{n_traces} traces merged …')

    print(f'Done -> {merged_path}  ({merged_path.stat().st_size/1e6:.1f} MB)')

    if remove_files:
        for i in range(1, n_traces + 1):
            Path(f'{base}{i}.out').unlink(missing_ok=True)
        print(f'Removed {n_traces} individual trace files.')

    return merged_path


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Merge gprMax B-scan outputs')
    parser.add_argument('base_path',  help='Base path without index or .out')
    parser.add_argument('n_traces',   type=int, help='Number of traces')
    parser.add_argument('--remove-files', action='store_true')
    args = parser.parse_args()
    merge_bscan(args.base_path, args.n_traces, args.remove_files)
