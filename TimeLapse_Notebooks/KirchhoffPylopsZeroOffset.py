"""
KirchhoffPylopsZeroOffset.py
Zero-offset Kirchhoff demigration / migration operator (PyLops-compatible).

The standard pylops.waveeqprocessing.Kirchhoff has dimsd = (ns, nr, nt).
For zero-offset B-scan with ns = nr = n_traces this becomes an n^2 x nt data
cube — completely infeasible for typical survey sizes.

This module restricts to the diagonal case (source i co-located with receiver i)
so the data shape collapses from (ns, nr, nt) to (n, nt).

Data shape  (dimsd) : (n, nt)   — one trace per co-located Tx/Rx pair
Model shape (dims)  : (nx, nz)  — reflectivity image on x × z grid

Typical usage
-------------
K = Kirchhoff(z, x, t, srcs, recs, vel, wav, wavcenter, mode='analytic')
image = (K.H @ data.flatten()).reshape(len(x), len(z)).T   # → (nz, nx)
"""

import warnings
import numpy as np
from pylops import LinearOperator
from pylops.signalprocessing import Convolve1D
from pylops.waveeqprocessing.kirchhoff import Kirchhoff as _PyLopsKirchhoff


class Kirchhoff(LinearOperator):
    """
    Zero-offset Kirchhoff demigration / migration operator.

    Restricts the full PyLops Kirchhoff to zero-offset geometry:
    source i is co-located with receiver i for all i = 0..n-1.
    This reduces the data dimension from (ns, nr, nt) → (n, nt).

    Parameters
    ----------
    z           : (nz,)   depth axis [m]
    x           : (nx,)   image spatial axis (= trace x-positions) [m]
    t           : (nt,)   time axis [ns]
    srcs        : (2, n)  source positions [m];  row 0 = x, row 1 = z
    recs        : (2, n)  receiver positions [m]; must satisfy srcs==recs for ZO
    vel         : float   propagation velocity [m/ns]
    wav         : (nw,)   wavelet array
    wavcenter   : int     index of wavelet centre sample (zero-delay)
    mode        : str     'analytic' (constant vel) or 'eikonal'
    dynamic     : bool    ignored; always False for zero-offset
    aperture    : ignored
    angleaperture : ignored
    dtype       : str
    name        : str

    Attributes
    ----------
    dims  : (nx, nz)   model shape  — use .reshape(nx, nz).T to get (nz, nx) image
    dimsd : (n, nt)    data shape
    """

    def __init__(
        self,
        z, x, t, srcs, recs, vel, wav, wavcenter,
        mode='analytic', dynamic=False, aperture=None, angleaperture=None,
        dtype='float64', name='KZO',
    ):
        ns, nr = srcs.shape[1], recs.shape[1]
        if ns != nr:
            raise ValueError(
                f'Zero-offset requires ns == nr, got ns={ns}, nr={nr}'
            )
        n  = ns
        nz, nx = len(z), len(x)
        ni = nx * nz
        nt = len(t)
        dt = float(t[1] - t[0])

        # Compute one-way traveltime tables via the standard PyLops method.
        # trav_srcs : (ni, ns),   trav_recs : (ni, nr)
        # For zero-offset: total TWTT for pair i = trav_srcs[:,i] + trav_recs[:,i]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            trav_srcs, trav_recs, _, _, _, _ = _PyLopsKirchhoff._traveltime_table(
                z, x, srcs, recs, vel, mode=mode
            )

        # Two-way traveltime: (ni, n)
        trav_zo = trav_srcs + trav_recs

        itrav = (trav_zo / dt).astype('int32')
        travd = trav_zo / dt - itrav
        valid = (itrav >= 0) & (itrav < nt - 1)

        # Apply angle aperture if provided
        if angleaperture is not None:
            X, Z = np.meshgrid(x, z, indexing="ij")
            dx = srcs[0] - X.ravel()[:, np.newaxis]
            dz = srcs[1] - Z.ravel()[:, np.newaxis]
            # Use arctan2(dx, dz) where dz is usually positive (surface to depth)
            angle_deg = np.rad2deg(np.abs(np.arctan2(np.abs(dx), np.abs(dz))))
            
            if isinstance(angleaperture, (int, float)):
                valid_angle = angle_deg <= float(angleaperture)
            elif isinstance(angleaperture, (list, tuple, np.ndarray)) and len(angleaperture) >= 2:
                valid_angle = angle_deg <= float(angleaperture[1]) # Hard threshold on max angle
            else:
                valid_angle = True
                
            valid = valid & valid_angle

        # Store only what _matvec / _rmatvec need
        self.n     = n
        self.ni    = ni
        self.nt    = nt
        self.dt    = dt
        self.itrav = itrav    # (ni, n) int32
        self.travd = travd    # (ni, n) float64
        self.valid = valid    # (ni, n) bool
        self._tr   = np.arange(n, dtype=np.intp)   # trace indices 0..n-1

        # Pre-compute valid-entry index arrays for fast scatter in _matvec
        ii_v, i_v = np.where(valid)
        self._ii_v  = ii_v                     # valid image-point indices
        self._i_v   = i_v                      # corresponding trace indices
        self._i0_v  = itrav[ii_v, i_v]        # integer traveltime sample
        self._td_v  = travd[ii_v, i_v]        # fractional part

        # Wavelet convolution operator on (n, nt) data along axis=1
        self.cop = Convolve1D((n, nt), h=wav, offset=wavcenter, axis=1, dtype=dtype)

        super().__init__(
            dtype=np.dtype(dtype),
            dims=(nx, nz),
            dimsd=(n, nt),
            name=name,
        )

    # ------------------------------------------------------------------
    def _matvec(self, x):
        """Forward (demigration): reflectivity (ni,) → B-scan (n*nt,)."""
        m = x.ravel()
        d = np.zeros((self.n, self.nt), dtype=self.dtype)
        # Scatter-interpolate reflectivity into data at two-way traveltimes
        m_v = m[self._ii_v]
        np.add.at(d, (self._i_v, self._i0_v),     m_v * (1.0 - self._td_v))
        np.add.at(d, (self._i_v, self._i0_v + 1), m_v * self._td_v)
        return self.cop._matvec(d.ravel())

    # ------------------------------------------------------------------
    def _rmatvec(self, x):
        """Adjoint (migration): B-scan (n*nt,) → reflectivity (ni,)."""
        # Step 1: wavelet cross-correlation along time axis (each trace independently)
        d = self.cop._rmatvec(x.ravel()).reshape(self.n, self.nt)

        # Step 2: vectorised gather-and-sum over all image points and traces
        #
        # For each image point ii and trace i:
        #   samp[ii, i] = (1-td)*d[i, i0] + td*d[i, i0+1]
        # then  m[ii] = sum_i  samp[ii, i]  (sum valid traces)
        #
        # Implemented with broadcast advanced indexing:
        #   d[n_idx, i0]  where n_idx=(1,n) broadcasts against i0=(ni,n)
        i0  = np.clip(self.itrav, 0, self.nt - 2)          # (ni, n)
        td  = self.travd                                     # (ni, n)
        # n_idx[np.newaxis, :] has shape (1, n) → broadcasts to (ni, n)
        samp = ((1.0 - td) * d[self._tr[np.newaxis, :], i0    ]
                +       td  * d[self._tr[np.newaxis, :], i0 + 1])   # (ni, n)
        return (samp * self.valid).sum(axis=1)               # (ni,)
