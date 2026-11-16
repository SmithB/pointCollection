#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
linear_interp.py
Linear interpolation of three-dimensional grids, compiled with numba

pointCollection.grid.data.interp uses linearGridInterpolator for fields with
three dimensions (y, x, and time, in either axis order) when numba is
installed.  It gives the same results as the
scipy.interpolate.RegularGridInterpolator that interp uses otherwise
(method='linear', bounds_error=False, fill_value=nan), with the same
operations in the same order, but works point by point, without the
temporary arrays that the scipy N-dimensional path makes for each corner of
the grid cells.  The kernel is the one used by pointAdvection.advection for
time-variable velocity grids, for a single field and with NaN outside the
grid.

numba is optional: without it, interp uses RegularGridInterpolator.

UPDATE HISTORY:
    Written 11/2026
"""

import numpy as np

try:
    import numba
except ImportError:
    numba = None

if numba is not None:
    @numba.njit(cache=True)
    def _cell_index(axis, c, regular, a0, scale):
        """
        Index of the cell of an increasing axis that contains c, the same as
        ``np.searchsorted(axis, c, side='right') - 1`` clipped to [0, n-2]
        (and so as scipy's find_interval_ascending with extrapolation)

        For a regular axis, the index is estimated from the spacing (a0 is
        the first node and scale is the inverse of the spacing) and then
        corrected so that it is exact; otherwise it is found by binary search.
        """
        n = axis.size
        if regular:
            f = (c - a0)*scale
            i = int(min(max(f, 0.0), n - 1.0))
            while (i > 0) and (c < axis[i]):
                i -= 1
            while (i < n - 1) and (c >= axis[i+1]):
                i += 1
        else:
            lo, hi = 0, n
            while lo < hi:
                mid = (lo + hi) >> 1
                if c < axis[mid]:
                    hi = mid
                else:
                    lo = mid + 1
            i = lo - 1
        return min(max(i, 0), n - 2)

    @numba.njit(cache=True)
    def _linear_3d_kernel(a0, a1, a2, r0, r1, r2, values, c0, c1, c2):
        """
        Linear interpolation of values on the grid with axes a0, a1, a2 at the
        points (c0, c1, c2), with the same operations, in the same order, as
        scipy.interpolate.RegularGridInterpolator (method='linear'): the
        weights are multiplied in axis order, and the corners are summed
        starting from zero, with the last axis varying fastest.  Points that
        are outside the grid or have a NaN coordinate get NaN.

        r0, r1, r2 flag the axes that are evenly spaced.  values may have any
        memory layout.
        """
        n = c0.size
        out = np.empty(n)
        s0 = (a0.size - 1)/(a0[-1] - a0[0])
        s1 = (a1.size - 1)/(a1[-1] - a1[0])
        s2 = (a2.size - 1)/(a2[-1] - a2[0])
        for k in range(n):
            p0 = c0[k]
            p1 = c1[k]
            p2 = c2[k]
            # outside the grid, or NaN (for which the comparisons are all False)
            if not ((p0 >= a0[0]) and (p0 <= a0[-1]) and
                    (p1 >= a1[0]) and (p1 <= a1[-1]) and
                    (p2 >= a2[0]) and (p2 <= a2[-1])):
                out[k] = np.nan
                continue
            i0 = _cell_index(a0, p0, r0, a0[0], s0)
            d0 = (p0 - a0[i0])/(a0[i0+1] - a0[i0])
            i1 = _cell_index(a1, p1, r1, a1[0], s1)
            d1 = (p1 - a1[i1])/(a1[i1+1] - a1[i1])
            i2 = _cell_index(a2, p2, r2, a2[0], s2)
            d2 = (p2 - a2[i2])/(a2[i2+1] - a2[i2])
            # sum the contributions of the corners of the cell
            v = 0.0
            for e0 in range(2):
                w0 = (1 - d0) if (e0 == 0) else d0
                for e1 in range(2):
                    w01 = w0*((1 - d1) if (e1 == 0) else d1)
                    for e2 in range(2):
                        w = w01*((1 - d2) if (e2 == 0) else d2)
                        v = v + values[i0 + e0, i1 + e1, i2 + e2]*w
            out[k] = v
        return out


class linearGridInterpolator:
    """
    Linear interpolator for a field on a three-dimensional grid

    A replacement for scipy.interpolate.RegularGridInterpolator with
    method='linear', bounds_error=False and fill_value=nan, with identical
    results, that uses a numba kernel.  Make one with from_grid, which
    returns None when the kernel can't be used for the grid (numba is not
    installed, the grid is not three-dimensional, an axis has fewer than two
    points or is not strictly monotonic, or the values are complex), so the
    caller can fall back to RegularGridInterpolator.

    Attributes
    ----------
    grid: tuple
        the axes of the grid, increasing (as in RegularGridInterpolator)
    values: np.ndarray
        the values, oriented to match the increasing axes
    method: str
        'linear'
    bounds_error: bool
        False: points outside the grid get fill_value
    fill_value: float
        NaN
    """
    method = 'linear'
    bounds_error = False
    fill_value = np.nan

    def __init__(self, grid, values, regular):
        self.grid = grid
        self.values = values
        self.regular = regular

    @classmethod
    def from_grid(cls, points, values):
        """
        Set up an interpolator for values on the grid defined by points

        Parameters
        ----------
        points: tuple of np.ndarray
            the axes of the grid, in the order of the dimensions of values
        values: np.ndarray
            values on the grid

        Returns
        -------
        linearGridInterpolator or NoneType
            None if the kernel can't be used for this grid
        """
        if numba is None or len(points) != 3:
            return None
        values = np.asarray(values)
        if values.ndim != 3 or np.iscomplexobj(values):
            return None
        if values.dtype not in (np.float32, np.float64):
            # as in RegularGridInterpolator, which converts values that
            # aren't floating point to float
            values = values.astype(np.float64)
        axes = []
        for dim, axis in enumerate(points):
            axis = np.asarray(axis, dtype=np.float64).ravel()
            if axis.size < 2 or axis.size != values.shape[dim] or not np.all(np.isfinite(axis)):
                return None
            steps = np.diff(axis)
            if np.all(steps < 0):
                # flip descending axes, as RegularGridInterpolator does
                axis = axis[::-1]
                values = np.flip(values, axis=dim)
            elif not np.all(steps > 0):
                return None
            axes.append(np.ascontiguousarray(axis))
        # (nearly) evenly spaced axes, where the kernel finds the cells from
        # the spacing instead of by searching
        regular = []
        for axis in axes:
            steps = np.diff(axis)
            regular.append(bool(np.all(np.abs(steps - steps[0]) <= 1e-3*steps[0])))
        return cls(tuple(axes), values, tuple(regular))

    def __call__(self, xi):
        """
        Interpolate at coordinates

        Parameters
        ----------
        xi: tuple of array-like, or array-like
            coordinates in the order of the grid's axes: a tuple of three
            arrays (broadcast against each other), or an array whose last
            dimension has length 3

        Returns
        -------
        np.ndarray
            interpolated values, with the shape of the broadcast coordinates
            (or of xi without its last dimension)
        """
        if isinstance(xi, tuple):
            coords = np.broadcast_arrays(*[np.asarray(c, dtype=np.float64) for c in xi])
        else:
            xi = np.asarray(xi, dtype=np.float64)
            if xi.shape[-1] != 3:
                raise ValueError(f'the requested sample points xi have dimension {xi.shape[-1]}, '
                                 'but this linearGridInterpolator has dimension 3')
            coords = [xi[..., dim] for dim in range(3)]
        if len(coords) != 3:
            raise ValueError(f'the requested sample points xi have dimension {len(coords)}, '
                             'but this linearGridInterpolator has dimension 3')
        shape = coords[0].shape
        coords = [np.ascontiguousarray(c).ravel() for c in coords]
        result = _linear_3d_kernel(*self.grid, *self.regular, self.values, *coords)
        return result.reshape(shape)
