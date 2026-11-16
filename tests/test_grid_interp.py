"""
grid.data.interp of a three-dimensional field uses the numba kernel in
grid.linear_interp when numba is installed.  Its results must be identical to
those of scipy's RegularGridInterpolator (which interp used before, and still
uses without numba): inside the grid, on nodes and cell edges, at the grid's
edges, outside it (NaN), at NaN coordinates, next to NaN values, for float32
and integer values, with irregular and descending axes, for either time axis,
with scalar times, and on gridded output.
"""
import numpy as np
import pytest
import scipy.interpolate
import pointCollection as pc
from pointCollection.grid.linear_interp import linearGridInterpolator

numba = pytest.importorskip('numba')


def axes_and_values(rng, irregular=False, dtype=np.float64, nan_frac=0.02):
    y = 1000. + 250.*np.arange(13)
    x = -5000. + 250.*np.arange(17)
    t = 2019. + 0.25*np.arange(9)
    if irregular:
        x = np.cumsum(rng.uniform(100., 400., 17)) - 5000.
        t = np.sort(rng.uniform(2019., 2021., 9))
    values = rng.normal(size=(y.size, x.size, t.size)).astype(dtype)
    if nan_frac:
        values[rng.random(values.shape) < nan_frac] = np.nan
    return [y, x, t], values


def sample_points(rng, axes, n=20000):
    """random points, inside and out, plus nodes, cell edges, grid edges and NaNs"""
    pts = [rng.uniform(a.min() - 0.1*np.ptp(a), a.max() + 0.1*np.ptp(a), n) for a in axes]
    # nodes, and points on the cell edges of one dimension
    for dim, a in enumerate(axes):
        idx = rng.integers(0, a.size, n//10)
        pts[dim][:n//10] = a[idx]
    # the grid's edges (first and last nodes) and NaN coordinates
    pts[0][-30:-20] = axes[0].max()
    pts[1][-20:-10] = axes[1].min()
    pts[2][-10:-5] = axes[2].max()
    pts[1][-5:] = np.nan
    return pts


def rgi(axes, values, pts):
    return scipy.interpolate.RegularGridInterpolator(axes, values, bounds_error=False)(tuple(pts))


@pytest.mark.parametrize('irregular', [False, True])
@pytest.mark.parametrize('dtype', [np.float64, np.float32, np.int32])
def test_kernel_matches_regular_grid_interpolator(irregular, dtype):
    rng = np.random.default_rng(1)
    axes, values = axes_and_values(rng, irregular=irregular, dtype=np.float64,
                                   nan_frac=0.02 if dtype != np.int32 else 0)
    values = values.astype(dtype) if dtype != np.int32 else (100*values).astype(dtype)
    pts = sample_points(rng, axes)
    interpolator = linearGridInterpolator.from_grid(tuple(axes), values)
    assert interpolator is not None
    result = interpolator(tuple(pts))
    expected = rgi(axes, values, pts)
    assert result.dtype == expected.dtype
    np.testing.assert_array_equal(result, expected)
    # and the same for an (N, 3) array of points
    np.testing.assert_array_equal(interpolator(np.c_[tuple(pts)]), expected)


def test_kernel_matches_with_descending_axes():
    rng = np.random.default_rng(2)
    axes, values = axes_and_values(rng)
    axes = [axes[0][::-1], axes[1], axes[2][::-1]]
    pts = sample_points(rng, axes)
    interpolator = linearGridInterpolator.from_grid(tuple(axes), values)
    np.testing.assert_array_equal(interpolator(tuple(pts)), rgi(axes, values, pts))


@pytest.mark.parametrize('axes_change', ['one_point', 'not_monotonic', 'two_d'])
def test_grids_the_kernel_cannot_take_give_none(axes_change):
    rng = np.random.default_rng(3)
    axes, values = axes_and_values(rng)
    if axes_change == 'one_point':
        axes[2] = axes[2][:1]
        values = values[:, :, :1]
    elif axes_change == 'not_monotonic':
        axes[1] = axes[1].copy()
        axes[1][3] = axes[1][5]
    else:
        axes = axes[:2]
        values = values[:, :, 0]
    assert linearGridInterpolator.from_grid(tuple(axes), values) is None


@pytest.mark.parametrize('t_axis', [0, 2])
def test_grid_interp_matches_regular_grid_interpolator(t_axis):
    rng = np.random.default_rng(4)
    (y, x, t), values = axes_and_values(rng)
    if t_axis == 0:
        values = np.moveaxis(values, 2, 0).copy()
        grid = pc.grid.data(t_axis=0).from_dict({'x': x, 'y': y, 't': t, 'z': values})
        axes = [t, y, x]
    else:
        grid = pc.grid.data().from_dict({'x': x, 'y': y, 't': t, 'z': values})
        axes = [y, x, t]
    py, px, pt = sample_points(rng, [y, x, t])
    result = grid.interp(px, py, pt)
    assert isinstance(grid.interpolator['z'], linearGridInterpolator)
    pts = [pt, py, px] if t_axis == 0 else [py, px, pt]
    np.testing.assert_array_equal(result, rgi(axes, values, pts))
    # one time for all of the points
    pts = [np.full_like(px, t[3] + 0.1), py, px] if t_axis == 0 else [py, px, np.full_like(px, t[3] + 0.1)]
    np.testing.assert_array_equal(grid.interp(px, py, t[3] + 0.1), rgi(axes, values, pts))


def test_gridded_interp_matches_regular_grid_interpolator():
    rng = np.random.default_rng(5)
    (y, x, t), values = axes_and_values(rng)
    grid = pc.grid.data().from_dict({'x': x, 'y': y, 't': t, 'z': values})
    xg = np.linspace(x[0] - 100, x[-1] + 100, 31)
    yg = np.linspace(y[0], y[-1], 17)
    tg = np.array([2019.1, 2020.3])
    result = grid.interp(xg, yg, tg, gridded=True)
    YY, XX, TT = np.meshgrid(yg, xg, tg, indexing='ij')
    np.testing.assert_array_equal(result, rgi([y, x, t], values, [YY, XX, TT]))
    assert result.shape == (yg.size, xg.size, tg.size)


def test_two_d_fields_and_bands_still_use_scipy():
    rng = np.random.default_rng(6)
    (y, x, t), values = axes_and_values(rng)
    grid = pc.grid.data().from_dict({'x': x, 'y': y, 't': t, 'z': values, 'z2': values[:, :, 0].copy()})
    grid.interp(x[:3], y[:3], field='z2')
    grid.interp(x[:3], y[:3], band=1, replace=True)
    assert isinstance(grid.interpolator['z2'], scipy.interpolate.RegularGridInterpolator)
    assert isinstance(grid.interpolator['z'], scipy.interpolate.RegularGridInterpolator)


def test_without_numba_interp_uses_scipy(monkeypatch):
    import pointCollection.grid.linear_interp as linear_interp
    monkeypatch.setattr(linear_interp, 'numba', None)
    rng = np.random.default_rng(7)
    (y, x, t), values = axes_and_values(rng)
    grid = pc.grid.data().from_dict({'x': x, 'y': y, 't': t, 'z': values})
    py, px, pt = sample_points(rng, [y, x, t], n=100)
    result = grid.interp(px, py, pt)
    assert isinstance(grid.interpolator['z'], scipy.interpolate.RegularGridInterpolator)
    np.testing.assert_array_equal(result, rgi([y, x, t], values, [py, px, pt]))
