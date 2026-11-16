"""
A grid cut down to a single row or column keeps the spacing of the grid or
file it was cut from: crop, copy_subset and index (the subsetting), and
from_h5, from_nc and from_geotif (reading for bounds).  mosaic.from_list
uses that spacing, so a mosaic whose bounds end exactly on the edge of a
tile, which then contributes a single row, gets both of its spacings.
"""
import numpy as np
import pytest
import pointCollection as pc

DX, DY = 2., 5.


def make_grid(x0=0., y0=0., nx=11, ny=9, nt=None):
    x = x0 + DX*np.arange(nx)
    y = y0 + DY*np.arange(ny)
    if nt is None:
        z = np.arange(ny*nx, dtype=float).reshape(ny, nx)
        return pc.grid.data().from_dict({'x': x, 'y': y, 'z': z})
    t = np.arange(nt, dtype=float)
    z = np.arange(ny*nx*nt, dtype=float).reshape(ny, nx, nt)
    return pc.grid.data().from_dict({'x': x, 'y': y, 't': t, 'z': z})


def test_full_grid_spacing_comes_from_its_coordinates():
    assert make_grid().spacing == [DX, DY]


@pytest.mark.parametrize('bounds, expected_shape', [
    ([[0., 20.], [10., 10.]], (1, 11)),   # one row
    ([[4., 4.], [0., 40.]], (9, 1)),      # one column
    ([[4., 4.], [10., 10.]], (1, 1)),     # one cell
])
def test_crop_to_a_single_row_or_column_keeps_the_spacing(bounds, expected_shape):
    cropped = make_grid().cropped(*bounds)
    assert cropped.z.shape == expected_shape
    assert cropped.spacing == [DX, DY]


def test_copy_subset_and_copies_of_subsets_keep_the_spacing():
    grid = make_grid(nt=3)
    one_row = grid.copy_subset([slice(2, 3), slice(None), slice(None)])
    assert one_row.z.shape == (1, 11, 3)
    assert one_row.spacing == [DX, DY]
    # a copy (or a crop of a copy) of the single-row grid keeps it too
    assert one_row.copy().spacing == [DX, DY]
    assert one_row.copy_meta().spacing == [DX, DY]
    assert one_row.cropped([[0., 4.], [0., 40.]]).spacing == [DX, DY]


def test_a_two_value_dimension_still_uses_the_coordinates():
    # every other column of a single row: x spacing from the subset, y inherited
    sub = make_grid().copy_subset([slice(2, 3), slice(None, None, 2)])
    assert sub.spacing == [2*DX, DY]


def test_scalar_stored_spacing_is_unchanged_without_coordinates():
    grid = pc.grid.data()
    grid.spacing = 250.
    assert grid.spacing == 250.


@pytest.mark.parametrize('fmt', ['h5', 'nc', 'tif'])
def test_reading_bounds_that_hold_a_single_row_keeps_the_file_spacing(tmp_path, fmt):
    grid = make_grid()
    filename = str(tmp_path/f'grid.{fmt}')
    if fmt == 'h5':
        grid.to_h5(filename, replace=True)
    elif fmt == 'nc':
        grid.to_nc(filename, replace=True)
    else:
        grid.to_geotif(filename, srs_epsg=3031)
    one_row = pc.grid.data().from_file(filename, bounds=[[0., 20.], [10., 10.]])
    assert one_row.z.shape == (1, 11)
    assert one_row.spacing == [DX, DY]
    one_col = pc.grid.data().from_file(filename, bounds=[[4., 4.], [0., 40.]])
    assert one_col.z.shape == (9, 1)
    assert one_col.spacing == [DX, DY]


def test_reading_with_skip_gives_the_skipped_spacing(tmp_path):
    filename = str(tmp_path/'grid.h5')
    make_grid().to_h5(filename, replace=True)
    one_row = pc.grid.data().from_h5(filename, bounds=[[0., 20.], [10., 10.]], skip=2)
    assert one_row.spacing == [2*DX, 2*DY]


@pytest.mark.parametrize('order', [1, -1])
def test_mosaic_bounds_on_a_tile_edge(tmp_path, order):
    """
    Two tiles, one above the other, sharing no rows; the mosaic's bounds end
    exactly on the bottom row of the upper tile, which then gives one row.
    With the upper tile first, the mosaic used to take its x spacing from that
    row and never find a y spacing.
    """
    lower = make_grid(y0=0., nt=2)      # y 0 ... 40
    upper = make_grid(y0=45., nt=2)     # y 45 ... 85
    upper.z += 1000.
    files = []
    for name, grid in (('lower', lower), ('upper', upper)):
        files.append(str(tmp_path/f'{name}.h5'))
        grid.to_h5(files[-1], replace=True)
    bounds = [[0., 20.], [20., 45.]]
    mos = pc.grid.mosaic().from_list(files[::-order], bounds=bounds)
    assert mos.spacing == [DX, DY]
    assert np.array_equal(mos.y, np.arange(20., 45.1, DY))
    assert np.array_equal(mos.x, lower.x)
    # the rows come from the right tiles
    assert np.array_equal(mos.z[:-1], lower.cropped(*bounds).z)
    assert np.array_equal(mos.z[-1:], upper.cropped(*bounds).z)
