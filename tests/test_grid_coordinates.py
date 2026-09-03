#!/usr/bin/env python3
"""
Pytest suite for the coordinate generalization in pointCollection.grid.data.

The refactor stores coordinate names in self.coordinates (default ['y','x','time'])
and exposes three helper properties: _row_coord, _col_coord, _band_coord.  All
spatial methods use these names to look up coordinate arrays rather than
hard-coding 'x', 'y', 't'.  Tests are grouped by the layer of functionality
they exercise, from construction through I/O.
"""

import numpy as np
import pytest
import pointCollection as pc


# ---------------------------------------------------------------------------
# Fixtures — small grids reused across multiple tests
# ---------------------------------------------------------------------------

@pytest.fixture
def grid_2d():
    """2-D grid with default coordinates ['y', 'x']."""
    d = pc.grid.data(coordinates=['y', 'x'])
    d.y = np.arange(4.0)
    d.x = np.arange(5.0)
    d.assign({'z': np.arange(20.0).reshape(4, 5)})
    d.__update_size_and_shape__()
    d.__update_extent__()
    return d


@pytest.fixture
def grid_3d():
    """3-D grid using the default coordinates ['y', 'x', 'time'] (t_axis=2)."""
    d = pc.grid.data()   # defaults to ['y', 'x', 'time'], t_axis=2
    d.y    = np.arange(4.0)
    d.x    = np.arange(5.0)
    d.time = np.array([1.0, 2.0, 3.0])
    d.assign({'z': np.arange(60.0).reshape(4, 5, 3)})
    d.__update_size_and_shape__()
    d.__update_extent__()
    return d


@pytest.fixture
def grid_3d_t0():
    """3-D grid with t_axis=0 and coordinates ['time', 'y', 'x']."""
    d = pc.grid.data(t_axis=0)
    d.y    = np.arange(4.0)
    d.x    = np.arange(5.0)
    d.time = np.array([1.0, 2.0, 3.0])
    d.assign({'z': np.arange(60.0).reshape(3, 4, 5)})
    d.__update_size_and_shape__()
    d.__update_extent__()
    return d


@pytest.fixture
def grid_3d_custom():
    """3-D grid with custom coordinates ['longitude', 'latitude', 'elevation']."""
    d = pc.grid.data(coordinates=['longitude', 'latitude', 'elevation'])
    d.longitude = np.arange(4.0)
    d.latitude  = np.arange(5.0)
    d.elevation = np.array([100.0, 200.0, 300.0])
    d.assign({'z': np.arange(60.0).reshape(4, 5, 3)})
    d.__update_size_and_shape__()
    d.__update_extent__()
    return d


# ---------------------------------------------------------------------------
# Construction: self.coordinates and self.t_axis
# ---------------------------------------------------------------------------

def test_default_coordinates_are_y_x_time():
    # __init__ with no arguments produces the canonical default coordinate list.
    d = pc.grid.data()
    assert d.coordinates == ['y', 'x', 'time']


def test_default_t_axis_is_2():
    # __init__ with no arguments leaves t_axis at 2.
    assert pc.grid.data().t_axis == 2


def test_t_axis_0_sets_time_first_coordinates():
    # Passing t_axis=0 reorders coordinates so time is the leading axis.
    d = pc.grid.data(t_axis=0)
    assert d.coordinates == ['time', 'y', 'x']
    assert d.t_axis == 0


def test_custom_coordinates_stored_verbatim():
    # An explicit coordinates list is kept exactly as supplied.
    coords = ['longitude', 'latitude', 'elevation']
    assert pc.grid.data(coordinates=coords).coordinates == coords


def test_2d_coordinates_have_no_band():
    # A two-element coordinates list signals a 2-D (no-band) grid.
    d = pc.grid.data(coordinates=['y', 'x'])
    assert d.coordinates == ['y', 'x']
    assert d._band_coord is None


def test_t_axis_is_settable_after_construction():
    # t_axis can be updated via its property setter.
    d = pc.grid.data()
    d.t_axis = 0
    assert d.t_axis == 0
    assert d._t_axis == 0


# ---------------------------------------------------------------------------
# Helper properties: _row_coord, _col_coord, _band_coord
# ---------------------------------------------------------------------------

def test_coord_helpers_default_t_axis_2():
    # Default coordinates ['y','x','time'] map row→'y', col→'x', band→'time'.
    d = pc.grid.data()
    assert d._row_coord  == 'y'
    assert d._col_coord  == 'x'
    assert d._band_coord == 'time'


def test_coord_helpers_t_axis_0():
    # With t_axis=0 and ['time','y','x'], spatial helpers still return 'y'/'x'.
    d = pc.grid.data(t_axis=0)
    assert d._row_coord  == 'y'
    assert d._col_coord  == 'x'
    assert d._band_coord == 'time'


def test_coord_helpers_custom_3d():
    # Custom 3-D coordinates map positions 0→row, 1→col, 2→band.
    d = pc.grid.data(coordinates=['longitude', 'latitude', 'elevation'])
    assert d._row_coord  == 'longitude'
    assert d._col_coord  == 'latitude'
    assert d._band_coord == 'elevation'


def test_coord_helpers_2d_band_is_none():
    # A 2-D coordinate list gives _band_coord of None.
    assert pc.grid.data(coordinates=['y', 'x'])._band_coord is None


# ---------------------------------------------------------------------------
# spacing property
# ---------------------------------------------------------------------------

def test_spacing_default_coordinates():
    # spacing reads self.x and self.y when coordinates are the default.
    d = pc.grid.data()
    d.x = np.array([0.0, 2.0, 4.0])
    d.y = np.array([0.0, 3.0, 6.0])
    assert d.spacing == [2.0, 3.0]


def test_spacing_custom_coordinates():
    # spacing reads the named column and row arrays for custom coordinates.
    d = pc.grid.data(coordinates=['longitude', 'latitude', 'elevation'])
    d.longitude = np.array([0.0, 1.5, 3.0])   # row coord → y-like
    d.latitude  = np.array([0.0, 2.5, 5.0])   # col coord → x-like
    assert d.spacing == [2.5, 1.5]


# ---------------------------------------------------------------------------
# __update_extent__
# ---------------------------------------------------------------------------

def test_extent_default_coordinates():
    # __update_extent__ builds [xmin, xmax, ymin, ymax] from self.x and self.y.
    d = pc.grid.data()
    d.x = np.array([1.0, 2.0, 3.0])
    d.y = np.array([10.0, 20.0])
    d.__update_extent__()
    assert d.extent == [1.0, 3.0, 10.0, 20.0]


def test_extent_custom_coordinates():
    # __update_extent__ uses the named col/row arrays for custom coordinates.
    d = pc.grid.data(coordinates=['longitude', 'latitude', 'elevation'])
    d.longitude = np.array([10.0, 20.0])        # row coord
    d.latitude  = np.array([1.0, 2.0, 3.0])     # col coord
    d.__update_extent__()
    assert d.extent == [1.0, 3.0, 10.0, 20.0]


# ---------------------------------------------------------------------------
# __update_size_and_shape__
# ---------------------------------------------------------------------------

def test_shape_default_2d(grid_2d):
    # 2-D default grid reports shape [ny, nx].
    assert grid_2d.shape == [4, 5]


def test_shape_default_3d(grid_3d):
    # 3-D default grid (t_axis=2) reports shape [ny, nx, nt].
    assert grid_3d.shape == [4, 5, 3]


def test_shape_t_axis_0(grid_3d_t0):
    # t_axis=0 grid reports shape [nt, ny, nx].
    assert grid_3d_t0.shape == [3, 4, 5]


def test_shape_custom_3d(grid_3d_custom):
    # Custom coordinate grid reports shape following the coordinates list order.
    assert grid_3d_custom.shape == [4, 5, 3]


# ---------------------------------------------------------------------------
# from_dict and assign: coordinate names must not appear in self.fields
# ---------------------------------------------------------------------------

def test_from_dict_default_coords_excluded_from_fields():
    # from_dict does not add 'x', 'y', or 'time' to self.fields.
    d = pc.grid.data()
    d.from_dict({'x': np.arange(5.0), 'y': np.arange(4.0),
                 'time': np.array([1.0, 2.0, 3.0]),
                 'z': np.ones((4, 5, 3))})
    for coord in ('x', 'y', 'time'):
        assert coord not in d.fields
    assert 'z' in d.fields


def test_from_dict_custom_coords_excluded_from_fields():
    # from_dict does not add custom coordinate names to self.fields.
    d = pc.grid.data(coordinates=['longitude', 'latitude', 'elevation'])
    d.from_dict({'longitude': np.arange(4.0),
                 'latitude':  np.arange(5.0),
                 'elevation': np.array([100.0, 200.0, 300.0]),
                 'z': np.ones((4, 5, 3))})
    for coord in ('longitude', 'latitude', 'elevation'):
        assert coord not in d.fields
    assert 'z' in d.fields


def test_assign_custom_coord_excluded_from_fields():
    # assign does not add a custom coordinate name to self.fields.
    d = pc.grid.data(coordinates=['longitude', 'latitude', 'elevation'])
    d.assign({'longitude': np.arange(4.0), 'z': np.ones((4, 5))})
    assert 'longitude' not in d.fields
    assert 'z' in d.fields


def test_assign_non_coord_field_added_to_fields():
    # assign adds an ordinary data field to self.fields.
    d = pc.grid.data()
    d.assign({'w': np.ones((4, 5))})
    assert 'w' in d.fields


# ---------------------------------------------------------------------------
# copy and copy_meta
# ---------------------------------------------------------------------------

def test_copy_preserves_coordinates(grid_3d_custom):
    # copy() transfers the coordinates list to the new object.
    assert grid_3d_custom.copy().coordinates == grid_3d_custom.coordinates


def test_copy_preserves_t_axis(grid_3d_t0):
    # copy() transfers _t_axis to the new object.
    assert grid_3d_t0.copy()._t_axis == grid_3d_t0._t_axis


def test_copy_preserves_custom_coord_arrays(grid_3d_custom):
    # copy() reproduces all named coordinate arrays for custom coordinates.
    c = grid_3d_custom.copy()
    np.testing.assert_array_equal(c.longitude, grid_3d_custom.longitude)
    np.testing.assert_array_equal(c.latitude,  grid_3d_custom.latitude)
    np.testing.assert_array_equal(c.elevation, grid_3d_custom.elevation)


def test_copy_preserves_default_coord_arrays(grid_3d):
    # copy() reproduces x, y, and time arrays for default coordinates.
    c = grid_3d.copy()
    np.testing.assert_array_equal(c.x,    grid_3d.x)
    np.testing.assert_array_equal(c.y,    grid_3d.y)
    np.testing.assert_array_equal(c.time, grid_3d.time)


def test_copy_preserves_data_field(grid_3d_custom):
    # copy() reproduces the data arrays alongside the coordinate arrays.
    np.testing.assert_array_equal(grid_3d_custom.copy().z, grid_3d_custom.z)


def test_copy_meta_carries_coordinates_not_fields(grid_3d_custom):
    # copy_meta() copies coordinates and their arrays but leaves fields empty.
    m = grid_3d_custom.copy_meta()
    assert m.coordinates == grid_3d_custom.coordinates
    np.testing.assert_array_equal(m.longitude, grid_3d_custom.longitude)
    assert m.fields == []


def test_copy_meta_default_carries_coord_arrays(grid_3d):
    # copy_meta() copies x, y, and time arrays for default coordinates.
    m = grid_3d.copy_meta()
    np.testing.assert_array_equal(m.x,    grid_3d.x)
    np.testing.assert_array_equal(m.y,    grid_3d.y)
    np.testing.assert_array_equal(m.time, grid_3d.time)


# ---------------------------------------------------------------------------
# select_slices
# ---------------------------------------------------------------------------

def test_select_slices_keys_match_default_coords():
    # select_slices keys the returned dict by 'x', 'y', 'time' for default coords.
    d = pc.grid.data()
    slices, _ = d.select_slices(None, np.arange(10.0), np.arange(8.0), None)
    assert {'x', 'y', 'time'} <= slices.keys()


def test_select_slices_keys_match_custom_coords():
    # select_slices uses _col_coord, _row_coord, _band_coord as dict keys.
    d = pc.grid.data(coordinates=['longitude', 'latitude', 'elevation'])
    slices, _ = d.select_slices(None, np.arange(10.0), np.arange(8.0), None)
    assert {'latitude', 'longitude', 'elevation'} <= slices.keys()


def test_select_slices_bounds_restrict_x_and_y():
    # select_slices returns slices that keep only values within the given bounds.
    d = pc.grid.data()
    x = np.arange(10.0)
    y = np.arange(8.0)
    slices, _ = d.select_slices([[2.0, 5.0], [1.0, 4.0]], x, y, None)
    assert x[slices['x']][0]  >= 2.0 and x[slices['x']][-1] <= 5.0
    assert y[slices['y']][0]  >= 1.0 and y[slices['y']][-1] <= 4.0


# ---------------------------------------------------------------------------
# index
# ---------------------------------------------------------------------------

def test_index_default_slices_spatial_coords(grid_3d):
    # index() subsets x and y by the supplied row/col indices.
    g = grid_3d.copy().index(np.array([0, 1]), np.array([0, 1, 2]))
    np.testing.assert_array_equal(g.y, grid_3d.y[[0, 1]])
    np.testing.assert_array_equal(g.x, grid_3d.x[[0, 1, 2]])


def test_index_custom_slices_named_spatial_coords(grid_3d_custom):
    # index() subsets the named row/col coordinate arrays for custom coordinates.
    g = grid_3d_custom.copy().index(np.array([0, 1, 2]), np.array([0, 1]))
    np.testing.assert_array_equal(g.longitude, grid_3d_custom.longitude[[0, 1, 2]])
    np.testing.assert_array_equal(g.latitude,  grid_3d_custom.latitude[[0, 1]])


def test_index_band_default_slices_time(grid_3d):
    # index() subsets the time array when band_ind is supplied.
    g = grid_3d.copy().index(slice(None), slice(None), band_ind=np.array([0, 2]))
    np.testing.assert_array_equal(g.time, grid_3d.time[[0, 2]])


def test_index_band_custom_slices_named_band_coord(grid_3d_custom):
    # index() subsets the named band coordinate array for custom coordinates.
    g = grid_3d_custom.copy().index(slice(None), slice(None),
                                    band_ind=np.array([0, 2]))
    np.testing.assert_array_equal(g.elevation, grid_3d_custom.elevation[[0, 2]])


def test_index_updates_shape(grid_3d):
    # index() recomputes self.shape after slicing all three dimensions.
    g = grid_3d.copy().index(np.array([0, 1]), np.array([0, 1, 2]),
                              band_ind=np.array([0, 2]))
    assert g.shape == [2, 3, 2]


# ---------------------------------------------------------------------------
# interp
# ---------------------------------------------------------------------------

def test_interp_2d_default_interior_point():
    # interp() evaluates a 2-D field at an interior point using self.x / self.y.
    d = pc.grid.data(coordinates=['y', 'x'])
    d.y = np.array([0.0, 1.0, 2.0])
    d.x = np.array([0.0, 1.0, 2.0])
    yg, xg = np.meshgrid(d.y, d.x, indexing='ij')
    d.assign({'z': yg + xg})       # z[i,j] = y[i] + x[j]
    result = d.interp(np.array([0.5]), np.array([0.5]), field='z')
    np.testing.assert_allclose(result, [1.0], atol=1e-10)


def test_interp_2d_custom_interior_point():
    # interp() uses the named row/col arrays to build the interpolant for custom coords.
    d = pc.grid.data(coordinates=['longitude', 'latitude'])
    d.longitude = np.array([0.0, 1.0, 2.0])   # row coord
    d.latitude  = np.array([0.0, 1.0, 2.0])   # col coord
    rg, cg = np.meshgrid(d.longitude, d.latitude, indexing='ij')
    d.assign({'z': rg + cg})
    result = d.interp(np.array([0.5]), np.array([0.5]), field='z')
    np.testing.assert_allclose(result, [1.0], atol=1e-10)


def test_interp_3d_default_on_grid_node(grid_3d):
    # interp() evaluates a 3-D field at a grid node and returns the exact value.
    # z[0,0,0] = 0.0 for the fixture (arange reshaped)
    result = grid_3d.interp(np.array([0.0]), np.array([0.0]),
                            t=np.array([1.0]), field='z')
    np.testing.assert_allclose(result, [0.0], atol=1e-10)


# ---------------------------------------------------------------------------
# as_points
# ---------------------------------------------------------------------------

def test_as_points_default_2d_count(grid_2d):
    # as_points() returns ny*nx points for a 2-D grid.
    pts = grid_2d.as_points(keep_all=True)
    assert len(pts.x) == 4 * 5


def test_as_points_default_2d_has_x_and_y(grid_2d):
    # as_points() output has 'x' and 'y' attributes for default coordinates.
    pts = grid_2d.as_points(keep_all=True)
    assert hasattr(pts, 'x') and hasattr(pts, 'y')


def test_as_points_default_3d_has_time(grid_3d):
    # as_points() output carries a time (or t) attribute for a 3-D default grid.
    pts = grid_3d.as_points(keep_all=True)
    assert hasattr(pts, 'time') or hasattr(pts, 't')


def test_as_points_custom_3d_uses_coord_names(grid_3d_custom):
    # as_points() uses custom coordinate names as output attribute names.
    pts = grid_3d_custom.as_points(keep_all=True)
    assert hasattr(pts, 'longitude')
    assert hasattr(pts, 'latitude')
    assert hasattr(pts, 'elevation')


def test_as_points_custom_3d_count(grid_3d_custom):
    # as_points() returns nlon*nlat*nelev points for a custom 3-D grid.
    pts = grid_3d_custom.as_points(keep_all=True)
    assert len(pts.longitude) == 4 * 5 * 3


# ---------------------------------------------------------------------------
# HDF5 round-trip
# ---------------------------------------------------------------------------

def test_h5_roundtrip_default(tmp_path, grid_3d):
    # to_h5 / from_h5 preserves x, y, time, and z for default coordinates.
    fname = str(tmp_path / 'default.h5')
    grid_3d.to_h5(fname, replace=True)
    d2 = pc.grid.data().from_h5(fname)
    np.testing.assert_array_equal(d2.x,    grid_3d.x)
    np.testing.assert_array_equal(d2.y,    grid_3d.y)
    np.testing.assert_array_equal(d2.time, grid_3d.time)
    np.testing.assert_allclose(d2.z, grid_3d.z)


def test_h5_roundtrip_t_axis_0(tmp_path, grid_3d_t0):
    # to_h5 / from_h5 preserves data layout for t_axis=0 grids.
    fname = str(tmp_path / 't0.h5')
    grid_3d_t0.to_h5(fname, replace=True)
    d2 = pc.grid.data(t_axis=0).from_h5(fname, t_axis=0)
    np.testing.assert_array_equal(d2.x,    grid_3d_t0.x)
    np.testing.assert_array_equal(d2.y,    grid_3d_t0.y)
    np.testing.assert_array_equal(d2.time, grid_3d_t0.time)
    np.testing.assert_allclose(d2.z, grid_3d_t0.z)


def test_h5_roundtrip_custom(tmp_path, grid_3d_custom):
    # to_h5 / from_h5 preserves named coordinate arrays and data for custom coords.
    fname = str(tmp_path / 'custom.h5')
    grid_3d_custom.to_h5(fname, replace=True)
    d2 = pc.grid.data(coordinates=['longitude', 'latitude', 'elevation']).from_h5(
        fname,
        xname='latitude', yname='longitude', timename='elevation')
    np.testing.assert_array_equal(d2.longitude, grid_3d_custom.longitude)
    np.testing.assert_array_equal(d2.latitude,  grid_3d_custom.latitude)
    np.testing.assert_array_equal(d2.elevation, grid_3d_custom.elevation)
    np.testing.assert_allclose(d2.z, grid_3d_custom.z)


# ---------------------------------------------------------------------------
# netCDF round-trip
# ---------------------------------------------------------------------------

def test_nc_roundtrip_default(tmp_path, grid_3d):
    # to_nc / from_nc preserves x, y, time, and z for default coordinates.
    fname = str(tmp_path / 'default.nc')
    grid_3d.to_nc(fname, replace=True)
    d2 = pc.grid.data().from_nc(fname)
    np.testing.assert_array_equal(d2.x,    grid_3d.x)
    np.testing.assert_array_equal(d2.y,    grid_3d.y)
    np.testing.assert_array_equal(d2.time, grid_3d.time)
    np.testing.assert_allclose(d2.z, grid_3d.z)


def test_nc_roundtrip_2d(tmp_path, grid_2d):
    # to_nc / from_nc preserves x, y, and z for a 2-D grid with no time axis.
    fname = str(tmp_path / '2d.nc')
    grid_2d.to_nc(fname, replace=True)
    d2 = pc.grid.data(coordinates=['y', 'x']).from_nc(fname)
    np.testing.assert_array_equal(d2.x, grid_2d.x)
    np.testing.assert_array_equal(d2.y, grid_2d.y)
    np.testing.assert_allclose(d2.z, grid_2d.z)


def test_nc_roundtrip_custom(tmp_path, grid_3d_custom):
    # to_nc / from_nc preserves named coordinate arrays and data for custom coords.
    fname = str(tmp_path / 'custom.nc')
    grid_3d_custom.to_nc(fname, replace=True)
    d2 = pc.grid.data(coordinates=['longitude', 'latitude', 'elevation']).from_nc(
        fname,
        xname='latitude', yname='longitude', timename='elevation')
    np.testing.assert_array_equal(d2.longitude, grid_3d_custom.longitude)
    np.testing.assert_array_equal(d2.latitude,  grid_3d_custom.latitude)
    np.testing.assert_array_equal(d2.elevation, grid_3d_custom.elevation)
    np.testing.assert_allclose(d2.z, grid_3d_custom.z)
