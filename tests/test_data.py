"""
Tests for pointCollection.data
"""
import os
import numpy as np
import pytest
import pointCollection as pc

TEST_H5 = os.path.join(os.path.dirname(__file__), '..', 'test_data',
                       'ATL06_20190205041106_05910210_209_01.h5')


# ---------------------------------------------------------------------------
# __init__ / assign
# ---------------------------------------------------------------------------

def test_fields_dict_assigns_data():
    x = np.array([1.0, 2.0, 3.0])
    y = np.array([4.0, 5.0, 6.0])
    D = pc.data(fields={'x': x, 'y': y})
    assert isinstance(D.fields, list)
    assert set(D.fields) == {'x', 'y'}
    np.testing.assert_array_equal(D.x, x)
    np.testing.assert_array_equal(D.y, y)


def test_fields_dict_is_list_not_dict():
    D = pc.data(fields={'x': np.zeros(3), 'y': np.zeros(3)})
    _ = D.fields[0]          # would fail if fields were a dict
    assert len(D.fields) == 2


def test_fields_list_unchanged():
    fields = ['x', 'y', 'z']
    D = pc.data(fields=fields)
    assert D.fields == fields


def test_fields_none_default():
    D = pc.data()
    assert D.fields == []


def test_assign_sets_shape_from_array():
    D = pc.data()
    D.assign({'x': np.arange(5.0), 'y': np.arange(5.0)})
    assert D.shape == (5,)
    assert D.size == 5


def test_assign_scalar_with_array_broadcasts():
    D = pc.data()
    D.assign({'x': np.arange(4.0), 'scale': 2.0})
    assert D.shape == (4,)
    np.testing.assert_array_equal(D.scale, np.full(4, 2.0))


def test_assign_all_scalars_no_shape_raises():
    D = pc.data()
    with pytest.raises(ValueError, match="cannot determine shape"):
        D.assign({'scale': 1.0})


def test_assign_kwargs_does_not_mutate_caller_dict():
    # |= on a dict is in-place; using | creates a new dict so the caller's dict is unchanged
    D = pc.data(fields={'x': np.arange(4.0)})
    original = {'y': np.arange(4.0)}
    D.assign(original, z=np.ones(4))
    assert 'z' not in original


# ---------------------------------------------------------------------------
# __update_size_and_shape__
# ---------------------------------------------------------------------------

def test_update_size_and_shape_explicit():
    D = pc.data()
    D.__update_size_and_shape__(shape=(10, 3))
    assert D.shape == (10, 3)
    assert D.size == 30


# ---------------------------------------------------------------------------
# from_h5 with multiple groups
# ---------------------------------------------------------------------------

def test_from_h5_multiple_groups_reads_all():
    groups = ['gt1l/land_ice_segments/dem', 'gt1l/land_ice_segments/fit_statistics']
    D = pc.data().from_h5(TEST_H5, group=groups)
    assert 'dem_h' in D.fields
    assert 'h_mean' in D.fields
    assert hasattr(D, 'dem_h')
    assert hasattr(D, 'h_mean')
    assert D.dem_h.shape == D.h_mean.shape


# ---------------------------------------------------------------------------
# copy_subset by_row behaviour
# ---------------------------------------------------------------------------

def test_copy_subset_by_row_true_slices_rows():
    D = pc.data(fields={'x': np.arange(6.0).reshape(3, 2)}, columns=2)
    sub = D.copy_subset(np.array([0, 2]), by_row=True)
    np.testing.assert_array_equal(sub.x, np.array([[0., 1.], [4., 5.]]))


def test_copy_subset_columns_auto_sets_by_row():
    # columns >= 1 forces row-wise slicing regardless of the by_row default
    D = pc.data(fields={'x': np.arange(6.0).reshape(3, 2)}, columns=2)
    sub = D.copy_subset(np.array([0, 2]))
    np.testing.assert_array_equal(sub.x, np.array([[0., 1.], [4., 5.]]))


def test_copy_subset_1d_ravels():
    # 1-D fields use element-wise indexing (no by_row auto-override)
    D = pc.data(fields={'x': np.arange(6.0)})
    sub = D.copy_subset(np.array([0, 2, 4]))
    np.testing.assert_array_equal(sub.x, np.array([0., 2., 4.]))


# ---------------------------------------------------------------------------
# get_xy / get_latlon roundtrip
# ---------------------------------------------------------------------------

def test_blockmedian_raises_on_multicolumn_field():
    D = pc.data(fields={
        'x': np.array([0., 0., 1., 1.]),
        'y': np.array([0., 0., 0., 0.]),
        'z': np.ones((4, 2)),          # multi-column — must not be used as median field
    })
    with pytest.raises(ValueError, match="must be 1-D"):
        D.blockmedian(1.0, field='z')


def test_blockmedian_applies_index_to_multicolumn_fields():
    # Four points in two spatial bins; each bin contains two identical rows.
    # After blockmedian the averaged rows should equal the original rows.
    x = np.array([0.1, 0.2, 1.1, 1.2])
    y = np.zeros(4)
    z = np.array([1.0, 2.0, 3.0, 4.0])
    v = np.array([[1., 10.], [2., 20.], [3., 30.], [4., 40.]])  # multi-column
    D = pc.data(fields={'x': x, 'y': y, 'z': z, 'v': v}, columns=2)
    D.blockmedian(1.0, field='z')
    assert D.v.ndim == 2
    assert D.v.shape[1] == 2


def test_to_h5_attrs_written_to_group(tmp_path):
    import h5py
    D = pc.data(fields={'x': np.arange(5.0), 'y': np.arange(5.0)})
    D.attrs = {'source': 'test', 'version': 1}
    outfile = str(tmp_path / 'out.h5')
    D.to_h5(outfile, group='/data')
    with h5py.File(outfile, 'r') as f:
        assert 'source' in f['/data'].attrs
        assert 'version' in f['/data'].attrs
        # attrs must NOT appear on a dataset
        assert 'source' not in f['/data/x'].attrs


def test_coordinates_roundtrip_h5(tmp_path):
    import h5py
    x = np.arange(5.0)
    y = np.arange(5.0) * 2
    z = np.ones(5)
    D = pc.data(fields={'x': x, 'y': y, 'z': z})
    outfile = str(tmp_path / 'coords.h5')
    D.to_h5(outfile, group='/data')
    # verify coordinates attr written on non-coordinate field, not on coordinate fields
    with h5py.File(outfile, 'r') as f:
        assert f['/data/z'].attrs.get('coordinates') == 'x y'
        assert 'coordinates' not in f['/data/x'].attrs
        assert 'coordinates' not in f['/data/y'].attrs
    # verify coordinates recovered on read
    D2 = pc.data().from_h5(outfile, group='/data')
    assert D2.coordinates == ['x', 'y']


def test_coordinates_roundtrip_h5_custom(tmp_path):
    lon = np.linspace(-10.0, 10.0, 6)
    lat = np.linspace(60.0, 70.0, 6)
    z = np.ones(6)
    D = pc.data(fields={'lon': lon, 'lat': lat, 'z': z}, coordinates=['lon', 'lat'])
    outfile = str(tmp_path / 'custom_coords.h5')
    D.to_h5(outfile, group='/data')
    D2 = pc.data().from_h5(outfile, group='/data')
    assert D2.coordinates == ['lon', 'lat']


def test_get_xy_get_latlon_roundtrip():
    lat = np.array([-75.0, -80.0])
    lon = np.array([0.0, 45.0])
    D = pc.data(fields={'latitude': lat, 'longitude': lon})
    D.get_xy(EPSG=3031)
    assert 'x' in D.fields and 'y' in D.fields
    D.get_latlon(EPSG=3031)
    np.testing.assert_allclose(D.latitude, lat, atol=1e-6)
    np.testing.assert_allclose(D.longitude, lon, atol=1e-6)


# ---------------------------------------------------------------------------
# configurable coordinates
# ---------------------------------------------------------------------------

def test_custom_coordinates_bounds():
    cx = np.array([1.0, 2.0, 3.0])
    cy = np.array([10.0, 20.0, 30.0])
    D = pc.data(fields={'cx': cx, 'cy': cy}, coordinates=['cx', 'cy'])
    xr, yr = D.bounds()
    np.testing.assert_array_equal(xr, np.array([1.0, 3.0]))
    np.testing.assert_array_equal(yr, np.array([10.0, 30.0]))


def test_custom_coordinates_coords():
    cx = np.array([1.0, 2.0])
    cy = np.array([3.0, 4.0])
    D = pc.data(fields={'cx': cx, 'cy': cy}, coordinates=['cx', 'cy'])
    result = D.coords()
    # coords() returns fields in self.coordinates order
    np.testing.assert_array_equal(result[0], cx)
    np.testing.assert_array_equal(result[1], cy)


def test_coords_custom_order():
    cx = np.array([1.0, 2.0])
    cy = np.array([3.0, 4.0])
    D = pc.data(fields={'cx': cx, 'cy': cy}, coordinates=['cx', 'cy'])
    result = D.coords(order=['cy', 'cx'])
    np.testing.assert_array_equal(result[0], cy)
    np.testing.assert_array_equal(result[1], cx)


def test_single_coordinate_coords():
    d = np.array([0.0, 1.0, 2.0])
    D = pc.data(fields={'dist': d}, coordinates=['dist'])
    result = D.coords()
    assert len(result) == 1
    np.testing.assert_array_equal(result[0], d)


def test_single_coordinate_bounds():
    d = np.array([1.0, 3.0, 5.0])
    D = pc.data(fields={'dist': d}, coordinates=['dist'])
    result = D.bounds()
    assert len(result) == 1
    np.testing.assert_array_equal(result[0], np.array([1.0, 5.0]))


def test_single_coordinate_crop():
    d = np.array([0.0, 1.0, 2.0, 3.0])
    z = np.array([10.0, 20.0, 30.0, 40.0])
    D = pc.data(fields={'dist': d, 'z': z}, coordinates=['dist'])
    D.crop([[0.5, 2.5]])
    np.testing.assert_array_equal(D.dist, np.array([1.0, 2.0]))
    np.testing.assert_array_equal(D.z, np.array([20.0, 30.0]))


def test_single_coordinate_blockmedian():
    d = np.array([0.1, 0.2, 1.1, 1.2])
    z = np.array([1.0, 3.0, 5.0, 7.0])
    D = pc.data(fields={'dist': d, 'z': z}, coordinates=['dist'])
    D.blockmedian(1.0, field='z')
    assert D.size == 2


def test_custom_coordinates_crop():
    cx = np.array([0.0, 1.0, 2.0, 3.0])
    cy = np.array([0.0, 0.0, 0.0, 0.0])
    D = pc.data(fields={'cx': cx, 'cy': cy}, coordinates=['cx', 'cy'])
    D.crop([[0.5, 2.5], [-1.0, 1.0]])
    np.testing.assert_array_equal(D.cx, np.array([1.0, 2.0]))


def test_coordinates_copied_in_copy_subset():
    cx = np.array([0.0, 1.0, 2.0])
    cy = np.array([5.0, 6.0, 7.0])
    D = pc.data(fields={'cx': cx, 'cy': cy}, coordinates=['cx', 'cy'])
    sub = D.copy_subset(np.array([0, 2]))
    assert sub.coordinates == ['cx', 'cy']
    xr, yr = sub.bounds()
    np.testing.assert_array_equal(xr, np.array([0.0, 2.0]))


def test_default_coordinates_are_xy():
    D = pc.data()
    assert D.coordinates == ['x', 'y']
    assert D._x_coord == 'x'
    assert D._y_coord == 'y'
