"""
Tests for writing boolean fields from pointCollection.grid.data

Booleans are written as bytes (uint8) with a nodata value of 255.
"""
import numpy as np
import pytest
import h5py
import pointCollection as pc


def _make_grid(three_d=False):
    x = np.arange(4.)
    y = np.arange(3.)
    if three_d:
        t = np.arange(2.)
        mask = np.zeros((3, 4, 2), dtype=bool)
        mask[:, 2:, 1] = True
        return pc.grid.data().from_dict({'x': x, 'y': y, 't': t,
                                         'z': np.ones(mask.shape), 'mask': mask})
    mask = np.zeros((3, 4), dtype=bool)
    mask[1:, 1:] = True
    return pc.grid.data().from_dict({'x': x, 'y': y, 'z': np.arange(12.).reshape(3, 4),
                                     'mask': mask})


@pytest.mark.parametrize('nocompression', [False, True])
@pytest.mark.parametrize('three_d', [False, True])
def test_to_h5_writes_bool_as_bytes(tmp_path, nocompression, three_d):
    grid = _make_grid(three_d)
    out_file = str(tmp_path / 'bool.h5')
    grid.to_h5(out_file, group='g', replace=True, nocompression=nocompression)
    with h5py.File(out_file, 'r') as h5f:
        assert 'mask' in h5f['g']
        assert h5f['g/mask'].dtype == np.uint8
        assert h5f['g/mask'].fillvalue == 255
        np.testing.assert_array_equal(np.array(h5f['g/mask']), grid.mask.astype(np.uint8))
        # other fields are unchanged
        assert h5f['g/z'].dtype == np.float64
    # reading back gives the mask values, with zeros intact
    temp = pc.grid.data().from_h5(out_file, group='g')
    np.testing.assert_array_equal(temp.mask, grid.mask.astype(float))
    # the in-memory field is still boolean
    assert grid.mask.dtype == bool


def test_to_h5_overwrites_bool_field(tmp_path):
    grid = _make_grid()
    out_file = str(tmp_path / 'bool.h5')
    grid.to_h5(out_file, group='g', replace=True)
    grid.mask = ~grid.mask
    grid.to_h5(out_file, group='g')
    with h5py.File(out_file, 'r') as h5f:
        assert h5f['g/mask'].dtype == np.uint8
        np.testing.assert_array_equal(np.array(h5f['g/mask']), grid.mask.astype(np.uint8))


@pytest.mark.parametrize('three_d', [False, True])
def test_to_nc_writes_bool_as_bytes(tmp_path, three_d):
    netCDF4 = pytest.importorskip('netCDF4')
    grid = _make_grid(three_d)
    out_file = str(tmp_path / 'bool.nc')
    grid.to_nc(out_file, replace=True)
    with netCDF4.Dataset(out_file, 'r') as ncf:
        assert 'mask' in ncf.variables
        var = ncf.variables['mask']
        assert var.dtype == np.uint8
        assert var._FillValue == 255
        var.set_auto_mask(False)
        np.testing.assert_array_equal(np.array(var[...]), grid.mask.astype(np.uint8))
    temp = pc.grid.data().from_nc(out_file)
    np.testing.assert_array_equal(temp.mask, grid.mask.astype(float))
