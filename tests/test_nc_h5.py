"""
Tests for reading netCDF4 files through h5py (grid.nc_h5), the path a remote
from_nc() takes so that a windowed read costs the chunks the window touches
rather than the whole granule.

No network access is required.  engine='h5py' exercises the same reader
against a local file, and a stand-in filesystem maps s3:// URIs onto local
files for the remote branch, so what is under test is the plumbing and the
netCDF4 API the adapter presents, not S3 itself.
"""
import io
import numpy as np
import pytest
import pointCollection as pc

netCDF4 = pytest.importorskip('netCDF4')


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

def write_packed_nc(path, ydesc=True, group=None, nt=4, ny=7, nx=9):
    """
    A file to_nc() cannot write: int16 data with scale_factor/add_offset and
    a _FillValue, a grid_mapping variable, and a dimension with no coordinate
    variable.  netCDF4 applies the scaling on the way out even with masking
    off, and h5py applies none of it, so this is the fixture that catches a
    reader that quietly returns different numbers.
    """
    with netCDF4.Dataset(path, 'w') as ds:
        root = ds.createGroup(group) if group else ds
        root.createDimension('x', nx)
        root.createDimension('y', ny)
        root.createDimension('time', nt)
        root.createDimension('nv', 2)     # netCDF dimension, not a variable
        root.createVariable('x', 'f8', ('x',))[:] = np.arange(nx) * 100.
        yv = np.arange(ny) * 100.
        root.createVariable('y', 'f8', ('y',))[:] = yv[::-1] if ydesc else yv
        root.createVariable('time', 'f8', ('time',))[:] = 2019. + np.arange(nt)
        z = root.createVariable('z', 'i2', ('time', 'y', 'x'),
                                fill_value=np.int16(-9999), zlib=True)
        z.scale_factor = np.float32(0.01)
        z.add_offset = np.float32(5.)
        z.grid_mapping = 'crs'
        # write storage values directly: assigning through netCDF4's packing
        # would turn a fill value into (fill - add_offset)/scale_factor
        z.set_auto_maskandscale(False)
        raw = np.arange(nt * ny * nx, dtype='i2').reshape(nt, ny, nx)
        raw[0, 1, 1] = -9999
        raw[2, 3, 4] = -9999
        z[:] = raw
        crs = root.createVariable('crs', 'S1', ())
        crs.spatial_epsg = '3413'
        crs.false_easting = 0.
    return path


@pytest.fixture
def nc_2d(tmp_path):
    path = str(tmp_path / 'grid_2d.nc')
    pc.grid.data().from_dict({
        'x': np.arange(9) * 100.,
        'y': np.arange(7) * 100.,
        'z': np.arange(63, dtype=float).reshape((7, 9))}).to_nc(path, replace=True)
    return path


@pytest.fixture
def nc_3d(tmp_path):
    path = str(tmp_path / 'grid_3d.nc')
    pc.grid.data().from_dict({
        'x': np.arange(9) * 100.,
        'y': np.arange(7) * 100.,
        't': 2019. + np.arange(4),
        'z': np.arange(7 * 9 * 4, dtype=float).reshape((7, 9, 4))}).to_nc(path, replace=True)
    return path


class FakeS3FS:
    """
    Stand-in for s3fs.S3FileSystem that maps s3:// URIs onto local files and
    records what was opened, and with what block size.
    """
    def __init__(self, mapping, wrapper=None):
        self.mapping = dict(mapping)
        self.opened = []
        self.block_sizes = []
        self.wrapper = wrapper

    def open(self, path, mode='rb', block_size=None):
        self.opened.append((path, mode))
        self.block_sizes.append(block_size)
        fd = open(self.mapping[path], mode)
        return self.wrapper(fd) if self.wrapper else fd

    def exists(self, path):
        return path in self.mapping


class CountingFile(io.RawIOBase):
    """file object that records how many bytes were actually read from it"""
    def __init__(self, fd):
        self.fd = fd
        self.bytes_read = 0

    def read(self, size=-1):
        data = self.fd.read(size)
        self.bytes_read += len(data)
        return data

    def readinto(self, buffer):
        n = self.fd.readinto(buffer)
        self.bytes_read += n or 0
        return n

    def seek(self, offset, whence=0):
        return self.fd.seek(offset, whence)

    def tell(self):
        return self.fd.tell()

    def seekable(self):
        return True

    def readable(self):
        return True

    def close(self):
        self.fd.close()


def assert_same_read(path, **kwargs):
    """the h5py reader and the netCDF4 reader must return the same grid"""
    expected = pc.grid.data().from_nc(path, engine='netcdf4', **kwargs)
    got = pc.grid.data().from_nc(path, engine='h5py', **kwargs)
    assert sorted(got.fields) == sorted(expected.fields)
    for field in expected.fields:
        a, b = getattr(expected, field), getattr(got, field)
        assert a.shape == b.shape, field
        assert a.dtype == b.dtype, field
        assert np.array_equal(a, b, equal_nan=True), field
    for coord in ('x', 'y', 't', 'time'):
        a, b = getattr(expected, coord, None), getattr(got, coord, None)
        assert (a is None) == (b is None), coord
        if a is not None:
            assert np.array_equal(np.asarray(a), np.asarray(b)), coord
    assert getattr(got, 'projection', None) == getattr(expected, 'projection', None)
    assert got.extent == expected.extent
    return got


# ---------------------------------------------------------------------------
# the h5py reader returns what the netCDF4 reader returns
# ---------------------------------------------------------------------------

READS = {
    'whole': {},
    'bounds': {'bounds': [[100, 500], [200, 600]]},
    'skip': {'skip': 2},
    'bounds_and_skip': {'bounds': [[0, 800], [0, 600]], 'skip': 2},
    'meta_only': {'meta_only': True},
    'field': {'field': 'z'},
    'fill_value': {'fill_value': -999.},
}


@pytest.mark.parametrize('case', sorted(READS))
def test_h5py_matches_netcdf4_2d(nc_2d, case):
    assert_same_read(nc_2d, **READS[case])


@pytest.mark.parametrize('case', sorted(READS))
def test_h5py_matches_netcdf4_3d(nc_3d, case):
    assert_same_read(nc_3d, **READS[case])


def test_bounds_outside_the_raster_behaves_the_same(nc_2d):
    """bounds that miss the raster raise out of select_slices(), either way"""
    for engine in ('netcdf4', 'h5py'):
        with pytest.raises(IndexError):
            pc.grid.data().from_nc(nc_2d, engine=engine,
                                   bounds=[[1e5, 2e5], [1e5, 2e5]])


@pytest.mark.parametrize('kwargs', [
    {'bands': [1, 3]},                    # scattered bands: not a plain slice
    {'bands': [0, 1, 2]},                 # contiguous bands
    {'bands': [2]},
    {'t_range': [2020, 2021]},
    {'t_axis': 0},
    {'t_axis': 0, 'bounds': [[100, 500], [200, 600]]},
    {'bands': [0, 3], 'skip': 2},
])
def test_h5py_matches_netcdf4_bands(nc_3d, kwargs):
    assert_same_read(nc_3d, **kwargs)


@pytest.mark.parametrize('ydesc', [True, False])
@pytest.mark.parametrize('group', [None, 'grp'])
@pytest.mark.parametrize('kwargs', [
    {}, {'bounds': [[100, 500], [200, 600]]}, {'skip': 2}, {'bands': [0, 2, 3]}])
def test_h5py_matches_netcdf4_packed(tmp_path, ydesc, group, kwargs):
    """
    scale_factor/add_offset, _FillValue, grid_mapping, a descending y axis
    (which makes select_slices() emit the negative-step slices h5py rejects)
    and a dimension with no coordinate variable.
    """
    path = write_packed_nc(str(tmp_path / 'packed.nc'), ydesc=ydesc, group=group)
    if group:
        kwargs = dict(kwargs, group=group)
    out = assert_same_read(path, **kwargs)
    # scaling really is being applied, and 'nv'/'crs' are not fields
    assert out.z.dtype == np.float32
    assert sorted(out.fields) == ['z']
    assert out.projection['spatial_epsg'] == '3413'


def test_packed_fill_value_is_converted(tmp_path):
    """
    A packed variable's data come back scaled while its _FillValue is stored
    raw, so the two only match once the fill value is put on the same scale.
    """
    path = write_packed_nc(str(tmp_path / 'packed.nc'), ydesc=False)
    for engine in ('netcdf4', 'h5py'):
        out = pc.grid.data().from_nc(path, engine=engine)
        # two fill values were written, at (time, y, x) = (0,1,1) and (2,3,4)
        assert np.isnan(out.z[1, 1, 0])
        assert np.isnan(out.z[3, 4, 2])
        assert np.count_nonzero(np.isnan(out.z)) == 2
        # and nothing else was swept up: -9999 unscaled is not a value here
        assert np.isfinite(out.z[0, 0, 0])


def test_packed_fill_matches_netcdf4_masking(tmp_path):
    """
    netCDF4's own masking (which from_nc turns off) knows which values are
    fills, so its mask is an independent check on which values from_nc
    converted.
    """
    path = write_packed_nc(str(tmp_path / 'packed.nc'), ydesc=False)
    with netCDF4.Dataset(path) as ds:      # auto-mask on: netCDF4 finds fills
        masked = np.transpose(ds.variables['z'][:], [1, 2, 0])
    for engine in ('netcdf4', 'h5py'):
        out = pc.grid.data().from_nc(path, engine=engine)
        assert np.array_equal(np.isnan(out.z), np.ma.getmaskarray(masked))


def test_packed_fill_value_respects_fill_value_argument(tmp_path):
    path = write_packed_nc(str(tmp_path / 'packed.nc'), ydesc=False)
    for engine in ('netcdf4', 'h5py'):
        out = pc.grid.data().from_nc(path, engine=engine, fill_value=-999.)
        assert out.z[1, 1, 0] == -999.
        assert np.count_nonzero(out.z == -999.) == 2


def test_unpacked_fill_value_is_unchanged(tmp_path):
    """a variable with no scaling must behave exactly as it did before"""
    path = str(tmp_path / 'plain.nc')
    with netCDF4.Dataset(path, 'w') as ds:
        ds.createDimension('x', 4)
        ds.createDimension('y', 3)
        ds.createVariable('x', 'f8', ('x',))[:] = np.arange(4) * 100.
        ds.createVariable('y', 'f8', ('y',))[:] = np.arange(3) * 100.
        z = ds.createVariable('z', 'f4', ('y', 'x'), fill_value=np.float32(-9999.))
        z.set_auto_maskandscale(False)
        values = np.arange(12, dtype='f4').reshape(3, 4)
        values[1, 2] = -9999.
        z[:] = values
    for engine in ('netcdf4', 'h5py'):
        out = pc.grid.data().from_nc(path, engine=engine)
        assert np.isnan(out.z[1, 2])
        assert np.count_nonzero(np.isnan(out.z)) == 1


def test_h5py_applies_scale_factor(tmp_path):
    """guard against the silent-wrong-numbers failure directly"""
    path = write_packed_nc(str(tmp_path / 'packed.nc'), ydesc=False)
    out = pc.grid.data().from_nc(path, engine='h5py')
    with netCDF4.Dataset(path) as ds:
        ds.set_auto_mask(False)
        expected = np.transpose(ds.variables['z'][:], [1, 2, 0])
    # everything but the fill values, which from_nc turns into fill_value
    valid = ~np.isnan(out.z)
    assert np.count_nonzero(valid) == out.z.size - 2
    assert np.array_equal(out.z[valid], expected[valid])
    # unscaled data would be off by a factor of 100 and an offset
    with pc.io_utils.open_h5(path) as h5f:
        raw = np.transpose(h5f['z'][:], [1, 2, 0])
    assert not np.array_equal(out.z[valid], raw[valid])


# ---------------------------------------------------------------------------
# remote reads go through h5py, and only fetch what they need
# ---------------------------------------------------------------------------

def test_remote_read_matches_local(nc_3d):
    fs = FakeS3FS({'s3://bucket/grid.nc': nc_3d})
    remote = pc.grid.data().from_nc('s3://bucket/grid.nc', fs=fs,
                                    bounds=[[100, 500], [200, 600]])
    local = pc.grid.data().from_nc(nc_3d, bounds=[[100, 500], [200, 600]])
    assert np.array_equal(remote.z, local.z)
    assert fs.opened == [('s3://bucket/grid.nc', 'rb')]


def test_remote_read_with_netcdf4_engine(nc_3d):
    """engine='netcdf4' keeps the old behaviour: read the file whole"""
    fs = FakeS3FS({'s3://bucket/grid.nc': nc_3d})
    remote = pc.grid.data().from_nc('s3://bucket/grid.nc', fs=fs, engine='netcdf4')
    assert np.array_equal(remote.z, pc.grid.data().from_nc(nc_3d).z)
    assert fs.block_sizes == [None]      # a whole-file read wants big blocks


def test_remote_read_passes_block_size(nc_3d):
    fs = FakeS3FS({'s3://bucket/grid.nc': nc_3d})
    pc.grid.data().from_nc('s3://bucket/grid.nc', fs=fs)
    # FakeS3FS.open() accepts block_size, so the default must have reached it
    assert fs.block_sizes == [pc.io_utils.DEFAULT_REMOTE_BLOCK_SIZE]

    fs = FakeS3FS({'s3://bucket/grid.nc': nc_3d})
    pc.grid.data().from_nc('s3://bucket/grid.nc', fs=fs, block_size=64 * 1024)
    assert fs.block_sizes == [64 * 1024]


def test_remote_read_survives_fs_without_block_size(nc_3d):
    """an fs whose open() takes no block_size must still work"""
    class PlainFS(FakeS3FS):
        def open(self, path, mode='rb'):
            self.opened.append((path, mode))
            return open(self.mapping[path], mode)

    fs = PlainFS({'s3://bucket/grid.nc': nc_3d})
    out = pc.grid.data().from_nc('s3://bucket/grid.nc', fs=fs)
    assert out.z is not None
    assert fs.opened == [('s3://bucket/grid.nc', 'rb')]


def test_windowed_remote_read_touches_a_fraction_of_the_file(tmp_path):
    """
    The point of the change: a bounded read must not pull the whole file.
    """
    import os
    path = str(tmp_path / 'chunked.nc')
    n, chunk = 400, 50
    rng = np.random.default_rng(0)
    with netCDF4.Dataset(path, 'w') as ds:
        ds.createDimension('x', n)
        ds.createDimension('y', n)
        ds.createVariable('x', 'f8', ('x',))[:] = np.arange(n) * 100.
        ds.createVariable('y', 'f8', ('y',))[:] = np.arange(n) * 100.
        z = ds.createVariable('z', 'f4', ('y', 'x'), zlib=True,
                              chunksizes=(chunk, chunk))
        z[:] = rng.random((n, n), dtype='f4')
    size = os.path.getsize(path)

    counters = []

    def wrapper(fd):
        counters.append(CountingFile(fd))
        return counters[-1]

    fs = FakeS3FS({'s3://bucket/chunked.nc': path}, wrapper=wrapper)
    window = pc.grid.data().from_nc('s3://bucket/chunked.nc', fs=fs,
                                    bounds=[[10e3, 15e3], [10e3, 15e3]])
    whole = pc.grid.data().from_nc(path)
    read = counters[0].bytes_read
    assert window.z.shape == (51, 51)
    assert np.array_equal(window.z, whole.z[100:151, 100:151])
    # the window spans 2x2 chunks out of 64; allow generous slack for metadata
    assert read < size / 4, f'read {read} of {size} bytes'


# ---------------------------------------------------------------------------
# files h5py cannot open
# ---------------------------------------------------------------------------

def make_netcdf3(path):
    with netCDF4.Dataset(path, 'w', format='NETCDF3_CLASSIC') as ds:
        ds.createDimension('x', 9)
        ds.createDimension('y', 7)
        ds.createVariable('x', 'f8', ('x',))[:] = np.arange(9) * 100.
        ds.createVariable('y', 'f8', ('y',))[:] = np.arange(7) * 100.
        ds.createVariable('z', 'f8', ('y', 'x'))[:] = np.arange(63.).reshape(7, 9)
    return path


def test_netcdf3_remote_falls_back_to_netcdf4(tmp_path):
    path = make_netcdf3(str(tmp_path / 'classic.nc'))
    fs = FakeS3FS({'s3://bucket/classic.nc': path})
    remote = pc.grid.data().from_nc('s3://bucket/classic.nc', fs=fs)
    assert np.array_equal(remote.z, pc.grid.data().from_nc(path).z)
    # h5py is tried first, then the file is read whole for netCDF4
    assert fs.opened == [('s3://bucket/classic.nc', 'rb')] * 2
    assert fs.block_sizes == [pc.io_utils.DEFAULT_REMOTE_BLOCK_SIZE, None]


def test_netcdf3_explicit_h5py_engine_raises(tmp_path):
    path = make_netcdf3(str(tmp_path / 'classic.nc'))
    with pytest.raises(OSError):
        pc.grid.data().from_nc(path, engine='h5py')


def test_unknown_engine_rejected(nc_2d):
    with pytest.raises(ValueError):
        pc.grid.data().from_nc(nc_2d, engine='zarr')


# ---------------------------------------------------------------------------
# the netCDF4 API the adapter presents
# ---------------------------------------------------------------------------

def test_adapter_presents_netcdf4_api(tmp_path):
    from pointCollection.grid import nc_h5
    path = write_packed_nc(str(tmp_path / 'packed.nc'), ydesc=False, group='grp')
    with pc.grid.data().nc_open(path, engine='h5py') as ds:
        grp = ds.groups['grp']
        assert set(grp.variables.keys()) == {'x', 'y', 'time', 'z', 'crs'}
        assert 'nv' not in grp.variables      # a dimension, not a variable
        z = grp.variables['z']
        assert z.dimensions == ('time', 'y', 'x')
        assert z.shape == (4, 7, 9)
        assert hasattr(z, '_FillValue')
        assert z.getncattr('_FillValue') == -9999
        assert z.getncattr('grid_mapping') == 'crs'      # decoded, not bytes
        assert not hasattr(z, 'units')
        # HDF5 bookkeeping attributes stay hidden, as in netCDF4
        assert 'DIMENSION_LIST' not in z.ncattrs()
        assert '_Netcdf4Dimid' not in z.ncattrs()
        assert set(z.ncattrs()) == {'_FillValue', 'scale_factor', 'add_offset',
                                    'grid_mapping'}
        assert grp.variables['x'].dimensions == ('x',)
        assert isinstance(ds, nc_h5.H5Dataset)
        ds.set_auto_mask(False)      # no-op, as netCDF4's is here


def test_adapter_slicing(tmp_path):
    path = write_packed_nc(str(tmp_path / 'packed.nc'), ydesc=False)
    with pc.grid.data().nc_open(path, engine='h5py') as ds:
        z = ds.variables['z']
        full = z[:]
        # negative steps, which h5py itself rejects
        assert np.array_equal(z[:, ::-1, :], full[:, ::-1, :])
        assert np.array_equal(z[::-1, ::-2, ::-1], full[::-1, ::-2, ::-1])
        assert np.array_equal(z[:, 5:1:-1, :], full[:, 5:1:-1, :])
        # index lists, contiguous and scattered
        assert np.array_equal(z[[0, 1, 2], :, :], full[[0, 1, 2], :, :])
        assert np.array_equal(z[[0, 3], :, :], full[[0, 3], :, :])
        assert np.array_equal(z[np.array([3, 0]), :, :], full[[3, 0], :, :])
        # integer index drops its axis; short tuples and Ellipsis
        assert np.array_equal(z[2], full[2])
        assert np.array_equal(z[2, ::-1], full[2, ::-1])
        assert np.array_equal(z[..., ::-1], full[..., ::-1])
        assert np.array_equal(z[0:2], full[0:2])
        # an empty selection
        assert z[:, 0:0, :].shape == (4, 0, 9)
        # a scalar variable, which netCDF4 lets you index with [:]
        crs = ds.variables['crs']
        assert crs.shape == ()
        assert np.shape(crs[:]) == () and np.shape(crs[...]) == ()


def test_adapter_closes_the_byte_source(nc_2d):
    fs = FakeS3FS({'s3://bucket/grid.nc': nc_2d})
    ds = pc.grid.data().nc_open('s3://bucket/grid.nc', fs=fs)
    source = ds._source
    ds.close()
    assert source.closed, 'the remote file object outlived the dataset'
