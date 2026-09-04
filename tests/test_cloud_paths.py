"""
Tests for reading ancillary grids and tiling schemas from a remote (s3://)
location rather than local disk.

No network access is required.  A stand-in filesystem maps s3:// URIs onto
local files, so the remote branch of each reader is exercised for real while
the bytes come off the local disk -- what is under test is the plumbing
(is_remote_path -> fs.open -> h5py/netCDF4/json), not S3 itself.
"""
import json
import numpy as np
import pytest
import pointCollection as pc


class FakeS3FS:
    """
    Minimal stand-in for s3fs.S3FileSystem.

    Maps s3:// URIs onto local paths and records what was opened, so a test
    can assert the remote branch was taken rather than a local fallback.
    """
    def __init__(self, mapping):
        self.mapping = dict(mapping)
        self.opened = []

    def open(self, path, mode='rb'):
        self.opened.append((path, mode))
        return open(self.mapping[path], mode)

    def exists(self, path):
        return path in self.mapping


@pytest.fixture
def grid():
    """a small 2-band grid with identifiable values"""
    return pc.grid.data().from_dict({
        'x': np.arange(4, dtype=float) * 100.,
        'y': np.arange(3, dtype=float) * 100.,
        't': np.array([2019., 2020.]),
        'z': np.arange(24, dtype=float).reshape((3, 4, 2))})


# ---------------------------------------------------------------------------
# io_utils.as_gdal_path
# ---------------------------------------------------------------------------

def test_as_gdal_path_translates_uris():
    assert pc.io_utils.as_gdal_path('s3://bucket/key.tif') == '/vsis3/bucket/key.tif'
    assert pc.io_utils.as_gdal_path('gs://bucket/key.tif') == '/vsigs/bucket/key.tif'
    assert pc.io_utils.as_gdal_path('https://host/key.tif') == '/vsicurl/https://host/key.tif'


def test_as_gdal_path_leaves_local_and_vsi_paths_alone():
    for path in ('/local/key.tif', 'relative/key.tif', '/vsis3/bucket/key.tif'):
        assert pc.io_utils.as_gdal_path(path) == path


# ---------------------------------------------------------------------------
# io_utils.get_s3fs(daac=None)
# ---------------------------------------------------------------------------

def test_get_s3fs_none_daac_uses_default_credentials():
    """
    daac=None must build a plain s3fs session rather than an earthaccess one:
    the buckets we own are not DAAC holdings and earthaccess credentials do
    not reach them.
    """
    s3fs = pytest.importorskip('s3fs')
    fs = pc.io_utils.get_s3fs(daac=None)
    assert isinstance(fs, s3fs.S3FileSystem)
    # cached, and keyed separately from the DAAC sessions
    assert pc.io_utils.get_s3fs(daac=None) is fs


# ---------------------------------------------------------------------------
# grid.data.h5_open / from_h5
# ---------------------------------------------------------------------------

def test_h5_open_remote_reads_through_fs(tmp_path, grid):
    local = str(tmp_path / 'grid.h5')
    grid.to_h5(local, group='/')
    fs = FakeS3FS({'s3://bucket/grid.h5': local})

    with pc.grid.data().h5_open('s3://bucket/grid.h5', fs=fs) as h5f:
        assert 'z' in h5f
    assert fs.opened == [('s3://bucket/grid.h5', 'rb')]


def test_from_h5_remote_matches_local(tmp_path, grid):
    local = str(tmp_path / 'grid.h5')
    grid.to_h5(local, group='/')
    fs = FakeS3FS({'s3://bucket/grid.h5': local})

    remote_read = pc.grid.data().from_h5('s3://bucket/grid.h5', fs=fs)
    local_read = pc.grid.data().from_h5(local)

    assert np.allclose(remote_read.z, local_read.z)
    assert fs.opened, 'remote branch was not taken'


def test_from_h5_local_path_ignores_fs(tmp_path, grid):
    """a local path must not be routed through the filesystem object"""
    local = str(tmp_path / 'grid.h5')
    grid.to_h5(local, group='/')
    fs = FakeS3FS({})

    assert pc.grid.data().from_h5(local, fs=fs).z is not None
    assert fs.opened == []


# ---------------------------------------------------------------------------
# grid.data.nc_open / from_nc
# ---------------------------------------------------------------------------

def test_from_nc_remote_matches_local(tmp_path, grid):
    pytest.importorskip('netCDF4')
    local = str(tmp_path / 'grid.nc')
    grid.to_nc(local)
    fs = FakeS3FS({'s3://bucket/grid.nc': local})

    remote_read = pc.grid.data().from_nc('s3://bucket/grid.nc', fs=fs)
    local_read = pc.grid.data().from_nc(local)

    assert np.allclose(remote_read.z, local_read.z)
    assert fs.opened == [('s3://bucket/grid.nc', 'rb')]


# ---------------------------------------------------------------------------
# tilingSchema.from_file
# ---------------------------------------------------------------------------

def test_tiling_schema_from_remote_json(tmp_path):
    local = str(tmp_path / 'schema.json')
    pc.tilingSchema(tile_spacing=2.e5, format_str='FOO_E%d_N%d',
                    scale=1.e3).to_json(local)
    fs = FakeS3FS({'s3://bucket/tiles/schema.json': local})

    ts = pc.tilingSchema().from_file('s3://bucket/tiles/schema.json', fs=fs)
    assert ts.tile_spacing == 2.e5
    # the directory default is derived from the schema's own location, so a
    # remote schema must point at the remote directory holding its tiles
    assert ts.directory == 's3://bucket/tiles'
    assert fs.opened == [('s3://bucket/tiles/schema.json', 'r')]


# ---------------------------------------------------------------------------
# geoIndex.from_file / query_ATL11_cloud with a remote index
# ---------------------------------------------------------------------------
# The per-granule geoIndex lives in our own bucket even when the granules it
# indexes are DAAC holdings, so it is read with the default AWS credentials
# rather than the earthaccess session used for the granule itself.

import os
import shutil
from pointCollection.scripts.query_ATL11_cloud import read_ATL11_granule_cloud_items

TEST_H5 = os.path.join(os.path.dirname(__file__), '..', 'test_data',
                       'ATL06_20190205041106_05910210_209_01.h5')
SRS_PROJ4 = ('+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 '
             '+datum=WGS84 +units=m +no_defs')


def _write_index(tmp_path):
    """build a geoIndex over the ATL06 fixture and return its local path"""
    index_file = str(tmp_path / 'index.h5')
    pc.geoIndex(delta=[1.e4, 1.e4], SRS_proj4=SRS_PROJ4).for_file(
        TEST_H5, 'ATL06', number=0).to_file(index_file)
    return index_file


def test_geoindex_from_remote_file(tmp_path):
    index_file = _write_index(tmp_path)
    fs = FakeS3FS({'s3://bucket/index.h5': index_file})

    gI = pc.geoIndex().from_file('s3://bucket/index.h5', fs=fs)
    assert gI.attrs is not None
    assert fs.opened == [('s3://bucket/index.h5', 'rb')]


def test_remote_index_is_not_reported_missing(tmp_path):
    """
    Regression: os.path.isfile() is False for any URI, so a remote index used
    to be reported missing and skipped -- which returned an empty tile instead
    of failing, and would have silently emptied every tile of a DPS run.
    """
    granule = tmp_path / os.path.basename(TEST_H5)
    shutil.copy(TEST_H5, granule)
    index_file = _write_index(tmp_path)
    index_fs = FakeS3FS({'s3://bucket/index.h5': index_file})

    import warnings
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        items = read_ATL11_granule_cloud_items(
            str(granule), 's3://bucket/index.h5', [-1.e8, 1.e8], [-1.e8, 1.e8],
            index_fs=index_fs)
    assert not [w for w in caught if 'missing geoIndex' in str(w.message)]
    assert items is not None and len(items) > 0


def test_missing_remote_index_still_skips(tmp_path):
    """a remote index that really is absent must still warn and skip"""
    index_fs = FakeS3FS({})
    import warnings
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        items = read_ATL11_granule_cloud_items(
            'ATL11_044110_0331_007_04.h5', 's3://bucket/no_such_index.h5',
            [-1, 1], [-1, 1], index_fs=index_fs)
    assert items is None
    assert any('missing geoIndex' in str(w.message) for w in caught)
