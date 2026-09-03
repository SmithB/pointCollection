"""
Tests for pointCollection.io_utils and the cloud-related geoIndex additions
(remote_file= override, resolve_path() remote short-circuit).

No network access is required -- these tests exercise the substitution
logic using local files only.
"""
import os
import shutil
import numpy as np
import pointCollection as pc
from pointCollection.scripts.query_ATL11_cloud import (
    read_ATL11_granule_cloud, read_ATL11_granule_cloud_items)

TEST_H5 = os.path.join(os.path.dirname(__file__), '..', 'test_data',
                       'ATL06_20190205041106_05910210_209_01.h5')

SRS_PROJ4 = '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'


# ---------------------------------------------------------------------------
# io_utils
# ---------------------------------------------------------------------------

def test_is_remote_path():
    assert pc.io_utils.is_remote_path('s3://bucket/key.h5')
    assert pc.io_utils.is_remote_path('https://example.com/key.h5')
    assert not pc.io_utils.is_remote_path('/local/path/file.h5')
    assert not pc.io_utils.is_remote_path('relative/path/file.h5')
    assert not pc.io_utils.is_remote_path(None)


def test_strip_pair_suffix():
    assert pc.io_utils.strip_pair_suffix('foo.h5:pair1') == 'foo.h5'
    assert pc.io_utils.strip_pair_suffix('foo.h5') == 'foo.h5'
    assert pc.io_utils.strip_pair_suffix(None) is None


# ---------------------------------------------------------------------------
# geoIndex.resolve_path() remote short-circuit
# ---------------------------------------------------------------------------

def test_resolve_path_remote_short_circuit():
    gI = pc.geoIndex()
    gI.attrs['dir_root'] = '/some/dir'
    gI.filename = '/some/dir/index.h5'
    remote = 's3://bucket/key.h5:pair1'
    assert gI.resolve_path(remote, dir_root='/some/other/dir') == remote


# ---------------------------------------------------------------------------
# geoIndex query_xy_box(remote_file=...) end-to-end substitution
# ---------------------------------------------------------------------------

def test_query_xy_box_remote_file_override(tmp_path):
    second_path = tmp_path / 'copy_of_fixture.h5'
    shutil.copy(TEST_H5, second_path)

    gI = pc.geoIndex(delta=[1.e4, 1.e4], SRS_proj4=SRS_PROJ4).for_file(TEST_H5, 'ATL06', number=0)
    xr, yr = [-1.e8, 1.e8], [-1.e8, 1.e8]

    D = gI.query_xy_box(xr, yr, remote_file=str(second_path))
    assert D is not None and len(D) > 0
    for Di in D:
        assert Di.filename == str(second_path)
        assert Di.latitude.size > 0


# ---------------------------------------------------------------------------
# query_ATL11_cloud.read_ATL11_granule_cloud[_items]
# ---------------------------------------------------------------------------
# 'remote_file' here is a local path -- is_remote_path() sees it isn't a
# URI and opens it as an ordinary local file, so these exercise the
# index-lookup/version-check/substitution logic without needing S3 access.

def test_read_ATL11_granule_cloud_items_returns_raw_list(tmp_path):
    # basename must match what the index indexes -- the version-mismatch
    # check in read_ATL11_granule_cloud_items() compares basenames
    second_path = tmp_path / os.path.basename(TEST_H5)
    shutil.copy(TEST_H5, second_path)

    index_file = str(tmp_path / 'index.h5')
    pc.geoIndex(delta=[1.e4, 1.e4], SRS_proj4=SRS_PROJ4).for_file(
        TEST_H5, 'ATL06', number=0).to_file(index_file)
    xr, yr = [-1.e8, 1.e8], [-1.e8, 1.e8]

    items = read_ATL11_granule_cloud_items(str(second_path), index_file, xr, yr)
    assert isinstance(items, list) and len(items) > 0
    for Di in items:
        assert Di.filename == str(second_path)
        assert Di.latitude.size > 0


def test_read_ATL11_granule_cloud_merges_items(tmp_path):
    second_path = tmp_path / os.path.basename(TEST_H5)
    shutil.copy(TEST_H5, second_path)

    index_file = str(tmp_path / 'index.h5')
    pc.geoIndex(delta=[1.e4, 1.e4], SRS_proj4=SRS_PROJ4).for_file(
        TEST_H5, 'ATL06', number=0).to_file(index_file)
    xr, yr = [-1.e8, 1.e8], [-1.e8, 1.e8]

    items = read_ATL11_granule_cloud_items(str(second_path), index_file, xr, yr)
    merged = read_ATL11_granule_cloud(str(second_path), index_file, xr, yr)
    assert merged is not None
    assert merged.size == sum(Di.size for Di in items)


def test_read_ATL11_granule_cloud_items_missing_index_returns_none(tmp_path):
    items = read_ATL11_granule_cloud_items(
        'not-a-real-granule.h5', str(tmp_path / 'no_such_index.h5'), [-1, 1], [-1, 1])
    assert items is None


def test_read_ATL11_granule_cloud_items_version_mismatch(tmp_path):
    second_path = tmp_path / 'copy_of_fixture.h5'
    shutil.copy(TEST_H5, second_path)

    index_file = str(tmp_path / 'index.h5')
    pc.geoIndex(delta=[1.e4, 1.e4], SRS_proj4=SRS_PROJ4).for_file(
        TEST_H5, 'ATL06', number=0).to_file(index_file)
    xr, yr = [-1.e8, 1.e8], [-1.e8, 1.e8]

    # 'mismatched.h5' has a different basename than what the index actually
    # indexes (the original TEST_H5 fixture) -- should be flagged
    import pytest
    with pytest.raises(ValueError):
        read_ATL11_granule_cloud_items('mismatched.h5', index_file, xr, yr)

    import warnings
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        items = read_ATL11_granule_cloud_items('mismatched.h5', index_file, xr, yr,
                                                version_mismatch='skip')
    assert items is None
    assert len(w) == 1
