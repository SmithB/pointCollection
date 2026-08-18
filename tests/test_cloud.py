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
