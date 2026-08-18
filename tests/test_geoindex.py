"""
Tests for pointCollection.geoIndex
"""
import os
import numpy as np
import pytest
import pointCollection as pc

# NOTE: from_xy() stores each bin's offset_start/offset_end as *inclusive*
# first/last row indices, but readers (data.py, ATL06/data.py) treat
# index_range as an exclusive-end Python slice. get_data()'s trim_last_point
# parameter (default False) compensates by adding 1 to offset_end before
# handing it to the 'h5', 'ATL06', 'ATL11', and 'ATM_Qfit' branches, so the
# last point of a multi-point bin is included by default; trim_last_point=True
# reproduces the old (last-point-dropped) behavior. Not applied to
# 'indexed_h5'/'indexed_h5_from_matlab' (offsets there can be a -1 sentinel,
# or come from an externally-built index of unverified convention) or a
# user-supplied `function` (always gets the raw, unmodified offsets).


def test_query_xy_relative_dir_root(tmp_path, monkeypatch):
    # geoIndex.query_xy() (full_path=True, the default) used to resolve each
    # file's path via resolve_path(), then get_data() resolved it again,
    # doubling any *relative* dir_root/index-path prefix and making the file
    # unreadable. Use a dir_root with more than one path component so the
    # doubled prefix can't accidentally self-cancel.
    data_dir = tmp_path / 'data' / 'subdir'
    data_dir.mkdir(parents=True)
    sub_dir = tmp_path / 'sub'
    sub_dir.mkdir()

    x = np.array([0., 1., 20.])
    y = np.array([0., 0., 0.])
    D = pc.data(fields={'x': x, 'y': y})
    D.to_h5(str(data_dir / 'file0.h5'))

    gi = pc.geoIndex(delta=[10, 10]).from_xy((x, y), filename='file0.h5',
                                              file_type='h5', number=0)
    gi.filename = str(sub_dir / 'index.h5')
    gi.attrs['dir_root'] = '../data/subdir/'

    monkeypatch.chdir(sub_dir)
    result = gi.query_xy((np.array([0.]), np.array([0.])), full_path=True,
                          get_data=True, fields=['x', 'y'])

    assert result
    np.testing.assert_array_equal(np.sort(result[0].x), [0., 1.])
    assert result[0].filename == os.path.join('..', 'data', 'subdir', 'file0.h5')


# ---------------------------------------------------------------------------
# self-contained (index + data in one file) via for_file(self_contained=True)
# ---------------------------------------------------------------------------

def test_self_contained_index_and_data(tmp_path):
    x = np.array([0., 1., 2., 20.])
    y = np.array([0., 0., 0., 0.])
    D = pc.data(fields={'x': x, 'y': y})
    path = str(tmp_path / 'combined.h5')
    D.to_h5(path, group='mydata')

    gi = pc.geoIndex(delta=[10, 10]).for_file(path, 'h5', group='mydata',
                                               self_contained=True)
    gi.to_file(path)  # writes the 'index' group into the SAME file

    gi2 = pc.geoIndex().from_file(path)
    result = gi2.query_xy((np.array([0.]), np.array([0.])), full_path=True,
                           get_data=True)

    assert result
    np.testing.assert_array_equal(np.sort(result[0].x), [0., 1., 2.])


def test_self_contained_relative_index_path(tmp_path, monkeypatch):
    # regression test: the ':group' self-reference used to double-prepend
    # directories when the index's own file path was relative *and included
    # a directory component* -- a bare filename (no '/') doesn't trigger the
    # bug, since prepending an empty dirname is a no-op, so this deliberately
    # references the file as 'sub/combined.h5' rather than chdir-ing into
    # 'sub' and using a bare name.
    sub_dir = tmp_path / 'sub'
    sub_dir.mkdir()
    x = np.array([0., 1., 2., 20.])
    y = np.array([0., 0., 0., 0.])
    D = pc.data(fields={'x': x, 'y': y})
    D.to_h5(str(sub_dir / 'combined.h5'), group='mydata')

    monkeypatch.chdir(tmp_path)
    rel_path = os.path.join('sub', 'combined.h5')
    gi = pc.geoIndex(delta=[10, 10]).for_file(rel_path, 'h5', group='mydata',
                                               self_contained=True)
    gi.to_file(rel_path)

    gi2 = pc.geoIndex().from_file(rel_path)  # self.filename is relative, with a directory component
    # full_path=True: query_xy() builds the ':group' reference itself and
    # never calls resolve_path() on it.
    result = gi2.query_xy((np.array([0.]), np.array([0.])), full_path=True,
                           get_data=True)
    assert result
    np.testing.assert_array_equal(np.sort(result[0].x), [0., 1., 2.])

    # full_path=False: get_data() calls resolve_path() on the raw entry
    # itself (already_resolved=False) -- this is what resolve_path()'s
    # self-referential-filename guard specifically protects.
    result = gi2.query_xy((np.array([0.]), np.array([0.])), full_path=False,
                           get_data=True)
    assert result
    np.testing.assert_array_equal(np.sort(result[0].x), [0., 1., 2.])


def test_trim_last_point_legacy_behavior(tmp_path):
    x = np.array([0., 1., 2., 20.])
    y = np.array([0., 0., 0., 0.])
    D = pc.data(fields={'x': x, 'y': y})
    path = str(tmp_path / 'combined.h5')
    D.to_h5(path, group='mydata')

    gi = pc.geoIndex(delta=[10, 10]).for_file(path, 'h5', group='mydata',
                                               self_contained=True)
    gi.to_file(path)

    gi2 = pc.geoIndex().from_file(path)

    # default (trim_last_point=False): the bin's last point is included
    result = gi2.query_xy((np.array([0.]), np.array([0.])), full_path=True,
                           get_data=True)
    assert result
    np.testing.assert_array_equal(np.sort(result[0].x), [0., 1., 2.])

    # trim_last_point=True: reproduces the old (last-point-dropped) behavior
    result = gi2.query_xy((np.array([0.]), np.array([0.])), full_path=True,
                           get_data=True, trim_last_point=True)
    assert result
    np.testing.assert_array_equal(np.sort(result[0].x), [0., 1.])


def test_adjacent_bins_no_gap_no_duplicate(tmp_path):
    # points 0 and 1 fall in bin x=0 (rows 0-1), points 2 and 3 fall in the
    # adjacent bin x=10 (rows 2-3) -- the two bins' offset ranges are
    # touching (row 1's bin ends where row 2's bin starts), which is exactly
    # the case query_xy()'s merge/cleanup step collapses into a single
    # contiguous read. Querying both bins together should return all 4
    # points exactly once each -- no gap at the boundary, no duplicate.
    x = np.array([0., 3., 10., 13.])
    y = np.array([0., 0., 0., 0.])
    D = pc.data(fields={'x': x, 'y': y})
    path = str(tmp_path / 'combined.h5')
    D.to_h5(path, group='mydata')

    gi = pc.geoIndex(delta=[10, 10]).for_file(path, 'h5', group='mydata',
                                               self_contained=True)
    gi.to_file(path)

    gi2 = pc.geoIndex().from_file(path)
    result = gi2.query_xy((np.array([0., 10.]), np.array([0., 0.])),
                           full_path=True, get_data=True)

    assert result
    assert len(result) == 1  # touching bins merge into a single read
    np.testing.assert_array_equal(np.sort(result[0].x), [0., 3., 10., 13.])


def test_self_contained_requires_group(tmp_path):
    x = np.array([0., 1., 2.])
    y = np.array([0., 0., 0.])
    D = pc.data(fields={'x': x, 'y': y})
    path = str(tmp_path / 'combined.h5')
    D.to_h5(path, group='mydata')

    with pytest.raises(ValueError):
        pc.geoIndex(delta=[10, 10]).for_file(path, 'h5', self_contained=True)


def test_self_contained_requires_h5_type(tmp_path):
    path = str(tmp_path / 'combined.h5')
    with pytest.raises(ValueError):
        pc.geoIndex(delta=[10, 10]).for_file(path, 'ATL06', group='mydata',
                                              self_contained=True)


def test_self_contained_nc_extension(tmp_path):
    # combined index+data file, but named '.nc' -- still real HDF5 bytes
    # written via to_h5()/to_file(), proving nothing in the new code paths
    # assumes a '.h5' extension.
    x = np.array([0., 1., 2., 20.])
    y = np.array([0., 0., 0., 0.])
    D = pc.data(fields={'x': x, 'y': y})
    path = str(tmp_path / 'combined.nc')
    D.to_h5(path, group='mydata')

    gi = pc.geoIndex(delta=[10, 10]).for_file(path, 'h5', group='mydata',
                                               self_contained=True)
    gi.to_file(path)

    gi2 = pc.geoIndex().from_file(path)
    result = gi2.query_xy((np.array([0.]), np.array([0.])), full_path=True,
                           get_data=True)

    assert result
    np.testing.assert_array_equal(np.sort(result[0].x), [0., 1., 2.])


# ---------------------------------------------------------------------------
# netCDF4 read support (via h5py, since netCDF4 is an HDF5 container format)
# ---------------------------------------------------------------------------

def test_read_netcdf_via_h5_type(tmp_path):
    netCDF4 = pytest.importorskip('netCDF4')

    # build the fixture with the netCDF4 library directly (not pc.data), to
    # genuinely test cross-tool compatibility rather than self-consistency.
    nc_path = str(tmp_path / 'external.nc')
    with netCDF4.Dataset(nc_path, 'w') as nc:
        grp = nc.createGroup('mydata')
        grp.createDimension('n', 4)
        xvar = grp.createVariable('x', 'f8', ('n',))
        yvar = grp.createVariable('y', 'f8', ('n',))
        xvar[:] = [0., 1., 2., 20.]
        yvar[:] = [0., 0., 0., 0.]

    D = pc.data().from_h5(nc_path, group='mydata')
    np.testing.assert_array_equal(np.sort(D.x), [0., 1., 2., 20.])

    gi = pc.geoIndex(delta=[10, 10]).for_file(nc_path, 'h5', group='mydata')
    result = gi.query_xy((np.array([0.]), np.array([0.])), full_path=True,
                          get_data=True, fields={'mydata': ['x', 'y']})

    assert result
    np.testing.assert_array_equal(np.sort(result[0].x), [0., 1., 2.])
