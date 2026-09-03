"""
Tests for pointCollection.indexedH5
"""
import numpy as np
import h5py
import pytest
import pointCollection as pc


@pytest.fixture
def indexed_file(tmp_path):
    """an indexedH5 file in the bin-major layout that to_file() writes"""
    rng = np.random.default_rng(0)
    n = 60
    D = pc.data().from_dict({'x': rng.uniform(-2.e4, 2.e4, n),
                             'y': rng.uniform(-2.e4, 2.e4, n),
                             'z': np.arange(n, dtype=float),
                             'time': np.zeros(n)})
    out_file = str(tmp_path / 'tile.h5')
    pc.indexedH5.data(bin_W=(1.e4, 1.e4)).to_file(D, out_file)
    return out_file


# ---------------------------------------------------------------------------
# fields=None.  read() did `fields.copy()` unconditionally, so passing None --
# which is geoIndex.get_data()'s default, and means 'read everything' for the
# other file types -- raised AttributeError. geoIndex swallowed it, so an
# indexedH5 query with no explicit fields silently returned nothing.
# ---------------------------------------------------------------------------

def test_read_fields_none_reads_every_field(indexed_file):
    xy = [np.array([0.]), np.array([0.])]
    D_none = pc.indexedH5.data(filename=indexed_file).read(xy, fields=None)
    D_explicit = pc.indexedH5.data(filename=indexed_file).read(
        xy, fields=['x', 'y', 'z', 'time'])
    assert sorted(D_none.fields) == ['time', 'x', 'y', 'z']
    assert D_none.size == D_explicit.size > 0
    np.testing.assert_array_equal(np.sort(D_none.z), np.sort(D_explicit.z))


def test_read_default_fields_unchanged(indexed_file):
    # only an explicit None means 'everything'; the signature default is
    # still the ['x','y','time'] subset
    D = pc.indexedH5.data(filename=indexed_file).read(
        [np.array([0.]), np.array([0.])])
    assert sorted(D.fields) == ['time', 'x', 'y']


def test_fields_in_file_bin_major(indexed_file):
    with h5py.File(indexed_file, 'r') as h5f:
        # top-level keys here are bins ('<x>E_<y>N'), not fields
        assert all('E_' in key for key in h5f.keys())
        assert sorted(pc.indexedH5.data.fields_in_file(h5f)) == ['time', 'x', 'y', 'z']


def test_fields_in_file_with_top_level_index(tmp_path):
    # the other layout: one dataset per field at the top level, alongside an
    # 'INDEX' group that is not itself a field
    flat_file = str(tmp_path / 'flat.h5')
    with h5py.File(flat_file, 'w') as h5f:
        for field in ['x', 'y', 'z', 'time']:
            h5f.create_dataset(field, data=np.arange(5, dtype=float))
        index = h5f.create_group('INDEX')
        index.create_dataset('bin_x', data=np.array([0.]))
        index.create_dataset('bin_y', data=np.array([0.]))
    with h5py.File(flat_file, 'r') as h5f:
        assert sorted(pc.indexedH5.data.fields_in_file(h5f)) == ['time', 'x', 'y', 'z']
