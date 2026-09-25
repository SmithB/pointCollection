"""
open_remote's block cache: a windowed read does not refetch blocks.

No network access is required.  CountingFS is a real fsspec filesystem whose
files are local bytes, and it counts the bytes each range request fetches, so
fsspec's own caches -- the thing under test -- run exactly as they do over S3.

The HDF5 file has many compressed fields in several groups, like an ATL11
granule.  Reading ONE window of every field (what a geoIndex query does, one
merged range per granule and beam pair) with fsspec's default one-block
'readahead' cache fetches some blocks again and again: on real ATL11 a
~1,400-point range read pulled ~30 MB, against ~7 MB with the block cache,
and a whole IS tile 1.71 GB against 0.26 GB (2026-09-25).  This file shows the
same effect at a smaller ratio.
"""
import h5py
import numpy as np
import pytest
from fsspec.spec import AbstractBufferedFile, AbstractFileSystem

import pointCollection as pc

N = 100_000
NAMES = [f'group{g}/field{i}' for g in range(4) for i in range(10)]
WINDOW = slice(50_000, 51_400)


class CountingFile(AbstractBufferedFile):
    def _fetch_range(self, start, end):
        data = self.fs.blobs[self.path][start:end]
        self.fs.fetched += len(data)
        return data


class CountingFS(AbstractFileSystem):
    protocol = 'counting'
    # fsspec reuses a filesystem instance built with the same arguments; each
    # test needs its own byte count
    cachable = False

    def __init__(self, blobs):
        super().__init__()
        self.blobs = blobs
        self.fetched = 0

    def _open(self, path, mode='rb', block_size=None, cache_type='readahead',
              cache_options=None, **kwargs):
        return CountingFile(self, path, mode, block_size=block_size or 5 * 2**20,
                            cache_type=cache_type, cache_options=cache_options,
                            size=len(self.blobs[path]))


@pytest.fixture(scope='module')
def granule(tmp_path_factory):
    path = tmp_path_factory.mktemp('remote_cache') / 'granule.h5'
    rng = np.random.default_rng(0)
    with h5py.File(path, 'w') as h5f:
        # create every field first and fill them afterwards, so the chunks
        # are laid out apart from the metadata, as in a written granule
        dsets = [h5f.create_dataset(name, shape=(N,), dtype='f8', chunks=(10_000,),
                                    compression='gzip') for name in NAMES]
        for ds in dsets:
            ds[:] = rng.normal(size=N)
    return path.read_bytes()


def read_window(fd):
    with h5py.File(fd, 'r') as h5f:
        return {name: h5f[name][WINDOW] for name in NAMES}


def test_windowed_read_does_not_refetch(granule):
    fs = CountingFS({'granule.h5': granule})
    with pc.io_utils.open_remote('granule.h5', fs=fs,
                                 block_size=pc.io_utils.DEFAULT_REMOTE_BLOCK_SIZE) as fd:
        cached = read_window(fd)

    # the same read on the old one-block cache
    old = CountingFS({'granule.h5': granule})
    with old.open('granule.h5', 'rb', block_size=pc.io_utils.DEFAULT_REMOTE_BLOCK_SIZE,
                  cache_type='readahead') as fd:
        plain = read_window(fd)

    assert old.fetched > 1.5 * fs.fetched
    for name in NAMES:
        assert np.array_equal(cached[name], plain[name])


def test_block_cache_requested_only_with_a_block_size():
    seen = []

    class RecordingFS:
        def open(self, path, mode='rb', **kwargs):
            seen.append(kwargs)

    pc.io_utils.open_remote('s3://b/k.h5', fs=RecordingFS(), block_size=256 * 1024)
    pc.io_utils.open_remote('s3://b/k.h5', fs=RecordingFS())
    assert seen[0]['cache_type'] == 'blockcache'
    assert (seen[0]['cache_options']['maxblocks'] * 256 * 1024
            == pc.io_utils.DEFAULT_REMOTE_CACHE_BYTES)
    assert seen[1] == {}      # a whole-file read keeps the filesystem's defaults
