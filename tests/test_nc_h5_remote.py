"""
Opt-in tests that read real remote files, to confirm that a windowed
from_nc() moves a small fraction of a large granule.

Skipped by default: set PC_TEST_REMOTE=1 for the https test (no credentials
needed) and PC_TEST_S3=1 for the NSIDC direct-S3 tests, which need earthaccess
credentials and only work in us-west-2.  Outside that region, earthaccess.open()
over https supports range requests and works the same way.
"""
import os
import time
import numpy as np
import pytest
import pointCollection as pc


def fetched_bytes(file_object):
    """bytes an fsspec file object actually pulled, if it tracks that"""
    cache = getattr(file_object, 'cache', None)
    return getattr(cache, 'total_requested_bytes', None)


class RecordedFile:
    """
    Proxy that snapshots the byte counter on close: fsspec drops the cache
    (and the counter with it) when a file is closed, and from_nc() closes the
    file it read from.
    """
    def __init__(self, fd):
        self.fd = fd
        self.fetched = None

    def close(self):
        if self.fetched is None:
            self.fetched = fetched_bytes(self.fd)
        self.fd.close()

    def __getattr__(self, name):
        return getattr(self.fd, name)


class RecordingFS:
    """wraps a filesystem so a test can look at the file objects afterwards"""
    def __init__(self, fs):
        self.fs = fs
        self.files = []

    def open(self, path, mode='rb', **kwargs):
        fd = RecordedFile(self.fs.open(path, mode, **kwargs))
        self.files.append(fd)
        return fd

    def __getattr__(self, name):
        return getattr(self.fs, name)


@pytest.mark.skipif(not os.environ.get('PC_TEST_REMOTE'),
                    reason='set PC_TEST_REMOTE=1 to read over the network')
def test_windowed_https_read_is_chunkwise():
    """
    A public, gzip-chunked netCDF4 file read over plain https range requests:
    no credentials, so this covers the remote plumbing anywhere.  Reading it
    at two block sizes also shows the block size doing what it is there for.
    """
    fsspec = pytest.importorskip('fsspec')
    url = ('https://noaa-goes16.s3.amazonaws.com/ABI-L2-SSTF/2020/001/00/'
           'OR_ABI-L2-SSTF-M6_G16_s20200010000216_e20200010059524_'
           'c20200010106082.nc')
    fs = RecordingFS(fsspec.filesystem('https'))
    size = fs.size(url)

    # x and y are packed int16 with a scale_factor, so take them as from_nc
    # sees them rather than raw out of h5py
    grid = pc.grid.data().from_nc(url, fs=fs, meta_only=True)
    x, y = np.sort(grid.x), np.sort(grid.y)
    bounds = [[x[2000], x[2500]], [y[2000], y[2500]]]

    reads = {}
    for block_size in (5 * 2**20, pc.io_utils.DEFAULT_REMOTE_BLOCK_SIZE,
                       64 * 1024):
        start = time.time()
        window = pc.grid.data().from_nc(url, fields=['SST'], fs=fs,
                                        bounds=bounds, block_size=block_size)
        elapsed = time.time() - start
        reads[block_size] = fs.files[-1].fetched
        assert window.SST.shape == (501, 501)
        print(f'\n{window.SST.shape} window, {block_size//1024} KiB blocks: '
              f'{reads[block_size]/2**20:.1f} MiB of a {size/2**20:.0f} MiB '
              f'file in {elapsed:.1f} s')

    # the whole point: a window costs a fraction of the granule.  This file
    # is only 29 MiB, so the fraction is modest -- the win grows with the
    # granule, which is what the ATL14 test below measures.
    tuned = reads[pc.io_utils.DEFAULT_REMOTE_BLOCK_SIZE]
    assert tuned is not None and tuned < size / 3
    # ... but only once the block size is sane.  fsspec caches one block, so
    # scattered chunk reads at 5 MiB a block re-fetch enough to beat
    # downloading the file outright -- which is why block_size is a parameter.
    assert tuned < reads[5 * 2**20] / 4


def centered_window(grid, half_width=30e3):
    """a 60 km window in the middle of a grid read with meta_only=True"""
    x0 = 0.5 * (grid.x[0] + grid.x[-1])
    y0 = 0.5 * (grid.y[0] + grid.y[-1])
    return [[x0 - half_width, x0 + half_width],
            [y0 - half_width, y0 + half_width]]


def granule_url(short_name, version='005', must_contain=()):
    """
    Direct-S3 URL of one granule.  ATL14/ATL15 collections carry the older
    0328 cycle range under the same short_name, so the caller filters on
    '_0329_': picking up the wrong one is a silent wrong-data bug, not an
    error.
    """
    earthaccess = pytest.importorskip('earthaccess')
    earthaccess.login()
    urls = []
    for result in earthaccess.search_data(short_name=short_name, version=version):
        links = result.data_links(access='direct')
        if links and all(token in links[0] for token in must_contain):
            urls.append(links[0])
    if not urls:
        pytest.skip(f'no {short_name} granule matching {must_contain}')
    return urls[0]


@pytest.mark.skipif(not os.environ.get('PC_TEST_S3'),
                    reason='set PC_TEST_S3=1 to read from NSIDC (us-west-2 only)')
def test_atl14_window_is_chunkwise():
    """
    ~1.4 GiB granule, 2-D float32, chunks (2491, 1401) gzip.  A 60 km window
    should move single-digit-to-tens of MiB and take about a second.
    """
    url = granule_url('ATL14', must_contain=('_0329_', '_GL_'))
    fs = RecordingFS(pc.io_utils.get_s3fs(daac='NSIDC'))
    size = fs.size(url)

    bounds = centered_window(pc.grid.data().from_nc(url, fs=fs, meta_only=True))
    start = time.time()
    window = pc.grid.data().from_nc(url, fields=['h'], fs=fs, bounds=bounds)
    elapsed = time.time() - start

    read = fs.files[-1].fetched
    print(f'\nATL14 {window.h.shape} window: {read/2**20:.1f} MiB of a '
          f'{size/2**20:.0f} MiB file in {elapsed:.1f} s')
    assert np.isfinite(window.h).any()
    assert read is not None and read < size / 10


@pytest.mark.skipif(not os.environ.get('PC_TEST_S3'),
                    reason='set PC_TEST_S3=1 to read from NSIDC (us-west-2 only)')
def test_atl15_banded_window_is_chunkwise():
    """
    ATL15 delta_h is the awkward case: chunks (8, 686, 386) gzip, so the band
    slicing has to be right or the read pulls whole time chunks it never uses.
    """
    url = granule_url('ATL15', must_contain=('_0329_', 'ATL15_A1_01km'))
    fs = RecordingFS(pc.io_utils.get_s3fs(daac='NSIDC'))
    size = fs.size(url)

    meta = pc.grid.data().from_nc(url, group='delta_h', fs=fs, meta_only=True)
    bounds = centered_window(meta)
    start = time.time()
    window = pc.grid.data().from_nc(url, group='delta_h', fields=['delta_h'],
                                    fs=fs, bounds=bounds,
                                    rdcc_nbytes=64 * 2**20)
    elapsed = time.time() - start

    read = fs.files[-1].fetched
    print(f'\nATL15 {window.delta_h.shape} window: {read/2**20:.1f} MiB of a '
          f'{size/2**20:.0f} MiB file in {elapsed:.1f} s')
    assert window.delta_h.shape[:2] == (window.y.size, window.x.size)
    assert read is not None and read < size / 10

    # a band subset must not cost more than the whole time series
    banded = pc.grid.data().from_nc(url, group='delta_h', fields=['delta_h'],
                                    fs=fs, bounds=bounds, bands=[0, 1])
    assert banded.delta_h.shape[2] == 2
    assert fs.files[-1].fetched <= read
