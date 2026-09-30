"""
grid.mosaic.from_list and make_mosaic.py: reading the tiles with a block_size,
in a process pool (workers), and from a remote directory must give exactly the
mosaic the original serial, local read gives.

No network: remote files live in fsspec's in-memory filesystem.  Only a
FORKED worker sees it (and the patched session getter), so the remote tests
with workers switch the pool to fork; the local ones run the default,
forkserver.
"""
import importlib
import os
import sys
import fsspec
import numpy as np
import pytest
import pointCollection as pc

CENTERS = [0., 40., 80.]
T = np.arange(4.)


def write_tile(filename, xc, yc, rng):
    """a 60 x 60 tile at 2-unit spacing, a 2-D and a 3-D group, some NaNs"""
    x = xc + np.arange(-30., 30.1, 2.)
    y = yc + np.arange(-30., 30.1, 2.)
    z0 = rng.normal(size=(y.size, x.size))
    dz = rng.normal(size=(y.size, x.size, T.size))
    dz[rng.random(dz.shape) < 0.05] = np.nan
    pc.grid.data().from_dict({'x': x, 'y': y, 'z0': z0,
                              'cell_area': np.ones_like(z0)}).to_h5(filename, group='z0', replace=True)
    pc.grid.data().from_dict({'x': x, 'y': y, 't': T, 'dz': dz}).to_h5(filename, group='dz', replace=False)


@pytest.fixture(scope='module')
def tiles(tmp_path_factory):
    d = tmp_path_factory.mktemp('tiles')
    rng = np.random.default_rng(3)
    files = []
    for xc in CENTERS:
        for yc in CENTERS:
            files.append(str(d / f'E{int(xc)}_N{int(yc)}.h5'))
            write_tile(files[-1], xc, yc, rng)
    return sorted(files)


@pytest.fixture
def memory_tiles(tiles, monkeypatch):
    """the same tiles in a memory filesystem, returned as memory:// URIs"""
    fs = fsspec.filesystem('memory')
    uris = []
    for f in tiles:
        uri = 'memory://bucket/region/matched/' + os.path.basename(f)
        with open(f, 'rb') as fh:
            fs.pipe(uri, fh.read())
        uris.append(uri)
    monkeypatch.setattr(pc.io_utils, 'get_s3fs', lambda daac=None, **kw: fs)
    # pc.grid.mosaic is the class; the module is shadowed by it
    monkeypatch.setattr(importlib.import_module('pointCollection.grid.mosaic'), '_START_METHOD', 'fork')
    yield uris
    fs.rm('memory://bucket', recursive=True)


CASES = {
    'weighted by band': dict(group='dz', fields=['dz'], pad=4, feather=8, by_band=True),
    'weighted all bands': dict(group='dz', fields=['dz'], pad=4, feather=8, by_band=False),
    'weighted 2-D': dict(group='z0', fields=['z0', 'cell_area'], pad=4, feather=8, by_band=False),
    'replace': dict(group='z0', fields=['z0'], pad=None, feather=None),
    'selected bands': dict(group='dz', fields=['dz'], pad=4, feather=8, by_band=True, bands=[1, 3]),
    'all fields': dict(group='z0', fields=None, pad=4, feather=8, by_band=False),
}


def mosaic(files, **kwargs):
    return pc.grid.mosaic().from_list(list(files), **kwargs)


def assert_same(a, b):
    assert a.fields == b.fields
    assert np.array_equal(a.x, b.x) and np.array_equal(a.y, b.y)
    for field in a.fields:
        assert np.array_equal(getattr(a, field), getattr(b, field), equal_nan=True), field


@pytest.mark.parametrize('case', CASES)
@pytest.mark.parametrize('options', [dict(workers=3), dict(block_size=4096), dict(workers=2, block_size=4096)],
                         ids=['workers', 'block_size', 'both'])
def test_local_read_options_change_nothing(tiles, case, options):
    assert_same(mosaic(tiles, **CASES[case]), mosaic(tiles, **CASES[case], **options))


# fork warns once the test process has threads (fsspec's loop); that is the
# reason the default is forkserver, and harmless for an in-memory filesystem
FORK_WARNING = pytest.mark.filterwarnings('ignore:This process .* is multi-threaded:DeprecationWarning')


@FORK_WARNING
@pytest.mark.parametrize('case', CASES)
def test_remote_tiles_give_the_local_mosaic(tiles, memory_tiles, case):
    reference = mosaic(tiles, **CASES[case])
    assert_same(reference, mosaic(memory_tiles, **CASES[case], block_size=4096))
    assert_same(reference, mosaic(memory_tiles, **CASES[case], block_size=4096, workers=3))


def test_default_start_method_is_forkserver():
    # plain fork deadlocks or fails ("not fork-safe") after s3fs has started
    assert importlib.import_module('pointCollection.grid.mosaic')._START_METHOD == 'forkserver'


def test_block_size_reaches_the_remote_open(memory_tiles, monkeypatch):
    seen = []
    open_remote = pc.io_utils.open_remote
    def spy(filename, *args, block_size=None, **kwargs):
        seen.append(block_size)
        return open_remote(filename, *args, block_size=block_size, **kwargs)
    monkeypatch.setattr(pc.io_utils, 'open_remote', spy)
    mosaic(memory_tiles, **CASES['weighted 2-D'], block_size=12345)
    assert seen and set(seen) == {12345}


def test_unreadable_tile_is_dropped_the_same_way(tiles, tmp_path):
    bad = str(tmp_path / 'E999_N999.h5')
    with open(bad, 'w') as fh:
        fh.write('not an hdf5 file')
    files = tiles[:4] + [bad] + tiles[4:]
    reference = mosaic(files, **CASES['weighted by band'])
    assert_same(reference, mosaic(files, **CASES['weighted by band'], workers=3))
    assert_same(reference, mosaic(tiles, **CASES['weighted by band']))


@pytest.mark.parametrize('kwargs', [dict(by_band=True), dict(by_band=True, bands=[0, 2]),
                                    dict(by_band=False), dict(by_band=False, bands=[0, 2])],
                         ids=['by band', 'by band, selected bands', 'all bands', 'selected bands'])
def test_grids_in_memory_pass_through(tiles, kwargs):
    grids = [pc.grid.mosaic().from_file(f, group='dz', fields=['dz']) for f in tiles[:3]]
    mixed = grids + tiles[3:]
    kwargs = dict(group='dz', fields=['dz'], pad=4, feather=8, **kwargs)
    assert_same(mosaic(mixed, **kwargs), mosaic(mixed, **kwargs, workers=2, block_size=4096))


@pytest.mark.parametrize('bands', [None, [1, 3]])
def test_in_memory_mosaic_by_band_matches_its_file(tiles, bands):
    # add_to_band used to slice an in-memory mosaic to a plain grid.data and
    # fail on update_spacing
    grids = [pc.grid.mosaic().from_file(f, group='dz', fields=['dz']) for f in tiles]
    kwargs = dict(group='dz', fields=['dz'], pad=4, feather=8, by_band=True, bands=bands)
    assert_same(mosaic(tiles, **kwargs), mosaic(grids, **kwargs))


def test_glob_remote_is_sorted_and_keeps_the_scheme(memory_tiles):
    found = pc.io_utils.glob_remote('memory://bucket/region/matched/E*_N0.h5')
    assert found == sorted(found)
    assert [os.path.basename(f) for f in found] == ['E0_N0.h5', 'E40_N0.h5', 'E80_N0.h5']
    assert all(f.startswith('memory://') for f in found)


def run_make_mosaic(monkeypatch, argv):
    from pointCollection.scripts import make_mosaic
    monkeypatch.setattr(sys, 'argv', ['make_mosaic.py'] + argv)
    make_mosaic.main()


@FORK_WARNING
def test_make_mosaic_remote_directory(tiles, memory_tiles, tmp_path, monkeypatch):
    common = ['-g', 'E*.h5', '-w', '-p', '4', '-f', '8', '--in_group', 'dz/', '-F', 'dz', '-R']
    run_make_mosaic(monkeypatch, ['-d', os.path.dirname(tiles[0]), '-O', str(tmp_path / 'local.h5')] + common)
    run_make_mosaic(monkeypatch, ['-d', 'memory://bucket/region/matched', '-O', str(tmp_path / 'remote.h5'),
                                  '-j', '2'] + common)
    local = pc.grid.data().from_h5(str(tmp_path / 'local.h5'), group='dz')
    remote = pc.grid.data().from_h5(str(tmp_path / 'remote.h5'), group='dz')
    assert np.array_equal(local.dz, remote.dz, equal_nan=True)


@pytest.mark.parametrize('output', ['mosaic.h5', 's3://bucket/mosaic.h5'])
def test_make_mosaic_remote_directory_needs_a_local_absolute_output(memory_tiles, monkeypatch, output):
    with pytest.raises(SystemExit) as exit_info:
        run_make_mosaic(monkeypatch, ['-d', 'memory://bucket/region/matched', '-g', 'E*.h5',
                                      '--in_group', 'z0/', '-F', 'z0', '-O', output])
    assert exit_info.value.code == 2
