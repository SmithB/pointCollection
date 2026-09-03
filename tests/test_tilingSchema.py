"""
Tests for pointCollection.tilingSchema.

Covers the 'source' addition that lets a schema point at a remote (e.g.
EarthAccess) collection instead of a local directory (no network access is
required), the mapping_function / mapping_function_name resolution, and
write_tiles().
"""
import os
import json
import numpy as np
import pytest
import pointCollection as pc


def test_source_round_trips_through_json(tmp_path):
    ts = pc.tilingSchema(tile_spacing=2.e5, format_str='FOO_E%d_N%d', scale=1.e3,
                          source={'type': 'EarthAccess', 'short_name': 'ATL11XO'})
    json_file = str(tmp_path / 'schema.json')
    ts.to_json(json_file)

    with open(json_file) as fh:
        raw = json.load(fh)
    assert raw['source'] == {'type': 'EarthAccess', 'short_name': 'ATL11XO'}

    ts2 = pc.tilingSchema().from_file(json_file)
    assert ts2.source == {'type': 'EarthAccess', 'short_name': 'ATL11XO'}
    # a remote source must not get a local directory default
    assert ts2.directory is None


def test_local_directory_default_unaffected(tmp_path):
    ts = pc.tilingSchema(tile_spacing=2.e5, format_str='FOO_E%d_N%d', scale=1.e3)
    json_file = str(tmp_path / 'sub' / 'schema.json')
    os.makedirs(os.path.dirname(json_file))
    ts.to_json(json_file)

    ts2 = pc.tilingSchema().from_file(json_file)
    assert ts2.source is None
    assert ts2.directory == os.path.dirname(json_file)


def test_tile_filename_bare_name_when_source_set():
    ts = pc.tilingSchema(tile_spacing=2.e5, format_str='ATL11XO_AR_E%d_N%d_c01_007_03',
                          format_variables=['x', 'y'], scale=1.e3, extension='.h5',
                          source={'type': 'EarthAccess', 'short_name': 'ATL11XO'})
    name = ts.tile_filename([-2.2e6, -1.4e6])
    assert name == 'ATL11XO_AR_E-2200_N-1400_c01_007_03.h5'
    # bare basename, no directory join
    assert os.path.sep not in name


def test_tile_filename_joins_directory_when_local():
    ts = pc.tilingSchema(tile_spacing=2.e5, format_str='E%d_N%d', format_variables=['x', 'y'],
                          scale=1.e3, extension='.h5', directory='/some/dir')
    name = ts.tile_filename([2.e5, 4.e5])
    assert name == os.path.join('/some/dir', 'E200_N400.h5')


def test_resolve_files_for_box_local_mode(tmp_path):
    # build two small tile files matching the naming convention, plus leave
    # one candidate tile unbuilt to confirm it's silently dropped
    ts = pc.tilingSchema(tile_spacing=2.e5, format_str='E%d_N%d', format_variables=['x', 'y'],
                          scale=1.e3, extension='.h5', bin_size=1.e4, directory=str(tmp_path))

    for xy in ([0., 0.], [2.e5, 0.]):
        path = ts.tile_filename(xy)
        pc.data(fields={'x': np.array([xy[0]]), 'y': np.array([xy[1]])}).to_h5(path)

    resolved, fs = ts.resolve_files_for_box([[-1.e5, 3.e5], [-1.e5, 1.e5]])
    assert fs is None
    assert set(resolved.keys()) == {'E0_N0.h5', 'E200_N0.h5'}
    assert resolved['E0_N0.h5'] == ts.tile_filename([0., 0.])

# ---------------------------------------------------------------------------
# mapping_function_name / mapping_function parsing
# ---------------------------------------------------------------------------
# tilingSchema.__init__ used to default `mapping_function=np.round` (instead
# of None), so `if mapping_function is not None:` was always true unless the
# caller *also* explicitly passed mapping_function=, silently discarding
# whatever `mapping_function_name` string was passed (e.g. 'floor') and
# forcing 'round'. Introduced in commit 5fe55c6 (2025-11-12); fixed by
# reverting the default back to None, restoring the original (pre-5fe55c6)
# behavior where the string alone is honored and lazily resolved by
# set_mapping_function() on first use.

def test_default_mapping_function_is_round():
    tS = pc.tilingSchema()
    assert tS.mapping_function_name == 'round'
    tS.tile_xy(xy=[np.array([0.]), np.array([0.])])
    assert tS.mapping_function is np.round


def test_mapping_function_name_only_selects_floor():
    # this is exactly the call pattern used by ATL11's make_ATL11xo_tiles.py:
    # only mapping_function_name is passed, mapping_function is left at its
    # default. Before the fix this silently stayed 'round'.
    tS = pc.tilingSchema(mapping_function_name='floor', tile_spacing=200000.)
    assert tS.mapping_function_name == 'floor'
    tS.tile_xy(xy=[np.array([0.]), np.array([0.])])
    assert tS.mapping_function is np.floor


def test_mapping_function_name_only_selects_round():
    tS = pc.tilingSchema(mapping_function_name='round', tile_spacing=200000.)
    tS.tile_xy(xy=[np.array([0.]), np.array([0.])])
    assert tS.mapping_function is np.round


def test_explicit_mapping_function_object_still_works():
    # explicitly passing the function object (bypassing the name) must
    # keep working, and should derive a matching mapping_function_name.
    tS = pc.tilingSchema(mapping_function=np.floor, tile_spacing=200000.)
    assert tS.mapping_function is np.floor
    assert tS.mapping_function_name == 'floor'


def test_explicit_mapping_function_overrides_conflicting_name():
    # when both are given and disagree, the explicit function object wins.
    tS = pc.tilingSchema(mapping_function_name='floor', mapping_function=np.round,
                          tile_spacing=200000.)
    assert tS.mapping_function is np.round
    assert tS.mapping_function_name == 'round'


# ---------------------------------------------------------------------------
# tile_xy(): floor and round must actually produce different tile
# assignments for the same data, once mapping_function_name is honored
# ---------------------------------------------------------------------------

def test_tile_xy_floor_vs_round_return_dict():
    tile_spacing = 200000.
    x = np.array([310000., 350000., 390000.])
    y = np.array([310000., 350000., 390000.])

    tS_floor = pc.tilingSchema(mapping_function_name='floor', tile_spacing=tile_spacing)
    tS_round = pc.tilingSchema(mapping_function_name='round', tile_spacing=tile_spacing)

    floor_keys = list(tS_floor.tile_xy(xy=[x.copy(), y.copy()], return_dict=True).keys())
    round_keys = list(tS_round.tile_xy(xy=[x.copy(), y.copy()], return_dict=True).keys())

    # floor: all points fall in [200000, 400000) -> corner (200000, 200000)
    assert floor_keys == [(200000.0, 200000.0)]
    # round: all points are nearest to 400000 -> center (400000, 400000)
    assert round_keys == [(400000.0, 400000.0)]


def test_tile_xy_return_dict_recovers_all_points():
    tile_spacing = 200000.
    x = np.array([310000., 350000., 390000.])
    y = np.array([310000., 350000., 390000.])

    tS = pc.tilingSchema(mapping_function_name='floor', tile_spacing=tile_spacing)
    bin_dict = tS.tile_xy(xy=[x, y], return_dict=True)
    (key, ii), = bin_dict.items()
    np.testing.assert_array_equal(np.sort(ii), [0, 1, 2])


# ---------------------------------------------------------------------------
# tile_bounds(): must not crash when called before any tile_xy() call has
# lazily resolved self.mapping_function, and must give the correct box for
# both conventions (round -> centered, floor -> corner-anchored).
# ---------------------------------------------------------------------------

def test_tile_bounds_floor_lazy_resolution():
    tile_spacing = 200000.
    tS = pc.tilingSchema(mapping_function_name='floor', tile_spacing=tile_spacing)
    # tile_xy() has never been called yet, so self.mapping_function is still
    # None -- tile_bounds() must resolve it itself rather than crashing.
    bounds = tS.tile_bounds(xy=[300000., 300000.])
    np.testing.assert_array_equal(bounds[0], [200000., 400000.])
    np.testing.assert_array_equal(bounds[1], [200000., 400000.])


def test_tile_bounds_round_lazy_resolution():
    tile_spacing = 200000.
    tS = pc.tilingSchema(mapping_function_name='round', tile_spacing=tile_spacing)
    # 350000 is unambiguously nearest the round-tile centered on 400000
    # (300000 would sit exactly on a tile boundary between two centers).
    bounds = tS.tile_bounds(xy=[350000., 350000.])
    np.testing.assert_array_equal(bounds[0], [300000., 500000.])
    np.testing.assert_array_equal(bounds[1], [300000., 500000.])


# ---------------------------------------------------------------------------
# scheme round-trip: from_file() updates mapping_function_name but not the
# already-resolved mapping_function; keeping mapping_function=None until
# first use (rather than eagerly resolving it in __init__) means a freshly
# constructed schema correctly picks up a loaded scheme's mapping function.
# ---------------------------------------------------------------------------

def test_scheme_json_roundtrip_preserves_floor(tmp_path):
    tile_spacing = 200000.
    tS = pc.tilingSchema(mapping_function_name='floor', tile_spacing=tile_spacing)
    json_file = str(tmp_path / 'scheme.json')
    tS.to_json(json_file)

    tS2 = pc.tilingSchema().from_file(json_file)
    assert tS2.mapping_function_name == 'floor'
    tS2.tile_xy(xy=[np.array([0.]), np.array([0.])])
    assert tS2.mapping_function is np.floor


# ---------------------------------------------------------------------------
# write_tiles(): __init__ sets data_format='indexedH5' (the pc.indexedH5 class
# name), but write_tiles used to test for 'indexed_h5', so neither branch ever
# matched and the call silently wrote nothing.  Both spellings now work, and
# bin_size falls back to the schema's own value instead of defaulting to None.
# ---------------------------------------------------------------------------

def _test_data(n=400):
    rng = np.random.default_rng(1)
    return pc.data().from_dict({'x': rng.uniform(-1.5e5, 1.5e5, n),
                                'y': rng.uniform(-1.5e5, 1.5e5, n),
                                'z': np.arange(n, dtype=float),
                                'time': np.zeros(n)})


def _tile_contents(directory):
    """(n_files, n_points, n_bins) over the indexedH5 tiles in a directory"""
    import h5py
    n_points = 0
    n_bins = 0
    files = sorted(os.listdir(directory))
    for f in files:
        with h5py.File(os.path.join(directory, f), 'r') as h5f:
            n_bins += len(h5f.keys())
            n_points += sum(h5f[g]['x'].size for g in h5f.keys())
    return len(files), n_points, n_bins


@pytest.mark.parametrize('data_format', ['indexedH5', 'indexed_h5', 'indexedh5'])
def test_write_tiles_writes_data(tmp_path, data_format):
    D = _test_data()
    tS = pc.tilingSchema(tile_spacing=1.e5, bin_size=1.e4, directory=str(tmp_path))
    tS.data_format = data_format
    tS.write_tiles(D)
    n_files, n_points, n_bins = _tile_contents(str(tmp_path))
    assert n_files > 1
    assert n_points == D.size      # every point landed in exactly one tile


def test_write_tiles_bin_size_defaults_to_schema(tmp_path):
    # write_tiles() took bin_size=None and passed it straight to
    # indexedH5.data(bin_W=(None, None)); it now falls back to self.bin_size.
    D = _test_data()
    fine, coarse = tmp_path / 'fine', tmp_path / 'coarse'
    for d, b in [(fine, 1.e4), (coarse, 5.e4)]:
        d.mkdir()
        pc.tilingSchema(tile_spacing=1.e5, bin_size=b, directory=str(d)).write_tiles(D)
    assert _tile_contents(str(fine))[2] > _tile_contents(str(coarse))[2]
    assert _tile_contents(str(fine))[1] == _tile_contents(str(coarse))[1] == D.size


def test_write_tiles_explicit_bin_size_overrides(tmp_path):
    D = _test_data()
    override_dir, coarse_dir = tmp_path / 'override', tmp_path / 'coarse'
    override_dir.mkdir()
    coarse_dir.mkdir()
    # a fine-binned schema told to write at 5e4 must match a 5e4 schema
    pc.tilingSchema(tile_spacing=1.e5, bin_size=1.e4,
                    directory=str(override_dir)).write_tiles(D, bin_size=5.e4)
    pc.tilingSchema(tile_spacing=1.e5, bin_size=5.e4,
                    directory=str(coarse_dir)).write_tiles(D)
    assert _tile_contents(str(override_dir))[2] == _tile_contents(str(coarse_dir))[2]


def test_write_tiles_rejects_unknown_data_format(tmp_path):
    tS = pc.tilingSchema(tile_spacing=1.e5, bin_size=1.e4, directory=str(tmp_path))
    tS.data_format = 'netcdf'
    with pytest.raises(ValueError, match='not understood'):
        tS.write_tiles(_test_data())


def test_write_tiles_rejects_remote_source(tmp_path):
    # tile_filename() returns bare granule names for a remote schema, so
    # writing would scatter files into the working directory
    tS = pc.tilingSchema(tile_spacing=1.e5, bin_size=1.e4,
                         source={'type': 'EarthAccess', 'short_name': 'ATL11XO'})
    with pytest.raises(ValueError, match='remote source'):
        tS.write_tiles(_test_data())
