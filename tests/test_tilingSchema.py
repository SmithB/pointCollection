"""
Tests for pointCollection.tilingSchema, focused on the 'source' addition
that lets a schema point at a remote (e.g. EarthAccess) collection instead
of a local directory. No network access is required.
"""
import os
import json
import numpy as np
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
