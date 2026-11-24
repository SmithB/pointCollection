"""
grid.mosaic.from_list(tile_diff_fields=...) and make_mosaic.py --tile_rms_diff:
the RMS difference between overlapping tiles must equal a brute-force RMS over
every pair of tiles that contribute to a cell (nonzero weight, finite value).
"""
import itertools
import os
import sys
import numpy as np
import pytest
import pointCollection as pc

CENTERS = [0., 40., 80.]
T = np.arange(4.)
PAD, FEATHER = 4, 8


def write_tile(filename, xc, yc, rng, offset=0.):
    """a 60 x 60 tile at 2-unit spacing (tiles overlap by 20), 2-D and 3-D groups, some NaNs"""
    x = xc + np.arange(-30., 30.1, 2.)
    y = yc + np.arange(-30., 30.1, 2.)
    z0 = offset + rng.normal(size=(y.size, x.size))
    dz = rng.normal(size=(y.size, x.size, T.size))
    dz[rng.random(dz.shape) < 0.05] = np.nan
    pc.grid.data().from_dict({'x': x, 'y': y, 'z0': z0,
                              'cell_area': np.ones_like(z0)}).to_h5(filename, group='z0', replace=True)
    pc.grid.data().from_dict({'x': x, 'y': y, 't': T, 'dz': dz}).to_h5(filename, group='dz', replace=False)


def make_tiles(d, offset=0.):
    rng = np.random.default_rng(5)
    files = []
    for xc in CENTERS:
        for yc in CENTERS:
            files.append(str(d / f'E{int(xc)}_N{int(yc)}.h5'))
            write_tile(files[-1], xc, yc, rng, offset=offset)
    return sorted(files)


@pytest.fixture(scope='module')
def tiles(tmp_path_factory):
    return make_tiles(tmp_path_factory.mktemp('tiles'))


def brute_force(files, group, field, weighted):
    """count and pairwise RMS difference of the contributing tile values, cell by cell"""
    grids = [pc.grid.mosaic().from_file(f, group=group, fields=[field]) for f in files]
    x = np.unique(np.concatenate([g.x for g in grids]))
    y = np.unique(np.concatenate([g.y for g in grids]))
    shape = [y.size, x.size] + list(getattr(grids[0], field).shape[2:])
    values = []
    for g in grids:
        v = np.full(shape, np.nan)
        data = getattr(g, field).reshape([g.y.size, g.x.size] + shape[2:])
        if weighted:
            g.weights(pad=PAD, feather=FEATHER)
            contributes = g.weight > 0
        else:
            contributes = np.ones((g.y.size, g.x.size), dtype=bool)
        data = np.where(contributes.reshape(contributes.shape + (1,)*(len(shape)-2)), data, np.nan)
        v[np.ix_(np.searchsorted(y, g.y), np.searchsorted(x, g.x))] = data
        values.append(v)
    values = np.array(values)
    count = np.sum(np.isfinite(values), axis=0)
    ss = np.zeros(shape)
    n_pairs = np.zeros(shape)
    for a, b in itertools.combinations(values, 2):
        both = np.isfinite(a) & np.isfinite(b)
        ss[both] += (a[both]-b[both])**2
        n_pairs += both
    rms = np.full(shape, np.nan)
    rms[n_pairs > 0] = np.sqrt(ss[n_pairs > 0]/n_pairs[n_pairs > 0])
    return count, rms


CASES = {
    'weighted 2-D': dict(group='z0', field='z0', weighted=True, by_band=False),
    'weighted all bands': dict(group='dz', field='dz', weighted=True, by_band=False),
    'weighted by band': dict(group='dz', field='dz', weighted=True, by_band=True),
    'replace 2-D': dict(group='z0', field='z0', weighted=False, by_band=False),
    'replace 3-D': dict(group='dz', field='dz', weighted=False, by_band=False),
}


@pytest.mark.parametrize('case', CASES)
def test_tile_rms_diff_matches_brute_force(tiles, case):
    c = CASES[case]
    kwargs = dict(pad=PAD, feather=FEATHER) if c['weighted'] else dict(pad=None, feather=None)
    m = pc.grid.mosaic().from_list(list(tiles), group=c['group'], fields=[c['field']],
                                   by_band=c['by_band'], tile_diff_fields=[c['field']], **kwargs)
    count, rms = brute_force(tiles, c['group'], c['field'], c['weighted'])
    out_count = getattr(m, c['field']+'_tile_count')
    out_rms = getattr(m, c['field']+'_tile_rms_diff')
    assert out_count.shape == count.shape
    assert np.array_equal(out_count, count)
    assert np.allclose(out_rms, rms, equal_nan=True, rtol=1e-12, atol=1e-12)
    # some cells have two or more tiles, and the rest have no RMS
    assert np.any(count >= 2) and np.all(np.isnan(out_rms[count < 2]))
    assert m.field_dims[c['field']+'_tile_rms_diff'] == m.field_dims[c['field']]


def test_tile_rms_diff_leaves_the_mosaic_unchanged(tiles):
    kwargs = dict(group='dz', fields=['dz'], pad=PAD, feather=FEATHER, by_band=False)
    plain = pc.grid.mosaic().from_list(list(tiles), **kwargs)
    with_diff = pc.grid.mosaic().from_list(list(tiles), tile_diff_fields=['dz'], **kwargs)
    assert plain.fields == ['dz']
    assert with_diff.fields == ['dz', 'dz_tile_rms_diff', 'dz_tile_count']
    assert np.array_equal(plain.dz, with_diff.dz, equal_nan=True)


def test_tile_rms_diff_is_stable_for_large_values(tiles, tmp_path):
    # elevations of thousands of meters must give the same differences
    offset_tiles = make_tiles(tmp_path, offset=4000.)
    kwargs = dict(group='z0', fields=['z0'], pad=PAD, feather=FEATHER, by_band=False, tile_diff_fields=['z0'])
    a = pc.grid.mosaic().from_list(list(tiles), **kwargs)
    b = pc.grid.mosaic().from_list(list(offset_tiles), **kwargs)
    assert np.allclose(a.z0_tile_rms_diff, b.z0_tile_rms_diff, equal_nan=True, rtol=0, atol=1e-10)


def test_tile_rms_diff_needs_a_mosaic_field(tiles):
    with pytest.raises(ValueError, match='not one of the mosaic fields'):
        pc.grid.mosaic().from_list(list(tiles), group='z0', fields=['z0'], pad=PAD, feather=FEATHER,
                                   tile_diff_fields=['cell_area'])


def test_make_mosaic_writes_tile_rms_diff(tiles, tmp_path, monkeypatch):
    from pointCollection.scripts import make_mosaic
    out = str(tmp_path / 'mosaic.h5')
    monkeypatch.setattr(sys, 'argv', ['make_mosaic.py', '-d', os.path.dirname(tiles[0]), '-g', 'E*.h5',
                                      '-w', '-p', str(PAD), '-f', str(FEATHER), '--in_group', 'dz/',
                                      '-F', 'dz', '--tile_rms_diff', 'dz', '-O', out, '-R'])
    make_mosaic.main()
    written = pc.grid.data().from_h5(out, group='dz', fields=['dz', 'dz_tile_rms_diff', 'dz_tile_count'])
    count, rms = brute_force(tiles, 'dz', 'dz', True)
    assert np.array_equal(written.dz_tile_count, count)
    assert np.allclose(written.dz_tile_rms_diff, rms, equal_nan=True, rtol=1e-12, atol=1e-12)
