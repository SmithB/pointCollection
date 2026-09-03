"""
Tests for pointCollection.grid.data raster reading
"""
import numpy as np
import pytest
import pointCollection as pc

gdal = pytest.importorskip('osgeo.gdal')
osr = pytest.importorskip('osgeo.osr')


def _make_geotif(path, n_bands=3, with_time=True, t0=2000.):
    """
    write a small multi-band geotiff

    Band n is filled with the value n, so the band a read returns is
    identifiable from its data, and (optionally) tagged with a 'time'
    metadata value of t0 + n - 1.
    """
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(str(path), 8, 6, n_bands, gdal.GDT_Float32)
    ds.SetGeoTransform((0., 100., 0., 600., 0., -100.))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(3031)
    ds.SetProjection(srs.ExportToWkt())
    for ii in range(n_bands):
        band = ds.GetRasterBand(ii + 1)
        band.WriteArray(np.full((6, 8), float(ii + 1), dtype=np.float32))
        if with_time:
            band.SetMetadata({'time': str(t0 + ii)})
    ds.FlushCache()
    ds = None
    return str(path)


@pytest.fixture
def timed_tif(tmp_path):
    """a 3-band geotiff whose bands carry 'time' metadata"""
    return _make_geotif(tmp_path / 'timed.tif', with_time=True)


@pytest.fixture
def plain_tif(tmp_path):
    """a 3-band geotiff with no time metadata"""
    return _make_geotif(tmp_path / 'plain.tif', with_time=False)


# ---------------------------------------------------------------------------
# from_gdal(bands=...).  Band numbers are 1-based (the GDAL convention), and
# are used both to index the time array (t[bands-1]) and to iterate over when
# reading. The t[bands-1] subtraction, added with t_range support, only works
# on an array, so a list or tuple of bands raised TypeError for any raster
# carrying 'time' band metadata -- and from_geotif() swallows exceptions from
# from_gdal() unless verbose=True, so the read silently returned an object
# with no 'z'. Every in-repo caller (geoIndex's DEM/geotif branches) passes a
# list. Coercing with np.atleast_1d().astype(int) covers a scalar, list,
# tuple or array; np.array() alone would leave a scalar as a 0-d array, which
# then fails to iterate.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('bands,expected_values,expected_t', [
    ([1], [1.], [2000.]),
    ([1, 3], [1., 3.], [2000., 2002.]),
    ((1, 2), [1., 2.], [2000., 2001.]),
    (1, [1.], [2000.]),
    (np.array([1, 2]), [1., 2.], [2000., 2001.]),
])
def test_from_geotif_band_forms(timed_tif, bands, expected_values, expected_t):
    D = pc.grid.data().from_geotif(timed_tif, bands=bands)
    assert hasattr(D, 'z'), 'from_gdal failed and from_geotif swallowed it'
    np.testing.assert_array_equal(np.unique(D.z), expected_values)
    np.testing.assert_array_equal(np.atleast_1d(D.t), expected_t)
    # a single band is returned 2-d, several are stacked on the third axis
    n_bands = len(expected_values)
    assert D.z.shape == ((6, 8) if n_bands == 1 else (6, 8, n_bands))


@pytest.mark.parametrize('bands', [[1], [1, 3], (1, 2), 1, np.array([1, 2])])
def test_from_geotif_bands_without_time_metadata(plain_tif, bands):
    # the same band forms must work when there is no time metadata, which is
    # the path that skips the t[bands-1] indexing entirely
    D = pc.grid.data().from_geotif(plain_tif, bands=bands)
    n_bands = len(np.atleast_1d(bands))
    np.testing.assert_array_equal(np.unique(D.z), np.atleast_1d(bands).astype(float))
    assert D.z.shape == ((6, 8) if n_bands == 1 else (6, 8, n_bands))


def test_from_geotif_all_bands(timed_tif):
    D = pc.grid.data().from_geotif(timed_tif)
    assert D.z.shape == (6, 8, 3)
    np.testing.assert_array_equal(np.unique(D.z), [1., 2., 3.])
    np.testing.assert_array_equal(D.t, [2000., 2001., 2002.])


def test_from_geotif_t_range_selects_bands(timed_tif):
    # t_range picks the bands itself, as 1-based numbers, then indexes t with
    # them the same way
    D = pc.grid.data().from_geotif(timed_tif, t_range=[2000., 2001.])
    assert D.z.shape == (6, 8, 2)
    np.testing.assert_array_equal(np.unique(D.z), [1., 2.])
    np.testing.assert_array_equal(D.t, [2000., 2001.])
