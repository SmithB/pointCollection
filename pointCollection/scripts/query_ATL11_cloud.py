#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Query ATL11 granules directly from NASA Earthdata Cloud (S3) for a region
of interest, using a set of pre-built per-granule geoIndex files to locate
the rows to read.

Workflow:
    1. (offline, on another machine) build one geoIndex .h5 file per ATL11
       granule across the archive:
           index_glob.py --type ATL11 -g <granule.h5> --index_file <name> --Relative
       and stage them under a shared, FUSE-mounted bucket, one subdirectory
       per (cycles, release, version) combination -- see
       index_path_for_granule() for the exact layout.
    2. search Earthdata Cloud for ATL11 granules intersecting a bounding box
    3. for each found granule, open its pre-built geoIndex (a normal local
       read) to get the row-offset ranges for the ROI, then read exactly
       those rows from the granule's real S3 location.

@author: ben
"""

import os
import re
import warnings
import argparse
import pointCollection as pc

_ATL11_RE = re.compile(
    r'ATL11_(?P<trackregion>\d{6})_(?P<cycles>\d{4})_(?P<release>\d{3})_(?P<version>\d{2})\.h5$')


def release_spec(release_string):
    """
    Parse a release string into the (cycles, release, version) triple that
    identifies one ATL11-family product generation.

    ATL1415 spells a generation as e.g. '007_cycle_03_31_v04', which is a
    DIRECTORY-NAMING convention on the processing side, not a CMR selector:
    make_ATL11_index.py takes only split('_')[0] from it.  The same three
    numbers appear in a granule name (ATL11_000103_0331_007_04.h5) and in the
    per-granule index subtree (ATL11_index_0331_007_04), which is what makes
    the triple the right thing to filter on.

    THE ATL11 AND ATL11XO GENERATIONS ARE NORMALLY DIFFERENT -- rel_006_0331
    pairs --ATL11_release=007_cycle_03_31_v04 with
    --ATL11xo_version=007_cycle_03_30_v03, i.e. different cycle ranges AND
    different versions under the same release 007.  So a caller must pass the
    spec for the product it is actually searching; there is deliberately no
    module-level "current release" to fall back on, because a single global
    would silently filter an ATL11XO search by ATL11's cycles.

    Parameters
    ----------
    release_string : str
        '007_cycle_03_31_v04', optionally prefixed 'ATL11_'; or an already
        parsed '0331_007_04'.

    Returns
    -------
    (cycles, release, version) : tuple of str, e.g. ('0331', '007', '04')
    """
    text = str(release_string)
    if text.startswith('ATL11_'):
        text = text[len('ATL11_'):]

    m = re.match(r'^(?P<release>\d{3})_cycle_(?P<c0>\d{2})_(?P<c1>\d{2})_v(?P<version>\d{2})$',
                 text)
    if m is not None:
        return (m['c0'] + m['c1'], m['release'], m['version'])

    m = re.match(r'^(?P<cycles>\d{4})_(?P<release>\d{3})_(?P<version>\d{2})$', text)
    if m is not None:
        return (m['cycles'], m['release'], m['version'])

    raise ValueError(
        f'cannot parse release string {release_string!r}; expected '
        "'<release>_cycle_<c0>_<c1>_v<vv>' (e.g. 007_cycle_03_31_v04) or "
        "'<cycles>_<release>_<version>' (e.g. 0331_007_04)")


def granule_matches_release(granule_basename, spec):
    """True if an ATL11-family granule name carries the given release triple."""
    m = _ATL11_RE.search(granule_basename)
    if m is None:
        return False
    return (m['cycles'], m['release'], m['version']) == tuple(spec)


def find_ATL11_granules(bounding_box, short_name='ATL11', granule_release=None,
                        **search_kwargs):
    """
    Search Earthdata Cloud for ATL11 granules intersecting bounding_box.

    Parameters
    ----------
    bounding_box : iterable
        (lon_min, lat_min, lon_max, lat_max), as expected by
        earthaccess.search_data().
    short_name : str, optional
    granule_release : str or tuple, optional
        Restrict results to one product generation, as a release string
        ('007_cycle_03_31_v04') or an already parsed (cycles, release,
        version) triple.  CMR's own `version` filter selects the COLLECTION
        version (007) and cannot separate cycle ranges or granule versions
        within it, so the filtering is done on the granule names.

        Left None the search returns every generation CMR holds for
        short_name.  That is only safe while exactly one generation has a
        staged per-granule index -- a second staged index subtree would make
        every tile read both and double-count its data.  Pass this whenever
        the caller knows which generation it indexed.
    **search_kwargs :
        additional keyword arguments passed to earthaccess.search_data()
        (e.g. version, cycle, temporal).

    Returns
    -------
    list of earthaccess.DataGranule
    """
    import earthaccess
    _try_earthaccess_login()
    granules = earthaccess.search_data(short_name=short_name,
                                       bounding_box=bounding_box, **search_kwargs)
    if granule_release is None:
        return granules

    spec = (release_spec(granule_release) if isinstance(granule_release, str)
            else tuple(granule_release))
    kept = []
    for granule in granules:
        try:
            url = granule.data_links(access='direct')[0]
        except (IndexError, KeyError):
            continue
        if granule_matches_release(os.path.basename(url), spec):
            kept.append(granule)
    return kept


def _try_earthaccess_login():
    """
    Best-effort Earthdata login; see pc.io_utils.try_earthaccess_login().

    Kept as a name here because callers outside this module import it from
    here, but the implementation lives in io_utils with get_s3fs(), so the
    login policy and the credential policy stay in one place.
    """
    return pc.io_utils.try_earthaccess_login()


def index_path_for_granule(granule_basename, index_root):
    """
    Naming convention for a granule's pre-built geoIndex file, matching the
    layout already in use in the shared bucket: one subdirectory per
    (cycles, release, version) combination, named
    'ATL11_index_<cycles>_<release>_<version>', containing one index file
    per granule using the granule's own basename (e.g.
    '<index_root>/ATL11_index_0331_007_04/ATL11_044110_0331_007_04.h5').
    """
    m = _ATL11_RE.search(granule_basename)
    if m is None:
        raise ValueError(f'query_ATL11_cloud: cannot parse ATL11 basename {granule_basename!r}')
    subdir = f"ATL11_index_{m.group('cycles')}_{m.group('release')}_{m.group('version')}"
    return os.path.join(index_root, subdir, granule_basename)


def read_ATL11_granule_cloud_items(s3_url, index_file, xr, yr, fields=None, fs=None,
                                   index_fs=None, version_mismatch='error',
                                   missing_index='error'):
    """
    Find the rows of a cloud ATL11 granule falling within [xr, yr], using a
    pre-built per-granule geoIndex to locate them, and return them as the
    raw list of pointCollection.data objects geoIndex.query_xy_box() itself
    returns (one item per matched beam pair/offset segment) -- i.e. without
    merging them into a single object. Callers that need per-pair
    granularity (e.g. ATL1415.read_ATL11_at(), which computes per-pair
    slope-derived sigma_corr) should use this directly; callers that just
    want the concatenated data should use read_ATL11_granule_cloud().

    Parameters
    ----------
    s3_url : str
        the granule's direct-access S3 URI (from earthaccess).
    index_file : str
        path to the granule's pre-built geoIndex .h5 file (read locally).
    xr, yr : 2-element iterables
        bounding box, in the geoIndex's projected coordinates.
    fields : list or dict, optional
    fs : s3fs.S3FileSystem, optional
        reused across calls to avoid re-deriving S3 credentials per granule.
        This is the session for the GRANULE, i.e. a DAAC (earthaccess) one.
    index_fs : s3fs.S3FileSystem, optional
        session for a remote index_file.  The index is ours and lives in our
        own bucket, so it needs the default AWS credentials rather than the
        DAAC session in `fs`; if None, one is derived on demand.
    version_mismatch : {'error', 'skip'}, optional
        what to do when the granule indexed by index_file does not match
        the cloud-found granule s3_url (e.g. the index was built from an
        older release). 'error' (default) raises a ValueError; 'skip'
        issues a warning and returns None instead of reading the granule.

    Returns
    -------
    list of pointCollection.data, or None if the granule was skipped or
    had no data within [xr, yr]

    Raises
    ------
    ValueError
        if version_mismatch=='error' and the indexed and cloud-found
        granules differ, or if version_mismatch is not 'error' or 'skip'.
    """
    if version_mismatch not in ('error', 'skip'):
        raise ValueError(f"version_mismatch must be 'error' or 'skip', got {version_mismatch!r}")
    if missing_index not in ('error', 'skip'):
        raise ValueError(f"missing_index must be 'error' or 'skip', got {missing_index!r}")

    if pc.io_utils.is_remote_path(index_file) and index_fs is None:
        index_fs = pc.io_utils.get_s3fs(daac=None)
    # os.path.isfile() is False for any URI, so a remote index has to be
    # checked against the bucket -- otherwise every granule would be reported
    # missing and skipped, and the tile would come back empty rather than failing
    if not pc.io_utils.path_exists(index_file, fs=index_fs, assume_remote_exists=False):
        # A MISSING INDEX IS A STAGING ERROR, NOT A DATA CONDITION.  This used
        # to warn and skip, which silently produced a partial fit: a tile whose
        # index tree was half-uploaded (there are ~8100 files) came back with
        # fewer granules than it should have, and nothing recorded that.  With
        # find_ATL11_granules(granule_release=...) filtering the search to the
        # generation we actually indexed, every granule that reaches here is
        # one we expect to have an index for, so its absence is a real fault.
        #
        # missing_index='skip' restores the old behaviour for a caller that
        # genuinely wants best-effort reads.
        msg = (f'query_ATL11_cloud: missing geoIndex {index_file} for granule '
               f'{s3_url}.  The granule matches the requested release, so its '
               f'index should be staged; reading on would silently return a '
               f'partial fit for this tile.')
        if missing_index == 'error':
            raise FileNotFoundError(msg)
        warnings.warn(msg + '  Skipping (missing_index="skip").')
        return None

    s3_basename = os.path.basename(s3_url)
    gI = pc.geoIndex().from_file(index_file, fs=index_fs)
    indexed_name = pc.io_utils.strip_pair_suffix(gI.attrs.get('file_0'))
    if indexed_name is not None:
        indexed_basename = os.path.basename(indexed_name)
        if indexed_basename != s3_basename:
            msg = (f'query_ATL11_cloud: index {index_file} indexes {indexed_basename}, '
                   f'which does not match granule {s3_basename}')
            if version_mismatch == 'error':
                raise ValueError(msg)
            warnings.warn(msg + ', skipping')
            return None

    D = gI.query_xy_box(xr, yr, remote_file=s3_url, fs=fs, fields=fields)
    if D is None or len(D) == 0:
        return None
    return D


def read_ATL11_granule_cloud(s3_url, index_file, xr, yr, fields=None, fs=None,
                             index_fs=None, version_mismatch='error',
                             missing_index='error'):
    """
    Read the rows of a cloud ATL11 granule falling within [xr, yr], using a
    pre-built per-granule geoIndex to locate them, concatenated into a
    single pointCollection.data object. See read_ATL11_granule_cloud_items()
    for the parameters (identical) and for a version that returns the
    unmerged per-pair list instead.

    Returns
    -------
    pointCollection.data, or None if the granule was skipped
    """
    D = read_ATL11_granule_cloud_items(s3_url, index_file, xr, yr, fields=fields,
                                        fs=fs, index_fs=index_fs,
                                        version_mismatch=version_mismatch,
                                        missing_index=missing_index)
    if D is None:
        return None
    return pc.data().from_list(D)


def query_ATL11_cloud(bounding_box, index_root, xr, yr, fields=None,
                       version_mismatch='error', missing_index='error',
                       granule_release=None, verbose=False):
    """
    Search Earthdata Cloud for ATL11 granules in bounding_box, and read the
    rows within [xr, yr] from each, using pre-built per-granule geoIndex
    files under index_root (see index_path_for_granule() for the expected
    layout).

    version_mismatch and missing_index are passed through to
    read_ATL11_granule_cloud(); see its docstring.

    granule_release restricts the search to one product generation -- see
    find_ATL11_granules().  Pass it whenever more than one generation could be
    in CMR, which is normally: index_root holds one subtree per generation, and
    an unfiltered search would read every subtree that happens to be staged.

    Returns
    -------
    pointCollection.data containing the concatenated results
    """
    granules = find_ATL11_granules(bounding_box, granule_release=granule_release)
    fs = pc.io_utils.get_s3fs()
    if verbose:
        print(f'query_ATL11_cloud: found {len(granules)} candidate granules')

    results = []
    for granule in granules:
        s3_url = granule.data_links(access='direct')[0]
        basename = os.path.basename(s3_url)
        index_file = index_path_for_granule(basename, index_root)
        if verbose:
            print(f'query_ATL11_cloud: reading {basename}')
        D = read_ATL11_granule_cloud(s3_url, index_file, xr, yr, fields=fields,
                                      fs=fs, version_mismatch=version_mismatch,
                                      missing_index=missing_index)
        if D is not None:
            results.append(D)

    return pc.data().from_list(results)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bounding_box', '-b', type=float, nargs=4, required=True,
                         metavar=('lon_min', 'lat_min', 'lon_max', 'lat_max'),
                         help='search bounding box, in geographic coordinates')
    parser.add_argument('--xr', type=float, nargs=2, required=True, help='x range, in the index projected coordinates')
    parser.add_argument('--yr', type=float, nargs=2, required=True, help='y range, in the index projected coordinates')
    parser.add_argument('--index_root', '-i', type=str, required=True,
                         help="root directory containing per-(cycles,release,version) "
                              "'ATL11_index_<cycles>_<release>_<version>' subdirectories "
                              "of per-granule geoIndex files")
    parser.add_argument('--out_file', '-o', type=str, required=True, help='output h5 file')
    parser.add_argument('--version_mismatch', choices=['error', 'skip'], default='error',
                         help="what to do when a pre-built index's granule doesn't match the "
                              "cloud-found granule: 'error' (default) raises, 'skip' warns and skips")
    parser.add_argument('--verbose', '-v', action='store_true')
    args = parser.parse_args()

    D = query_ATL11_cloud(args.bounding_box, args.index_root, args.xr, args.yr,
                           version_mismatch=args.version_mismatch, verbose=args.verbose)
    D.to_h5(args.out_file)


if __name__ == '__main__':
    main()
