#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
I/O helpers shared by geoIndex, tilingSchema and data.from_h5(): cloud-aware
access, letting HDF5 files be read from either local disk or a remote (e.g.
s3://) location, and normalization of the data-format names used to tag the
sources in an index or a tiling schema.
"""
import re

# key -> (filesystem, expires_at or None).  Brokered DAAC credentials expire
# (MAAP/NSIDC issue roughly four hours), so the session cannot simply be cached
# for the life of the process: a long read phase would start 403-ing partway
# through, which at a per-tile fan-out looks like a random data error rather
# than an expiry.  A session with no expiry (daac=None, the worker's own
# identity) caches as before.
_S3FS_CACHE = {}

# Re-derive this many seconds BEFORE the stated expiry, so a read that starts
# just under the wire does not expire mid-request.
_CREDENTIAL_SAFETY_MARGIN = 600

# DAAC -> the DAAC's own S3-credentials endpoint, for the MAAP path in
# get_s3fs().  MAAP brokers these: maap.aws.earthdata_s3_credentials(<uri>)
# returns short-lived accessKeyId/secretAccessKey/sessionToken using MAAP's own
# auth, so a MAAP DPS worker needs no Earthdata credentials of its own -- no
# .netrc, nothing at rest.  A DAAC that is not listed simply falls through to
# earthaccess, which is still the right answer off MAAP.
MAAP_S3_CREDENTIALS_ENDPOINTS = {
    'NSIDC': 'https://data.nsidc.earthdatacloud.nasa.gov/s3credentials',
}

# Block size for remote reads that pull windows out of a large file.  fsspec's
# 5 MiB default is sized for reading a file end to end; a windowed read of a
# chunked HDF5 file touches scattered chunks, and the read-ahead is then mostly
# waste -- a 60 km window out of ATL14 fetches 35 MiB in 5 MiB blocks and
# 6.6 MiB in 256 KiB ones.
DEFAULT_REMOTE_BLOCK_SIZE = 256 * 1024

# pc.indexedH5 is the class that reads and writes this format, so 'indexedH5'
# is its canonical name.  geoIndex files written before that spelling was
# settled on, and calling code following the geoIndex file_type convention,
# spell the same format 'indexed_h5' (and 'indexedh5' turns up by hand), so
# accept any of them wherever a file_type / data_format is given.  Keys are
# lowercased with underscores removed; values are the canonical spelling.
FILE_TYPE_ALIASES = {'indexedh5': 'indexedH5'}


def canonical_file_type(file_type):
    """
    Map alternate spellings of a data-format name onto the canonical one.

    Parameters
    ----------
    file_type : str, bytes, or None
        Format name, as passed to geoIndex.for_file() or stored in a
        geoIndex 'type_N' attribute.  bytes (as older HDF5 files may
        return) are decoded to str.

    Returns
    -------
    str or None
        The canonical name if `file_type` is a recognized alias (e.g.
        'indexed_h5' -> 'indexedH5'), otherwise `file_type` unchanged.
        Names that are not aliases keep their case, so unrelated types
        ('ATL06', 'indexed_h5_from_matlab', ...) pass through untouched.
    """
    if isinstance(file_type, bytes):
        file_type = file_type.decode('utf-8')
    if not isinstance(file_type, str):
        return file_type
    return FILE_TYPE_ALIASES.get(file_type.lower().replace('_', ''), file_type)

def is_remote_path(filename):
    """
    True if filename looks like a URI (e.g. s3://..., https://...) rather
    than a local path.
    """
    return isinstance(filename, str) and re.match(r'^[a-zA-Z][a-zA-Z0-9+.\-]*://', filename) is not None

def strip_pair_suffix(filename):
    """
    Remove a trailing ':pairN' suffix (used by geoIndex for ATL06/ATL11
    per-beam-pair entries) from a filename, if present.
    """
    return None if filename is None else re.sub(r':pair\d+$', '', filename)

def get_s3fs(daac='NSIDC', **kwargs):
    """
    Return a cached s3fs.S3FileSystem.

    Parameters
    ----------
    daac : str or None, default 'NSIDC'
        If a DAAC name, the session carries the short-lived in-region
        credentials that DAAC's cloud buckets require.  On MAAP those come
        from MAAP's own credential broker (see _s3fs_from_maap), which needs
        no Earthdata credentials at all; everywhere else, and whenever the
        broker declines, they come from earthaccess.get_s3fs_session().
        If None, an ordinary s3fs.S3FileSystem() is returned, which picks up
        whatever the default AWS credential chain provides (environment,
        ~/.aws, or an instance/task role).  That is the right choice for
        buckets we own rather than read from a DAAC -- ancillary rasters,
        masks and tiling schemas on s3://maap-ops-workspace/... -- since
        earthaccess credentials do not grant access to them.

    Sessions are cached by (daac, kwargs), so repeated calls are a dict lookup
    -- but a cached session whose brokered credentials are within
    _CREDENTIAL_SAFETY_MARGIN of expiring is discarded and re-derived.  Callers
    doing a long sequence of reads should therefore call this per item rather
    than hoisting one filesystem out of the loop, which costs a dict lookup and
    is what makes the refresh actually take effect.
    """
    import time

    key = (daac, tuple(sorted(kwargs.items())))
    entry = _S3FS_CACHE.get(key)
    if entry is not None:
        fs, expires_at = entry
        if expires_at is None or time.time() < expires_at - _CREDENTIAL_SAFETY_MARGIN:
            return fs
        # Otherwise fall through and re-derive: the credentials this session
        # carries are about to stop working.

    if daac is None:
        import s3fs
        # The default credential chain refreshes itself; nothing to expire here.
        _S3FS_CACHE[key] = (s3fs.S3FileSystem(**kwargs), None)
    else:
        fs, expires_at = _s3fs_from_maap(daac, **kwargs)
        if fs is None:
            import earthaccess
            # earthaccess manages its own session; we do not know its expiry.
            fs, expires_at = earthaccess.get_s3fs_session(daac=daac, **kwargs), None
        _S3FS_CACHE[key] = (fs, expires_at)
    return _S3FS_CACHE[key][0]


def _s3fs_from_maap(daac, **kwargs):
    """
    Build an s3fs session from MAAP-brokered DAAC credentials.

    Returns (filesystem, expires_at) where expires_at is a POSIX timestamp, or
    (None, None) -- with a warning saying why -- whenever this is not a MAAP
    environment or the broker will not answer, so the caller falls back to
    earthaccess.  Off MAAP this costs one dict lookup and returns (None, None).

    This exists because a MAAP DPS worker has NO Earthdata credentials: it runs
    as root with no ~/.netrc, and earthaccess's netrc and environment
    strategies both come up empty there.  What it does have is MAAP's own auth
    ($MAAP_PGT, and a maap_token from config), which
    maap.aws.earthdata_s3_credentials() exchanges for the DAAC's temporary
    read credentials.  See docs.maap-project.org, science/NISAR/NISAR_access.html.

    The credentials are short-lived -- MAAP/NSIDC issue roughly four hours --
    so the response's `expiration` is parsed and handed back for get_s3fs() to
    cache against.  An Antarctic tile whose ATL11 read phase runs long would
    otherwise start 403-ing partway through, and at a per-tile fan-out that
    reads as a random data error rather than as an expiry.

    Every failure warns rather than passing silently: a session that quietly
    came from somewhere other than where you think is exactly the kind of
    failure that only shows up later, as a permission error with no obvious
    cause.
    """
    import os
    import warnings

    endpoint = MAAP_S3_CREDENTIALS_ENDPOINTS.get(str(daac).upper())
    if endpoint is None:
        return None, None
    if not os.environ.get('MAAP_PGT'):
        # The ADE and DPS workers both set it; its absence means this is not a
        # MAAP environment, which is not worth warning about.
        return None, None

    try:
        from maap.maap import MAAP
    except ImportError as exc:
        warnings.warn(f'MAAP_PGT is set but maap-py is not importable ({exc}); '
                      f'falling back to earthaccess for {daac}.')
        return None, None

    try:
        creds = MAAP(
            maap_host=os.environ.get('MAAP_API_HOST', 'api.maap-project.org')
        ).aws.earthdata_s3_credentials(endpoint)
        return _s3fs_with_credentials(creds, **kwargs), _expiry_timestamp(creds)
    except Exception as exc:
        warnings.warn(f'MAAP could not broker {daac} credentials from {endpoint} '
                      f'({type(exc).__name__}: {exc}); falling back to earthaccess.')
        return None, None


def try_earthaccess_login():
    """
    Log in to Earthdata if we can, and carry on if we cannot.

    A CMR metadata search needs no authentication -- only granule READS do, and
    those get their credentials from get_s3fs(), which on MAAP uses MAAP's
    broker rather than earthaccess.  Callers used to write
    earthaccess.login(strategy='netrc'), which hard-codes the ONE strategy a
    MAAP DPS worker cannot satisfy: it runs as root with no ~/.netrc, so the
    search raised LoginStrategyUnavailable before ever reaching CMR.

    Bare login() tries environment, then netrc, then interactive, so a local
    user's existing setup still works.  Failure warns rather than raising,
    because the caller very likely does not need it.
    """
    import warnings
    import earthaccess
    try:
        earthaccess.login()
    except Exception as exc:
        warnings.warn(f'earthaccess.login() failed ({type(exc).__name__}: {exc}); '
                      'continuing, since a CMR search needs no credentials.  '
                      'Granule reads get their credentials separately, via '
                      'pointCollection.io_utils.get_s3fs().')


def _expiry_timestamp(creds):
    """
    POSIX timestamp for a credential response's `expiration`, or None.

    The field arrives as e.g. '2026-09-08 22:19:41+00:00'.  An unparseable or
    absent value is treated as a SHORT lifetime rather than an unlimited one:
    guessing "no expiry" from a value we failed to read is how a session
    outlives its credentials.
    """
    import datetime
    import warnings

    raw = creds.get('expiration')
    if raw is None:
        warnings.warn('credential response carried no expiration; '
                      're-deriving conservatively.')
        return _conservative_expiry()
    try:
        when = datetime.datetime.fromisoformat(str(raw))
    except ValueError:
        warnings.warn(f'could not parse credential expiration {raw!r}; '
                      're-deriving conservatively.')
        return _conservative_expiry()
    if when.tzinfo is None:
        when = when.replace(tzinfo=datetime.timezone.utc)
    return when.timestamp()


def _conservative_expiry():
    """A short fallback lifetime for a credential response we could not read."""
    import time
    return time.time() + 1800 + _CREDENTIAL_SAFETY_MARGIN


def _s3fs_with_credentials(creds, **kwargs):
    """s3fs session from an earthdata_s3_credentials() response."""
    import s3fs
    missing = [k for k in ('accessKeyId', 'secretAccessKey', 'sessionToken')
               if k not in creds]
    if missing:
        raise KeyError(f'credential response is missing {missing}')
    return s3fs.S3FileSystem(anon=False,
                             key=creds['accessKeyId'],
                             secret=creds['secretAccessKey'],
                             token=creds['sessionToken'],
                             **kwargs)

def open_remote(filename, mode='rb', fs=None, block_size=None, daac='NSIDC'):
    """
    Open a remote (e.g. s3://) file as a file object.

    Parameters
    ----------
    filename : str
    mode : str, default 'rb'
    fs : s3fs.S3FileSystem or NoneType, default None
        Filesystem to open with.  If None, a cached session is obtained via
        get_s3fs(daac=daac).
    block_size : int or NoneType, default None
        Bytes fetched per range request.  None leaves the filesystem's own
        default (5 MiB for s3fs) in place; see DEFAULT_REMOTE_BLOCK_SIZE for
        why a windowed read wants a smaller one.  Passed per file rather than
        to the session, so it applies to a caller-supplied fs too -- including
        an earthaccess DAAC session, whose constructor takes no such argument.
    daac : str or NoneType, default 'NSIDC'
        DAAC whose credentials are needed, if fs is None.  None selects the
        default AWS credential chain; see get_s3fs().

    Returns
    -------
    file object
    """
    if fs is None:
        fs = get_s3fs(daac=daac)
    if block_size is None:
        return fs.open(filename, mode)
    try:
        return fs.open(filename, mode, block_size=block_size)
    except TypeError:
        # a filesystem (or a stand-in) whose open() takes no block_size
        return fs.open(filename, mode)

def as_gdal_path(filename):
    """
    Translate a URI into the /vsi... path GDAL uses for the same object.

    GDAL cannot open an 's3://bucket/key' URI directly, but it reads the same
    object through its /vsis3/ virtual filesystem, which uses the ordinary AWS
    credential chain (AWS_* environment variables or ~/.aws).  Local paths and
    paths that are already /vsi... are returned unchanged, so callers can pass
    everything through this on the way to gdal.Open().
    """
    if not is_remote_path(filename):
        return filename
    scheme, _, rest = filename.partition('://')
    scheme = scheme.lower()
    if scheme == 's3':
        return '/vsis3/' + rest
    if scheme == 'gs':
        return '/vsigs/' + rest
    if scheme in ('http', 'https', 'ftp'):
        return '/vsicurl/' + filename
    # unknown scheme: hand it to GDAL as-is and let GDAL report the problem
    return filename

def path_exists(filename, fs=None, assume_remote_exists=True):
    """
    Check whether a local or remote file exists.

    Parameters
    ----------
    filename : str
    fs : s3fs.S3FileSystem, optional
        filesystem to use for a remote existence check. If None, a cached
        session is obtained via get_s3fs().
    assume_remote_exists : bool, optional
        if True (default), remote paths are assumed to exist without a
        network round-trip (appropriate when the path came from a search
        that already confirmed the granule exists, e.g. earthaccess/CMR).
    """
    import os
    if is_remote_path(filename):
        if assume_remote_exists:
            return True
        return (fs or get_s3fs()).exists(filename)
    return os.path.isfile(filename)

def open_h5(filename, mode='r', fs=None, block_size=None):
    """
    Open an HDF5 file for reading, whether local or remote.

    Parameters
    ----------
    filename : str
    mode : str, optional
    fs : s3fs.S3FileSystem, optional
        filesystem to use for a remote open. If None, a cached session is
        obtained via get_s3fs().
    block_size : int, optional
        bytes fetched per range request for a remote file.  None leaves the
        filesystem default; see open_remote().

    Returns
    -------
    h5py.File
    """
    import h5py
    if is_remote_path(filename):
        return h5py.File(open_remote(filename, fs=fs, block_size=block_size), mode)
    return h5py.File(filename, mode)
