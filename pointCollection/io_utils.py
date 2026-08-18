#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cloud-aware I/O helpers, letting geoIndex and data.from_h5() read HDF5 files
from either local disk or a remote (e.g. s3://) location.
"""
import re

_S3FS_CACHE = {}

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
    Return a cached, authenticated s3fs.S3FileSystem for the given DAAC,
    created via earthaccess.get_s3fs_session(). Sessions are cached by
    (daac, kwargs) so repeated calls don't re-derive credentials.
    """
    key = (daac, tuple(sorted(kwargs.items())))
    if key not in _S3FS_CACHE:
        import earthaccess
        _S3FS_CACHE[key] = earthaccess.get_s3fs_session(daac=daac, **kwargs)
    return _S3FS_CACHE[key]

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

def open_h5(filename, mode='r', fs=None):
    """
    Open an HDF5 file for reading, whether local or remote.

    Parameters
    ----------
    filename : str
    mode : str, optional
    fs : s3fs.S3FileSystem, optional
        filesystem to use for a remote open. If None, a cached session is
        obtained via get_s3fs().

    Returns
    -------
    h5py.File
    """
    import h5py
    if is_remote_path(filename):
        return h5py.File((fs or get_s3fs()).open(filename, 'rb'), mode)
    return h5py.File(filename, mode)
