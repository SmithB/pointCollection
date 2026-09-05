#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read a netCDF4 file through h5py, presenting the netCDF4 API that
grid.data.from_nc() uses.

netCDF4.Dataset accepts a filename or an in-memory buffer, but never a file
object, so a remote netCDF4 file can only be read whole.  netCDF4 files are
HDF5 underneath, though, and h5py does take a file object -- so pointing h5py
at an fsspec/s3fs file gives chunk-wise reads over range requests, and a
windowed read costs the chunks the window touches rather than the whole
granule.

The classes here wrap an h5py.File so that from_nc() cannot tell the
difference: .groups, .variables, var.shape, var.dimensions, var.ncattrs(),
var.getncattr(), attribute access (hasattr(var, '_FillValue')) and the
set_auto_*() no-ops.  Three differences between the libraries are papered
over here rather than in from_nc():

  * netCDF4 applies scale_factor/add_offset unless set_auto_scale(False) is
    called, and from_nc() only turns off masking.  h5py applies neither, so
    _scale() replicates what netCDF4 does -- otherwise packed variables come
    back with different values and no error.
  * h5py rejects negative-step slices, which select_slices() emits whenever
    the y axis descends, and it broadcasts multiple index arrays the way
    numpy does.  _plan() rewrites an index tuple into a read h5py accepts
    plus per-axis fixups applied afterwards.
  * HDF5 datasets that are netCDF dimensions without being netCDF variables
    (dimensions with no coordinate variable) are hidden, as are the HDF5
    bookkeeping attributes that netCDF4 keeps out of ncattrs().
"""
import posixpath

import numpy as np

# attributes netCDF4 does not report through ncattrs()
HIDDEN_ATTRIBUTES = {'CLASS', 'NAME', 'DIMENSION_LIST', 'REFERENCE_LIST',
                     'DIMENSION_LABELS', '_Netcdf4Dimid', '_Netcdf4Coordinates',
                     '_nc3_strict'}

# NAME attribute of an HDF5 dimension scale that is not a netCDF variable
PHONY_DIMENSION_PREFIX = 'This is a netCDF dimension but not a netCDF variable'


def _as_str(value):
    """decode an HDF5 string (bytes) the way netCDF4 hands it back"""
    if isinstance(value, bytes):
        return value.decode('utf-8', errors='replace')
    return value


def _attribute_value(value):
    """
    Convert an h5py attribute value into what netCDF4 would have returned:
    a scalar for a single-valued attribute, str for text.
    """
    if isinstance(value, np.ndarray):
        if value.size == 1:
            return _attribute_value(value.reshape(-1)[0])
        if value.dtype.kind in ('S', 'O'):
            return [_as_str(item) for item in value]
        return value
    return _as_str(value)


def _is_phony_dimension(dset):
    """True if dset is a netCDF dimension that is not also a netCDF variable"""
    name = _as_str(dset.attrs.get('NAME'))
    return isinstance(name, str) and name.startswith(PHONY_DIMENSION_PREFIX)


def _is_dimension_scale(dset):
    return _as_str(dset.attrs.get('CLASS')) == 'DIMENSION_SCALE'


def apply_scaling(data, scale=None, offset=None):
    """
    Apply scale_factor/add_offset to data the way netCDF4 does.

    Mirrors netCDF4.Variable.__getitem__, including its 1.0/0.0 special
    cases, so that reading through h5py returns the same numbers.
    """
    if scale is None and offset is None:
        return data
    if scale is not None and offset is not None:
        if offset != 0.0 or scale != 1.0:
            return data * scale + offset
        return data.astype(np.asarray(scale).dtype)
    if scale is not None:
        return data * scale if scale != 1.0 else data
    return data + offset if offset != 0.0 else data


def scaled_fill_value(variable):
    """
    A variable's _FillValue on the scale of the data that comes back from it.

    netCDF4 scales packed data on the way out but reports _FillValue as it is
    stored, so comparing the two directly never matches: a packed variable's
    fill values survive into the output unconverted.  Putting the fill value
    through the same arithmetic as the data makes the comparison meaningful.
    Works for a netCDF4.Variable as well as for an H5Variable.

    Parameters
    ----------
    variable : netCDF4.Variable or H5Variable

    Returns
    -------
    The fill value, scaled if the variable is packed, or None if the variable
    has no _FillValue.
    """
    attributes = variable.ncattrs()
    if '_FillValue' not in attributes:
        return None
    fill = variable.getncattr('_FillValue')
    scale = variable.getncattr('scale_factor') if 'scale_factor' in attributes else None
    offset = variable.getncattr('add_offset') if 'add_offset' in attributes else None
    if scale is None and offset is None:
        return fill
    # as stored, so that the arithmetic matches what the data went through
    return apply_scaling(np.asarray(fill, dtype=variable.dtype), scale, offset)[()]


class H5Variable:
    """an h5py.Dataset presented as a netCDF4.Variable"""

    def __init__(self, dset):
        self._dset = dset

    # -- identity ---------------------------------------------------------
    @property
    def name(self):
        return posixpath.basename(self._dset.name)

    @property
    def shape(self):
        return self._dset.shape

    @property
    def dtype(self):
        return self._dset.dtype

    @property
    def ndim(self):
        return self._dset.ndim

    def __len__(self):
        return len(self._dset)

    @property
    def dimensions(self):
        """
        Names of the variable's dimensions, as netCDF4 reports them.

        In HDF5 these are the dimension scales attached to each axis.  A
        coordinate variable is its own dimension scale, and a scale is not
        attached to itself, so that case falls back to the dataset's own name.
        """
        dset = self._dset
        names = []
        for axis, dim in enumerate(dset.dims):
            if len(dim):
                names.append(posixpath.basename(dim[0].name))
            elif dim.label:
                names.append(posixpath.basename(_as_str(dim.label)))
            elif dset.ndim == 1 and _is_dimension_scale(dset):
                names.append(self.name)
            else:
                names.append(f'phony_dim_{axis}')
        return tuple(names)

    # -- attributes -------------------------------------------------------
    def ncattrs(self):
        return [name for name in self._dset.attrs
                if name not in HIDDEN_ATTRIBUTES]

    def getncattr(self, name):
        if name in HIDDEN_ATTRIBUTES:
            raise AttributeError(name)
        try:
            return _attribute_value(self._dset.attrs[name])
        except KeyError:
            raise AttributeError(name)

    def __getattr__(self, name):
        # reached only when normal lookup fails, so this is the netCDF4
        # spelling of attribute access: var._FillValue, var.grid_mapping
        if name.startswith('_') and name.endswith('__'):
            raise AttributeError(name)
        try:
            dset = self.__dict__['_dset']
        except KeyError:
            raise AttributeError(name)
        if name in HIDDEN_ATTRIBUTES or name not in dset.attrs:
            raise AttributeError(name)
        return _attribute_value(dset.attrs[name])

    # -- reading ----------------------------------------------------------
    def _plan(self, key):
        """
        Rewrite an index tuple into one h5py accepts, plus the per-axis
        fixups to apply to the array that comes back.

        h5py takes only forward slices, and index arrays on more than one
        axis would broadcast rather than select an outer product, so index
        arrays are read as the range they span and selected afterwards.
        """
        if any(index is Ellipsis for index in key):
            # expand to explicit slices so each entry maps to one axis
            where = [i for i, index in enumerate(key) if index is Ellipsis][0]
            fill = (slice(None),) * (self._dset.ndim - (len(key) - 1))
            key = key[:where] + fill + key[where + 1:]
        read = []
        fixups = []
        for index in key:
            if isinstance(index, (int, np.integer)):
                # an integer index drops the axis, so it needs no fixup
                read.append(int(index))
                continue
            size = self._dset.shape[len(read)] if len(read) < self._dset.ndim else None
            if isinstance(index, slice):
                step = 1 if index.step is None else index.step
                if step > 0:
                    read.append(index)
                    fixups.append(None)
                    continue
                # negative step: read forwards, reverse afterwards
                positions = range(*index.indices(size))
                if len(positions) == 0:
                    read.append(slice(0, 0))
                    fixups.append(None)
                    continue
                read.append(slice(positions[-1], positions[0] + 1, -step))
                fixups.append(slice(None, None, -1))
                continue
            # a sequence of indices (bands, typically)
            indices = np.asarray(index)
            if indices.dtype == bool:
                indices = np.flatnonzero(indices)
            indices = np.asarray(indices, dtype=np.int64)
            indices = np.where(indices < 0, indices + size, indices)
            if indices.size == 0:
                read.append(slice(0, 0))
                fixups.append(None)
                continue
            start, stop = int(indices.min()), int(indices.max()) + 1
            steps = np.unique(np.diff(indices)) if indices.size > 1 else np.array([1])
            if steps.size == 1 and steps[0] > 0:
                # evenly spaced and increasing: a slice reads it exactly
                read.append(slice(start, stop, int(steps[0])))
                fixups.append(None)
            else:
                # read the span the indices cover, then select within it
                read.append(slice(start, stop))
                fixups.append(indices - start)
        return tuple(read), fixups

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            key = (key,)
        if self._dset.ndim == 0:
            # netCDF4 lets a scalar variable be indexed with [:], [...] or
            # [()], and hands back a 0-d array; h5py accepts only the last two
            return self._scale(np.asarray(self._dset[()]))
        read, fixups = self._plan(key)
        data = self._dset[read]
        # one axis at a time: a tuple of index arrays would broadcast
        for axis, fixup in enumerate(fixups):
            if fixup is None:
                continue
            data = data[(slice(None),) * axis + (fixup,)]
        return self._scale(data)

    def _scale(self, data):
        """
        Apply scale_factor/add_offset the way netCDF4 does.

        from_nc() calls set_auto_mask(False) but not set_auto_scale(False),
        so netCDF4 scales packed data on the way out.
        """
        attrs = self._dset.attrs
        return apply_scaling(
            data,
            _attribute_value(attrs['scale_factor']) if 'scale_factor' in attrs else None,
            _attribute_value(attrs['add_offset']) if 'add_offset' in attrs else None)


class _Mapping:
    """
    dict-like view of the datasets or groups in an h5py group.

    Group lookup is more forgiving than netCDF4's: a path with separators or
    a leading '/' resolves too, since from_nc() takes the group name from its
    caller.
    """
    def __init__(self, group, wrap, want_group):
        self._group = group
        self._wrap = wrap
        self._want_group = want_group

    def _members(self):
        import h5py
        members = {}
        for name, item in self._group.items():
            if self._want_group:
                if isinstance(item, h5py.Group):
                    members[name] = item
            elif isinstance(item, h5py.Dataset) and not _is_phony_dimension(item):
                members[name] = item
        return members

    def __getitem__(self, name):
        members = self._members()
        if name in members:
            return self._wrap(members[name])
        if self._want_group:
            # nested or '/'-prefixed group path
            key = name.strip('/')
            if key in ('', '.'):
                return self._wrap(self._group)
            try:
                item = self._group[key]
            except KeyError:
                raise KeyError(name)
            return self._wrap(item)
        raise KeyError(name)

    def __contains__(self, name):
        try:
            self[name]
        except KeyError:
            return False
        return True

    def keys(self):
        return self._members().keys()

    def values(self):
        return [self._wrap(item) for item in self._members().values()]

    def items(self):
        return [(name, self._wrap(item)) for name, item in self._members().items()]

    def __iter__(self):
        return iter(self._members())

    def __len__(self):
        return len(self._members())


class H5Group:
    """an h5py.Group presented as a netCDF4.Group"""

    def __init__(self, group):
        self._group = group

    @property
    def variables(self):
        return _Mapping(self._group, H5Variable, want_group=False)

    @property
    def groups(self):
        return _Mapping(self._group, H5Group, want_group=True)

    def __getitem__(self, name):
        import h5py
        item = self._group[name]
        return H5Group(item) if isinstance(item, h5py.Group) else H5Variable(item)

    def __contains__(self, name):
        return name in self._group

    def ncattrs(self):
        return [name for name in self._group.attrs
                if name not in HIDDEN_ATTRIBUTES]

    def getncattr(self, name):
        try:
            return _attribute_value(self._group.attrs[name])
        except KeyError:
            raise AttributeError(name)

    def __getattr__(self, name):
        if name.startswith('_') and name.endswith('__'):
            raise AttributeError(name)
        try:
            group = self.__dict__['_group']
        except KeyError:
            raise AttributeError(name)
        if name in HIDDEN_ATTRIBUTES or name not in group.attrs:
            raise AttributeError(name)
        return _attribute_value(group.attrs[name])

    # netCDF4's masking and scaling switches.  Masking is off either way in
    # h5py; scaling is always applied, matching netCDF4's default.
    def set_auto_mask(self, value):
        pass

    def set_auto_scale(self, value):
        if not value:
            raise NotImplementedError(
                'the h5py netCDF4 reader always applies scale_factor/add_offset')

    def set_auto_maskandscale(self, value):
        self.set_auto_scale(value)

    def set_always_mask(self, value):
        pass


class H5Dataset(H5Group):
    """
    an h5py.File presented as a netCDF4.Dataset.

    Closing this closes the byte source underneath as well: closing an
    h5py.File does not close the fsspec file object it was opened from.
    """

    def __init__(self, h5f, source=None):
        super().__init__(h5f)
        self._h5f = h5f
        self._source = source

    @property
    def filepath(self):
        return self._h5f.filename

    def close(self):
        try:
            self._h5f.close()
        finally:
            if self._source is not None and hasattr(self._source, 'close'):
                self._source.close()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


def open_nc_as_h5(source, mode='r', rdcc_nbytes=None, rdcc_nslots=None):
    """
    Open a netCDF4 file with h5py and present it as a netCDF4.Dataset.

    Parameters
    ----------
    source : str or file-like
        Local path, or an open file object (e.g. from s3fs) to read through
        range requests.
    mode : str, default 'r'
    rdcc_nbytes : int or NoneType
        Size of the HDF5 chunk cache, in bytes.  Worth raising above the 1 MiB
        default for files with chunks larger than that -- ATL15's delta_h
        chunks are (8, 686, 386) float32, 8.5 MiB apiece.
    rdcc_nslots : int or NoneType
        Number of chunk slots in that cache.

    Returns
    -------
    H5Dataset
    """
    import h5py
    kwargs = {}
    if rdcc_nbytes is not None:
        kwargs['rdcc_nbytes'] = rdcc_nbytes
        # HDF5 wants roughly 10 slots per cacheable chunk, and a prime
        kwargs['rdcc_nslots'] = rdcc_nslots if rdcc_nslots is not None else 5003
    elif rdcc_nslots is not None:
        kwargs['rdcc_nslots'] = rdcc_nslots
    h5f = h5py.File(source, mode=mode, **kwargs)
    return H5Dataset(h5f, source=None if isinstance(source, str) else source)
