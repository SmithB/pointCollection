#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mosaic.py
Routines for creating a weighted mosaic from a series of tiles

UPDATE HISTORY:
    Updated 09/2026: from_list reads remote files with a block_size, and in a
        process pool (workers); tiles are still added in list order
    Updated 06/2023: calculate x and y arrays using np.arange and spacing
    updated 03/2021: change scheme for calculating weights, raised cosine as default
    Updated 03/2020: check number of dimensions of z if only a single band
    Written 03/2020
"""

import collections
import itertools
import os
import pickle
import numpy as np
from .data import data
import pointCollection as pc

# extensions whose readers (from_h5, from_nc) take a block_size for a remote file
_BLOCK_SIZE_FORMATS = ('h5', 'hdf', 'hdf5', 'nc', 'netcdf')

def _read_kwargs(item, block_size, **kwargs):
    """
    from_file keyword arguments for one item: block_size is added only for a
    remote HDF5 or netCDF file, the readers that take it (a local file has no
    blocks, and from_geotif has no such argument).
    """
    if block_size is not None and pc.io_utils.is_remote_path(item) and \
            os.path.splitext(item)[1][1:].lower() in _BLOCK_SIZE_FORMATS:
        kwargs['block_size'] = block_size
    return kwargs

def _read_item(job):
    """
    Read one file for a mosaic: (item, meta_only, kwargs) -> (grid, exception).

    Module-level so a process pool can pickle it.  An exception is returned,
    not raised, so the caller can raise it where the serial read would have
    -- inside the same try/except, with the same consequence for the tile.
    """
    item, meta_only, kwargs = job
    try:
        if meta_only:
            return pc.grid.data().from_file(item, meta_only=True, **kwargs), None
        return pc.grid.mosaic().from_file(item, **kwargs), None
    except Exception as exc:
        try:
            pickle.dumps(exc)
        except Exception:
            # an exception that cannot cross the process boundary would fail
            # the whole pool; carry its text instead
            exc = RuntimeError(f'{type(exc).__name__}: {exc}')
        return None, exc

# How a reading pool starts its workers.  NOT plain fork: by the time a mosaic
# reads remote tiles the parent has an s3fs event-loop thread (listing the
# tiles starts it), and a child forked from a multi-threaded process can
# deadlock -- s3fs itself refuses, "This class is not fork-safe".  forkserver
# forks each worker from a clean single-threaded server that has imported
# pointCollection once, so a worker still starts fast (40 GL tiles, 8
# workers: 9 s, after a one-off ~10 s server start; serial 32 s).
_START_METHOD = 'forkserver'

def _init_worker():
    """
    Drop any s3fs session a worker inherited (possible only under a fork
    start method): its event loop and connections belong to the parent.  The
    worker builds its own on first use.
    """
    pc.io_utils._S3FS_CACHE.clear()

def _ordered_reads(pool, jobs, window):
    """
    Yield _read_item(job) for each job, IN ORDER, with at most `window` reads
    in flight -- so results are consumed in list order (the summation order
    of the serial loop) and at most `window` tiles wait in memory.
    """
    jobs = iter(jobs)
    pending = collections.deque(pool.submit(_read_item, job)
                                for job in itertools.islice(jobs, window))
    while pending:
        result = pending.popleft().result()
        for job in itertools.islice(jobs, 1):
            pending.append(pool.submit(_read_item, job))
        yield result

class _TileReader:
    """
    Reads the string items of a mosaic's input list, serially or in a process
    pool, always yielding in list order.  Non-string items (grids already in
    memory) pass through untouched.

    Threads would not help: h5py holds its global lock for the whole of a
    read from a Python file object, so remote reads in threads run one at a
    time (measured: 8 threads = serial).  Processes do not share that lock.
    """
    def __init__(self, workers=1, block_size=None):
        self.workers = max(1, int(workers or 1))
        self.block_size = block_size
        self.pool = None

    def __enter__(self):
        if self.workers > 1:
            import concurrent.futures
            import multiprocessing
            method = _START_METHOD
            if method not in multiprocessing.get_all_start_methods():
                method = 'spawn'
            context = multiprocessing.get_context(method)
            if method == 'forkserver':
                # effective only before the server starts; it then lasts for
                # the life of this process, so later pools start at once
                context.set_forkserver_preload(['pointCollection'])
            self.pool = concurrent.futures.ProcessPoolExecutor(
                self.workers, mp_context=context, initializer=_init_worker)
        return self

    def __exit__(self, *exc_info):
        if self.pool is not None:
            self.pool.shutdown()
            self.pool = None

    def read(self, items, meta_only=False, **kwargs):
        """
        yield (item, grid, exception) for each item of `items`, in order;
        grid is the item itself (exception None) for a non-string item
        """
        items = list(items)
        jobs = [(item, meta_only, _read_kwargs(item, self.block_size, **kwargs))
                for item in items if isinstance(item, str)]
        if self.pool is None:
            results = map(_read_item, jobs)
        else:
            results = _ordered_reads(self.pool, jobs, 2*self.workers)
        for item in items:
            if isinstance(item, str):
                grid, exc = next(results)
                yield item, grid, exc
            else:
                yield item, item, None

class mosaic(data):
    def __init__(self, spacing=None, **kwargs):
        #self.x=None
        #self.y=None
        #self.t=None
        #self.z=None
        super().__init__(**kwargs)
        self.invalid=None
        self.weight=None
        self.extent=[np.inf,-np.inf,np.inf,-np.inf]
        self.dimensions=[None,None,None]
        self.field_dims={}
        # copy spacing so that mosaics never share a spacing list
        self.spacing=[None, None] if spacing is None else list(spacing)
        self.tile_weight=None
        self.fill_value=np.nan
        self.normalized=True
        self.use_time=False
        self.fields=[]

    def from_grid(self, source, copy=False):
        """make a mosaic object from a pc.data.grid object"""

        for field in ['x','y','projection','filename','extent','time', 't', 't_axis']:
            if hasattr(source, field):
                setattr(self, field, getattr(source, field))
        for field in source.fields:
            if copy:
                self.assign({field:getattr(source, field).copy()})
            else:
                self.assign({field:getattr(source, field)})
        self.__update_size_and_shape__()
        self.__update_extent__()
        return self

    def update_spacing(self, temp):
        """
        update the step size of mosaic
        """
        # try automatically getting spacing of tile
        try:
            dx = temp.x[1] - temp.x[0]
            if not dx == 0:
                self.spacing[0] = dx
            dy = temp.y[1] - temp.y[0]
            if not dy == 0:
                self.spacing[1] = dy
        except:
            pass
        return self

    def update_bounds(self, temp):
        """
        update the bounds of mosaic
        """
        if (temp.extent[0] < self.extent[0]):
            self.extent[0] = np.copy(temp.extent[0])
        if (temp.extent[1] > self.extent[1]):
            self.extent[1] = np.copy(temp.extent[1])
        if (temp.extent[2] < self.extent[2]):
            self.extent[2] = np.copy(temp.extent[2])
        if (temp.extent[3] > self.extent[3]):
            self.extent[3] = np.copy(temp.extent[3])
        return self

    def update_dimensions(self, temp):
        """
        update the dimensions of the mosaic with new extents
        """
        # get number of bands
        t_name=None
        for this_t_name in['t','time']:
            if hasattr(temp, this_t_name):
                t_attr=getattr(temp, this_t_name)
                if hasattr(t_attr, '__len__') and len(t_attr) > 0:
                    self.dimensions[2]=len(t_attr)
                    setattr(self, this_t_name, t_attr.copy())
                    t_name=this_t_name
        if t_name is None:
            self.dimensions[2] = 1
        # calculate y dimensions with new extents
        self.dimensions[0] = np.int64((self.extent[3] - self.extent[2])/self.spacing[1]) + 1
        # calculate x dimensions with new extents
        self.dimensions[1] = np.int64((self.extent[1] - self.extent[0])/self.spacing[0]) + 1
        # calculate x and y arrays
        self.x = self.extent[0] + self.spacing[0]*np.arange(self.dimensions[1])
        self.y = self.extent[2] + self.spacing[1]*np.arange(self.dimensions[0])
        return self

    def image_coordinates(self, temp, validate=False):
        """
        get the image coordinates
        """
        iy = np.rint((temp.y[:,None]-self.extent[2])/self.spacing[1]).astype(np.int64)
        ix = np.rint((temp.x[None,:]-self.extent[0])/self.spacing[0]).astype(np.int64)

        if validate:
            iy1 = np.flatnonzero((iy >= 0) & (iy < self.shape[0]))[:, None]
            iy0 = iy[iy1.ravel(),:]
            ix1 = np.flatnonzero((ix >= 0) & (ix < self.shape[1]))[None,:]
            ix0 = ix[:,ix1.ravel()]
            return(iy0, ix0, iy1, ix1)
        else:
            return (iy,ix)

    def setup_bounds_from_list(self, in_list,
                group=None,
                fields=None,
                bounds=None,
                bands=None,
                reader=None):
        """
        Set the mosaic's extent and spacing from its inputs, removing from
        in_list any that cannot be read or fall outside bounds.  reader, a
        _TileReader, reads the files' metadata (in a pool, if it has one);
        without one they are read here, one at a time.
        """
        if reader is None:
            reader = _TileReader()
        metadata = reader.read(in_list.copy(), meta_only=True, group=group, bands=bands)
        for item, meta, read_error in metadata:
            if isinstance(item, str):
                # read tile grid from file
                try:
                    if read_error is not None:
                        raise read_error
                    temp=meta
                    if bounds is not None:
                        temp=temp.cropped(*bounds)
                    if temp is not None and (len(temp.x)>0) and (len(temp.y) > 0):
                        #update the mosaic bounds to include this tile
                        self.update_bounds(temp)
                        if self.spacing[0] is None:
                            self.update_spacing(temp)
                    else:
                        in_list.remove(item)
                except Exception:
                    print(f"failed to read group {group} "+ str(item))
                    in_list.remove(item)
            else:
                if bands is not None:
                    item=item[:,:,bands]
                if bounds is not None:
                    item=item.cropped(*bounds)
                if item is not None and (len(item.x)>0) and (len(item.y) > 0):
                    self.update_bounds(item)
                    if self.spacing[0] is None:
                        self.update_spacing(item)
                    temp=item
                else:
                    in_list.remove(item)

        self.update_dimensions(temp)
        self.__update_extent__()
        self.__update_size_and_shape__()


    def setup_fields(self, item, group=None, fields=None, bands=None, reader=None):
        '''
        Set up fields based on an input data item
        '''

        # create output mosaic, weights, and mask
        # read data grid from the first tile HDF5, use it to set the field dimensions

        if isinstance(item, str):
            if reader is None:
                reader = _TileReader()
            _, prototype, read_error = next(reader.read([item], group=group, fields=fields, bands=bands))
            if read_error is not None:
                raise read_error
        else:
            if bands is not None:
                item=item[:,:,bands]
            prototype=item
        if len(prototype.fields) == 0:
            message = f"pointCollection.grid.mosaic.py: did not find fields {fields} in file {prototype.filename}"
            return message
        # if time is not specified, squeeze extra dimenstions out of inputs
        self.use_time=False
        for field in ['time','t']:
            if hasattr(prototype, field) and getattr(prototype, field) is not None:
                self.use_time=True
        if not self.use_time:
            for field in prototype.fields:
                setattr(prototype, field, np.squeeze(getattr(prototype, field)))
        if fields is None:
            fields=prototype.fields
        these_fields=[field for field in fields if field in prototype.fields]
        self.field_dims={field:getattr(prototype, field).ndim for field in these_fields}
        for field in these_fields:
            self.assign({field:np.zeros(self.dimensions[0:self.field_dims[field]])})
        self.invalid = np.ones(self.dimensions,dtype=bool)
        self.__update_size_and_shape__()


    def replace(self, item, group=None, fields=None, bands=None):
        """
        Overwrite a section of the mosaic with an input file or mosaic
        """

        if fields is None:
            fields=self.fields.copy()
        # read data grid from HDF5
        if isinstance(item, str):
            temp=pc.grid.mosaic().from_file(item, group=group, fields=fields, bands=bands)
        else:
            if bands is not None:
                item=item[:,:,bands]
            temp=item
        if not temp.overlaps(self):
            return

        these_fields=[field for field in fields if field in temp.fields]
        # get the image coordinates of the input file
        iy0, ix0, iy1, ix1 = self.image_coordinates(temp, validate=True)
        for field in these_fields:
            try:
                if self.field_dims[field]==3:
                    field_data=np.atleast_3d(getattr(temp, field))
                    for band in range(self.dimensions[2]):
                        getattr(self, field)[iy0,ix0,band] = field_data[iy1,ix1,band]
                    self.invalid[iy0,ix0] = False
                else:
                    field_data=getattr(temp, field)
                    getattr(self, field)[iy0,ix0] = field_data[iy1, ix1]
                    self.invalid[iy0, ix0] = False
            except IndexError as e:
                thestr = f"problem with field {field}"
                if isinstance(item, str):
                    thestr += f" in group {group} in file {item}"
                else:
                    thestr += f" in item {item}"
                print(thestr)
                raise(e)


    def add(self, item, fields, group=None, use_time=False,
            pad=0, feather=0, bands=None):
        """
        Add all bands from an item to a mosaic.

        This method propagates invalid values from the input to the output, and
        adds all bands at once

        Parameters
        ----------
        item : str or mosaic
            Item to be added to the mosaic.
        fields : iterable
            Fields to be added.
        group : str optional
            Group from which to read fields, if item is a file  The default is None.
        use_time : bool, optional
            If true, preserve singleton dimensions in inputs. The default is False.
        pad : float, optional
            pad the weights by this distance. The default is 0.
        feather : float optional
            feather weights by this distance. The default is 0.
        bands : iterable, optional
            Bands to read from the item. The default is None (read all bands).

        Returns
        -------
        None.

        """

        if fields is None:
            fields=self.fields.copy()

        self.normalized=False
        # read data grid from file
        if isinstance(item, str):
            temp=pc.grid.mosaic().from_file(item, group=group, fields=fields, bands=bands)
        else:
            if bands is not None:
                item=item[:,:,bands]
            if isinstance(item, self.__class__):
                temp=item
            else:
                temp=pc.grid.mosaic().from_grid(item)
        if not temp.overlaps(self):
            return
        temp.update_spacing(temp)
        if not use_time:
            for field in temp.fields:
                setattr(temp, field, np.squeeze(getattr(temp, field)))
        these_fields=[field for field in fields if field in temp.fields]
        # copy weights for tile
        if self.tile_weight is not None:
            temp.weight=self.tile_weight.copy()
        else:
            temp.weights(pad=pad, feather=feather)
            self.tile_weights=temp.weight.copy()
        # get the image coordinates of the input file
        iy0, ix0, iy1, ix1 = self.image_coordinates(temp, validate=True)
        for field in these_fields:
            try:
                if self.field_dims[field]==3:
                    field_data=np.atleast_3d(getattr(temp, field))
                    bands=range(self.dimensions[2])
                    for band in bands:
                        getattr(self, field)[iy0,ix0,band] += field_data[iy1,ix1,band]*temp.weight[iy1, ix1]
                    self.invalid[iy0,ix0] = False
                else:
                    field_data=getattr(temp, field).copy()
                    getattr(self, field)[iy0, ix0] += field_data*temp.weight[iy1,ix1]
                    self.invalid[iy0,ix0] = False
            except (IndexError, ValueError) as e:
                thestr = f"problem with field {field}"
                if isinstance(item, str):
                    thestr += f"in group {group} in file {item}"
                else:
                    thestr += f"in item {item}"
                print(thestr)
                raise(e)
        # add weights to total weight matrix
        self.weight[iy0,ix0] += temp.weight[iy1,ix1]


    def add_to_band(self, item, fields, group=None,
                        pad=0, feather=0, band=None, in_band=None, out_band=None):
        """
        Add an item to a band of a mosaic.

        This method treats invalids in inputs as providing no data, rather
        than propagating the invalid to the output

        Parameters
        ----------
        item : str or mosaic
            Item to be added to the mosaic.
        fields : iterable
            Fields to be added.
        group : str optional
            Group from which to read fields, if item is a file  The default is None.
        use_time : bool, optional
            If true, preserve singleton dimensions in inputs. The default is False.
        pad : float, optional
            pad the weights by this distance. The default is 0.
        feather : float optional
            feather weights by this distance. The default is 0.
        in_band: band to read from input object
        out_band: band in self to write

        Returns
        -------
        None.

        """
        if in_band is None:
            in_band=band

        if out_band is None and in_band is not None:
            out_band = in_band

        if fields is None:
            fields=self.fields.copy()

        self.normalized=False
        # read data grid from file
        if isinstance(item, str):
            if in_band is None:
                temp=pc.grid.mosaic().from_file(item, group=group, fields=fields)
            else:
                temp=pc.grid.mosaic().from_file(item, group=group, fields=fields, bands=[in_band])
        else:
            if isinstance(item, self.__class__):
                if in_band is None:
                    temp=item
                else:
                    # slicing returns a plain grid.data (grid.data.__copy__)
                    temp=pc.grid.mosaic().from_grid(item[:,:,in_band])
            else:
                if in_band is None:
                    temp=pc.grid.mosaic().from_grid(item)
                else:
                    temp=pc.grid.mosaic().from_grid(item[:,:,in_band])
        temp.update_spacing(temp)
        if not temp.overlaps(self):
            return

        for field in temp.fields:
            setattr(temp, field, np.squeeze(getattr(temp, field)))
        these_fields=[field for field in fields if field in temp.fields]
        # copy weights for tile
        if self.tile_weight is not None:
            temp.weight=self.tile_weight.copy()
        else:
            temp.weights(pad=pad, feather=feather)
            self.tile_weight=temp.weight.copy()
        # get the image coordinates of the input file
        iy0, ix0, iy1, ix1 = self.image_coordinates(temp, validate=True)
        for field in these_fields:
            try:
                field_data=getattr(temp, field).copy()
                valid_mask = np.isfinite(field_data)
                if not np.all(valid_mask):
                    temp.weight[valid_mask==0]=0
                    field_data[valid_mask==0]=0
                if getattr(self, field).ndim==3:
                    getattr(self, field)[iy0, ix0, out_band] += field_data[iy1, ix1] * temp.weight[iy1,ix1]
                else:
                    getattr(self, field)[iy0, ix0] += field_data[iy1, ix1] * temp.weight[iy1,ix1]
                self.invalid[iy0, ix0] = False
            except (IndexError, ValueError) as e:
                thestr = f"problem with field {field}"
                if isinstance(item, str):
                    thestr += f"in group {group} in file {item}"
                else:
                    thestr += f"in item {item}"
                print(thestr)
                raise(e)
        # add weights to total weight matrix

        self.weight[iy0,ix0] += temp.weight[iy1,ix1]

    def normalize(self, fields=None, band=None, by_weight=True):
        """
        Normalize the mosaic fields by the sum of the weights

        Returns
        -------
        None.

        """

        # find invalid points:
        if by_weight:
            i_zero= np.flatnonzero( (self.weight == 0) | self.invalid)
            i_nonzero = np.flatnonzero(self.weight)
        else:
            i_zero= np.flatnonzero( self.invalid )

        # find valid points
        for field in self.fields:
            if self.field_dims[field]==3 and band is None:
                band_list = range(self.dimensions[2])
            elif band is not None:
                band_list=[band]
            else:
                band_list=[None]
            for this_band in band_list:
                if this_band is None:
                    if by_weight:
                        getattr(self, field)[:, :, this_band].ravel()[i_nonzero] /= self.weight.ravel()[i_nonzero]
                    getattr(self, field)[:, :, this_band].ravel()[i_zero] = self.fill_value
                else:
                    if by_weight:
                        iy_nz, ix_nz = np.unravel_index(i_nonzero, self.shape[0:2])
                        i_out = np.ravel_multi_index((iy_nz, ix_nz, np.zeros_like(ix_nz)+this_band), self.shape)
                        #temp_out=getattr(self, field)[:, :, this_band]
                        #temp_out.ravel()[i_nonzero] /= self.weight.ravel()[i_nonzero]
                        #getattr(self, field)[:,:,this_band]=temp_out
                        getattr(self, field).ravel()[i_out] /= self.weight.ravel()[i_nonzero]
                    if len(i_zero) > 0:
                        iy_z, ix_z = np.unravel_index(i_zero, self.shape[0:2])
                        i_out = np.ravel_multi_index((iy_z, ix_z, np.zeros_like(ix_z)+this_band), self.shape)
                        getattr(self, field).ravel()[i_out] = self.fill_value

        if self.weight is not None:
            self.weight[:]=0
        self.normalized=True

    def raised_cosine_weights(self, pad, feather):
        """
        Create smoothed weighting function using a raised cosine function
        """
        weights=[]
        for dim, xy, delta in zip([0, 1], [self.y, self.x], self.spacing):
            xy0 = np.mean(xy)
            eps_grid=0.01*delta
            W = xy[-1]-xy[0]
            dist = np.abs(xy-xy0)
            wt = np.zeros_like(dist)
            i_feather = (dist >= W/2 - pad - feather - eps_grid) & ( dist <= W/2 -pad )
            wt_feather = 0.5 + 0.5 * np.sin( -np.pi*(dist[i_feather] - (W/2 - pad - feather/2)) / (feather+2*delta))
            wt[ i_feather ] = wt_feather
            wt[ dist < W/2 - pad - feather - eps_grid] = 1
            wt[ dist > W/2 - pad + eps_grid] = 0
            weights += [wt]
        self.weight *= weights[0][:,None].dot(weights[1][None,:])

    def gaussian_weights(self, pad, feather):
        """
        Create smoothed weighting function using a Gaussian function
        """
        weights=[]
        for dim, xy in zip([0, 1], [self.x, self.y]):
            xy0 = np.mean(xy)
            W = xy[-1]-xy[0]
            dist = np.abs(xy-xy0)
            wt = np.zeros_like(dist)
            i_feather = (dist >= W/2 - pad - feather) & ( dist <= W/2 -pad )
            wt_feather = np.exp(-((xy[i_feather]-xy0)/(feather/2.))**2)
            wt[ i_feather ] = wt_feather
            wt[ dist <= W/2 - pad - feather ] = 1
            wt[ dist >= W/2 - pad] = 0
            weights += [wt]
        self.weight *= weights[0][:,None].dot(weights[1][None,:])

    def pad_edges(self, pad):
        """
        Pad the edges of the weights with zeros
        """
        weights=[]
        for dim, xy in zip([0, 1], [self.x, self.y]):
            xy0 = np.mean(xy)
            W = xy[-1]-xy[0]
            dist = np.abs(xy-xy0)
            wt=np.ones_like(dist)
            wt[ dist >= W/2 - pad] = 0
            weights += [wt]
        self.weight *= weights[0][:,None].dot(weights[1][None,:])

    def weights(self, pad=0, feather=0, mode='raised cosine'):
        """
        Create a weight matrix for a given grid
        """
        # find dimensions of matrix
        sh = getattr(self, self.fields[0]).shape
        if len(sh)==3:
            ny, nx, nband = sh
        else:
            ny, nx = sh
        # allocate for weights matrix
        self.weight = np.ones((ny,nx), dtype=float)
        # feathering the weight matrix
        if feather:
            if mode == 'raised cosine':
                self.raised_cosine_weights(pad, feather)
            elif mode == 'gaussian':
                self.gaussian_weights(pad, feather)
        if pad:
            self.pad_edges(pad)

        return self

    def from_list(self, in_list,
                  bounds=None,
                  fields=None,
                  group='/',
                  pad=0,
                  feather=0,
                  by_band=True,
                  verbose=False,
                  spacing=[None, None],
                  bands=None,
                  block_size=None,
                  workers=1,
                  ):
        """
        Generate a mosaic from a list of inputs.

        Inputs can be strings (indicating files) or pointCollection.grid or
        pointCollection.mosaic objects.  A string may be a URI
        (e.g. s3://bucket/key.h5).

        Parameters
        ----------
        in_list : iterable
            Objects to mosaic.
        bounds : iterable of iterables, optional
            specifies bounds of output mosaic, [[xmin, xmax], [ymin, ymax]].
            if None, the bounds will be determined from the inputs
            The default is None.
        fields : iterable of str, optional
            Specifies fields from inputs to be mosaicked.  If not specified, all
            will be mosaicked.  The default is None.
        pad : float, optional
            elements within pad will be removed from inputs before mosacking.
            The default is 0.
        feather : float, optional
            Smooth blending length for inputs. The default is 0.
        group : str, optional
            group in hdf5 or netcdf4 files to read. The default is '/'.
        bands : iterable, optional
            Bands (e.g. time slices) to read from each input, in order. If
            not specified, all bands in each input are read. The default is None.
        block_size : int, optional
            Bytes per range request for a remote HDF5 or netCDF input.  None
            leaves the filesystem's default, which for s3fs is 50 MiB: a
            mosaic reads a few fields from each tile, so with the default a
            remote tile is read almost whole.  io_utils.DEFAULT_REMOTE_BLOCK_SIZE
            suits this read.  Ignored for local files.  The default is None.
        workers : int, optional
            Read the files in this many processes (a pool).  Tiles are still
            added in list order, so the result is the same as a serial read;
            at most 2 x workers read tiles wait in memory.  Remote reads are
            latency-bound, so this is where the speed-up is.  Each worker
            costs its own interpreter and imports, ~0.35 GiB resident
            (measured, 2026-09): budget workers x 0.35 GiB on top of the
            mosaic.  The default is 1: files are read one at a time, in this
            process.

        Returns
        -------
        self
            pointCollection.grid.mosaic object containing the mosaicked data.

        """
        weight = (pad is not None and pad > 0) or (feather is not None and feather>0)

        with _TileReader(workers=workers, block_size=block_size) as reader:
            # with neither option the loops below hand the file names to add,
            # add_to_band and replace, which read them: the original code path
            prefetch = reader.pool is not None or block_size is not None

            def items(**read_kwargs):
                """(item, what to add, bands for the add, read error) in order"""
                if not prefetch:
                    for item in in_list:
                        yield item, item, read_kwargs.get('bands'), None
                    return
                for item, grid, read_error in reader.read(in_list, **read_kwargs):
                    # a file read here is already band-selected; an in-memory
                    # grid is band-selected by the method it goes to
                    yield item, grid, (None if isinstance(item, str) else read_kwargs.get('bands')), read_error

            self.setup_bounds_from_list(in_list, group=group, fields=fields, bounds=bounds,
                                        bands=bands, reader=reader)
            message = self.setup_fields(in_list[0], group=group, fields=fields, bands=bands,
                                        reader=reader)
            if message is not None:
                return message
            # add, add_to_band and replace read self.fields when fields is None
            read_fields = fields if fields is not None else self.fields.copy()
            # check if using a weighted summation scheme for calculating mosaic
            if weight:
                if by_band:
                    if len(self.shape)>2:
                        band_list=range(self.shape[2])
                    else:
                        band_list=[None]
                    for band in band_list:
                        self.invalid = np.ones(self.dimensions[0:2],dtype=bool)
                        self.weight = np.zeros((self.dimensions[0],self.dimensions[1]))
                        # if specific input bands were requested, map the output
                        # band index back to the corresponding input band
                        in_band = bands[band] if (bands is not None and band is not None) else band
                        in_bands = None if in_band is None else [in_band]
                        for item, grid, grid_bands, read_error in items(group=group, fields=read_fields, bands=in_bands):
                            if read_error is not None:
                                raise read_error
                            if prefetch and isinstance(item, str):
                                self.add_to_band(grid, group=group, fields=fields, pad=pad, feather=feather, in_band=None, out_band=band)
                            else:
                                self.add_to_band(item, group=group, fields=fields, pad=pad, feather=feather, in_band=in_band, out_band=band)
                        self.normalize(band=band)
                else:
                    self.invalid = np.ones(self.dimensions[0:2],dtype=bool)
                    self.weight = np.zeros((self.dimensions[0],self.dimensions[1]))
                    # for each file in the list
                    for item, grid, grid_bands, read_error in items(group=group, fields=read_fields, bands=bands):
                        try:
                            if read_error is not None:
                                raise read_error
                            self.add(grid, group=group, fields=fields, pad=pad, feather=feather, bands=grid_bands)
                        except Exception as e:
                            print(f"mosaic.from_list : problem with {item} for group={group} and fields={fields}")
                            print(e)
                    self.normalize()
            else:
                # overwrite the mosaic with each subsequent tile
                # for each file in the list
                self.invalid = np.ones(self.dimensions[0:2],dtype=bool)
                for item, grid, grid_bands, read_error in items(group=group, fields=read_fields, bands=bands):
                    if read_error is not None:
                        raise read_error
                    self.replace(grid, group=group, fields=fields, bands=grid_bands)
                self.normalize(by_weight=False)

        return self
