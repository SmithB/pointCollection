#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 15:35:13 2018

@author: ben
"""

from osgeo import gdal, gdalconst, osr, ogr
import numpy as np

import re
import io
import os
import bz2
import gzip
import uuid
import h5py
import pyproj
import netCDF4
import posixpath
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import scoreatpercentile
import pointCollection as pc
from . import DEM_date
from .fill_edges import fill_edges, smooth_corrected
import shapely

class data(object):
    """Class that holds gridded data.

    A pointCollection.grid.data object contains one or more fields containing
    2 or 3 dimensional arrays of data, and coordinates specifying the location
    of those data.  Data can be read from and written to hdf5, netCDF, and geotif
    files.

    In the representation of a pointCollection.grid object, the x and y coodinates
    specify the centers of the grid cells, similar to what is assumed in numpy
    interpolation routines.  The rows and columns of the grid have the same
    order as the coordinates arrays, and the order of the coordinates is assumed
    to be increasing.  This means that pointCollection.grid objects may be seen
    as 'upside down' in some plotting routines.
    """

    def __init__(self, fields=None, fill_value=np.nan, t_axis=2):
        self.x=None
        self.y=None
        self.projection=None
        self.filename=None
        self.extent=None
        self.img_extent=None
        self.interpolator={}
        self.nan_interpolator={}
        self.fill_value=fill_value
        self.time=None
        self.size=None
        self.shape=None
        self.t_axis=t_axis
        self.xform=None

        if fields is None:
            self.fields=list()
        else:
            self.fields=fields
        for field in self.fields:
            setattr(self, field, None)

    def __copy__(self, fields=None):
        """Return a copy of a grid, optionally with a subset of fields."""
        if fields is None:
            fields=self.fields
        temp=pc.grid.data()
        for field in ['x','y','projection','filename','extent','time', 't', 't_axis']:
            if hasattr(self, field):
                setattr(temp, field, getattr(self, field))
        for field in fields:
            if hasattr(self, field):
                setattr(temp, field, getattr(self, field).copy())
        temp.fields=fields.copy()
        temp.__update_size_and_shape__()
        temp.__update_extent__()
        return temp

    def __repr__(self):
        """Return a printable representation of a grid."""
        out=f"{self.__class__} with shape {self.shape},"+"\n"
        out += "with fields:"+"\n"
        out += f"{self.fields}"
        return out

    def copy(self, fields=None):
        """Return a copy of the current dataset."""
        return self.__copy__()

    def copy_meta(self):
        """Return an empty dataset matching the current dataset."""
        temp=pc.grid.data()
        for field in ['x','y','projection','filename','extent','time', 't', 't_axis']:
            if hasattr(self, field):
                setattr(temp, field, getattr(self, field))
        temp.__update_size_and_shape__()
        temp.__update_extent__()
        return temp

    def __getitem__(self, *args, **kwargs):
        """Return a subset of the data."""
        return self.copy_subset(*args, **kwargs)

    def __update_extent__(self):
        """Update the extent of the data to match its x and y fields."""
        try:
            self.extent=[np.min(self.x), np.max(self.x), np.min(self.y), np.max(self.y)]
            if len(self.x)>1:
                hdx=np.abs(self.x[1]-self.x[0])/2
            else:
                hdx=0
            if len(self.y)>1:
                hdy=np.abs(self.y[1]-self.y[0])/2
            else:
                hdy=0
            self.img_extent=[np.min(self.x)-hdx, np.max(self.x)+hdx, np.min(self.y)-hdy, np.max(self.y)+hdy]
        except ValueError:
            # usually happens when self.x or self.y is empty
            self.extent=[None, None, None, None]

    def __update_size_and_shape__(self):
        """Update the size and shape parameters of the object to match that of its data fields."""
        for field in ['z']+self.fields:
            try:
                self.size=getattr(self, field).size
                self.shape=getattr(self, field).shape
                break
            except Exception:
                pass

    def summary(self, return_table=True, return_dict=False):
         """
         Summarize the dimensions and statistics of a data object

         Parameters
         ----------
         return_table : boolean, optional
             If True, return a tab-delimeted (printable) table. The default is True.
         return_dict : boolean, optional
             If True. return a dictionary of parameters. The default is False.

         Returns
         -------
         dict or string or both
             string or dictionary giving the shape, number of nonzero entries,
             standard deviation, and minmax for each field.  If no output is
             specified, the summary is printed.

         """
         summary={}
         table = ['field \tshape \tn_finite \tSTD \t minmax']
         for field in self.fields:
             ff=getattr(self, field)
             finfo={'shape':ff.shape,
                             'n_finite':np.sum(np.isfinite(ff)),
                             'std':np.nanstd(ff),
                             'minmax':[np.nanmin(ff), np.nanmax(ff)]}

             table += [f'{field}\t{finfo["shape"]}\t{finfo["n_finite"]}\t{finfo["std"]:0.2e}\t{finfo["minmax"][0]:0.2e} {finfo["minmax"][1]:0.2e}']
             summary[field] = finfo
         out=[]
         if return_dict:
             out += [summary]
         if return_table:
             out += ['\n'.join(table)]
         if len(out)==1:
             return(out[0])
         else:
             return tuple(out)

    def from_dict(self, thedict):
        """
        Build a grid object from a python dictionary.

        Parameters
        ----------
        thedict: dict
            Dictionary of spatial grid variables
        """
        for field in thedict:
                setattr(self, field, thedict[field])
                if field not in self.fields and field not in ['x','y','time', 't']:
                    self.fields.append(field)
        self.__update_extent__()
        self.__update_size_and_shape__()
        return self

    def assign(self, *args, **kwargs):
        """Set field values."""
        if len(args):
            newdata=args[0]
        else:
            newdata=dict()
        if len(kwargs) > 0:
            newdata |= kwargs
        for field in newdata.keys():
            setattr(self, field, newdata[field])
            if field not in self.fields:
                self.fields.append(field)
        return self

    def from_list(self, D_list, t_axis=None, sort=False):
        """
        Build a grid object from a list of other grid objects.

        Parameters
        ----------
        D_list: list
            List of pointCollection grid objects
        t_axis: int or NoneType, default None
            Time axis for output concatenated pointCollection grid object
        sort: bool, default False
            Sort the output pointCollection grid object by time
        """
        if t_axis is None:
            t_axis = self.t_axis
        # get name of each field
        # get total number of time slices
        # get merged extent and spacing
        if len(self.fields)==0:
            fields=set()
            # number of time slices
            nt = 0
            # calculate merged extent
            xmin,xmax,ymin,ymax = [np.inf,-np.inf,np.inf,-np.inf]
            spacing = [None]*2
            for D in D_list:
                if hasattr(D,'fields'):
                    fields=fields.union(D.fields)
                # get length of time dimension
                try:
                    ntime = D.shape[D.t_axis]
                except:
                    nt += 1
                else:
                    nt += ntime
                # get spacing
                spacing[0] = D.x[1] - D.x[0]
                spacing[1] = D.y[1] - D.y[0]
                # get new extents of merged grid
                if (D.extent[0] < xmin):
                    xmin = np.copy(D.extent[0])
                if (D.extent[1] > xmax):
                    xmax = np.copy(D.extent[1])
                if (D.extent[2] < ymin):
                    ymin = np.copy(D.extent[2])
                if (D.extent[3] > ymax):
                    ymax = np.copy(D.extent[3])
            # convert unique fields to list
            self.fields=list(fields)
        # calculate x and y dimensions with new extents
        nx = np.int64((xmax - xmin)/spacing[0]) + 1
        ny = np.int64((ymax - ymin)/spacing[1]) + 1
        # calculate x and y arrays
        self.x = np.linspace(xmin,xmax,nx)
        self.y = np.linspace(ymin,ymax,ny)
        # try to extract times
        time = np.zeros((nt))
        i = 0
        for D in D_list:
            try:
                ntime = D.shape[D.t_axis]
            except:
                ntime = 1
            # try each of the possible time attributes
            if hasattr(D,'time') and D.time is not None:
                time[i:i+ntime] = D.time
                time_field='time'
            elif hasattr(D,'t'):
                time[i:i+ntime] = D.t
                time_field='t'
            i += ntime
        setattr(self, time_field, time)
        # if sorting by time
        if sort:
            isort = np.argsort(self.time)
            self.time = self.time[isort]
        # for each field
        for field in self.fields:
            if (t_axis == 0):
                data_field = np.zeros((nt,ny,nx))
            elif (t_axis == 2):
                data_field = np.zeros((ny,nx,nt))
            # counter for filling variables
            i = 0
            # for each object
            for D in D_list:
                try:
                    # calculate grid coordinates for merging fields
                    iy = np.array((D.y[:,None]-ymin)/spacing[1],dtype=int)
                    ix = np.array((D.x[None,:]-xmin)/spacing[0],dtype=int)
                    # number of time steps in field
                    try:
                        ntime = D.shape[D.t_axis]
                    except:
                        ntime = 1
                    # create temporary field with merged x/y dimensions
                    if (D.t_axis == 0):
                        temp = np.zeros((ntime,ny,nx))
                        temp[:,iy,ix] = getattr(D,field)
                    elif (D.t_axis == 2) or (len(D.shape) == 2):
                        temp = np.zeros((ny,nx,ntime))
                        temp[iy,ix,:] = np.atleast_3d(getattr(D,field))
                    # merge fields
                    if (t_axis == 0) and (D.t_axis == 0):
                        data_field[i:i+ntime,:,:] = temp[:]
                    elif (t_axis == 2) and (D.t_axis == 2):
                        data_field[:,:,i:i+ntime] = temp[:]
                    elif (t_axis == 0) and (D.t_axis == 2):
                        data_field[i:i+ntime,:,:] = np.transpose(temp[:],axes=(2,0,1))
                    elif (t_axis == 2) and (D.t_axis == 0):
                        data_field[:,:,i:i+ntime] = np.transpose(temp[:],axes=(1,2,0))
                    # add to counter
                    i += ntime
                except AttributeError:
                    print(f"Problem with field {field}")
            # sort data in field
            if sort and (t_axis == 0):
                data_field = data_field[isort,:,:]
            elif sort and (t_axis == 2):
                data_field = data_field[:,:,isort]
            # get the attribute for field
            setattr(self, field, data_field)
        # update the size and extent of the merged grid
        self.__update_extent__()
        self.__update_size_and_shape__()
        return self

    def from_file(self, raster_file, file_format=None, **kwargs):
        """
        Read a raster file.

        Parameters
        ----------
        raster_file: str
            geotif file
        file_format, str or NoneType, default None
            Raster file format

            - ``'geotif'``: geotiff
            - ``'h5'``: HDF5
            - ``'nc'``: netCDF4
        kwargs, dict, default {}
            Keyword arguments for the input reader
        """
        if file_format is None:
            file_format=os.path.splitext(raster_file)[1][1:]
        self.filename = raster_file
        if file_format.lower() in ('tif','tiff','geotif','geotiff','vrt'):
            return self.from_geotif(raster_file, **kwargs)
        elif file_format.lower() in ('h5','hdf','hdf5'):
            return self.from_h5(raster_file, **kwargs)
        elif file_format.lower() in ('nc','netcdf'):
            return self.from_nc(raster_file, **kwargs)

    def from_geotif(self, file, date_format=None, **kwargs):
        """
        Read a raster from a geotiff file.

        Parameters
        ----------
        file: str
            geotif file
        date_format, str or NoneType, default None
            Date format for geotiff file

            - ``'year'``: WorldView year from filename
            - ``'matlab'``: WorldView matlab date from filename
        """
        self.filename=file
        if date_format is not None:
            try:
                self.get_date(date_format)
            except Exception as e:
                print(f'pc.grid.data.from_geotif:\n\t while getting the date, encountered exception:\n\t{e}')
                print(f'\n\tfor filename {file}')
                print('\t\tsetting time to NaN')
                self.time=np.NaN
                return self
        try:
            ds=gdal.Open(file, gdalconst.GA_ReadOnly)
            self.from_gdal(ds, **kwargs)
        except Exception as e:
            if 'verbose' in kwargs and kwargs['verbose']:
                print(f'pc.grid.data.from_geotif:\n\t while reading the file, encountered exception:\n\t{e}')
                print(f'\n\tfor filename {file}')
        return self

    def from_gdal(self, ds, field='z', bands=None, bounds=None, extent=None,
                  skip=1, fill_value=np.nan, min_res=None, meta_only=False,
                  t_range=None,
                  verbose=False):
        """
        Make a pointCollection.grid.data from a gdal dataset.

        Parameters
        ----------
        ds : gdal dataset
            Can be a dataset from gdal.Open, or a memory dataset
        field : str, optional
            Fieldname for the read data. The default is 'z'.
        bands : list, optional
            Bands to read. The default is None.
        bounds : list-like, optional
            boundaries to read, [[xmin, xmax], [ymin, ymax]]. If not specified,
            read the whole file.  The default is None.
        extent : list-like, optional
            Extent of the file to read, [xmin, xmax, ymin, ymax].
            The default is None.
        skip : Integer, optional
            Specifies that every skip'th value should be read. The default is 1.
        min_res : float, optional
            Attempt to read with a skip value chosen to match min_res.
            The default is None.
        t_range : iterable, optional
            read bands that have metadata 'time' values between t_range[0] and t_range[1]
        meta_only : return raster extent without reading data
        verbose:  if true, report errors, etc
        Raises
        ------
        AttributeError
            If too many bands are requested, throws an error

        Returns
        -------
        TYPE
            pc.grid.data object containing the map data.
        """
        GT=ds.GetGeoTransform()

        if min_res is not None:
            skip=np.max([1, np.ceil(min_res/np.abs(GT[1]))]).astype(int)

        if fill_value is not None:
            self.fill_value = fill_value

        # check if there is time information
        bb=ds.GetRasterBand(1)
        if 'time' in bb.GetMetadata():
            t=np.zeros(ds.RasterCount)
            for this_band in range(ds.RasterCount):
                bb=ds.GetRasterBand(this_band+1)
                t[this_band]=float(bb.GetMetadata()['time'])
            if bands is not None:
                t=t[bands-1]
            elif t_range is not None:
                bands=np.flatnonzero((t>=t_range[0]) & (t<=t_range[1]))+1
                t=t[bands-1]
            self.t=t
        proj=ds.GetProjection()
        if bands is None:
            n_bands=ds.RasterCount
            bands=np.arange(n_bands, dtype=int)+1
        if not isinstance(bands, (list, tuple, np.ndarray)):
            bands=[bands]
        # get geolocation info, allocate outputs
        band=ds.GetRasterBand(1)
        nodataValue=band.GetNoDataValue()
        # ii and jj are the pixel center coordinates.  0,0 in GDAL is the upper-left
        # corner of the first pixel.
        ii=np.arange(0, band.XSize)+0.5
        jj=np.arange(0, band.YSize)+0.5
        x=GT[0]+GT[1]*ii
        y=GT[3]+GT[5]*jj

        if extent is not None:
            bounds=[[extent[0], extent[1]], [extent[2], extent[3]]]

        if bounds is not None:
            cols = np.where(( x>=bounds[0][0] ) & ( x<= bounds[0][1] ))[0]
            rows = np.where(( y>=bounds[1][0] ) & ( y<= bounds[1][1] ))[0]
        else:
            rows=np.arange(band.YSize, dtype=int)
            cols=np.arange(band.XSize, dtype=int)
        if skip >1:
            cols=cols[::skip]
            rows=rows[::skip]
        x=x[cols]
        y=y[rows]
        self.x=x
        self.y=y[::-1]
        self.projection=proj
        self.__update_extent__()
        if meta_only:
            return self

        z=list()
        for band_num in bands:
            if band_num > ds.RasterCount:
                raise AttributeError()
            band=ds.GetRasterBand(int(band_num))
            try:
                z.append(band.ReadAsArray(int(cols[0]), int(rows[0]), int(cols[-1]-cols[0]+1), int(rows[-1]-rows[0]+1))[::-1,:])
            except IndexError as e:
                    raise e
            if skip > 1:
                z[-1]=z[-1][::skip, ::skip]
        if len(bands)==1:
            z=z[0]
        else:
            z=np.stack(z, axis=2)
        ds=None

        if nodataValue is not None and np.isfinite(nodataValue):
            bad = z==np.array(nodataValue).astype(z.dtype)
            z = np.float64(z)
            z[bad] = self.fill_value
        else:
            z = np.float64(z)

        self.assign({field: z})
        self.__update_size_and_shape__()
        return self

    def h5_open(self, h5_file, mode='r', compression=None):
        """
        Open an HDF5 file with or without external compression.

        Parameters
        ----------
        h5_file: str
            HDF5 file
        mode: str, default 'r'
            Mode of opening the HDF5 file
        compression: str or NoneType, default None
            Compression format for the HDF5 file

            - ``None``: file is not externally compressed
            - ``'bzip'``
            - ``'gzip'``
        """
        if (compression is None):
            return h5py.File(h5_file,mode=mode)
        elif (compression == 'bzip'):
            # read bytes from bzip compressed file
            with bz2.BZ2File(h5_file) as fd:
                fid = io.BytesIO(fd.read())
                fid.seek(0)
                return h5py.File(fid, 'r')
        elif (compression == 'gzip'):
            # read gzip compressed file and extract into in-memory file object
            with gzip.open(h5_file,'r') as fd:
                fid = io.BytesIO(fd.read())
                fid.seek(0)
                return h5py.File(fid, 'r')

    def select_slices(self, bounds, x, y, bands):
        # get orientation of y-axis
        yorient = np.sign(y[1] - y[0])

        # reduce raster to bounds and orient to lower
        slices={}
        if (bounds is not None) and (yorient > 0):
            # indices to read
            xind, = np.nonzero((x >= bounds[0][0]) & (x <= bounds[0][1]))
            slices['x'] = slice(xind[0],xind[-1],1)
            yind, = np.nonzero((y >= bounds[1][0]) & (y <= bounds[1][1]))
            slices['y'] = slice(yind[0],yind[-1],1)
        elif (bounds is not None) and (yorient < 0):
            # indices to read with reversed y
            xind, = np.nonzero((x >= bounds[0][0]) & (x <= bounds[0][1]))
            slices['x'] = slice(xind[0],xind[-1],1)
            yind, = np.nonzero((y >= bounds[1][0]) & (y <= bounds[1][1]))
            slices['y'] = slice(yind[-1],yind[0],-1)
        elif (yorient < 0):
            # indices to read (all) with reversed y
            slices['x'] = slice(None,None,1)
            slices['y'] = slice(None,None,-1)
        else:
            # indices to read (all)
            slices['x'] = slice(None,None,1)
            slices['y'] = slice(None,None,1)

        if bands is None:
            slices['time'] = slice(None)
        else:
            slices['time'] = list(bands)
        return slices, yorient

    def choose_bands_by_time(self, t=None, bounds=None, t_range=None):
        """
        Select bands based on time stamps

        Parameters
        ----------
        t : iterable or float, optional
            time range. The default is None.
        bounds : iterable of iterables, optional
            bounds to be read, if it contains three values, the appropriate one will
            be used as the time bounds. The default is None.
        t_range : iterable, optional
            The time range to be read. The default is None.

        Returns
        -------
        bands : iterable
            Bands to be read from the file.
        t_range : iterable
            range of time values to be read.

        """
        bands=None
        if bounds is not None and len(bounds)==3:
            if self.t_axis==0:
                t_range=bounds[0]
                bounds=bounds[1:]
            else:
                t_range=bounds[2]
                bounds=bounds[:2]

        if t_range is not None:
            bands = np.flatnonzero((t>=t_range[0]) & (t<=t_range[1]))

        return bands, t_range


    def read_data(self, src, i0, i1, bands):
        """
        Read data from a data source (array or h5py dataset).

        Parameters
        ----------
        src : np.array, or h5py dataset
            source from which to read the data.
        i0, i1 : slice
            which values to read from first and second x-y dimension
        bands : slice
            which time slices to read

        Returns
        -------
        z : np.array
            data subsetted to match i0, i1, and bands

        """
        if len(src.shape) == 2:
            z = np.array(src[i0, i1])
        else:
            if bands is None:
                if len(src.shape) > 2:
                    if self.t_axis==2:
                        z = np.array(src[i0,i1,:])
                    elif self.t_axis==0:
                        z =  np.array(src[:,i0,i1])
                elif len(src.shape) == 2:
                    z = np.array(src[i0,i1])
            else:
                if self.t_axis==2:
                    z = np.array(src[i0,i1,bands])
                elif self.t_axis==0:
                    z = np.array(src[bands,i0,i1])
        return z

    def from_h5(self, h5_file, field_mapping=None, group='/', fields=None,
        xname='x', yname='y', timename='t',
        bounds=None,  skip=1, fill_value=None,
        t_axis=None, t_range=None, bands=None,
        compression=None, swap_xy=False, source_fillvalue=None):
        """
        Read a raster from an HDF5 file.

        Parameters
        ----------
        h5_file: str
            HDF5 file
        field_mapping: dict or NoneType, default None
            Field mapping for input and output variables
        group : str, default '/'
            HDF5 group to read variables
        fields: list or NoneType, default None
            Fields to read from the HDF5 file
        xname: str, default 'x'
            x-coordinate variable to read from the HDF5 file
        yname: str, default 'y'
            y-coordinate variable to read from the HDF5 file
        timename: str, default 't'
            time-coordinate variable to read from the HDF5 file
        bounds: list or NoneType, default None
            boundaries to read, [[xmin, xmax], [ymin, ymax]]. If not specified,
            read the whole file.
        t_range: range of time stamped bands to read.
        bands: list or NoneType, default None
            Bands to read. If not specified, reads all contained bands.
        skip: int, default 1
            Specifies that every skip'th value should be read.
        fill_value: float or NoneType, default None
            Value for invalid data points. If not specified, is the class default
        t_axis: int or NoneType, default None
            Time axis for output pointCollection grid object
        compression: str or NoneType, default None
            Compression format for the HDF5 file

            - ``None``: file is not externally compressed
            - ``'bzip'``
            - ``'gzip'``
        swap_xy: bool, default False
            Swap the orientation of x and y variables in the grid

        Returns
        -------
        self
            pc.grid.data object containing the map data.
        """
        if t_axis is not None:
            self.t_axis=t_axis
        if fill_value is not None:
            self.fill_value = fill_value

        if field_mapping is None:
            field_mapping={}
        self.filename=h5_file
        dims=[xname,yname,'t','time']
        if not group.startswith('/'):
            group='/'+group
        t=None
        src_t=None

        grid_mapping_name=None

        #default
        yorient=1
        with self.h5_open(h5_file, mode='r', compression=compression) as h5f:
            x=np.array(h5f[group][xname]).ravel()
            y=np.array(h5f[group][yname]).ravel()
            for time_var_name in set(['time','t', timename]):
                if time_var_name in h5f[group].keys():
                    t=h5f[group][time_var_name][:].copy()
                    timename=time_var_name
                    break
            if t is not None:
                src_t = t.copy()
            if bands is None:
                bands, t_range = self.choose_bands_by_time(t=t, bounds=bounds, t_range=t_range)
            if bands is not None and len(bands)==0:
                self.__update_extent__()
                self.__update_size_and_shape__()
                return self


            if t is not None and bands is not None:
                t=t[bands]
            # get orientation of y-axis
            if len(y) > 1:
                yorient = np.sign(y[1] - y[0])

            # if no field mapping provided, add everything in the group
            if len(field_mapping.keys())==0:
                for key in h5f[group].keys():
                    if key in dims or key in ['CRS','crs']:
                        continue
                    if fields is not None and key not in fields:
                        continue
                    if key not in field_mapping:
                        if hasattr(h5f[group][key],'shape'):
                            field_mapping.update({key:key})
            # reduce raster to bounds and orient to lower
            if bounds is not None:
                # indices to read
                xind, = np.nonzero((x >= bounds[0][0]) & (x <= bounds[0][1]))
                cols = slice(xind[0],xind[-1], skip)
                yind, = np.nonzero((y >= bounds[1][0]) & (y <= bounds[1][1]))
                rows = slice(yind[0],yind[-1], skip)
            else:
                # indices to read (all)
                rows = slice(None,None, skip)
                cols = slice(None,None, skip)

            if swap_xy:
                i1=rows
                i0=cols
            else:
                i0=rows
                i1=cols
            if not swap_xy:
                default_shape_2d = [len(y), len(x)]
            else:
                default_shape_2d = [len(x), len(y)]

            nT=1
            if src_t is not None:
                nT=len(src_t)

            if t_axis==0:
                default_shape_3d = [nT] + default_shape_2d
            else:
                default_shape_3d = default_shape_2d + [nT]


            # check that raster can be sliced
            if len(x[cols]) > 0 and len(y[rows]) > 0:
                for self_field in field_mapping:
                    f_field_name = posixpath.join(group,field_mapping[self_field])
                    if f_field_name not in h5f:
                        continue
                    f_field=h5f[f_field_name]
                    if f_field.shape in [tuple(default_shape_3d), tuple(default_shape_2d) ]:
                        z=self.read_data(f_field, i0, i1, bands)
                    elif src_t is not None and len(f_field)==np.prod(default_shape_3d):
                        z=self.read_data(np.array(f_field).reshape(default_shape_3d), i0, i1, bands)
                    elif len(f_field)==np.prod(default_shape_2d):
                        z=self.read_data(np.array(f_field).reshape(default_shape_3d), i0, i1, bands)
                    else:
                        raise(IndexError(f'from filename {h5_file}, field {f_field_name} has shape:{f_field.shape} incompatible with data shape:{default_shape_3d}.'))

                    # replace invalid values with nan
                    this_fillvalue=source_fillvalue
                    if this_fillvalue is None and hasattr(f_field, 'fillvalue'):
                        this_fillvalue=f_field.fillvalue
                    if this_fillvalue is not None:
                        try:
                            z[z == this_fillvalue] = self.fill_value
                        except ValueError:
                            # try setting the data type to the data type of the fill value
                            z=z.astype(self.fill_value.__class__)
                            z[z == this_fillvalue] = self.fill_value
                    # try to find grid mapping name from variable
                    try:
                        grid_mapping_name = f_field.attrs['grid_mapping']
                    except (KeyError,AttributeError):
                        pass
                    if swap_xy:
                        if len(z.shape) > 2:
                            if self.t_axis==0:
                                z=np.transpose(z, [0, 2, 1])
                            else:
                                z=np.transpose(z, [1, 0, 2])
                        else:
                            z=np.transpose(z, [1, 0])
                    # orient to lower:
                    if yorient == -1:
                        if len(z.shape) >2:
                            if self.t_axis ==0:
                                z=z[:,::-1,:]
                            else:
                                z=z[::-1,:,:]
                        else:
                            z=z[::-1,:]
                    # set output field
                    setattr(self, self_field, z)
                    if self_field not in self.fields:
                        self.fields.append(self_field)

                self.x=x[cols]
                self.y=y[rows]
                if yorient==-1:
                    y=y[::-1]

                if t is not None:
                    self.t=t
            # try to retrieve grid mapping and add to projection
            if grid_mapping_name is not None:
                self.projection = {}
                map_field_name = posixpath.join(group,grid_mapping_name)
                try:
                    for att_name,att_val in h5f[map_field_name].attrs.items():
                        self.projection[att_name] = att_val
                except:
                    pass
        self.__update_extent__()
        self.__update_size_and_shape__()
        return self

    def nc_open(self, nc_file, mode='r', compression=None):
        """
        Open a netCDF4 file with or without external compression.

        Parameters
        ----------
        nc_file: str
            netCDF4 file
        mode: str, default 'r'
            Mode of opening the netCDF4 file
        compression, str or NoneType, default None
            Compression format for the netCDF4 file

            - ``None``: file is not externally compressed
            - ``'bzip'``
            - ``'gzip'``
        """
        if (compression is None):
            return netCDF4.Dataset(nc_file, mode=mode)
        elif (compression == 'bzip'):
            # read bytes from bzipcompressed file
            with bz2.BZ2File(nc_file) as fd:
                return netCDF4.Dataset(uuid.uuid4().hex, mode=mode, memory=fd.read())
        elif (compression == 'gzip'):
            # read bytes from gzip compressed file
            with gzip.open(nc_file) as fd:
                return netCDF4.Dataset(uuid.uuid4().hex, mode=mode, memory=fd.read())

    def from_nc(self, nc_file, field_mapping=None, group='', fields=None,
        xname='x', yname='y', timename='time', bounds=None, t_range=None,
        bands=None, skip=1,
        fill_value=None,
        t_axis=None, compression=None):
        """
        Read a raster from a netCDF4 file.

        Parameters
        ----------
        nc_file: str
            netCDF4 file
        field_mapping: dict or NoneType, default None
            Field mapping for input and output variables
        group : str, default '/'
            netCDF4 group to read variables
        fields: list or NoneType, default None
            Fields to read from the netCDF4 file
        xname: str, default 'x'
            x-coordinate variable to read from the netCDF4 file
        yname: str, default 'y'
            y-coordinate variable to read from the netCDF4 file
        timename:str, default 'time'
            time-coordinate variable name to read from the netCDF4 file
        bounds: list or NoneType, default None
            boundaries to read, [[xmin, xmax], [ymin, ymax]]. If not specified,
            read the whole file.
        t_range: read time slices between these values (inclusive).  If None,
            read either values specified by bands, values specified by bounds,
            or the whole file
        bands: list or NoneType, default None
            Bands to read. If not specified, reads all contained bands.
        skip: int, default 1
            Specifies that every skip'th value should be read.
        fill_value: float or NoneType, default None
            Value for invalid data points. If not specified, is the class default
        t_axis: int or NoneType, default None
            Time axis for output pointCollection grid object
        compression: str or NoneType, default None
            Compression format for the netCDF4 file

            - ``None``: file is not externally compressed
            - ``'bzip'``
            - ``'gzip'``

        Returns
        -------
        self
            pc.grid.data object containing the map data.
        """
        dim_names={xname:'x',yname:'y', timename:'time'}

        if t_axis is not None:
            self.t_axis=t_axis

        if fill_value is not None:
            self.fill_value = fill_value

        if field_mapping is None:
            field_mapping={}
        self.filename=nc_file
        dims=[xname,yname,'t','time']
        t=None
        grid_mapping_name = None
        with self.nc_open(nc_file,mode='r',compression=compression) as fileID:
            # set automasking
            fileID.set_auto_mask(False)
            # check if reading from root group or sub-group
            ncf=fileID.groups[group] if group else fileID
            x=ncf.variables[xname][:].copy()
            y=ncf.variables[yname][:].copy()
            for this_time_var_name in set([timename, 'time','t']):
                if this_time_var_name in ncf.variables.keys():
                    t=ncf.variables[this_time_var_name][:].copy()
                    timename=this_time_var_name
                    break
            dim_names={xname:'x',yname:'y', timename:'time'}

            if bands is None:
                bands, t_range = self.choose_bands_by_time(t=t, bounds=bounds, t_range=t_range)
            if bands is not None and len(bands)==0:
                self.__update_extent__()
                self.__update_size_and_shape__()
                return self

            if t is not None and bands is not None:
                t=t[bands]

            # if no field mapping provided, add everything in the group
            if len(field_mapping.keys())==0:
                for key in ncf.variables.keys():
                    if key in dims:
                        continue
                    if fields is not None and key not in fields:
                        continue
                    if key not in field_mapping:
                        var = ncf.variables[key]
                        if hasattr(var,'shape') and var.shape:
                            field_mapping.update({key:key})

            slices=self.select_slices(bounds, x, y, bands)[0]
            # check that raster can be sliced
            if len(x[slices['x']]) == 0 or len(y[slices['y']]) == 0:
                self.__update_extent__()
                self.__update_size_and_shape__()
                return self

            for self_field in field_mapping:
                f_field_name=field_mapping[self_field]
                if f_field_name not in ncf.variables:
                    continue
                f_field=ncf.variables[f_field_name]

                # establish the order of output dimensions for this field
                # based on the time-axis dimension and the number of dimensions
                if len(f_field.shape)==2:
                    out_dims=['y','x']
                else:
                    if self.t_axis==2:
                        out_dims=['y','x','time']
                    else:
                        out_dims=['time','y','x']
                if f_field.dimensions is None:
                    f_dims = out_dims
                else:
                    f_dims = f_field.dimensions
                # order in which dimensions appear in the file
                this_dim_order = [dim_names[name] for name in f_dims if name in dim_names]
                # use the dimensions to make a tuple of slices with which to index the input field
                this_slice = tuple([slices[dim] for dim in this_dim_order])
                z = np.array(f_field[this_slice])

                # replace invalid values with nan
                if hasattr(f_field, '_FillValue'):
                    fill_value = f_field.getncattr('_FillValue')
                    try:
                        z[z == fill_value] = self.fill_value
                    except ValueError:
                        z=z.astype(float)
                        z[z == fill_value] = self.fill_value
                # find grid mapping name from variable
                if hasattr(f_field, 'grid_mapping'):
                    grid_mapping_name = f_field.getncattr('grid_mapping')

                # decide how to reorient the output product
                output_order = [this_dim_order.index(dim) for dim in out_dims if dim in this_dim_order]
                if not output_order==[0, 1, 2]:
                    z=z.transpose(output_order)

                # set output field
                self.assign({self_field:z})

            # reduce x and y to bounds
            self.x=x[slices['x']]
            self.y=y[slices['y']]
            if t is not None:
                self.t=t
            # try to retrieve grid mapping and add to projection
            if grid_mapping_name is not None:
                self.projection = {}
                try:
                    for att_name in ncf[grid_mapping_name].ncattrs():
                        self.projection[att_name] = ncf.variables[grid_mapping_name].getncattr(att_name)
                except:
                    pass
        self.__update_extent__()
        self.__update_size_and_shape__()
        return self

    def to_h5(self, out_file, fields=None, group='/', replace=False, nocompression=False, attributes={}, fill_value=None, overwrite_coords=False,  **kwargs):
        """Write a grid data object to an hdf5 file."""
        kwargs.setdefault('srs_proj4', None)
        kwargs.setdefault('srs_wkt', None)
        kwargs.setdefault('srs_epsg', None)
        kwargs.setdefault('dimensions', {})
        kwargs.setdefault('grid_mapping_name', 'crs')
        # check whether overwriting existing files
        # append to existing files as default
        mode = 'w' if replace else 'a'

        if fields is None:
            fields=self.fields
        if not group.startswith('/'):
            group='/'+group

        # update fill values in fields
        if fill_value is not None:
            self.replace_invalid(fields=fields, fill_value=fill_value)

        # get crs attributes
        self.crs_attributes(**kwargs)
        with h5py.File(out_file,mode) as h5f:
            # try adding file level attributes
            try:
                for att_name,att_val in attributes['ROOT'].items():
                    h5f.attrs[att_name] = att_val
            except Exception:
                pass

            # try creating the output group
            try:
                h5f.create_group(group)
            except Exception:
                pass

            # try adding group attributes
            try:
                for att_name,att_val in attributes[group].items():
                    h5f[group].attrs[att_name] = att_val
            except Exception:
                pass

            for field in ['x','y','time', 't'] + fields:
                f_field_name = posixpath.join(group,field)
                # if field exists, and overwrite_coords is True, overwrite it
                if field in h5f[group] and overwrite_coords:
                    if hasattr(self, field):
                        h5f[f_field_name][...] = getattr(self, field)
                else:
                    #Otherwise, try to create the dataset
                    try:
                        if nocompression or field in ['x','y','time']:
                            h5f.create_dataset(f_field_name, data=getattr(self, field))
                        else:
                            h5f.create_dataset(f_field_name, data=getattr(self, field),
                                chunks=True, compression="gzip", fillvalue=self.fill_value)
                    except Exception:
                         pass
                # try adding field attributes
                try:
                    for att_name,att_val in attributes[field].items():
                        h5f[f_field_name].attrs[att_name] = att_val
                except Exception:
                    pass
                # try adding dimensions
                try:
                    if field in ['x','y','time', 't'] and any(kwargs['dimensions']):
                        h5f[f_field_name].make_scale(field)
                    elif any(kwargs['dimensions']):
                        for i,dim in enumerate(kwargs['dimensions'][f_field_name]):
                            h5f[f_field_name].dims[i].attach_scale(h5f[dim])
                except Exception as exc:
                    pass
            # add crs attributes if applicable
            if self.crs:
                # add grid mapping attribute to each grid field
                for field in fields:
                    h5f[f_field_name].attrs['grid_mapping'] = kwargs['grid_mapping_name']
                # add grid mapping variable with projection attributes
                h5crs = h5f.create_dataset(kwargs['grid_mapping_name'], (), dtype=np.byte)
                for att_name,att_val in self.crs.items():
                    h5crs.attrs[att_name] = att_val
            if self.xform:
                h5f[group].attrs['xform_origin']=self.xform['origin']
                h5f[group].attrs['xform_basis_vectors']=self.xform['basis_vectors']


    def to_nc(self, out_file, fields=None, group='', replace=False, nocompression=False, attributes={}, fill_value=None, **kwargs):
        """Write a grid data object to a netCDF4 file."""
        kwargs.setdefault('srs_proj4', None)
        kwargs.setdefault('srs_wkt', None)
        kwargs.setdefault('srs_epsg', None)
        kwargs.setdefault('grid_mapping_name', 'crs')
        # check whether overwriting existing files
        # append to existing files as default
        mode = 'w' if replace else 'a'

        if fields is None:
            fields=self.fields
        # try getting the time variable name
        try:
            t_name = [field for field in ('t','time') if hasattr(self, field)
                and np.any(getattr(self,field))].pop()
        except Exception:
            t_name = None

        # update fill values in fields
        if fill_value is not None:
            self.replace_invalid(fields=fields, fill_value=fill_value)

        # get crs attributes
        self.crs_attributes(**kwargs)
        with netCDF4.Dataset(out_file,mode) as fileID:
            # try adding file level attributes
            try:
                for att_name,att_val in attributes['ROOT'].items():
                    fileID.setncattr(att_name, att_val)
            except Exception:
                pass
            # check if writing to root group or sub-group
            try:
                fileID.createGroup(group)
            except Exception:
                pass
            ncf=fileID.groups[group] if group else fileID
            # try adding group level attributes
            if group:
                try:
                    for att_name,att_val in attributes[group].items():
                        ncf.setncattr(att_name, att_val)
                except Exception:
                    pass

            # for each dimension variable
            for field in ['x','y','time','t']:
                # if field exists, overwrite it
                if field in ncf.variables.keys() and hasattr(self, field):
                    var = getattr(self, field)
                    if var is not None:
                        ncf.variables[field][:] = var
                elif hasattr(self, field):
                    var = getattr(self, field)
                    if var is not None:
                        ncf.createDimension(field, len(np.atleast_1d(var)))
                        _ = ncf.createVariable(field, var.dtype, (field,))
                        ncf.variables[field][:] = var
            for field in fields:
                if (self.t_axis == 2) and t_name is not None:
                    nc_dim=('y','x',t_name,)
                elif (self.t_axis == 0) and t_name is not None:
                    nc_dim=(t_name,'y','x',)
                else:
                    nc_dim=('y','x',)
                # if field exists, overwrite it
                if field in ncf.variables.keys() and hasattr(self, field):
                    var = getattr(self, field)
                    if var is not None:
                        ncf.variables[field][...] = var
                elif hasattr(self, field):
                    #Otherwise, try to create the dataset
                    var = getattr(self, field)
                    try:
                        ncf.createVariable(field, var.dtype, nc_dim,
                            zlib=np.logical_not(nocompression),
                            fill_value=self.fill_value)
                    except Exception:
                        pass
                    else:
                        ncf.variables[field][...] = var
                # try adding the field attributes
                try:
                    for att_name,att_val in attributes[field].items():
                        ncf[field].setncattr(att_name, att_val)
                except Exception:
                    pass
            # add crs attributes if applicable
            if self.crs:
                # add grid mapping attribute to each grid field
                for field in fields:
                    ncf[field].setncattr('grid_mapping', kwargs['grid_mapping_name'])
                # add grid mapping variable with projection attributes
                nccrs = ncf.createVariable(kwargs['grid_mapping_name'],np.byte,())
                for att_name,att_val in self.crs.items():
                    nccrs.setncattr(att_name, att_val)

    def to_geotif(self, out_file, **kwargs):
        """
        Write a grid object to a geotif.

        Parameters
        ----------
        out_file : str
            file name to write
        **kwargs :
            keywords to be passed to the to_gdal() method

        Returns
        -------
        None
        """
        if 'options' not in kwargs:
            options=["compress=LZW", "tiled=YES", "bigtiff=IF_SAFER"]
            kwargs['options'] = options
        out_ds=self.to_gdal(out_file=out_file, driver='GTiff',**kwargs)
        return out_ds

    def to_gdal(self, driver='MEM', out_file='', field='z', srs_proj4=None,
                srs_wkt=None, srs_epsg=None, EPSG=None,dtype=gdal.GDT_Float32,
                options=["compress=LZW"]):
        """
        Write a grid object to a gdal memory object.

        Parameters
        ----------
        driver: str, default 'MEM'
            GDAL driver for output raster
        out_file: str, default ''
            Output filename
        field: str, default 'z'
            Output field to write to raster
        srs_proj4: str or NoneType, default None
            PROJ4 projection string
        srs_wkt: str or NoneType, default None
            Well-Known Text (WKT) projection string
        srs_epsg: int or NoneType, default None
            EPSG projection code
        EPSG: int or Nonetype, default is None
            synonym for srs_epsg
        dtype: obj, default gdal.GDT_Float32
            GDAL data type for output raster
        options: list, default ["compress=LZW"]
            GDAL creation options for output raster
        """

        z=np.atleast_3d(getattr(self, field))
        ny,nx,nband = z.shape
        dx=np.abs(np.diff(self.x[0:2]))[0]
        dy=np.abs(np.diff(self.y[0:2]))[0]

        # no supported creation options with in memory rasters
        if driver=='MEM':
            options=[]

        # set up the dataset with creation options
        out_ds=gdal.GetDriverByName(driver).Create(out_file, nx, ny, nband, dtype, options=options)

        # top left x, w-e pixel resolution, rotation
        # top left y, rotation, n-s pixel resolution
        out_ds.SetGeoTransform((self.x.min()-dx/2, dx, 0, self.y.max()+dy/2, 0., -dy))

        if EPSG is not None:
            srs_epsg=EPSG
        # set the spatial projection reference information
        sr=osr.SpatialReference()
        if srs_proj4 is not None:
            sr.ImportFromProj4(srs_proj4)
        elif srs_wkt is not None:
            sr.ImportFromWkt(srs_wkt)
        elif srs_epsg is not None:
            sr.ImportFromEPSG(srs_epsg)
        else:
            raise ValueError("must specify at least one of srs_proj4, srs_wkt, srs_epsg")
        # export the spatial projection reference information to file
        out_ds.SetProjection(sr.ExportToWkt())
        # for each output band
        for band in range(nband):
            # change orientation to upper left corner
            out_ds.GetRasterBand(band+1).WriteArray(z[::-1,:,band])
            # set fill value for band
            try:
                fill_value = getattr(self,'fill_value')
                out_ds.GetRasterBand(band+1).SetNoDataValue(fill_value)
            except:
                pass
        if driver not in ('MEM',):
            out_ds.FlushCache()
            out_ds = None
        return out_ds

    def as_points(self, fields=None, keep_all=False):
        """Return a pointCollection.data object containing the points in the grid."""
        if fields is None:
            fields=self.fields


        if len(self.shape) == 2:
            x,y=np.meshgrid(self.x, self.y)
        else:
            if self.t_axis==0:
                if self.time is not None:
                    t, y, x = np.meshgrid(self.time, self.y, self.x, indexing='ij')
                else:
                    t, y, x = np.meshgrid(self.t, self.y, self.x, indexing='ij')
            elif self.t_axis==2:
                if self.time is not None:
                    y, x, t = np.meshgrid(self.y, self.x, self.time, indexing='ij')
                else:
                    y, x, t = np.meshgrid(self.y, self.x, self.t, indexing='ij')

        if keep_all:
            result =  pc.data(filename=self.filename).\
                from_dict({'x':x.ravel(),'y':y.ravel()})
            for field in fields:
                result.assign({field:getattr(self, field).ravel()})
        else:
            good=np.isfinite(getattr(self, fields[0])).ravel()
            result = pc.data(filename=self.filename).\
                from_dict({'x':x.ravel()[good],'y':y.ravel()[good]})
            for field in fields:
                result.assign({field:getattr(self, field).ravel()[good]})
        if len(self.shape)==2:
            if self.time is not None:
                result.assign({'time':self.time+np.zeros_like(getattr(result, field))})
        else:
            if self.time is not None:
                time_var='time'
            else:
                time_var='t'
            if keep_all:
                result.assign({time_var:t.ravel()})
            else:
                    result.assign({time_var:t.ravel()[good]})
        return result

    def get_latlon(self, srs_proj4=None, srs_wkt=None, srs_epsg=None):
        """
        Get the latitude and longitude of grid cells.

        Parameters
        ----------
        srs_proj4: str or NoneType, default None
            PROJ4 projection string
        srs_wkt: str or NoneType, default None
            Well-Known Text (WKT) projection string
        srs_epsg: int or NoneType, default None
            EPSG projection code

        Returns
        -------
        longitude: float
            longitude coordinates of grid cells
        latitude: float
            latitude coordinates of grid cells
        """
        # set the spatial projection reference information
        if srs_proj4 is not None:
            source = pyproj.CRS.from_proj4(srs_proj4)
        elif srs_wkt is not None:
            source = pyproj.CRS.from_wkt(srs_wkt)
        elif srs_epsg is not None:
            source = pyproj.CRS.from_epsg(srs_epsg)
        else:
            source = pyproj.CRS.from_string(self.projection)
        # target spatial reference (WGS84 latitude and longitude)
        target = pyproj.CRS.from_epsg(4326)
        # create transformation
        transformer = pyproj.Transformer.from_crs(source, target, always_xy=True)
        # create meshgrid of points in original projection
        x,y = np.meshgrid(self.x, self.y)
        # convert coordinates to latitude and longitude
        self.longitude,self.latitude = transformer.transform(x,y)
        return self

    def add_alpha_band(self, alpha=None, field='z', nodata_vals=None):
        """Add a transparencey band to a field."""
        if alpha is None:
            if nodata_vals is not None:
                alpha=np.ones_like(getattr(self, field)[:,:,0])
                if hasattr(nodata_vals, 'len') and len(nodata_vals)==3:
                    for ii in range(3):
                        alpha[~np.isfinite(getattr(self, field)[:,:,ii]) | (getattr(self, field)[:,:,ii]==nodata_vals[ii])]=0
                elif nodata_vals is not None:
                    alpha[np.all(~np.isfinite(getattr(self, field)) | (getattr(self, field)==nodata_vals), axis=2)]=0
            else:
                alpha=np.any(~np.isfinite(getattr(self, field)), axis=2)
        if len(getattr(self, field).shape)==3 and getattr(self, field).shape[2]==4:
            getattr(self, field)[:,:,-1]=alpha
        else:
            if len(alpha.shape)<3:
                alpha.shape=(alpha.shape[0], alpha.shape[1], 1)
            setattr(self, field, np.concatenate([getattr(self, field), alpha], axis=2))
        return self

    def get_date(self, date_format=None):
        """Get the date from the filename of a Worldview or tanDEM-x file."""
        if date_format is None or date_format == 'year':
            self.time = DEM_date.DEM_year(self.filename)
        elif date_format == 'matlab':
            self.time = DEM_date.DEM_MatlabDate(self.filename)
        return self

    def normalize(self, field='z', z0=[0., 255.], z1=[0., 1.], truncate=True, dtype=np.float64):
        """Normalize the z range of a grid object."""
        getattr(self, field)[:] = (getattr(self, field).astype(np.float64)-z0[0])/(z0[1]-z0[0])*(z1[1]-z1[0])+z1[0]
        if truncate:
            getattr(self, field)[getattr(self, field) < z1[0]] = z1[0]
            getattr(self, field)[getattr(self, field) > z1[1]] = z1[1]
        setattr(self, field, getattr(self, field).astype(dtype))
        return self

    def normalized(self, **kwargs):
        """Return a normalized copy."""
        return self.copy().normalize(**kwargs)

    def calc_gradient(self, field='z', band=0):
        """
        Calculate the gradient of a field.

        Parameters
        ----------
        field : str, optional
            DESCRIPTION. The default is 'z'.

        Returns
        -------
        None.
        """
        if len(getattr(self, field).shape) > 2:
            if self.t_axis==0:
                zz=getattr(self, field)[band,:,:]
            else:
                zz=getattr(self, field)[:,:,band]
        else:
            zz=getattr(self, field)

        gy, gx=np.gradient(zz, self.y, self.x)
        self.assign({field+'_x':gx, field+'_y':gy})

    def toRGB(self, cmap=None, field='z', bands=None, caxis=None, alpha=None):
        """Convert a field to RGB."""
        if bands is not None:
            if self.t_axis==0:
                setattr(self, field, getattr(self, field)[bands,:,:])
            elif self.t_axis==2:
                setattr(self, field, getattr(self, field)[:,:,bands])

        if caxis is None:
            caxis=[getattr(self, field).min(), getattr(self, field).max()]

        self.normalize(z0=caxis, field=field)
        if cmap is not None:
            setattr(self, field, cmap(getattr(self, field)))
        if alpha is not None:
            self.add_alpha_band(alpha)
        return self

    def index(self, row_ind, col_ind, band_ind=None, fields=None ):
        """
        subset the rows, columns, and bands of self

        Parameters
        ----------
        row_ind, col_ind, band_ind : bool, iterable of ints, or slice
            variables with which to slide self, band_ind is optional
        fields : iterable, optional
            fields to include in the output. The default is None.


        Returns
        -------
        pointCollection.grid.data
            indexed array

        """


        """Slice a grid by row or column."""
        if fields is None:
            fields=self.fields
        self.x=self.x[col_ind]
        self.y=self.y[row_ind]
        if band_ind is not None:
            for field in ['t','time']:
                try:
                    setattr(self, field, getattr(self, field)[band_ind])
                except:
                    pass
        for field in fields:
            if len(getattr(self, field).shape) == 2:
                setattr(self, field, getattr(self, field)[row_ind,:][:, col_ind])
            else:
                if self.t_axis==2:
                    if band_ind is None:
                        setattr(self, field, getattr(self, field)[row_ind,:, :][:, col_ind,:])
                    else:
                        setattr(self, field, getattr(self, field)[row_ind,:, :][:, col_ind,:][:, :, band_ind])
                elif self.t_axis==0:
                    if band_ind is None:
                        setattr(self, field, getattr(self, field)[:, row_ind,:][:, :, col_ind])
                    else:
                        setattr(self, field, getattr(self, field)[band_ind,:,:][:, row_ind, :][:, :, col_ind])
        self.__update_extent__()
        self.__update_size_and_shape__()
        return self

    def copy_subset(self, rc_ind, band_ind=None, fields=None):
        """
        return a subset of a copy of self

        Parameters
        ----------
        rc_ind : list
            slicing objects (boolean arrays, iterables of ints, or slices) for each
            dimension of self (following numpy rules)
        band_ind : boolean array, iterable of ints, or slice, optional
            slicing object for dimension 2. The default is None.
        fields : list, optional
            fields to include in the output. The default is None.

        Returns
        -------
        pointCollection.grid.data
            sliced copy of self.

        """
        if fields is None:
            fields=self.fields
        if len(rc_ind) > 2:
            if self.t_axis==2:
                return self.copy(fields=fields).index(rc_ind[0], rc_ind[1], band_ind=rc_ind[2])
            else:
                return self.copy(fields=fields).index(rc_ind[1], rc_ind[2], band_ind=rc_ind[0])
        return self.copy(fields=fields).index(rc_ind[0], rc_ind[1], band_ind=band_ind)

    def crop(self, XR, YR, TR=None, fields=None):
        '''
        Crop self to specified bounds


        Parameters
        ----------
        XR, YR, TR : iterables
            two-element iterables specifying the range in each dimension. TR is optional
        fields : iterable, optional
            strings specifying fields to include in the output. The default is None.

        Raises
        ------
        IndexError
            If TR is specified and neither self.time nor self.t is defined, the
            subset is not defined

        Returns
        -------
        pointCollection.grid.data
            cropped version of self

        '''

        col_ind = np.flatnonzero((self.x >= XR[0]) & (self.x <= XR[1]))
        row_ind = np.flatnonzero((self.y >= YR[0]) & (self.y <= YR[1]))
        time_ind = None
        if TR is not None:
            if self.time is not None:
                time_ind = np.flatnonzero((self.time >= TR[0]) & (self.time <= TR[1]))
            elif self.t is not None:
                time_ind = np.flatnonzero((self.t >= TR[0]) & (self.t <= TR[1]))
            else:
                raise IndexError('neither t nor time defined')
        try:
           self.index(row_ind, col_ind, fields=fields, band_ind=time_ind)
           return self
        except Exception as e:
           print("grid: self extent is: ", self.extent)
           print("XR is " + str(XR))
           print("YR is " + str(YR))
           print("Error is" )
           print(e)

    def cropped(self, *args, **kwargs):
        """Return a cropped copy of a grid."""
        return self.copy().crop(*args, **kwargs)

    def show(self, field='z', band=None, ax=None, xy_scale=1, gradient=False, ddt=None, stretch_pct=None, **kwargs):
        """Make a matplotlib image of a grid."""
        import matplotlib.pyplot as plt
        kwargs['extent']=np.array(self.img_extent)*xy_scale
        kwargs['origin']='lower'
        if 'interpolation' not in kwargs:
            kwargs['interpolation']='nearest'
        if ax is None:
            ax = plt.gca()
        if band is None:
            zz=getattr(self, field)
        elif (band is not None) and (self.t_axis==2):
            zz=getattr(self, field)[:,:,band]
        elif (band is not None) and (self.t_axis==0):
            zz=getattr(self, field)[band,:,:]

        if ddt is not None:
            if self.t is not None:
                t=self.t
            elif self.time is not None:
                t=self.time
            else:
                t=range(self.shape[self.t_axis])
            if self.t_axis==0:
                zz = getattr(self, field)[ddt[1],:,:]-getattr(self, field)[ddt[0],:,:]
            else:
                zz = getattr(self, field)[:,:,ddt[1]]-getattr(self, field)[:,:,ddt[0]]
            zz /= (t[ddt[1]]-t[ddt[0]])

        if gradient:
            zz=np.gradient(zz.squeeze(), self.x[1]-self.x[0], self.y[1]-self.y[0])[0]
            if 'stretch_pct' not in kwargs and 'clim' not in kwargs:
                stretch_pct=[5, 95]
            if 'cmap' not in kwargs:
                kwargs['cmap']='gray'

        if stretch_pct is not None:
            LH=scoreatpercentile(zz.ravel()[np.isfinite(zz.ravel())], stretch_pct)
            kwargs['vmin']=LH[0]
            kwargs['vmax']=LH[1]

        h_im = ax.imshow(zz, **kwargs)
        return h_im

    def interp(self, x, y, t=None, gridded=False, band=None, field='z', replace=False):
        '''
        interpolat a grid at a set of points

        Parameters
        ----------
        x, y, t: iterables
            Coordinates at which to interpolate self.  t is optional.
        gridded : bool, optional
            If true, x, y, and t are treated as the coordinates of a grid. The default is False.
        band : int, optional
           Band to be sampled. The default is None.
        field : str, optional
            field to sample. The default is 'z'.
        replace : bool, optional
            If true, the existing interpolator object defined for field is replaced.
            The default is False.

        Returns
        -------
        result : numpy.array
            interpolated values.

        '''

        """Interpolate a 2-D grid to a set of x and y points."""
        if (field not in self.interpolator) or replace:
            if band is not None:
                if (len(getattr(self, field).shape) > 2) and (self.t_axis==2):
                    z0 = getattr(self, field)[:,:,band].copy()
                elif (len(getattr(self, field).shape) > 2) and (self.t_axis==0):
                    z0 = getattr(self, field)[band,:,:].copy()
                else:
                    z0 = getattr(self, field).copy()
            else:
                z0 = getattr(self, field).copy()
            if len(getattr(self, field).shape)==2 or band is not None:
                grid_vars=(self.y, self.x)
            else:
                if self.t_axis==0:
                    grid_vars=(self.t, self.y, self.x)
                else:
                    grid_vars=(self.y, self.x, self.t)
            self.interpolator[field] = RegularGridInterpolator(
                                                    grid_vars,
                                                    z0, bounds_error=False)
        if gridded:
            if len(self.interpolator[field].grid)==3:
                if self.t_axis==0:
                    tg, yg, xg = np.meshgrid(t, y, x, indexing='ij')
                    result=self.interpolator[field]((tg, yg, xg))
                elif self.t_axis==2:
                    yg, xg, tg = np.meshgrid(y, x, t, indexing='ij')
                    result=self.interpolator[field]((yg, xg, tg))
            else:
                xg, yg=np.meshgrid(x, y)
                result=self.interpolator[field]((yg, xg))
        else:
            if len(self.interpolator[field].grid)==3:
                if self.t_axis==0:
                    result=self.interpolator[field]((t, y, x))
                else:
                    result=self.interpolator[field]((y, x, t))
            else:
                result=self.interpolator[field]((y,x))
        return result

    def bounds(self, pad=0):
        """
        Return the bounds of the x and y coordinates in self.

        Parameters
        ----------
        pad : float, int, optional
            amount by which to pad the returned bounds. The default is 0.

        Returns
        -------
        XR, YR: minimum and maximum of x and y

        """

        return [[np.min(self.x)-pad, np.max(self.x)+pad], [np.min(self.y)-pad, np.max(self.y)+pad]]

    def replace_invalid(self, fields=None, fill_value=np.nan):
        """Replace invalid values within a field with a new fill_value."""
        if fields is None:
            fields=self.fields
        # for each field
        for field in fields:
            the_field = getattr(self, field)
            indices = np.nonzero(np.isnan(the_field) | (the_field == self.fill_value))
            the_field[indices] = fill_value
            setattr(self, field, the_field)
        # update class fill value
        self.fill_value = fill_value
        return self

    def fill_smoothed(self, w_smooth=1, fields=None):
        if fields is None:
            fields=self.fields
        for field in fields:
            dim = self.t_axis if len(getattr(self,field).shape)>2 else None
            getattr(self, field)[:]=fill_edges(getattr(self, field), w_smooth=w_smooth, dim=dim)

    def rasterize_poly(self, in_geom, field='z', burn_value=1, raster_epsg=None, poly_epsg=None):
        """Rasterize a shapely polygon"""
        from osgeo import ogr
        mask_ds=self.to_gdal(srs_epsg=raster_epsg, field=field)
        shpDriver = ogr.GetDriverByName("memory")

        # create the spatial reference system, WGS84
        srs =  osr.SpatialReference()
        srs.ImportFromEPSG(poly_epsg)

        outDataSource = shpDriver.CreateDataSource('temp')
        if isinstance(in_geom, shapely.Polygon):
            outLayer = outDataSource.CreateLayer('temp', srs, geom_type=ogr.wkbPolygon)
        else:
            outLayer = outDataSource.CreateLayer('temp', srs, geom_type=ogr.wkbMultiPolygon)
        # Add an ID field
        idField = ogr.FieldDefn("id", ogr.OFTInteger)
        outLayer.CreateField(idField)

        geom=ogr.CreateGeometryFromWkb(in_geom.wkb)

        # Create the feature and set values
        featureDefn = outLayer.GetLayerDefn()
        feature = ogr.Feature(featureDefn)
        feature.SetGeometry(geom)
        feature.SetField("id", 1)
        outLayer.CreateFeature(feature)
        gdal.RasterizeLayer(mask_ds, [1], outLayer, burn_values=[burn_value])
        mask_ds=pc.grid.data().from_gdal(mask_ds)
        # write the field from the mask ds to self
        setattr(self, field, mask_ds.z)

    def crs_attributes(self, srs_proj4=None, srs_wkt=None, srs_epsg=None, **kwargs):
        """
        Return a dictionary of attributes for a projection.

        Parameters
        ----------
        srs_proj4: str or NoneType, default None
            PROJ4 projection string
        srs_wkt: str or NoneType, default None
            Well-Known Text (WKT) projection string
        srs_epsg: int or NoneType, default None
            EPSG projection code
        """
        # output projection attributes dictionary
        self.crs = {}
        # set the spatial projection reference information
        sr = osr.SpatialReference()
        if srs_proj4 is not None:
            sr.ImportFromProj4(srs_proj4)
        elif srs_wkt is not None:
            sr.ImportFromWkt(srs_wkt)
        elif srs_epsg is not None:
            sr.ImportFromEPSG(srs_epsg)
        else:
            return
        # convert proj4 string to dictionary
        proj4_dict = {p[0]:p[1] for p in re.findall(r'\+(.*?)\=([^ ]+)',sr.ExportToProj4())}
        # get projection attributes
        try:
            self.crs['spatial_epsg'] = int(sr.GetAttrValue('AUTHORITY',1))
        except Exception:
            pass
        self.crs['crs_wkt'] = sr.ExportToWkt()
        self.crs['spatial_ref'] = sr.ExportToWkt()
        self.crs['proj4_params'] = sr.ExportToProj4()
        # get datum attributes
        self.crs['semi_major_axis'] = sr.GetSemiMajor()
        self.crs['semi_minor_axis'] = sr.GetSemiMinor()
        self.crs['inverse_flattening'] = sr.GetInvFlattening()
        self.crs['reference_ellipsoid_name'] = sr.GetAttrValue('DATUM',0)
        self.crs['geographic_crs_name'] = sr.GetAttrValue('GEOGCS',0)
        # get projection attributes
        self.crs['projected_crs_name'] = sr.GetName()
        self.crs['grid_mapping_name'] = sr.GetAttrValue('PROJECTION',0).lower()
        self.crs['standard_name'] = sr.GetAttrValue('PROJECTION',0)
        self.crs['prime_meridian_name'] = sr.GetAttrValue('PRIMEM',0)
        self.crs['longitude_of_prime_meridian'] = float(sr.GetAttrValue('PRIMEM',1))
        try:
            self.crs['latitude_of_projection_origin'] = float(proj4_dict['lat_0'])
        except Exception:
            pass
        self.crs['standard_parallel'] = float(sr.GetProjParm(osr.SRS_PP_LATITUDE_OF_ORIGIN,1))
        self.crs['straight_vertical_longitude_from_pole'] = float(sr.GetProjParm(osr.SRS_PP_CENTRAL_MERIDIAN,1))
        self.crs['false_northing'] = float(sr.GetProjParm(osr.SRS_PP_FALSE_NORTHING,1))
        self.crs['false_easting'] = float(sr.GetProjParm(osr.SRS_PP_FALSE_EASTING,1))
        self.crs['scale_factor'] = float(sr.GetProjParm(osr.SRS_PP_SCALE_FACTOR,1))
        return self
