#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 15:35:13 2018

@author: ben
"""

from osgeo import gdal, gdalconst, osr
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
from scipy.interpolate import RectBivariateSpline
from scipy.stats import scoreatpercentile
import pointCollection as pc
from . import WV_date
#import os

class data(object):
    def __init__(self, fields=None, fill_value=np.nan, t_axis=2):
        self.x=None
        self.y=None
        self.projection=None
        self.filename=None
        self.extent=None
        self.interpolator={}
        self.nan_interpolator={}
        self.fill_value=fill_value
        self.time=None
        self.size=None
        self.shape=None
        self.t_axis=t_axis

        if fields is None:
            self.fields=list()
        else:
            self.fields=fields
        for field in self.fields:
            setattr(self, field, None)

    def __copy__(self, fields=None):
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
        return temp

    def __repr__(self):
        out=f"{self.__class__} with shape {self.shape},"+"\n"
        out += "with fields:"+"\n"
        out += f"{self.fields}"
        return out

    def copy(self, fields=None):
        return self.__copy__()

    def copy_meta(self):
        temp=pc.grid.data()
        for field in ['x','y','projection','filename','extent','time', 't', 't_axis']:
            if hasattr(self, field):
                setattr(temp, field, getattr(self, field))
        return temp
    
    def __getitem__(self, *args, **kwargs):
        """
        wrapper for the copy_subset() method
        """
        return self.copy_subset(*args, **kwargs)

    def __update_extent__(self):
        try:
            self.extent=[np.min(self.x), np.max(self.x), np.min(self.y), np.max(self.y)]
        except ValueError:
            # usually happens when self.x or self.y is empty
            self.extent=[None, None, None, None]

    def __update_size_and_shape__(self):
        """
        update the size and shape parameters of the object to match that of its data fields
        """
        for field in ['z']+self.fields:
            try:
                self.size=getattr(self, field).size
                self.shape=getattr(self, field).shape
                break
            except Exception:
                pass

    def from_dict(self, thedict):
        """
        build a grid object from a python dictionary

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

    def assign(self, newdata):
        for field in newdata.keys():
            setattr(self, field, newdata[field])
            if field not in self.fields:
                self.fields.append(field)
        return self

    def from_list(self, D_list, t_axis=None, sort=False):
        """
        build a grid object from a list of other grid objects

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
        nx = np.int((xmax - xmin)/spacing[0]) + 1
        ny = np.int((ymax - ymin)/spacing[1]) + 1
        # calculate x and y arrays
        self.x = np.linspace(xmin,xmax,nx)
        self.y = np.linspace(ymin,ymax,ny)
        # try to extract times
        self.time = np.zeros((nt))
        i = 0
        for D in D_list:
            try:
                ntime = D.shape[D.t_axis]
            except:
                ntime = 1
            # try each of the possible time attributes
            if getattr(D,'time'):
                self.time[i:i+ntime] = D.time
            elif getattr(D,'t'):
                self.time[i:i+ntime] = D.t
            i += ntime
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
                    iy = np.array((D.y[:,None]-xmin)/spacing[1],dtype=int)
                    ix = np.array((D.x[None,:]-ymin)/spacing[0],dtype=int)
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
        Wrapper function for reading a raster file

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
        if file_format.lower() in ('tif','tiff','geotif','geotiff'):
            return self.from_geotif(raster_file, **kwargs)
        elif file_format.lower() in ('h5','hdf','hdf5'):
            return self.from_h5(raster_file, **kwargs)
        elif file_format.lower() in ('nc','netcdf'):
            return self.from_nc(raster_file, **kwargs)

    def from_geotif(self, file, date_format=None, **kwargs):
        """
        Read a raster from a geotiff file

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
            self.get_date(date_format)

        ds=gdal.Open(file, gdalconst.GA_ReadOnly)
        self.from_gdal(ds, **kwargs)
        return self

    def from_gdal(self, ds, field='z', bands=None, bounds=None, extent=None, skip=1, fill_value=np.nan, min_res=None):
        """
        make a pointCollection.grid.data from a gdal dataset

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
        min_res : TYPE, optional
            Attempt to read with a skip value chosen to match min_res.
            The default is None.

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

        if skip >1:
            cols=cols[::skip]
            rows=rows[::skip]
        if nodataValue is not None and np.isfinite(nodataValue):
            bad = z==np.array(nodataValue).astype(z.dtype)
            z = np.float64(z)
            z[bad] = self.fill_value
        else:
            z = np.float64(z)
        x=x[cols]
        y=y[rows]
        self.x=x
        self.y=y[::-1]
        self.assign({field: z})
        self.projection=proj
        self.__update_extent__()
        self.__update_size_and_shape__()
        return self

    def h5_open(self, h5_file, mode='r', compression=None):
        """
        Open an HDF5 file with or without external compression

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


    def read_data(self, src, i0, i1, bands):
        '''
        read data from a data source (array or h5py dataset)

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

        '''
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
        xname='x', yname='y', bounds=None, bands=None, skip=1, fill_value=None,
        t_axis=None, compression=None, swap_xy=False, source_fillvalue=None):
        """
        Read a raster from an HDF5 file

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
        bounds: list or NoneType, default None
            boundaries to read, [[xmin, xmax], [ymin, ymax]]. If not specified,
            read the whole file.
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
        if group[0] != '/':
            group='/'+group
        t=None
        src_t=None

        grid_mapping_name=None

        #default
        yorient=1
        with self.h5_open(h5_file, mode='r', compression=compression) as h5f:
            x=np.array(h5f[group+'/'+xname]).ravel()
            y=np.array(h5f[group+'/'+yname]).ravel()
            if 't' in h5f[group]:
                t=np.array(h5f[group]['t'])
            elif 'time' in h5f[group]:
                t=np.array(h5f[group]['time'])
            if t is not None:
                src_t = t.copy()
            if t is not None and bands is not None:
                t=t[bands]
            # get orientation of y-axis
            if len(y) > 1:
                yorient = np.sign(y[1] - y[0])

            # if no field mapping provided, add everything in the group
            if len(field_mapping.keys())==0:
                for key in h5f[group].keys():
                    if key in dims:
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
                    f_field_name=group+'/'+field_mapping[self_field]
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
                        raise(IndexError(f'field {f_field_name} has shape:{f_field.shape} incompatible with data shape:{default_shape_3d}.'))

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
                try:
                    for att_name,att_val in h5f[group+'/'+grid_mapping_name].attrs.items():
                        self.projection[att_name] = att_val
                except:
                    pass
        self.__update_extent__()
        self.__update_size_and_shape__()
        return self

    def nc_open(self, nc_file, mode='r', compression=None):
        """
        Open a netCDF4 file with or without external compression

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
        xname='x', yname='y', bounds=None, bands=None, skip=1, fill_value=None,
        t_axis=None, compression=None):
        """
        Read a raster from a netCDF4 file

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
        bounds: list or NoneType, default None
            boundaries to read, [[xmin, xmax], [ymin, ymax]]. If not specified,
            read the whole file.
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
            if 't' in ncf.variables.keys():
                t=ncf.variables['t'][:].copy()
            elif 'time' in ncf.variables.keys():
                t=ncf.variables['time'][:].copy()

            if t is not None and bands is not None:
                t=t[bands]
            # get orientation of y-axis
            yorient = np.sign(y[1] - y[0])

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
            # reduce raster to bounds and orient to lower
            if (bounds is not None) and (yorient > 0):
                # indices to read
                xind, = np.nonzero((x >= bounds[0][0]) & (x <= bounds[0][1]))
                cols = slice(xind[0],xind[-1],1)
                yind, = np.nonzero((y >= bounds[1][0]) & (y <= bounds[1][1]))
                rows = slice(yind[0],yind[-1],1)
            elif (bounds is not None) and (yorient < 0):
                # indices to read with reversed y
                xind, = np.nonzero((x >= bounds[0][0]) & (x <= bounds[0][1]))
                cols = slice(xind[0],xind[-1],1)
                yind, = np.nonzero((y >= bounds[1][0]) & (y <= bounds[1][1]))
                rows = slice(yind[-1],yind[0],-1)
            elif (yorient < 0):
                # indices to read (all) with reversed y
                rows = slice(None,None,-1)
                cols = slice(None,None,1)
            else:
                # indices to read (all)
                rows = slice(None,None,1)
                cols = slice(None,None,1)

            # check that raster can be sliced
            if len(x[cols]) > 0 and len(y[rows]) > 0:
                for self_field in field_mapping:
                    f_field_name=field_mapping[self_field]
                    if f_field_name not in ncf.variables:
                        continue
                    f_field=ncf.variables[f_field_name]

                    if len(f_field.shape) == 2:
                        z = np.array(f_field[rows,cols])
                    else:
                        if len(f_field.shape) == 2:
                            z = np.array(f_field[rows,cols])
                        elif bands is None and len(f_field.shape) > 2:
                            if self.t_axis==2:
                                z = np.array(f_field[rows,cols,:])
                            elif self.t_axis==0:
                                z = np.array(f_field[:,rows,cols])
                        else:
                            if self.t_axis==2:
                                z = np.array(f_field[rows,cols,bands])
                            elif self.t_axis==0:
                                z = np.array(f_field[bands,rows,cols])
                    # replace invalid values with nan
                    if hasattr(f_field, '_FillValue'):
                        fill_value = f_field.getncattr('_FillValue')
                        try:
                            z[z == fill_value] = self.fill_value
                        except ValueError as e:
                            z=z.astype(float)
                            z[z == fill_value] = self.fill_value
                    # find grid mapping name from variable
                    if hasattr(f_field, 'grid_mapping'):
                        grid_mapping_name = f_field.getncattr('grid_mapping')
                    # set output field
                    setattr(self, self_field, z)
                    if self_field not in self.fields:
                        self.fields.append(self_field)
                # reduce x and y to bounds
                self.x=x[cols]
                self.y=y[rows]
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

    def to_h5(self, out_file, fields=None, group='/', replace=False, nocompression=False, attributes={}, fill_value=None, **kwargs):
        """
        write a grid data object to an hdf5 file
        """
        kwargs.setdefault('srs_proj4', None)
        kwargs.setdefault('srs_wkt', None)
        kwargs.setdefault('srs_epsg', None)
        kwargs.setdefault('grid_mapping_name', 'crs')
        # check whether overwriting existing files
        # append to existing files as default
        mode = 'w' if replace else 'a'

        if fields is None:
            fields=self.fields
        if group[0] != '/':
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
                # if field exists, overwrite it
                if field in h5f[group]:
                    if hasattr(self, field):
                        h5f[group+'/'+field][...] = getattr(self, field)
                else:
                    #Otherwise, try to create the dataset
                    try:
                        if nocompression or field in ['x','y','time']:
                            h5f.create_dataset(group+'/'+field, data=getattr(self, field))
                        else:

                            h5f.create_dataset(group+'/'+field, data=getattr(self, field),
                                chunks=True, compression="gzip", fillvalue=self.fill_value)
                    except Exception:
                         pass
                # try adding field attributes
                try:
                    for att_name,att_val in attributes[field].items():
                        h5f[group+'/'+field].attrs[att_name] = att_val
                except Exception:
                    pass
            # add crs attributes if applicable
            if self.crs:
                # add grid mapping attribute to each grid field
                for field in fields:
                    h5f[group+'/'+field].attrs['grid_mapping'] = kwargs['grid_mapping_name']
                # add grid mapping variable with projection attributes
                h5crs = h5f.create_dataset(kwargs['grid_mapping_name'], (), dtype=np.byte)
                for att_name,att_val in self.crs.items():
                    h5crs.attrs[att_name] = att_val

    def to_nc(self, out_file, fields=None, group='', replace=False, nocompression=False, attributes={}, fill_value=None, **kwargs):
        """
        write a grid data object to a netCDF4 file
        """
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
        except Exception as e:
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
                        ncv = ncf.createVariable(field, var.dtype, (field,))
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
                    except Exception as e:
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
        write a grid object to a geotif

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
        out_ds=self.to_gdal(out_file=out_file, driver='GTiff',**kwargs)
        return out_ds

    def to_gdal(self, driver='MEM', out_file='', field='z', srs_proj4=None, srs_wkt=None, srs_epsg=None, dtype=gdal.GDT_Float32, options=["compress=LZW"]):
        """
        Write a grid object to a gdal memory object

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
        """
        Return a pointCollection.data object containing the points in the grid
        """
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
        Get the latitude and longitude of grid cells

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
        """
        Get the date from the filename of a Worldview file
        """
        if date_format is None or date_format == 'year':
            self.time = WV_date.WV_year(self.filename)
        elif date_format == 'matlab':
            self.time = WV_date.WV_MatlabDate(self.filename)
        return self

    def normalize(self, field='z', z0=[0., 255.], z1=[0., 1.], truncate=True, dtype=np.float64):
        """
        Normalize the z range.
        """
        getattr(self, field)[:] = (getattr(self, field).astype(np.float64)-z0[0])/(z0[1]-z0[0])*(z1[1]-z1[0])+z1[0]
        if truncate:
            getattr(self, field)[getattr(self, field) < z1[0]] = z1[0]
            getattr(self, field)[getattr(self, field) > z1[1]] = z1[1]
        setattr(self, field, getattr(self, field).astype(dtype))
        return self

    def calc_gradient(self, field='z'):
        """
        calculate the gradient of a field

        Parameters
        ----------
        field : TYPE, optional
            DESCRIPTION. The default is 'z'.

        Returns
        -------
        None.
        """
        gy, gx=np.gradient(getattr(self, field), self.y, self.x)
        self.assign({field+'_x':gx, field+'_y':gy})

    def toRGB(self, cmap, field='z', caxis=None, alpha=None):
        """
        Convert a field to RGB
        """
        if caxis is None:
            caxis=[getattr(self, field).min(), getattr(self, field).max()]
        self.normalize(z0=caxis)
        setattr(self, field, cmap(getattr(self, field)))
        if alpha is not None:
            self.add_alpha_band(alpha)
        return self

    def index(self, row_ind, col_ind, fields=None, band_ind=None):
        """
        slice a grid by row or column
        """
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
        if fields is None:
            fields=self.fields
        if len(rc_ind) > 2:
            if self.t_axis==2:
                return self.copy(fields=fields).index(rc_ind[0], rc_ind[1], band_ind=rc_ind[2])
            else:
                return self.copy(fields=fields).index(rc_ind[1], rc_ind[2], band_ind=rc_ind[0])
        return self.copy(fields=fields).index(rc_ind[0], rc_ind[1], band_ind=band_ind)

    def crop(self, XR, YR, TR=None, fields=None):
        """
        Return a subset of a grid by x and y range
        """
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

    def show(self, field='z', band=None, ax=None, xy_scale=1, gradient=False, stretch_pct=None, **kwargs):
        import matplotlib.pyplot as plt
        kwargs['extent']=np.array(self.extent)*xy_scale
        kwargs['origin']='lower'
        if ax is None:
            ax = plt.gca()
        if band is None:
            zz=getattr(self, field)
        elif (band is not None) and (self.t_axis==2):
            zz=getattr(self, field)[:,:,band]
        elif (band is not None) and (self.t_axis==0):
            zz=getattr(self, field)[band,:,:]

        if gradient:
            zz=np.gradient(zz.squeeze(), self.x[1]-self.x[0], self.y[1]-self.y[0])[0]
            if 'stretch_pct' not in kwargs and 'clim' not in kwargs:
                stretch_pct=[5, 95]
                print(stretch_pct)
            if 'cmap' not in kwargs:
                kwargs['cmap']='gray'
        if stretch_pct is not None:
            LH=scoreatpercentile(zz.ravel()[np.isfinite(zz.ravel())], stretch_pct)
            kwargs['vmin']=LH[0]
            kwargs['vmax']=LH[1]

        h_im = ax.imshow(zz, **kwargs)
        return h_im

    def interp(self, x, y, gridded=False, band=0, field='z', replace=False):
        """
        interpolate a 2-D grid to a set of x and y points
        """
        if (field not in self.interpolator) or replace:
            if (len(getattr(self, field).shape) > 2) and (self.t_axis==2):
                z0 = getattr(self, field)[:,:,band].copy()
            elif (len(getattr(self, field).shape) > 2) and (self.t_axis==0):
                z0 = getattr(self, field)[band,:,:].copy()
            else:
                z0 = getattr(self, field).copy()
            # find where invalid
            NaN_mask =  (np.isfinite(z0)==0) | (z0 == self.fill_value)
            z0[NaN_mask] = 0

            if self.y[1]> self.y[0]:
                self.interpolator[field] = RectBivariateSpline(self.y, self.x, z0, kx=1, ky=1)
                if np.any(NaN_mask.ravel()):
                    self.nan_interpolator[field] = RectBivariateSpline(self.y, self.x, NaN_mask.astype(float), kx=1, ky=1)
            else:
                self.interpolator[field] = RectBivariateSpline(self.y[::-1], self.x, z0[::-1,:], kx=1, ky=1)
                if np.any(NaN_mask.ravel()):
                    self.nan_interpolator[field] = RectBivariateSpline(self.y[::-1], self.x, NaN_mask[::-1,:].astype(float), kx=1, ky=1)

        if gridded:
            result=np.zeros((len(y), len(x)))+np.NaN
            good_x = np.flatnonzero((x >= np.min(self.x)) & (x <= np.max(self.x)))
            good_y = np.flatnonzero((y >= np.min(self.y)) & (y <= np.max(self.y)))
            if (len(good_y)>0) and (len(good_x)>0):
                good_x = slice(good_x[0], good_x[-1]+1)
                good_y = slice(good_y[0], good_y[-1]+1)

                result[good_y, good_x] = self.interpolator[field](y[good_y], x[good_x])
                if field in self.nan_interpolator:
                    to_NaN=np.ones_like(result, dtype=bool)
                    to_NaN[good_y, good_x] = self.nan_interpolator[field](y[good_y], x[good_x])
                    result[to_NaN] = self.fill_value
        else:
            result = np.zeros_like(x) + self.fill_value
            good = (x >= np.min(self.x)) & (x <= np.max(self.x)) & \
                   (y >= np.min(self.y)) & (y <= np.max(self.y))

            result[good]=self.interpolator[field].ev(y[good], x[good])
            if field in self.nan_interpolator:
                to_NaN = good
                # nan_interpolator returns nonzero for NaN points in self.z
                to_NaN[good] = self.nan_interpolator[field].ev(y[good], x[good]) != 0
                result[to_NaN] = self.fill_value
        return result

    def bounds(self, pad=0):
        """
        Return the x and y bounds of a grid
        """
        return [[np.min(self.x)-pad, np.max(self.x)+pad], [np.min(self.y)-pad, np.max(self.y)+pad]]

    def replace_invalid(self, fields=None, fill_value=np.nan):
        """
        Replace invalid values within a field with a new fill_value
        """
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

    def crs_attributes(self, srs_proj4=None, srs_wkt=None, srs_epsg=None, **kwargs):
        """
        Return a dictionary of attributes for a projection

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
        except Exception as e:
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
        except Exception as e:
            pass
        self.crs['standard_parallel'] = float(sr.GetProjParm(osr.SRS_PP_LATITUDE_OF_ORIGIN,1))
        self.crs['straight_vertical_longitude_from_pole'] = float(sr.GetProjParm(osr.SRS_PP_CENTRAL_MERIDIAN,1))
        self.crs['false_northing'] = float(sr.GetProjParm(osr.SRS_PP_FALSE_NORTHING,1))
        self.crs['false_easting'] = float(sr.GetProjParm(osr.SRS_PP_FALSE_EASTING,1))
        self.crs['scale_factor'] = float(sr.GetProjParm(osr.SRS_PP_SCALE_FACTOR,1))
        return self
