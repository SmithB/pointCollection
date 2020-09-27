#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 15:35:13 2018

@author: ben
"""

from osgeo import gdal, gdalconst, osr
import numpy as np
import matplotlib.pyplot as plt
import h5py
import scipy.interpolate as si
from scipy.stats import scoreatpercentile
import pointCollection as pc
from . import WV_date
import os

class data(object):
    def __init__(self, fields=None):
        self.x=None
        self.y=None
        self.projection=None
        self.filename=None
        self.extent=None
        self.interpolator={}
        self.nan_interpolator={}
        self.time=None
        self.size=None
        self.shape=None
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
        for field in ['x','y','projection','filename','extent','time'] + fields:
            setattr(temp, field, getattr(self, field))
        temp.fields=fields.copy()
        temp.__update_size_and_shape__()
        return temp

    def __repr__(self):
        out=f"{self.__class__} with shape {self.shape},"+"\n"
        out += f"with fields:"+"\n"
        out += f"{self.fields}"
        return out

    def copy(self, fields=None):
        return self.__copy__()

    def __getitem__(self, *args, **kwargs):
        """
        wrapper for the copy_subset() method
        """
        return self.copy_subset(*args, **kwargs)

    def __update_extent__(self):
        self.extent=[np.min(self.x), np.max(self.x), np.min(self.y), np.max(self.y)]

    def __update_size_and_shape__(self):
        for field in ['z']+self.fields:
            try:
                self.size=getattr(self, field).size
                self.shape=getattr(self, field).shape
            except Exception:
                pass

    def from_dict(self, thedict):
        for field in thedict:
                setattr(self, field, thedict[field])
                if field not in self.fields and field not in ['x','y','time']:
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

    def from_geotif(self, file, field='z', bands=None, bounds=None, extent=None, skip=1, min_res=None, date_format=None):
        """
        Read a raster from a geotif
        """
        self.filename=file
        if date_format is not None:
            self.get_date(date_format)

        ds=gdal.Open(file, gdalconst.GA_ReadOnly)
        GT=ds.GetGeoTransform()
        
        if min_res is not None:
            skip=np.max([1, np.ceil(min_res/np.abs(GT[1]))]).astype(int)
        
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
            z.append(band.ReadAsArray(int(cols[0]), int(rows[0]), int(cols[-1]-cols[0]+1), int(rows[-1]-rows[0]+1))[::-1,:])
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
            z[bad] = np.NaN
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

    def from_h5(self, h5_file, field_mapping=None, group='/', fields=None, bounds=None, skip=1):
       """
       Read a raster from an hdf5 file
       """

       if field_mapping is None:
            field_mapping={}
       self.filename=h5_file
       dims=['x','y','t','time']
       if group[0] != '/':
            group='/'+group
       t=None
       with h5py.File(h5_file,'r') as h5f:
           x=np.array(h5f[group+'/x'])
           y=np.array(h5f[group+'/y'])
           if 't' in h5f[group]:
               t=np.array(h5f[group]['t'])
           elif 'time' in h5f[group]:
               t=np.array(h5f[group]['time'])

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
           if bounds is not None:
               cols = np.where(( x>=bounds[0][0] ) & ( x<= bounds[0][1] ))[0]
               rows = np.where(( y>=bounds[1][0] ) & ( y<= bounds[1][1] ))[0]
           else:
               rows=np.arange(y.size, dtype=int)
               cols=np.arange(x.size, dtype=int)
           if len(rows) > 0 and len(cols) > 0:
               for self_field in field_mapping:
                   f_field_name=group+'/'+field_mapping[self_field]
                   if f_field_name not in h5f:
                       continue
                   f_field=h5f[f_field_name]

                   if len(h5f[f_field_name].shape) == 2:
                       setattr(self, self_field,\
                               np.array(h5f[f_field_name][rows[0]:rows[-1]+1, cols[0]:cols[-1]+1]))
                   else:
                       setattr(self, self_field,\
                               np.array(f_field[rows[0]:rows[-1]+1, cols[0]:cols[-1]+1,:]))
                   if self_field not in self.fields:
                        self.fields.append(self_field)
               self.x=x[cols]
               self.y=y[rows]
               if t is not None:
                   self.t=t
       self.__update_extent__()
       self.__update_size_and_shape__()
       return self

    def to_h5(self, out_file, fields=None, group='/'):

        mode='w'
        if os.path.isfile(out_file):
            mode='r+'

        if fields is None:
            fields=self.fields
        if group[0] != '/':
            group='/'+group
        with h5py.File(out_file,mode) as h5f:
            try:
                h5f.create_group(group)
            except Exception:
                pass
            for field in ['x','y','time', 't'] + fields:
                try:
                    h5f.create_dataset(group+'/'+field, data=getattr(self, field))
                except Exception:
                    pass

    def to_geotif(self, out_file, field='z', srs_proj4=None, srs_wkt=None,  srs_epsg=None):
        """
        Write a grid object to a geotif.
        """

        z=getattr(self, field)

        nx=z.shape[1]
        ny=z.shape[0]
        if len(z.shape)>2:
            n_bands=z.shape[2]
        else:
            n_bands=1;
        dx=np.abs(np.diff(self.x[0:2]))[0]
        dy=np.abs(np.diff(self.y[0:2]))[0]

        out_ds=gdal.GetDriverByName('GTiff').Create(out_file, nx, ny, n_bands, gdal.GDT_Float32, options=["compress=LZW"])
        out_ds.SetGeoTransform((self.x.min()-dx/2, dx, 0, self.y.max()+dy/2, 0., -dy))
        sr=osr.SpatialReference()
        if srs_proj4 is not None:
            sr.ImportFromEPSG(srs_epsg)
        elif srs_wkt is not None:
            sr.ImportFromWKT(srs_wkt)
        elif srs_epsg is not None:
            sr.ImportFromEPSG(srs_epsg)

        else:
            raise ValueError("must specify at least one of srs_proj4, srs_wkt, srs_epsg")

        out_ds.SetProjection(sr.ExportToWkt())
        if n_bands == 1:
            out_ds.GetRasterBand(1).WriteArray(z[::-1,:])
        else:
            for band in range(n_bands):
                out_ds.GetRasterBand(band+1).WriteArray(z[::-1,:,band])
        out_ds.FlushCache()
        out_ds = None

    def as_points(self, field='z', keep_all=False):
        """
        Return a pointCollection.data object containing the points in the grid
        """
        x,y=np.meshgrid(self.x, self.y)
        if keep_all:
            result =  pc.data(filename=self.filename).\
                from_dict({'x':x.ravel(),'y':y.ravel(),'z':getattr(self, field).ravel()})
        else:
            good=np.isfinite(getattr(self, field)).ravel()
            result = pc.data(filename=self.filename).\
                from_dict({'x':x.ravel()[good],'y':y.ravel()[good],'z':getattr(self, field).ravel()[good]})
        if self.time is not None:
            result.assign({'time':self.time+np.zeros_like(getattr(result, field))})
        return result

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
        for field in fields:
            if len(getattr(self, field).shape) == 2:
                setattr(self, field, getattr(self, field)[row_ind,:][:, col_ind])
            else:
                if band_ind is None:
                    setattr(self, field, getattr(self, field)[row_ind,:, :][:, col_ind,:])
                else:
                    setattr(self, field, getattr(self, field)[row_ind,:, :][:, col_ind,band_ind])
                    self.t=self.t[band_ind]
        self.__update_extent__()
        self.__update_size_and_shape__()
        return self

    def copy_subset(self, rc_ind, band_ind=None, fields=None):
        return self.copy(fields=fields).index(rc_ind[0], rc_ind[1], band_ind=band_ind)

    def crop(self, XR, YR, fields=None):
        """
        Return a subset of a grid by x and y range
        """

        col_ind = np.where((self.x >= XR[0]) & (self.x <= XR[1]))[0]
        row_ind = np.where((self.y >= YR[0]) & (self.y <= YR[1]))[0]
        try:
           self.index(row_ind, col_ind, fields)
           return self
        except Exception as e:
           print("grid: self extent is: ", self.extent)
           print("XR is %s", XR)
           print("YR is %s", YR)
           print("Error is" )
           print(e)

    def show(self, field='z', band=None, ax=None, xy_scale=1, gradient=False, stretch_pct=None, **kwargs):
        kwargs['extent']=np.array(self.extent)*xy_scale
        kwargs['origin']='lower'
        if band is None:
            zz=getattr(self, field)
        else:
            zz=getattr(self, field)[:,:,band]

        if gradient:
            zz=np.gradient(zz.squeeze(), self.x[1]-self.x[0], self.y[1]-self.y[0])[0]
            if 'stretch_pct' not in kwargs:
                stretch_pct=[5, 95]
            if 'cmap' not in kwargs:
                kwargs['cmap']='gray'
        if stretch_pct is not None:
            LH=scoreatpercentile(zz.ravel()[np.isfinite(zz.ravel())], stretch_pct)
            kwargs['vmin']=LH[0]
            kwargs['vmax']=LH[1]

        if ax is None:
            h_im = plt.imshow(zz, **kwargs)
        else:
            h_im = ax.imshow(zz, **kwargs)
        return h_im

    def interp(self, x, y, gridded=False, band=0, field='z'):
        """
        interpolate a grid to a set of x and y points
        """
        if field not in self.interpolator:
            if len(getattr(self, field).shape) > 2:
                z0 = getattr(self, field)[:,:,band]
            else:
                z0 = getattr(self, field).copy()
            NaN_mask =  np.isfinite(z0)==0
            z0[NaN_mask] = 0

            if self.y[1]> self.y[0]:
                self.interpolator[field] = si.RectBivariateSpline(self.y, self.x, z0)
                if np.any(NaN_mask.ravel()):
                    self.nan_interpolator[field] = si.RectBivariateSpline(self.y, self.x, NaN_mask.astype(float), kx=1, ky=1)
            else:
                self.interpolator[field] = si.RectBivariateSpline(self.y[::-1], self.x, z0[::-1,:], kx=1, ky=1)
                if np.any(NaN_mask.ravel()):
                    self.nan_interpolator[field] = si.RectBivariateSpline(self.y[::-1], self.x, NaN_mask[::-1,:].astype(float), kx=1, ky=1)

        if gridded:
            result=np.zeros((len(y), len(x)))
            good_x = np.flatnonzero((x >= np.min(self.x)) & (x <= np.max(self.x)))
            good_y = np.flatnonzero((y >= np.min(self.y)) & (y <= np.max(self.y)))
            if (len(good_y)>0) and (len(good_x)>0):
                good_x = slice(good_x[0], good_x[-1]+1)
                good_y = slice(good_y[0], good_y[-1]+1)

                result[good_y, good_x] = self.interpolator[field](y[good_y], x[good_x])
                if field in self.nan_interpolator:
                    to_NaN=np.ones_like(result, dtype=bool)
                    to_NaN[good_y, good_x] = self.nan_interpoator(y[good_y], x[good_x])
                    result[to_NaN] = np.NaN
        else:
            result = np.zeros_like(x)+np.NaN
            good = (x >= np.min(self.x)) & (x <= np.max(self.x)) & \
                   (y >= np.min(self.y)) & (y <= np.max(self.y))

            result[good]=self.interpolator[field].ev(y[good], x[good])
        if field in self.nan_interpolator:
            to_NaN = good
            # nan_interpolator returns nonzero for NaN points in self.z
            to_NaN[good] = self.nan_interpolator[field].ev(y[good], x[good]) != 0
            result[to_NaN] = np.NaN
        return result

    def bounds(self, pad=0):
        """
        Return the x and y bounds of a grid
        """
        return [[np.min(self.x)-pad, np.max(self.x)+pad], [np.min(self.y)-pad, np.max(self.y)+pad]]
