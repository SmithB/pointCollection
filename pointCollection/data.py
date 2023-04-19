# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 14:28:30 2018

@author: ben
"""
import h5py
import numpy as np
import pointCollection as pc
import pyproj

class data(object):
    np.seterr(invalid='ignore')
    def __init__(self, fields=None, SRS_proj4=None, field_dict=None, columns=0, filename=None):

        if field_dict is None:
            self.field_dict=self.__default_field_dict__()
        else:
            self.field_dict=field_dict

        if fields is None:
            fields=list()
            if field_dict is not None:
                for group in self.field_dict.keys():
                    for field in self.field_dict[group]:
                        fields.append(field)
        if isinstance(fields, dict):
            self.fields=list(fields)
        self.fields=fields
        self.SRS_proj4=SRS_proj4
        self.columns=columns
        self.shape=None
        self.size=None
        self.filename=filename

    def __repr__(self):
        out=f"{self.__class__} with shape {self.shape},"+"\n"
        out += f"with fields:"+"\n"
        out += f"{self.fields}"
        return out

    def __default_field_dict__(self):
        """
        Define the default fields that get read from the h5 file
        """
        field_dict={None:('latitude','longitude','z')}
        return field_dict

    def __copy__(self):
        other=self.copy_attrs()
        for field in self.fields:
            setattr(other, field, getattr(self, field).copy())
        return other

    def __update_size_and_shape__(self):
        """
        When data size and shape may have changed, update the size and shape atttributes
        """
        self.size=0
        self.shape=[0]
        for field in self.fields:
            try:
                temp=getattr(self, field)
                if hasattr(temp, 'size'):
                    self.size=temp.size
                    self.shape=temp.shape
                    break
            except AttributeError:
                pass
        return self

    def __getitem__(self, *args, **kwargs):
        """
        wrapper for the copy_subset() method
        """
        return self.copy_subset(*args, **kwargs)

    def copy(self):
        """
        I do not know why this is necessary, but the copy method does not work
        without it
        """
        return self.__copy__()

    def copy_attrs(self):
        """
        copy attributes to a new data object
        """
        out=data()
        for field in ['fields', 'SRS_proj4','columns']:
            temp=getattr(self, field)
            if temp is not None:
                try:
                    setattr(out, field, temp.copy())
                except AttributeError:
                    setattr(out, field, temp)
        out.size=self.size
        out.shape=self.shape
        return out

    def from_h5(self, filename, group=None, field_dict=None, index_range=None):
        """
        read a data object from a HDF5 file
        """
        self.filename=filename

        with h5py.File(filename, 'r') as h5_f:
            nan_fields=list()
            if field_dict is None:
                if group is not None:
                    # build the field dict from the group
                    if not isinstance(group, (list, tuple)):
                        group=[group]
                    for this_group in group:
                        field_dict={this_group: [key for key in h5_f[this_group].keys() \
                                                 if isinstance(h5_f[this_group][key], h5py.Dataset)]}
                else:
                    field_dict=self.field_dict
            # make a slice out of whatever was provided in index_range
            if index_range is None:
                ind=slice(None)
            else:
                ind=slice(*index_range)

            for group in field_dict.keys():
                if group == '__calc_internal__':
                    self.fields += field_dict[group]
                    continue
                # loop over unique fields in the group
                for field in list(dict().fromkeys(field_dict[group])):
                    if field not in self.fields:
                        self.fields.append(field)
                    try:
                        # define the groups
                        if group is None:
                            ds=h5_f[field]
                        else:
                            ds=h5_f[group][field]
                        # read the data:
                        if self.columns==0 or self.columns is None or (ds.ndim==1):
                            setattr(self, field, np.array(ds[ind]))
                        else:
                            setattr(self, field, np.array(ds[ind,:]))
                        # check if there is a '_FillValue' defined in the dataset
                        if '_FillValue' in ds.attrs:
                            if ds.dtype in ['float32', 'float64']:
                                bad=getattr(self, field)==ds.attrs['_FillValue']
                                getattr(self, field)[bad]=np.NaN
                    except KeyError:
                        nan_fields.append(field)
                # find the first populated field
            if len(nan_fields) > 0:
                for field in self.fields:
                    if hasattr(self, field):
                        self.shape=getattr(self, field).shape
                        break
                if self.shape is not None:
                    for field in nan_fields:
                        setattr(self, field, np.zeros(self.shape)+np.NaN)
        if '__calc_internal__' in field_dict:
            try:
                self.__internal_field_calc__(field_dict)
            except Exception as e:
                print(f"pointCollection.data(): problem with __internal_field_calc__ for file {self.filename}, exception follows:")
                print(e)
                return None
        self.__update_size_and_shape__()
        return self

    def get_xy(self, proj4_string=None, EPSG=None):
        if proj4_string is not None:
            crs=proj4_string
        elif EPSG is not None:
            crs=EPSG
        if hasattr(pyproj, 'proj'):
            xy=np.array(pyproj.proj.Proj(crs)(self.longitude, self.latitude))
        else:
            if EPSG is not None:
                xy=np.array(pyproj.Proj("+init=epsg:"+str(EPSG))(self.longitude, self.latitude))
            else:
                xy=np.array(pyproj.Proj(crs)(self.longitude, self.latitude))
        self.x=xy[0,:].reshape(self.shape)
        self.y=xy[1,:].reshape(self.shape)
        if 'x' not in self.fields:
            self.fields += ['x','y']
        return self

    def get_latlon(self, proj4_string=None, EPSG=None):
        if proj4_string is not None:
            crs=proj4_string
        elif EPSG is not None:
            crs=EPSG
        if hasattr(pyproj, 'proj'):
            lonlat=np.array(pyproj.proj.Proj(crs)(self.x, self.y, inverse=True))
        else:
            if EPSG is not None:
                lonlat=np.array(pyproj.Proj("+init=epsg:"+str(EPSG))(self.x, self.y))
            else:
                lonlat=np.array(pyproj.Proj(crs)(self.x, self.y))
        self.assign({"longitude":lonlat[0,:].reshape(self.shape), \
                     "latitude":lonlat[1,:].reshape(self.shape)})
        return self

    def from_dict(self, dd, fields=None):
        """
        build a data object from a dictionary
        """
        if fields is not None:
            self.fields=fields
        else:
            self.fields=[key for key in dd.keys()]
        #work out a default size for arrays:
        try:
            default_shape=dd[next(iter(dd))].shape
        except StopIteration:
            print("HERE!")
        for field in self.fields:
            if field in dd:
                setattr(self, field, dd[field])
            else:
                setattr(self, field, np.zeros(default_shape)+np.NaN)
        self.__update_size_and_shape__()
        return self

    def from_list(self, D_list):
        """
        build a data object from a list of other data objects
        """
        if D_list is None:
            return
        if len(self.fields)==0:
            fields=set()
            for D in D_list:
                if hasattr(D,'fields'):
                    fields=fields.union(D.fields)
            self.fields=list(fields)

        for D in D_list:
            if D is not None:
                D.complete_fields(self.fields)
        try:
            for field in self.fields:
                data_list=[];
                for this_D in D_list:
                    if this_D is None:
                        continue
                    try:
                        if self.columns>1:
                            data_list.append(getattr(this_D,field))
                        else:
                            data_list.append(getattr(this_D, field).ravel())
                    except AttributeError:
                        print("Problem with field %s" % field)
                #data_list=[getattr(this_D, field).ravel() for this_D in D_list if this_D is not None]
                if len(data_list)>0:
                    setattr(self, field, np.concatenate(data_list, 0))
                else:
                    setattr(self, field, np.zeros((0)))
        except TypeError:
            for field in self.fields:
                setattr(self, field, getattr(D_list, field))
        self.__update_size_and_shape__()
        return self

    def index(self, index):
        """
        index the fields of a data object
        """
        for field in self.fields:
            try:
                setattr(self, field, getattr(self, field)[index])
            except IndexError:
                #print("IndexError for field %s, setting to NaN" % field)
                setattr(self, field, np.zeros(self.shape)[index]+np.NaN)
        self.__update_size_and_shape__()
        return self

    def blockmedian(self, scale, field='z'):
        """
        Apply a blockmedian subsetting operation to the fields
        """
        if self.size<2:
            return self
        ind = pc.pt_blockmedian(self.x, self.y, np.float64(getattr(self, field)), scale, return_index=True)[3]
        try:
            for field in self.fields:
                temp_field=getattr(self, field)
                setattr(self, field,  0.5*temp_field[ind[:,0]] + 0.5*temp_field[ind[:,1]])
        except IndexError:
            pass
        self.__update_size_and_shape__()
        return self

    def complete_fields(self, fields):
        for field in fields:
            if field not in self.fields:
                self.assign({field:np.zeros(self.shape)+np.NaN})

    def copy_subset(self, index, by_row=False, datasets=None):
        """
        return a copy of a subset of the object
        """
        dd=dict()
        if self.columns is not None and self.columns >=1 and by_row is not None or isinstance(index, slice):
            by_row=True
        if datasets is None:
            datasets=self.fields.copy()
        if (not isinstance(index, (slice, int, np.integer, float, np.float64))) and ((len(index) == 0) or ( (index.dtype == 'bool') and np.all(index==0))):
            dd={key:np.zeros([1,0]) for key in datasets}
        else:
            for field in datasets:
                temp_field=self.__dict__[field]
                try:
                    if temp_field.size==0:
                        dd[field]=temp_field.copy()
                        continue
                    if temp_field.ndim ==1:
                        dd[field]=temp_field[index]
                    else:
                        if by_row is not None and by_row:
                            dd[field]=temp_field[index,:]
                        else:
                            dd[field]=temp_field.ravel()[index.ravel()]
                except IndexError:
                    print("IndexError")
        return self.copy_attrs().from_dict(dd, fields=datasets)

    def cropped(self, bounds, return_index=False, **kwargs):
        """
        return a copy of a subset of the object falling within specified bounds

        the bounds are specified as (XR, YR)
        """
        ind = (self.x >= bounds[0][0]) & (self.x <= bounds[0][1]) &\
              (self.y >= bounds[1][0]) & (self.y <= bounds[1][1])
        if return_index:
            return ind
        return self.copy_subset(ind, **kwargs)

    def to_h5(self, fileOut, replace=True,  group='/', extensible=True, meta_dict=None):
        """
        write a data object to an hdf5 file
        """
        # check whether overwriting existing files
        # append to existing files as default
        mode = 'w' if replace else 'a'

        h5f_out=h5py.File(fileOut,mode)

        if group is not None:
            if not group in h5f_out:
                h5f_out.create_group(group)
        field_dict={}
        if meta_dict is None:
            field_dict = {field:field for field in self.fields}
        else:
            field_dict = {}
            for field, this_md in meta_dict.items():
                if this_md['source_field'] is not None and this_md['source_field'] in self.fields:
                    field_dict[field] = this_md['source_field']
                elif field in self.fields:
                    field_dict[field] = field

        for out_field, field in field_dict.items():

            this_data=getattr(self, field)
            maxshape=this_data.shape
            if extensible:
                maxshape=list(maxshape)
                maxshape[0]=None
            # try prepending the 'group' entry of meta_dict to the output field.
            # catch the exception thrown if meta_dict is None or if the 'group'
            # entry is not there
            out_field_name=out_field
            try:
                out_field_name = meta_dict[out_field]['group'] + out_field
            except (TypeError, KeyError) as e:
                pass
            out_field_name = group + '/' +out_field
            kwargs = dict( compression="gzip", maxshape=tuple(maxshape))
            if meta_dict is not None and 'precision' in meta_dict[out_field] and meta_dict[out_field]['precision'] is not None:
                kwargs['scaleoffset']=int(meta_dict[out_field]['precision'])
            h5f_out.create_dataset(out_field_name, data=this_data,  \
                                   **kwargs)
            if meta_dict is not None and out_field in meta_dict:
                for key, val in meta_dict[out_field].items():
                    if key not in ['group','source_field','precision']:
                        h5f_out[out_field].attrs[key]=str(val).encode('utf-8')
        h5f_out.close()

    def append_to_h5(self, file, group='/', ind_fields=['x','y','time']):

        """
        Append new data to a group in an existing file on disk.  Only those
        data that are not in the original dataset are appended.
        N.B.  NOT TESTED
        """

        with h5py.File(file,'r+') as h5f:
            if group in h5f:
                M_old = np.c_[tuple([np.array(h5f[group][field]).ravel() for field in ind_fields])]
                M_new=np.c_[tuple([getattr(self, field).ravel() for field in ind_fields])]
                new_ind = new_rows(M_new, M_old)
                old_N=M_old.shape[0]
            else:
                new_ind=np.ones(self.size, dtype=bool)
            new_N=new_ind.sum()
            for field in self.fields:
                if group in h5f and field in h5f[group]:
                    h5f[group][field].resize((old_N+new_N,))
                    h5f[group][field][-new_N:]=getattr(self, field)[new_ind]
                else:
                    this_data=getattr(self, field)
                    maxshape=list(this_data.shape)
                    maxshape=list(maxshape)
                    maxshape[0]=None
                    h5f.create_dataset(group+'/'+field,data=this_data,  \
                                   compression="gzip", maxshape=tuple(maxshape))


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

    def coords(self):
        """Return the x, y, and t coordinates"""
        if 'time' in self.fields:
            return (self.y, self.x, self.time)
        elif 't' in self.fields:
            return (self.y, self.x, self.t)
        else:
            return self.y, self.x

    def bounds(self, pad=0):
        """
        Parameters
        ----------
        pad : float, int, optional
            amount by which to pad the returned bounds. The default is 0.

        Returns
        -------
        XR, YR: minimum and maximum of x and y

        """
        return np.array([np.nanmin(self.x)-pad, np.nanmax(self.x)+pad]), \
                np.array([np.nanmin(self.y)-pad, np.nanmax(self.y)+pad])

    def ravel_fields(self):
        """
        ravel all the fields in self
        """
        for field in self.fields:
            setattr(self, field, getattr(self, field).ravel())
        self.__update_size_and_shape__()
        return self

def new_rows(A, B):
    """
    Return the rows in A that are not in B

    Helper function
    Solution from https://stackoverflow.com/questions/8317022/get-intersecting-rows-across-two-2d-numpy-arrays
    """
    _, ncols=A.shape
    dtype={'names':['f{}'.format(i) for i in range(ncols)],
       'formats':ncols * [A.dtype]}
    return ~np.in1d(A.view(dtype), B.view(dtype))
