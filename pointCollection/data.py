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
    """
    pointCollection.data objects are container objects holding point-like data.
    They can be read from hdf5 files, or constructed from dictionaries using the
    from_dict() method or using the data_dict argument to the __init__() method,
    and new data can be assigned to an existing object using the assign() mmethod.

    The data in a pointCollection.data object is stored as a set of fields that
    are listed in the object's fields attribute.  The object's methods will only
    be applied to these fields, which means that all data should be assigned
    during object creation or using the object's assign() method.  Each field
    is assumed to be a numpy.array object.

    Each object has size and shape attributes analagous to those of a numpy array.
    By default, they will match the first field listed in self.fields.

    Indexing a pointCollection.data with square brackets returns:
        -a field's data if the index is a string
        -a copy of the object sliced by the index otherwise
    """
    np.seterr(invalid='ignore')
    def __init__(self, data_dict=None, fields=None, SRS_proj4=None, EPSG=None, field_dict=None, columns=0, filename=None):
        """
        Initialize a new pointCollection.data object

        Parameters
        ----------
        data_dict : dict, optional
            Dictionary containing fields to be assigned to the new object. The default is None.
        fields : iterable, optional
            list of fields in the new object. The default is None.
        SRS_proj4 : str, optional
            proj4 string for the object's coordinate system. The default is None.
        EPSG : int, optional
            EPSG string for the object's coordinate system. The default is None.
        field_dict : dict, optional
            dictionary specifying how the object is to be read from an hdf5 file. The default is None.
        columns : int optional
            number of columns in each field in the object.  Fields are assumed to be
            of shape [Ndata, columns]. The default is 0, specifying a shape of [Ndata]
        filename : str, optional
            a filename indicating the source of the data. The default is None.

        Returns
        -------
        None.

        """
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
        self.EPSG=EPSG
        self.columns=columns
        self.shape=None
        self.size=None
        self.filename=filename
        if data_dict is not None:
            self.assign(data_dict)
        self.attrs={}
    def __repr__(self):
        out=f"{self.__class__} with shape {self.shape},"+"\n"
        out += "with fields:"+"\n"
        out += f"{self.fields}"
        return out

    def __default_field_dict__(self):
        """
        Define the default fields that get read from the h5 file
        """
        field_dict=None
        #old: {None:('latitude','longitude','z')}
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
        return a field or a subset of self

        If the input argument is a single string, return the field corresponding to that string.
        Otherwise, pass args and kwargs to the copy_subset() method
        """
        if isinstance(args[0], str):
            return getattr(self, args[0])
        return self.copy_subset(*args, **kwargs)

    def __setitem__(self, *args):
        self.assign({args[0]:args[1]})

    def __add__(self, x):
        return pc.data().from_list([self, x])

    def choose_crs(self, *args, proj4_string=None, EPSG=None, SRS_proj4=None, SRS_EPSG=None):
        '''
        Choose a coordinate system based on keyword arguments

        Parameters
        ----------
        *args: iterable, optional
            list of input arguments
        proj4_string: str, optional
            proj4 string for the destination coordinate system
        SRS_proj4:
            synonym for proj4_string
        EPSG: int, optional
            EPSG number for the destination coordinate system
        SRS_EPSG: int, optional
            synonym for EPSG

        Returns
        -------
        crs : str or int
            coordinate system to be passed to proj4.

        '''
        if len(args)>0 and args[0] is not None:
            return args[0]

        if SRS_proj4 is not None:
            proj4_string=SRS_proj4
        if SRS_EPSG is not None:
            EPSG=SRS_EPSG
        if proj4_string is not None:
            crs=proj4_string
        elif EPSG is not None:
            crs=EPSG
        else:
            if self.SRS_proj4 is None:
                crs=self.SRS_proj4
            else:
                crs=self.EPSG
        return crs

    def copy(self):
        """
        return a copy of self
        """
        return self.__copy__()

    def copy_attrs(self):
        """
        copy attributes to a new data object
        """
        out=data()
        for field in ['fields', 'SRS_proj4', 'EPSG', 'columns']:
            try:
                temp=getattr(self, field)
            except AttributeError:
                continue
            if temp is not None:
                try:
                    setattr(out, field, temp.copy())
                except AttributeError:
                    setattr(out, field, temp)
        out.size=self.size
        out.shape=self.shape
        return out

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

        strlen = len('field')
        for field in self.fields:
            strlen = np.maximum(len(field), strlen)
        strlen += 1
        fmt = '%'+str(strlen)+'s'

        table = [fmt %'field'+' \tshape \tn_finite \tSTD \t minmax']

        for field in self.fields:
            ff=getattr(self, field)
            finfo={'shape':ff.shape,
                            'n_finite':np.sum(np.isfinite(ff)),
                            'std':np.nanstd(ff),
                            'minmax':[np.nanmin(ff), np.nanmax(ff)]}

            table += [fmt % field +f' \t{finfo["shape"]}\t{finfo["n_finite"]}\t{finfo["std"]:0.2e}\t{finfo["minmax"][0]:0.2e} {finfo["minmax"][1]:0.2e}']
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

    def from_h5(self, filename, group=None, fields=None, field=None, field_dict=None, index_range=None):
        """
        read a data object from an HDF5 file

        Parameters:

        filename: str
            Filename to be read
        group: str, optional
            Group within the file to be read. If unspecified, will be treated as '/',
                unless the data object has a __default_field_dict__() function that does
                not return None
        fields: list, optional
            fields within the group to be read.  If None, all are read.
        field: str, optional
            field within the group to be read, equivalent to fields=['field']
        field_dict: dict, optional
            Dictionary specifying how fields are to be read.  Each key in field_dict specfies a
                group (or subgroup) within the file, and each value specifies a list of fields
                to be read.  If no such field appears in the file, the output will be filled with
                numpy.NaN values.  If a group is specifled as __calc_internal__, after the rest
                of the file has been read, the data object's __calc_internal__ function will be
                called to fill in the missing fields.

        """
        if filename is None:
            filename=self.filename
        else:
            self.filename=filename

        if fields is None and fields is not None:
            fields = [field]
        if fields is not None:
            field_dict = {group:fields}

        with h5py.File(filename, 'r') as h5_f:
            nan_fields=list()
            if field_dict is None:
                if group is None:
                    if self.field_dict is None:
                        group='/'
                if group is not None:
                    # build the field dict from the group
                    if not isinstance(group, (list, tuple)):
                        group=[group]
                    for this_group in group:
                        field_dict={this_group: [key for key in h5_f[this_group].keys() \
                                                 if isinstance(h5_f[this_group][key], h5py.Dataset)]}
                else:
                    if self.field_dict is not None:
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
                        if 'scale_factor' in ds.attrs:
                            scale_factor = float(ds.attrs['scale_factor'])
                        else:
                            scale_factor=1
                        if 'add_offset' in ds.attrs:
                            add_offset = float(ds.attrs['add_offset'])
                        else:
                            add_offset=0
                        if self.columns==0 or self.columns is None or (ds.ndim==1):
                            temp=np.array(ds[ind])
                        else:
                            temp=np.array(ds[ind,:])
                        # check if there is a '_FillValue' defined in the dataset
                        if '_FillValue' in ds.attrs:
                            if ds.dtype in ['float32', 'float64']:
                                bad = temp==ds.attrs['_FillValue']
                                temp[bad]=np.nan
                        if scale_factor != 1:
                            temp = temp.astype(float)*scale_factor
                        if add_offset != 0:
                            temp = temp.astype(float)+add_offset
                        setattr(self, field, temp)
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
                        setattr(self, field, np.zeros(self.shape)+np.nan)
        if '__calc_internal__' in field_dict:
            try:
                self.__internal_field_calc__(field_dict)
            except Exception as e:
                print(f"pointCollection.data(): problem with __internal_field_calc__ for file {self.filename}, exception follows:")
                print(e)
                return None
        self.__update_size_and_shape__()
        return self

    def get_xy(self, *args, **kwargs):

        '''
        Calculate projected coordinates from latitude and longitude fields

        x and y fields are calculated based on an object's latitude and
        longitude fields.

        Parametes
        *args: iterable
            list of arguments passed to self.choose_crs()
        **kwargs: dict
            keyword=pair arguments passed to self.choose_crs()
        Returns
        -------
        self: pointCollection.data
            Current object with updated x and y fields
        '''
        crs=self.choose_crs(*args,**kwargs)


        try:
            # this is compatible with pyproj 3.7.0
            # Note that EPSG:4326 takes latitude first
            xy=np.array(
                pyproj.Transformer.from_crs(pyproj.CRS(4326), pyproj.CRS(crs))\
                    .transform(self.latitude, self.longitude)
                )
        except Exception:
            if hasattr(pyproj, 'proj'):
                xy=np.array(pyproj.proj.Proj(crs)(self.longitude, self.latitude))
            else:
                if isinstance(crs, int):
                    xy=np.array(pyproj.Proj("+init=epsg:"+str(crs))(self.longitude, self.latitude))
                else:
                    xy=np.array(pyproj.Proj(crs)(self.longitude, self.latitude))
        self.x=xy[0,:].reshape(self.shape)
        self.y=xy[1,:].reshape(self.shape)
        if 'x' not in self.fields:
            self.fields += ['x','y']
        return self

    def get_latlon(self, *args, **kwargs):
        '''
        Calculate geographic coordinate fields from x and y fields

        latitude and longitude fields are calculated based on an object's x and
        y fields

        Parameters
        ----------
        *args: iterable
            list of arguments passed to self.choose_crs()
        **kwargs: dict
            keyword=pair arguments passed to self.choose_crs()

        Returns
        -------
        self: pointCollection.data
            Current object with updated latitude and longitude fields

        '''
        crs=self.choose_crs(*args,**kwargs)
        try:
            # Compatible with pyproj 3.7.0
            # Note that EPSG:4326 takes latitude first
            latlon = np.array(pyproj.Transformer.from_crs(pyproj.CRS(crs),
                                                 pyproj.CRS(4326)).transform(self.x, self.y))
            self.assign({"longitude":latlon[1,:].reshape(self.shape),\
                         "latitude":latlon[0,:].reshape(self.shape)})
        except Exception:
            if hasattr(pyproj, 'proj'):
                lonlat=np.array(pyproj.proj.Proj(crs)(self.x, self.y, inverse=True))
            else:
                if isinstance(crs, int):
                    lonlat=np.array(pyproj.Proj("+init=epsg:"+str(crs))(self.x, self.y))
                else:
                    lonlat=np.array(pyproj.Proj(crs)(self.x, self.y))
            self.assign({"longitude":lonlat[0,:].reshape(self.shape), \
                                 "latitude":lonlat[1,:].reshape(self.shape)})
        return self

    def from_dict(self, dd, fields=None):
        """
        create a pointCollection.data object from the entries in a dict

        Parameters
        ----------
        dd : dict
            Dictionary containing data values for fieldnames.  Each key specifies
            a fieldname, each value contains the data associated with that field
        fields : iterable, optional
            If specified, only these fields are read from dd. The default is None.

        Returns
        -------
        pointCollection.data
            Data object containing the values from the input dictionary

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
                setattr(self, field, np.zeros(default_shape)+np.nan)
        self.__update_size_and_shape__()
        return self

    def from_list(self, D_list):
        """
        Build a data object from a list of other data objects

        Constructs a pointCollection.data objects from an input list of objects.
        The constructed object will contain fields matching the union of all
        fields within the input objects.  If not every object contains every
        field, the missing values will be filled with numpy.NaN values.  Each
        field represents the concatenation of the input fields along its
        first axis.

        Parameters
        ----------
        D_list : iterable
        list of pointCollection.data objects.


        Returns
        -------
        pointCollection.data
            Object containing concatenated fields from input objects

        """

        if D_list is None:
            return
        if len(self.fields)==0:
            fields=[]
            if len(D_list) > 0:
                for D in D_list:
                    if hasattr(D,'fields'):
                        if not D.fields == fields:
                            for ff in D.fields:
                                if ff not in fields:
                                    fields += [ff]
            self.fields=fields
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

    def as_dict(self, fields=None):
        '''
        Return a dictionary containing the data from self.

        Parameters
        ----------
        fields : iterable, optional
            fields to include in the dictionary. The default is None.

        Returns
        -------
        dict
            dictionary containing the data from self

        '''
        if fields is None:
            fields=self.fields
        return {field:getattr(self, field) for field in fields}

    def index(self, index):
        """
        subset a data object in place

        Applies index to each field in self, assuming numpy rules for indexing
        apply.  If more control is needed over the indexing results, the
        copy_subset function contains more options.

        Parameters
        ----------
        index : bool, tuple, numpy array of ints, or slice
             Index used to subset self

        Returns
        -------
        pointCollection.data
            subsetted version of self.

        """
        for field in self.fields:
            try:
                setattr(self, field, getattr(self, field)[index])
            except IndexError:
                #print("IndexError for field %s, setting to NaN" % field)
                setattr(self, field, np.zeros(self.shape)[index]+np.nan)
        self.__update_size_and_shape__()
        return self

    def blockmedian(self, scale, field='z'):
        """
        Apply a blockmedian subsetting operation to an object's fields.


        Divides the points within an object into groups of size scale x scale,
        identifies the elements within each group representing the median of
        the specified field, and returns the values for all other fields in the
        object corresponding to these elements.  For groups containing an even
        number of elements, the average of the elements spanning the median is
        returned.

        Parameters
        ----------
        scale : float
            Scale over which the blockmedian is calculated.
        field : str, optional
            Points within the object are selected based on the values in this
            field. The default is 'z'.

        Returns
        -------
        pointCollection.data
            Object containing the meidan values

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
        """
        Fill in missing fields with np.nan


        Parameters
        ----------
        fields : iterable
            Strings containing fields to be filled

        Returns
        -------
        None.

        """

        for field in fields:
            if field not in self.fields:
                self.assign({field:np.zeros(self.shape)+np.nan})

    def copy_subset(self, index, by_row=False, datasets=None, fields=None):
        """
        return a subsetted version of self.

        Parameters
        ----------
        index : bool, tuple, numpy array, or slide
            Object used to subset the fields in self.
        by_row : bool, optional
            If True, slicing operation is applied to the first dimension of each field. The default is False.
        datasets : iterable, optional
            List of fields to be returned. If None, all fields are returned. The default is None.
        fields : TYPE, optional
            list of fields to be returned.  Synonym for datasets. The default is None.

        Returns
        -------
        pointCollection.data
            subsetted versino of self.

        """
        dd=dict()
        if self.columns is not None and self.columns >=1 and by_row is not None or isinstance(index, slice):
            by_row=True
        if fields is not None and datasets is None:
            datasets=fields
        if datasets is None:
            datasets=self.fields.copy()
        if isinstance(index, tuple):
            dd={key:getattr(self, key)[index] for key in datasets }
        elif (not isinstance(index, (slice, int, np.integer, float, np.float64))) and \
            ((len(index) == 0) or ( (index.dtype == 'bool') and np.all(index==0))):
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
        Return a cropped copy of an object

        Returns a copy self for which
            bounds[0][0] <= self.x <= bounds[0][1]
            and
            bounds[1][0] <= self.y <= bounds[1][1]

        Parameters
        ----------
        bounds : iterable
            A list of two iterables containing the range of x values and y values
            to include in the output.
        return_index : bool, optional
            If True, return only the indices of the points within the bounds.
            The default is False.
        **kwargs : dict
            keyword arguments that are passed to copy_subset.

        Returns
        -------
        pointCollection.data
            cropped version of the inptu data

        """

        ind = (self.x >= bounds[0][0]) & (self.x <= bounds[0][1]) &\
              (self.y >= bounds[1][0]) & (self.y <= bounds[1][1])
        if return_index:
            return ind
        return self.copy_subset(ind, **kwargs)

    def crop(self, bounds):
        """
        Crop current object in place.

        Parameters
        ----------
        bounds : iterable
            A list of two iterables containing the range of x values and y values
            to include in the output.
        return_index : bool, optional
            If True, return only the indices of the points within the bounds.
            The default is False.

        Returns
        -------
        None

        """

        self.index((self.x >= bounds[0][0]) & (self.x <= bounds[0][1]) &\
              (self.y >= bounds[1][0]) & (self.y <= bounds[1][1]))

    def to_h5(self, fileOut, replace=True,  group='/', extensible=True,
              dimension_fields = None,
              compression='gzip',
              chunks=True,
              meta_dict = None, DEBUG = False):
        """
        Write the contents of self to an hdf5 file.

        Saves the current data object to a file.  Optional parameters specify the
        group within the output file, and specify whether the file will be
        overwritten or if the data in self will be written to new fields and
        groups within the file.

        The meta_dict parameter allows specification of dataset attributes.  Special
        values include:
            'compression': overwrites the default gzip compression
            'source_field': alternative fieldname within self from which data
                will be copied
            'precision':
                value passed to the 'scaleoffset' keyword in h5py
        Other values in meta_dict are assigned as attributes to each dataset.

        Parameters
        ----------
        fileOut : str
            Filename to which the data are written.
        replace : bool, optional
            If true, any existing file is overwritten. The default is True.
        group : str, optional
            Group into which the data are written. The default is '/'.
        extensible : bool, optional
            If true, datasets are created as extensivble. The default is True.
        meta_dict : dict optional
            Dictionary containing attributes of output datasets. The default is None.
        DEBUG : bool, optional
            If True, report error messages from dataset creation

        Returns
        -------
        None.

        """
        # check whether overwriting existing files
        # append to existing files as default
        mode = 'w' if replace else 'a'


        h5f_out=h5py.File(fileOut,mode)

        if group is not None:
            if not group in h5f_out:
                h5f_out.create_group(group.encode('ascii'))
        field_dict={}
        if meta_dict is None:
            field_dict = {field:field for field in self.fields}
        else:
            field_dict = {}
            for out_field, this_md in meta_dict.items():
                if 'source_field' in this_md and this_md['source_field'] is not None and this_md['source_field'] in self.fields:
                    field_dict[out_field] = this_md['source_field']
                elif out_field in self.fields:
                    field_dict[out_field] = out_field

        if meta_dict is None:
            meta_dict = {out_field:{} for out_field in field_dict.keys()}

        # establish the coordinate fields first
        coordinate_fields=[]
        non_coordinate_fields=[]
        for out_field in field_dict.keys():
            if 'coordinate' in meta_dict[out_field] and out_field==meta_dict[out_field]['coordinate']:
                coordinate_fields += [out_field]
            else:
                non_coordinate_fields += [out_field]

        for out_field in coordinate_fields + non_coordinate_fields:
            this_data=getattr(self, field_dict[out_field])
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
                if DEBUG:
                    print(e)
            out_field_name = group + '/' +out_field
            kwargs = dict( compression=compression,
                          maxshape=tuple(maxshape))
            if meta_dict is not None and 'precision' in meta_dict[out_field] and meta_dict[out_field]['precision'] is not None:
                kwargs['scaleoffset']=int(meta_dict[out_field]['precision'])
            if 'datatype' in meta_dict[out_field]:
                dtype = meta_dict[out_field]['datatype'].lower()
                kwargs['dtype']=dtype
                if 'int' in dtype:
                    kwargs['fillvalue'] = np.iinfo(np.dtype(dtype)).max
                else:
                    kwargs['fillvalue'] = np.finfo(np.dtype(dtype)).max

            if 'fillvalue' in kwargs:
                this_data = np.nan_to_num(this_data,nan=kwargs['fillvalue'])
            if 'dtype' in kwargs:
                this_data = this_data.astype(kwargs['dtype'])

            # Create the dataset
            dset = h5f_out.create_dataset(out_field_name.encode('ASCII'),
                                data=this_data,  **kwargs)
            if out_field in meta_dict:
                for key, val in meta_dict[out_field].items():
                    if key.lower() not in ['group','source_field','precision','dimensions','coordinate']:
                        if isinstance(val, str):
                            h5f_out[out_field_name.encode('ASCII')].attrs[key] = str(val).encode('utf-8')
                        else:
                            h5f_out[out_field_name.encode('ASCII')].attrs[key] = val
            if 'coordinates' in meta_dict[out_field]:
                coords = meta_dict[out_field]['coordinates']
                if isinstance(coords, str):
                    coords = coords.split(',')
                for ind, coord in enumerate(coords):
                    dset.dims[ind].label=coord
            if 'fillvalue' in kwargs:
                dset.attrs['_FillValue'.encode('ASCII')] = kwargs['fillvalue']

        for key, val in self.attrs.items():
            if val is not None:
                h5f_out[group].attrs[key]=val

        for key in ['EPSG','SRS_proj4']:
            val=getattr(self, key)
            if val is not None:
                h5f_out[group].attrs[key] = val
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
        """
        Assign a field to self

        Method used to assign new values to an existing pointEollection.data
        object.  The inputs can either be a single dictionary containing
        fieldnames and values to be assigned, or a list of keyword=value
        pairs containing new fields and their values.

        Parameters
        ----------
        *args : list
            Dictionary containing fieldnames and values to be assigned.
        **kwargs : dict
            keyword=value pairs containing fields to be assigned.

        Returns
        -------
        pointCollection.data
            Updated data object.

        """

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
        """
        Return the coordinates of an object

        Returns the object's coordinates, in (y, x, time) order, or (y, x)
        if a 't' or 'time' field is not present.  Note that the coordinate
        order for this method is not consistent with that in the bounds()
        and crop() functions

        Returns
        -------
        tuple
            tuple of object coordinates

        """


        if 'time' in self.fields:
            return (self.y, self.x, self.time)
        elif 't' in self.fields:
            return (self.y, self.x, self.t)
        else:
            return self.y, self.x

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
        if len(self.x)==0:
            return None, None

        return np.array([np.nanmin(self.x)-pad, np.nanmax(self.x)+pad]), \
                np.array([np.nanmin(self.y)-pad, np.nanmax(self.y)+pad])

    def ravel_fields(self):
        """
        ravel all the fields in self

        Applies the numpy .ravel method to each field in self

        Parameters
        ---------
        None

        Returns
        -------
        None
        """
        for field in self.fields:
            setattr(self, field, getattr(self, field).ravel())
        self.__update_size_and_shape__()
        return self

    def map(self, fn, fields=None, copy=False, **kwargs):
        """
        apply a function to the fields of self

        Parameters
        ----------
        fn : function
        function to be applied to the fields of self.
        fields : iterable, optional
            fields to which the function is applied. If None, it is applied to
            all fields. The default is None.
        copy : bool, optional
            If True, operations are performed on a copy of self, otherwise
            they are performed in place.
        **kwargs:
            keyword args to pass to fn
        Returns
        -------
        pointCollection.data


        """
        if copy:
            this=self.copy()
        else:
            this=self

        if fields is None:
            fields=this.fields
        for field in fields:
            setattr(this, field, fn(getattr(this, field)))

        this.__update_size_and_shape__()
        return this

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
