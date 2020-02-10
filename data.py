# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 14:28:30 2018

@author: ben
"""
import h5py
import numpy as np
from osgeo import osr
from .pt_blockmedian import pt_blockmedian
import re
import os
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
        for field in self.fields:
            temp=getattr(self, field)
            if hasattr(temp, 'size'):
                self.size=temp.size
                self.shape=temp.shape
                break
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
                        field_dict={this_group: [key for key in h5_f[this_group].keys()]}
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
                for field in field_dict[group]:
                    if field not in self.fields:
                        self.fields.append(field)
                    try:
                        if group is None:
                            if self.columns==0 or self.columns is None:
                                setattr(self, field, np.array(h5_f[field][ind]).transpose())
                            else:
                                setattr(self, field, np.array(h5_f[field][:,ind]).transpose())
                        else:
                            if self.columns==0 or self.columns is None:
                                setattr(self, field, np.array(h5_f[group][field][ind]).transpose())
                            else:
                                setattr(self, field, np.array(h5_f[group][field][:,ind]).transpose())
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
            self.__internal_field_calc__(field_dict)
        self.__update_size_and_shape__()
        return self

    def from_qi(self, filename, field_dict=None, index_range=None):
        """
        read a data object from an ATM qi file
        """
        self.filename=filename
        # read the input file to get file information
        fd = os.open(os.path.expanduser(filename),os.O_RDONLY)
        file_info = os.fstat(fd)
        # open the filename in binary read mode
        fid = os.fdopen(fd, 'rb')
        # regular expression pattern for extracting parameters
        rx = re.compile(('(BLATM1B|ILATM1B|ILNSA1B)_((\d{4})|(\d{2}))'
        	'(\d{2})(\d{2})(.*?)\.qi$'),re.VERBOSE)
        # extract mission and other parameters from filename
        match_object = rx.match(os.path.basename(filename))
        # convert year, month and day to float variables
        year = np.float(match_object.group(2))
        month = np.float(match_object.group(5))
        day = np.float(match_object.group(6))
        # early ATM date strings omitted century and millenia (e.g. 93 for 1993)
        if match_object.group(4):
        	year = (year + 1900.0) if (year >= 90) else (year + 2000.0)

        # get the number of variables and the endianness of the file
        # assume initially big endian (all input data 32-bit integers)
        dtype = np.dtype('>i4')
        value, = np.fromfile(fid, dtype=dtype, count=1)
        fid.seek(0)
        # swap to little endian and reread first line
        if (value > 100):
            dtype = np.dtype('<i4')
            value, = np.fromfile(fid, dtype=dtype, count=1)
            fid.seek(0)
        # get the number of variables
        n_blocks = value//dtype.itemsize
        # read past first record
        np.fromfile(fid, dtype=dtype, count=n_blocks)

        # check that the number of blocks per record is less than MAXARG
        MAXARG = 14
        if (n_blocks > MAXARG):
        	raise Exception('ERROR: Unexpected number of variables')

        # read over header text
        header_count = 0
        header_text = ''
        value = np.full((n_blocks), -1, dtype=np.int32)
        while (value[0] < 0):
            # read past first record
            line = fid.read(n_blocks*dtype.itemsize)
            value = np.fromstring(line, dtype=dtype, count=n_blocks)
            header_text += str(line[dtype.itemsize:])
            header_count += dtype.itemsize*n_blocks
        # rewind file to previous record and remove last record from header text
        fid.seek(header_count)
        header_text = header_text[:-dtype.itemsize*n_blocks]
        setattr(self, 'header', header_text.replace('\x00','').rstrip())

        # make a slice out of whatever was provided in index_range
        # number of records to read with and without input subsetter
        if index_range is None:
            # number of records to read within file (file size - header size)
            n_records = (file_info.st_size - header_count)//n_blocks//dtype.itemsize
            subsetter = None
        else:
            # number of records in subsetter
            n_records = len(index_range)
            # convert from index_range indices into binary variable indices
            subsetter = header_count + dtype.itemsize*(np.array(*index_range)*n_blocks)

        # 10 word format = 0
        # 12 word format = 1
        # 14 word format = 2
        w = (n_blocks - 10)//2
        # scaling factors for each variable for the 3 word formats (14 max)
        scaling_table = [
            [1e3, 1e6, 1e6, 1e3, 1, 1, 1e3, 1e3, 1e3, 1e3],
            [1e3, 1e6, 1e6, 1e3, 1, 1, 1e3, 1e3, 1e3, 1.0e1, 1, 1e3],
            [1e3, 1e6, 1e6, 1e3, 1, 1, 1e3, 1e3, 1e3, 1, 1e6, 1e6, 1e3, 1e3]]
        # input variable names for the 3 word formats (14 max)
        variable_table = []
        # 10 word format
        variable_table.append(['rel_time','latitude','longitude','elevation',
            'xmt_sigstr','rcv_sigstr','azimuth','pitch','roll','time_hhmmss'])
        # 12 word format
        variable_table.append(['rel_time','latitude','longitude','elevation',
            'xmt_sigstr','rcv_sigstr','azimuth','pitch','roll',
            'gps_pdop','pulse_width','time_hhmmss'])
        # 14 word format
        variable_table.append(['rel_time','latitude','longitude','elevation',
            'xmt_sigstr','rcv_sigstr','azimuth','pitch','roll','passive_sig',
            'pass_foot_lat','pass_foot_long','pass_foot_synth_elev','time_hhmmss'])
        # input variable data types for the 3 word formats (14 max)
        dtype_table = []
        # 10 word format
        dtype_table.append(['f','f','f','f','i','i','f','f','f','f'])
        # 12 word format
        dtype_table.append(['f','f','f','f','i','i','f','f','f','f','i','f'])
        # 14 word format
        dtype_table.append(['f','f','f','f','i','i','f','f','f','i','f','f','f','f'])

        # dictionary with output variables
        ATM_L1b_input = {}
        for n,d in zip(variable_table[w],dtype_table[w]):
            ATM_L1b_input[n] = np.zeros((n_records), dtype=np.dtype(d))
        # extract hour, minute and second from time_hhmmss
        hour = np.zeros((n_records),dtype=np.float)
        minute = np.zeros((n_records),dtype=np.float)
        second = np.zeros((n_records),dtype=np.float)
        # for each record in the ATM Level-1b file
        for r in range(n_records):
            # set binary to point if using input subsetter from index_range
            if subsetter is not None:
                fid.seek(subsetter[r])
            # input data record r
            i = np.fromfile(fid,dtype=dtype,count=n_blocks)
            # read variable and scale to output format
            for v,n,d,s in zip(i,variable_table[w],dtype_table[w],scaling_table[w]):
                ATM_L1b_input[n][r] = v.astype(d)/s
            # unpack GPS time
            time_hhmmss = '{0:010.3f}'.format(ATM_L1b_input['time_hhmmss'][r])
            hour[r] = np.float(time_hhmmss[:2])
            minute[r] = np.float(time_hhmmss[2:4])
            second[r] = np.float(time_hhmmss[4:])
        # close the input file
        fid.close()

        # calculate dates differently than the HDF5 version 2 data as file names
        # may not contain times and version 1 data is not in UTC (no leap seconds)

        # calculate GPS time of ATM data (seconds since Jan 6, 1980 00:00:00)
        GPS_Time = self.calc_GPS_time(year,month,day,hour,minute,second)
        # leap seconds for converting from GPS time to UTC time
        leap_seconds = self.count_leap_seconds(GPS_Time)
        # calculate dates as J2000 days (UTC)
        ATM_L1b_input['days_J2k'] = (GPS_Time - leap_seconds)/86400.0 - 7300.0

        # if exporting all variable keys from the qi file
        if field_dict is None:
            field_dict = {key:None for key in ATM_L1b_input.keys()}
        # extract fields of interest using field dict keys
        for field in field_dict.keys():
            if field not in self.fields:
                self.fields.append(field)
            setattr(self, field, ATM_L1b_input[field])
        self.__update_size_and_shape__()
        # return the data and header text
        return self

    def calc_GPS_time(self, year, month, day, hour, minute, second):
        """
        Calculate the GPS time (seconds since Jan 6, 1980 00:00:00)
        """
        GPS = 367.*year - np.floor(7.*(year + np.floor((month+9.)/12.))/4.) - \
            np.floor(3.*(np.floor((year + (month - 9.)/7.)/100.) + 1.)/4.) + \
            np.floor(275.*month/9.) + day + 1721028.5 - 2444244.5
        GPS_Time = GPS*86400.0 + hour*3600.0 + minute*60.0 + second
        return GPS_Time

    def count_leap_seconds(self, GPS_Time):
        """
        Count number of leap seconds that have passed for given GPS times
        """
        # GPS times for leap seconds
        leaps = [46828800, 78364801, 109900802, 173059203, 252028804, 315187205,
            346723206, 393984007, 425520008, 457056009, 504489610, 551750411,
            599184012, 820108813, 914803214, 1025136015, 1119744016, 1167264017]
        # number of leap seconds prior to GPS_Time
        n_leaps = np.zeros_like(GPS_Time)
        for i,leap in enumerate(leaps):
            count = np.count_nonzero(GPS_Time >= leap)
            if (count > 0):
                indices, = np.nonzero(GPS_Time >= leap)
                n_leaps[indices] += 1.0
        return n_leaps

    # def get_xy(self, proj4_string=None, EPSG=None):
    #     """
    #     get projected coordinates for the data.  Adds 'x' and 'y' fields to the data, optionally returns 'self'
    #     """
    #     out_srs=osr.SpatialReference()
    #     if proj4_string is None and EPSG is not None:
    #         out_srs.ImportFromEPSG(EPSG)
    #     else:
    #         errCode=out_srs.ImportFromProj4(proj4_string)
    #         if errCode > 0:
    #             errCode=out_srs.ImportFromWkt(proj4_string)
    #     ll_srs=osr.SpatialReference()
    #     ll_srs.ImportFromEPSG(4326)
    #     ct=osr.CoordinateTransformation(ll_srs, out_srs)
    #     if self.latitude.size==0:
    #         self.x=np.zeros_like(self.latitude)
    #         self.y=np.zeros_like(self.latitude)
    #     else:
    #         xy=np.array(ct.TransformPoints(np.c_[self.longitude.ravel(), self.latitude.ravel()]))
    #         self.x=np.reshape(xy[:,0], self.latitude.shape)
    #         self.y=np.reshape(xy[:,1], self.longitude.shape)
    #     if 'x' not in self.fields:
    #         self.fields += ['x','y']
    #     return self

    def get_xy(self, proj4_string=None, EPSG=None):
        if proj4_string is not None:
            crs=proj4_string
        elif EPSG is not None:
            crs=EPSG
        xy=np.array(pyproj.proj.Proj(crs)(self.longitude, self.latitude))
        self.x=xy[0,:].reshape(self.shape)
        self.y=xy[1,:].reshape(self.shape)
        if 'x' not in self.fields:
            self.fields += ['x','y']
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
        default_shape=dd[next(iter(dd))].shape
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
                        data_list.append(getattr(this_D, field).ravel())
                    except AttributeError:
                        print("Problem with field %s" % field)
                #data_list=[getattr(this_D, field).ravel() for this_D in D_list if this_D is not None]
                setattr(self, field, np.concatenate(data_list, 0))
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
        ind=pt_blockmedian(self.x, self.y, np.float64(getattr(self, field)), scale, return_index=True)[3]
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
        if self.columns is not None and self.columns >=1 and by_row is not None:
            by_row=True
        if datasets is None:
            datasets=self.fields.copy()
        if (not isinstance(index, slice)) and ((len(index) == 0) or ( (index.dtype == 'bool') and np.all(index==0))):
            dd={key:np.zeros([1,0]) for key in datasets}
        else:
            for field in datasets:
                temp_field=self.__dict__[field]
                try:
                    if temp_field.size==0 :
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

    def to_h5(self, fileOut, replace=True, group='/'):
        """
        write a data object to an hdf5 file
        """
        if replace or not os.path.isfile(fileOut):
            if os.path.isfile(fileOut):
                os.remove(fileOut)
            h5f_out=h5py.File(fileOut,'w')
        else:
            h5f_out=h5py.File(fileOut,'r+')
        if group is not None:
            if not group in h5f_out:
                h5f_out.create_group(group)
        for field in self.fields:
            h5f_out.create_dataset(group+'/'+field,data=getattr(self,field),  compression="gzip")
        h5f_out.close()

    def assign(self,d):
        for key in d.keys():
            if key not in self.fields:
                self.fields.append(key)
            setattr(self, key, d[key])
        return self

    def coords(self):
        if 'time' in self.fields:
            return (self.y, self.x, self.time)
        else:
            return self.y, self.x

    def ravel_fields(self):
        """
        ravel all the fields in self
        """
        for field in self.fields:
            setattr(self, field, getattr(self, field).ravel())
        self.__update_size_and_shape__()
        return self
