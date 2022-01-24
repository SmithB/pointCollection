# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 10:46:21 2017

Class to read and manipulate ATL06 data.  Currently set up for Ben-style fake data, should be modified to work with the official ATL06 prodct foramt

@author: ben
"""
import numpy as np
from datetime import datetime, timedelta
import pointCollection as pc
import re
import os

class data(pc.data):
    np.seterr(invalid='ignore')

    def __default_field_dict__(self):
        return {None:['latitude','longitude','elevation'], \
                'instrument_parameters':['azimuth','rel_time'],\
                    '__calc_internal__':['days_J2k']}

    def __internal_field_calc__(self, field_dict):
        # find the date and time number in filename
        if 'days_J2k' in field_dict['__calc_internal__']:
            m=re.search(r"ATM1B.*_(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2}).*.h5",self.filename)
            # create list of date variables
            this_time=[int(m.group(ind+1)) for ind in range(6)]

            t0=datetime(*this_time[0:3]) + timedelta(hours=this_time[3], minutes=this_time[4], seconds=this_time[5])-datetime(2000, 1, 1, 0, 0, 0)
            t0=np.float64(t0.days)+np.float64(t0.seconds)/24./3600

            self.days_J2k = t0 + self.rel_time.astype(np.float64)/24./3600.

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
        rx = re.compile((r'(BLATM1B|ILATM1B|ILNSA1B)_((\d{4})|(\d{2}))(\d{2})(\d{2})(.*?)\.qi$'),re.VERBOSE)
        # extract mission and other parameters from filename
        match_object = rx.match(os.path.basename(filename))
        # convert year, month and day to float variables
        year = float(match_object.group(2))
        month = float(match_object.group(5))
        day = float(match_object.group(6))
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
        hour = np.zeros((n_records),dtype=float)
        minute = np.zeros((n_records),dtype=float)
        second = np.zeros((n_records),dtype=float)
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
            hour[r] = float(time_hhmmss[:2])
            minute[r] = float(time_hhmmss[2:4])
            second[r] = float(time_hhmmss[4:])
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
