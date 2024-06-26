# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 2020

Class to read and manipulate CryoSat-2 waveform data
Reads CryoSat Level-1b data products from baselines A, B and C
Reads CryoSat Level-1b netCDF4 data products from baseline D
Supported CryoSat Modes: LRM, SAR, SARin, FDM, SID, GDR

INPUTS:
    full_filename: full path of CryoSat .DBL or .nc file

PYTHON DEPENDENCIES:
    numpy: Scientific Computing Tools For Python
        http://www.numpy.org
        http://www.scipy.org/NumPy_for_Matlab_Users
    netCDF4: Python interface to the netCDF C library
        https://unidata.github.io/netcdf4-python/netCDF4/index.html

UPDATE HISTORY:
    Updated 08/2020: flake8 compatible binary regular expression strings
    Forked 02/2020 from read_cryosat_L1b.py
    Updated 11/2019: empty placeholder dictionary for baseline D DSD headers
    Updated 09/2019: added netCDF4 read function for baseline D
    Updated 04/2019: USO correction signed 32 bit int
    Updated 10/2018: updated header read functions for python3
    Updated 05/2016: using __future__ print and division functions
    Written 03/2016
"""
from __future__ import print_function
from __future__ import division

import numpy as np
import pointCollection as pc
import netCDF4
import re
import os

class data(pc.data):
    np.seterr(invalid='ignore')

    def __default_field_dict__(self):
        """
        Define the default fields that get read from the CryoSat-2 file
        """
        field_dict = {}
        field_dict['Location'] = ['days_J2k','Day','Second','Micsec','USO_Corr',
            'Mode_ID','SSC','Inst_config','Rec_Count','Lat','Lon','Alt','Alt_rate',
            'Sat_velocity','Real_beam','Baseline','ST_ID','Roll','Pitch','Yaw','MCD']
        field_dict['Data'] = ['TD', 'H_0','COR2','LAI','FAI','AGC_CH1','AGC_CH2',
            'TR_gain_CH1','TR_gain_CH2','TX_Power','Doppler_range','TR_inst_range',
            'R_inst_range','TR_inst_gain','R_inst_gain','Internal_phase',
            'External_phase','Noise_power','Phase_slope']
        field_dict['Geometry'] = ['dryTrop','wetTrop','InvBar','DAC','Iono_GIM',
            'Iono_model','ocTideElv','lpeTideElv','olTideElv','seTideElv','gpTideElv',
            'Surf_type','Corr_status','Corr_error']
        field_dict['Waveform_20Hz'] = ['Waveform','Linear_Wfm_Multiplier',
            'Power2_Wfm_Multiplier','N_avg_echoes']
        field_dict['METADATA'] = ['MPH','SPH']

        return field_dict

    def from_dbl(self, full_filename, field_dict=None, unpack=False, verbose=False):
        """
        Read CryoSat Level-1b data from binary formats
        """
        # file basename and file extension of input file
        fileBasename,fileExtension=os.path.splitext(os.path.basename(full_filename))

        # CryoSat file class
        # OFFL (Off Line Processing/Systematic)
        # NRT_ (Near Real Time)
        # RPRO (ReProcessing)
        # TEST (Testing)
        # TIxx (Stand alone IPF1 testing)
        # LTA_ (Long Term Archive)
        regex_class = 'OFFL|NRT_|RPRO|TEST|TIxx|LTA_'
        # CryoSat mission products
        # SIR1SAR_FR: Level 1 FBR SAR Mode (Rx1 Channel)
        # SIR2SAR_FR: Level 1 FBR SAR Mode (Rx2 Channel)
        # SIR_SIN_FR: Level 1 FBR SARin Mode
        # SIR_LRM_1B: Level-1 Product Low Rate Mode
        # SIR_FDM_1B: Level-1 Product Fast Delivery Marine Mode
        # SIR_SAR_1B: Level-1 SAR Mode
        # SIR_SIN_1B: Level-1 SARin Mode
        # SIR1LRC11B: Level-1 CAL1 Low Rate Mode (Rx1 Channel)
        # SIR2LRC11B: Level-1 CAL1 Low Rate Mode (Rx2 Channel)
        # SIR1SAC11B: Level-1 CAL1 SAR Mode (Rx1 Channel)
        # SIR2SAC11B: Level-1 CAL1 SAR Mode (Rx2 Channel)
        # SIR_SIC11B: Level-1 CAL1 SARin Mode
        # SIR_SICC1B: Level-1 CAL1 SARIN Exotic Data
        # SIR1SAC21B: Level-1 CAL2 SAR Mode (Rx1 Channel)
        # SIR2SAC21B: Level-1 CAL2 SAR Mode (Rx2 Channel)
        # SIR1SIC21B: Level-1 CAL2 SARin Mode (Rx1 Channel)
        # SIR2SIC21B: Level-1 CAL2 SARin Mode (Rx1 Channel)
        # SIR1LRM_0M: LRM and TRK Monitoring Data from Rx 1 Channel
        # SIR2LRM_0M: LRM and TRK Monitoring Data from Rx 2 Channel
        # SIR1SAR_0M: SAR Monitoring Data from Rx 1 Channel
        # SIR2SAR_0M: SAR Monitoring Data from Rx 1 Channel
        # SIR_SIN_0M: SARIN Monitoring Data
        # SIR_SIC40M: CAL4 Monitoring Data
        regex_products = ('SIR1SAR_FR|SIR2SAR_FR|SIR_SIN_FR|SIR_LRM_1B|SIR_FDM_1B|'
            'SIR_SAR_1B|SIR_SIN_1B|SIR1LRC11B|SIR2LRC11B|SIR1SAC11B|SIR2SAC11B|'
            'SIR_SIC11B|SIR_SICC1B|SIR1SAC21B|SIR2SAC21B|SIR1SIC21B|SIR2SIC21B|'
            'SIR1LRM_0M|SIR2LRM_0M|SIR1SAR_0M|SIR2SAR_0M|SIR_SIN_0M|SIR_SIC40M')
        # CRYOSAT LEVEL-1b PRODUCTS NAMING RULES
        # Mission Identifier
        # File Class
        # File Product
        # Validity Start Date and Time
        # Validity Stop Date and Time
        # Baseline Identifier
        # Version Number
        regex_pattern = r'(.*?)_({0})_({1})_(\d+T?\d+)_(\d+T?\d+)_(.*?)(\d+)'
        rx = re.compile(regex_pattern.format(regex_class,regex_products),re.VERBOSE)
        # extract file information from filename
        MI,CLASS,PRODUCT,START,STOP,BASELINE,VERSION=rx.findall(fileBasename).pop()

        # CryoSat-2 Mode record sizes
        i_size_timestamp = 12
        n_SARIN_BC_RW = 1024
        n_SARIN_RW = 512
        n_SAR_BC_RW = 256
        n_SAR_RW = 125
        n_LRM_RW = 128
        n_blocks = 20
        n_BeamBehaviourParams = 50
        # check baseline from file to set i_record_size and allocation function
        if (BASELINE == 'C'):
            # calculate total record sizes of each dataset group
            i_size_timegroup = i_size_timestamp + 4 + 2*2 + 6*4 + 3*3*4 + 3*2 + 4*4
            i_size_measuregroup = 8 + 4*17 + 8
            i_size_external_corr = 4*13 + 12
            i_size_1Hz_LRM = i_size_timestamp + 3*4 + 8 + n_LRM_RW*2 + 2*4 + 2*2
            i_size_1Hz_SAR = i_size_timestamp + 4*3 + 8 + n_SAR_RW*2 + 4 + 4 + 2 + 2
            i_size_1Hz_SARIN = i_size_timestamp + 4*3 + 8 + n_SARIN_RW*2 + 4 + 4 + 2 + 2
            i_size_LRM_waveform = n_LRM_RW*2 + 4 + 4 + 2 + 2
            i_size_SAR_waveform = n_SAR_BC_RW*2 + 4 + 4 + 2 + 2 + n_BeamBehaviourParams*2
            i_size_SARIN_waveform = n_SARIN_BC_RW*2 + 4 + 4 + 2 + 2 + n_SARIN_BC_RW*2 + \
                n_SARIN_BC_RW*4 + n_BeamBehaviourParams*2
            # Low-Resolution Mode Record Size
            i_record_size_LRM_L1b = n_blocks * (i_size_timegroup + \
                i_size_measuregroup + i_size_LRM_waveform) + i_size_external_corr + \
                i_size_1Hz_LRM
            # SAR Mode Record Size
            i_record_size_SAR_L1b = n_blocks * (i_size_timegroup + \
                i_size_measuregroup + i_size_SAR_waveform) + i_size_external_corr + \
                i_size_1Hz_SAR
            # SARIN Mode Record Size
            i_record_size_SARIN_L1b = n_blocks * (i_size_timegroup + \
                i_size_measuregroup + i_size_SARIN_waveform) + i_size_external_corr + \
                i_size_1Hz_SARIN
            # set read function for Baseline C
            read_cryosat_variables = self.cryosat_baseline_C
        else:
            # calculate total record sizes of each dataset group
            i_size_timegroup = i_size_timestamp + 4 + 2*2+ 6*4 + 3*3*4 + 4
            i_size_measuregroup = 8 + 4*17 + 8
            i_size_external_corr = 4*13 + 12
            i_size_1Hz_LRM = i_size_timestamp + 3*4 + 8 + n_LRM_RW*2 + 2*4 + 2*2
            i_size_1Hz_SAR = i_size_timestamp + 4*3 + 8 + n_SAR_RW*2 + 4 + 4 + 2 + 2
            i_size_1Hz_SARIN = i_size_timestamp + 4*3 + 8 + n_SARIN_RW*2 + 4 + 4 + 2 + 2
            i_size_LRM_waveform = n_LRM_RW*2 + 4 + 4 + 2 + 2
            i_size_SAR_waveform = n_SAR_RW*2 + 4 + 4 + 2 + 2 + n_BeamBehaviourParams*2
            i_size_SARIN_waveform = n_SARIN_RW*2 + 4 + 4 + 2 + 2 + n_SARIN_RW*2 + \
                n_SARIN_RW*4 + n_BeamBehaviourParams*2
            # Low-Resolution Mode Record Size
            i_record_size_LRM_L1b = n_blocks * (i_size_timegroup + \
                i_size_measuregroup + i_size_LRM_waveform) + i_size_external_corr + \
                i_size_1Hz_LRM
            # SAR Mode Record Size
            i_record_size_SAR_L1b = n_blocks * (i_size_timegroup + \
                i_size_measuregroup + i_size_SAR_waveform) + i_size_external_corr + \
                i_size_1Hz_SAR
            # SARIN Mode Record Size
            i_record_size_SARIN_L1b = n_blocks * (i_size_timegroup + \
                i_size_measuregroup + i_size_SARIN_waveform) + i_size_external_corr + \
                i_size_1Hz_SARIN
            # set read function for Baselines A and B
            read_cryosat_variables = self.cryosat_baseline_AB

        # get dataset MODE from PRODUCT portion of file name
        # set record sizes and DS_TYPE for read_DSD function
        self.MODE = re.findall('(LRM|SAR|SIN)', PRODUCT).pop()
        if (self.MODE == 'LRM'):
            i_record_size = i_record_size_LRM_L1b
            DS_TYPE = 'CS_L1B'
        elif (self.MODE == 'SAR'):
            i_record_size = i_record_size_SAR_L1b
            DS_TYPE = 'CS_L1B'
        elif (self.MODE == 'SIN'):
            i_record_size = i_record_size_SARIN_L1b
            DS_TYPE = 'CS_L1B'

        # read the input file to get file information
        fid = os.open(os.path.expanduser(full_filename),os.O_RDONLY)
        file_info = os.fstat(fid)
        os.close(fid)

        # num DSRs from SPH
        j_num_DSR = np.int32(file_info.st_size//i_record_size)
        # print file information
        if verbose:
            print(full_filename)
            print('{0:d} {1:d} {2:d}'.format(j_num_DSR,file_info.st_size,i_record_size))
            # Check if MPH/SPH/DSD headers
            if (j_num_DSR*i_record_size == file_info.st_size):
                print('No Header on file')
                print('The number of DSRs is: {0:d}'.format(j_num_DSR))
            else:
                print('Header on file')

        # Check if MPH/SPH/DSD headers
        if (j_num_DSR*i_record_size != file_info.st_size):
            # If there are MPH/SPH/DSD headers
            s_MPH_fields = self.read_MPH(full_filename)
            j_sph_size = np.int32(re.findall(r'[-+]?\d+',s_MPH_fields['SPH_SIZE']).pop())
            s_SPH_fields = self.read_SPH(full_filename, j_sph_size)
            # extract information from DSD fields
            s_DSD_fields = self.read_DSD(full_filename, DS_TYPE=DS_TYPE)
            # extract DS_OFFSET
            j_DS_start = np.int32(re.findall(r'[-+]?\d+',s_DSD_fields['DS_OFFSET']).pop())
            # extract number of DSR in the file
            j_num_DSR = np.int32(re.findall(r'[-+]?\d+',s_DSD_fields['NUM_DSR']).pop())
            # check the record size
            j_DSR_size = np.int32(re.findall(r'[-+]?\d+',s_DSD_fields['DSR_SIZE']).pop())
            #  minimum size is start of the read plus number of records to read
            j_check_size = j_DS_start + (j_DSR_size*j_num_DSR)
            if verbose:
                print('The offset of the DSD is: {0:d} bytes'.format(j_DS_start))
                print('The number of DSRs is {0:d}'.format(j_num_DSR))
                print('The size of the DSR is {0:d}'.format(j_DSR_size))
            # check if invalid file size
            if (j_check_size > file_info.st_size):
                raise IOError('File size error')
            # extract binary data from input CryoSat data file (skip headers)
            fid = open(os.path.expanduser(full_filename), 'rb')
            cryosat_header = fid.read(j_DS_start)
            # iterate through CryoSat file and fill output variables
            CS_L1b_mds = read_cryosat_variables(fid, j_num_DSR)
            # add headers to output dictionary as METADATA
            CS_L1b_mds['METADATA'] = {}
            CS_L1b_mds['METADATA']['MPH'] = s_MPH_fields
            CS_L1b_mds['METADATA']['SPH'] = s_SPH_fields
            CS_L1b_mds['METADATA']['DSD'] = s_DSD_fields
            # close the input CryoSat binary file
            fid.close()
        else:
            # If there are not MPH/SPH/DSD headers
            # extract binary data from input CryoSat data file
            fid = open(os.path.expanduser(full_filename), 'rb')
            # iterate through CryoSat file and fill output variables
            CS_L1b_mds = read_cryosat_variables(fid, j_num_DSR)
            # close the input CryoSat binary file
            fid.close()

        # if unpacking the units
        if unpack:
            CS_l1b_scale = self.cryosat_scaling_factors()
            # for each dictionary key
            for group in CS_l1b_scale.keys():
                # for each variable
                for key,val in CS_L1b_mds[group].items():
                    # check if val is the 20Hz waveform beam variables
                    if isinstance(val, dict):
                        # for each waveform beam variable
                        for k,v in val.items():
                            # scale variable
                            CS_L1b_mds[group][key][k] = CS_l1b_scale[group][key][k]*v.copy()
                    else:
                        # scale variable
                        CS_L1b_mds[group][key] = CS_l1b_scale[group][key]*val.copy()

        # calculate GPS time of CryoSat data (seconds since Jan 6, 1980 00:00:00)
        # from TAI time since Jan 1, 2000 00:00:00
        GPS_Time = self.calc_GPS_time(CS_L1b_mds['Location']['Day'],
            CS_L1b_mds['Location']['Second'], CS_L1b_mds['Location']['Micsec'])
        # leap seconds for converting from GPS time to UTC time
        leap_seconds = self.count_leap_seconds(GPS_Time)
        # calculate dates as J2000 days (UTC)
        CS_L1b_mds['Location']['days_J2k'] = (GPS_Time - leap_seconds)/86400.0 - 7300.0

        # parameters to extract
        if field_dict is None:
            field_dict = self.__default_field_dict__()
        # extract fields of interest using field dict keys
        for group,variables in field_dict.items():
            for field in variables:
                if field not in self.fields:
                    self.fields.append(field)
                setattr(self, field, CS_L1b_mds[group][field])

        # update size and shape of input data
        self.__update_size_and_shape__()
        # return the data and header text
        return self

    def from_nc(self, full_filename, field_dict=None, unpack=False, verbose=False):
        """
        Read CryoSat Level-1b data from netCDF4 format data
        """
        # file basename and file extension of input file
        fileBasename,fileExtension=os.path.splitext(os.path.basename(full_filename))

        # CryoSat file class
        # OFFL (Off Line Processing/Systematic)
        # NRT_ (Near Real Time)
        # RPRO (ReProcessing)
        # TEST (Testing)
        # TIxx (Stand alone IPF1 testing)
        # LTA_ (Long Term Archive)
        regex_class = 'OFFL|NRT_|RPRO|TEST|TIxx|LTA_'
        # CryoSat mission products
        # SIR1SAR_FR: Level 1 FBR SAR Mode (Rx1 Channel)
        # SIR2SAR_FR: Level 1 FBR SAR Mode (Rx2 Channel)
        # SIR_SIN_FR: Level 1 FBR SARin Mode
        # SIR_LRM_1B: Level-1 Product Low Rate Mode
        # SIR_FDM_1B: Level-1 Product Fast Delivery Marine Mode
        # SIR_SAR_1B: Level-1 SAR Mode
        # SIR_SIN_1B: Level-1 SARin Mode
        # SIR1LRC11B: Level-1 CAL1 Low Rate Mode (Rx1 Channel)
        # SIR2LRC11B: Level-1 CAL1 Low Rate Mode (Rx2 Channel)
        # SIR1SAC11B: Level-1 CAL1 SAR Mode (Rx1 Channel)
        # SIR2SAC11B: Level-1 CAL1 SAR Mode (Rx2 Channel)
        # SIR_SIC11B: Level-1 CAL1 SARin Mode
        # SIR_SICC1B: Level-1 CAL1 SARIN Exotic Data
        # SIR1SAC21B: Level-1 CAL2 SAR Mode (Rx1 Channel)
        # SIR2SAC21B: Level-1 CAL2 SAR Mode (Rx2 Channel)
        # SIR1SIC21B: Level-1 CAL2 SARin Mode (Rx1 Channel)
        # SIR2SIC21B: Level-1 CAL2 SARin Mode (Rx1 Channel)
        # SIR1LRM_0M: LRM and TRK Monitoring Data from Rx 1 Channel
        # SIR2LRM_0M: LRM and TRK Monitoring Data from Rx 2 Channel
        # SIR1SAR_0M: SAR Monitoring Data from Rx 1 Channel
        # SIR2SAR_0M: SAR Monitoring Data from Rx 1 Channel
        # SIR_SIN_0M: SARIN Monitoring Data
        # SIR_SIC40M: CAL4 Monitoring Data
        regex_products = ('SIR1SAR_FR|SIR2SAR_FR|SIR_SIN_FR|SIR_LRM_1B|SIR_FDM_1B|'
            'SIR_SAR_1B|SIR_SIN_1B|SIR1LRC11B|SIR2LRC11B|SIR1SAC11B|SIR2SAC11B|'
            'SIR_SIC11B|SIR_SICC1B|SIR1SAC21B|SIR2SAC21B|SIR1SIC21B|SIR2SIC21B|'
            'SIR1LRM_0M|SIR2LRM_0M|SIR1SAR_0M|SIR2SAR_0M|SIR_SIN_0M|SIR_SIC40M')
        # CRYOSAT LEVEL-1b PRODUCTS NAMING RULES
        # Mission Identifier
        # File Class
        # File Product
        # Validity Start Date and Time
        # Validity Stop Date and Time
        # Baseline Identifier
        # Version Number
        regex_pattern = r'(.*?)_({0})_({1})_(\d+T?\d+)_(\d+T?\d+)_(.*?)(\d+)'
        rx = re.compile(regex_pattern.format(regex_class,regex_products),re.VERBOSE)
        # extract file information from filename
        MI,CLASS,PRODUCT,START,STOP,BASELINE,VERSION=rx.findall(fileBasename).pop()
        print(full_filename) if verbose else None
        # get dataset MODE from PRODUCT portion of file name
        self.MODE = re.findall(r'(LRM|FDM|SAR|SIN)', PRODUCT).pop()
        # read level-2 CryoSat-2 data from netCDF4 file
        CS_L1b_mds = self.cryosat_baseline_D(full_filename, unpack=unpack)

        # calculate GPS time of CryoSat data (seconds since Jan 6, 1980 00:00:00)
        # from TAI time since Jan 1, 2000 00:00:00
        GPS_Time = self.calc_GPS_time(CS_L1b_mds['Location']['Day'],
            CS_L1b_mds['Location']['Second'], CS_L1b_mds['Location']['Micsec'])
        # leap seconds for converting from GPS time to UTC time
        leap_seconds = self.count_leap_seconds(GPS_Time)
        # calculate dates as J2000 days (UTC)
        CS_L1b_mds['Location']['days_J2k'] = (GPS_Time - leap_seconds)/86400.0 - 7300.0

        # parameters to extract
        if field_dict is None:
            field_dict = self.__default_field_dict__()
        # extract fields of interest using field dict keys
        for group,variables in field_dict.items():
            for field in variables:
                if field not in self.fields:
                    self.fields.append(field)
                setattr(self, field, CS_L1b_mds[group][field])

        # update size and shape of input data
        self.__update_size_and_shape__()
        # return the data and header text
        return self

    def calc_GPS_time(self, day, second, micsec):
        """
        Calculate the GPS time (seconds since Jan 6, 1980 00:00:00)
        """
        # TAI time is ahead of GPS by 19 seconds
        return (day + 7300.0)*86400.0 + second.astype('f') + micsec/1e6 - 19

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
                i_records,i_blocks = np.nonzero(GPS_Time >= leap)
                n_leaps[i_records,i_blocks] += 1.0
        return n_leaps

    def read_MPH(self, full_filename):
        """
        Read ASCII Main Product Header (MPH) block from an ESA PDS file
        """
        # read input data file
        with open(os.path.expanduser(full_filename), 'rb') as fid:
            file_contents = fid.read().splitlines()

        # Define constant values associated with PDS file formats
        # number of text lines in standard MPH
        n_MPH_lines = 41
        # check that first line of header matches PRODUCT
        if not bool(re.match(br'PRODUCT\=\"(.*)(?=\")',file_contents[0])):
            raise IOError('File does not start with a valid PDS MPH')
        # read MPH header text
        s_MPH_fields = {}
        for i in range(n_MPH_lines):
            # use regular expression operators to read headers
            if bool(re.match(br'(.*?)\=\"(.*)(?=\")',file_contents[i])):
                # data fields within quotes
                field,value=re.findall(br'(.*?)\=\"(.*)(?=\")',file_contents[i]).pop()
                s_MPH_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()
            elif bool(re.match(br'(.*?)\=(.*)',file_contents[i])):
                # data fields without quotes
                field,value=re.findall(br'(.*?)\=(.*)',file_contents[i]).pop()
                s_MPH_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()

        # Return block name array to calling function
        return s_MPH_fields

    def read_SPH(self, full_filename, j_sph_size):
        """
        Read ASCII Specific Product Header (SPH) block from a PDS file
        """
        # read input data file
        with open(os.path.expanduser(full_filename), 'rb') as fid:
            file_contents = fid.read().splitlines()

        # Define constant values associated with PDS file formats
        # number of text lines in standard MPH
        n_MPH_lines = 41
        # compile regular expression operator for reading headers
        rx = re.compile(br'(.*?)\=\"?(.*)',re.VERBOSE)
        # check first line of header matches SPH_DESCRIPTOR
        if not bool(re.match(br'SPH\_DESCRIPTOR\=',file_contents[n_MPH_lines+1])):
            raise IOError('File does not have a valid PDS DSD')
        # read SPH header text (no binary control characters)
        s_SPH_lines = [li for li in file_contents[n_MPH_lines+1:] if rx.match(li)
            and not re.search(br'[^\x20-\x7e]+',li)]

        # extract SPH header text
        s_SPH_fields = {}
        c = 0
        while (c < len(s_SPH_lines)):
            # check if line is within DS_NAME portion of SPH header
            if bool(re.match(br'DS_NAME',s_SPH_lines[c])):
                # add dictionary for DS_NAME
                field,value=re.findall(br'(.*?)\=\"(.*)(?=\")',s_SPH_lines[c]).pop()
                key = value.decode('utf-8').rstrip()
                s_SPH_fields[key] = {}
                for line in s_SPH_lines[c+1:c+7]:
                    if bool(re.match(br'(.*?)\=\"(.*)(?=\")',line)):
                        # data fields within quotes
                        dsfield,dsvalue=re.findall(br'(.*?)\=\"(.*)(?=\")',line).pop()
                        s_SPH_fields[key][dsfield.decode('utf-8')] = dsvalue.decode('utf-8').rstrip()
                    elif bool(re.match(br'(.*?)\=(.*)',line)):
                        # data fields without quotes
                        dsfield,dsvalue=re.findall(br'(.*?)\=(.*)',line).pop()
                        s_SPH_fields[key][dsfield.decode('utf-8')] = dsvalue.decode('utf-8').rstrip()
                # add 6 to counter to go to next entry
                c += 6
            # use regular expression operators to read headers
            elif bool(re.match(br'(.*?)\=\"(.*)(?=\")',s_SPH_lines[c])):
                # data fields within quotes
                field,value=re.findall(br'(.*?)\=\"(.*)(?=\")',s_SPH_lines[c]).pop()
                s_SPH_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()
            elif bool(re.match(br'(.*?)\=(.*)',s_SPH_lines[c])):
                # data fields without quotes
                field,value=re.findall(br'(.*?)\=(.*)',s_SPH_lines[c]).pop()
                s_SPH_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()
            # add 1 to counter to go to next line
            c += 1

        # Return block name array to calling function
        return s_SPH_fields

    def read_DSD(self, full_filename, DS_TYPE=None):
        """
        Read ASCII Data Set Descriptors (DSD) block from a PDS file
        """
        # read input data file
        with open(os.path.expanduser(full_filename), 'rb') as fid:
            file_contents = fid.read().splitlines()

        # Define constant values associated with PDS file formats
        # number of text lines in standard MPH
        n_MPH_lines = 41
        # number of text lines in a DSD header
        n_DSD_lines = 8

        # Level-1b CryoSat DS_NAMES within files
        regex_patterns = []
        if (DS_TYPE == 'CS_L1B'):
            regex_patterns.append(br'DS_NAME\="SIR_L1B_LRM[\s+]*"')
            regex_patterns.append(br'DS_NAME\="SIR_L1B_SAR[\s+]*"')
            regex_patterns.append(br'DS_NAME\="SIR_L1B_SARIN[\s+]*"')
        elif (DS_TYPE == 'SIR_L1B_FDM'):
            regex_patterns.append(br'DS_NAME\="SIR_L1B_FDM[\s+]*"')
        # find the DSD starting line within the SPH header
        c = 0
        Flag = False
        while ((Flag is False) and (c < len(regex_patterns))):
            # find indice within
            indice = [i for i,line in enumerate(file_contents[n_MPH_lines+1:]) if
                re.search(regex_patterns[c],line)]
            if indice:
                Flag = True
            else:
                c+=1
        # check that valid indice was found within header
        if not indice:
            raise IOError('Can not find correct DSD field')

        # extract s_DSD_fields info
        DSD_START = n_MPH_lines + indice[0] + 1
        s_DSD_fields = {}
        for i in range(DSD_START,DSD_START+n_DSD_lines):
            # use regular expression operators to read headers
            if bool(re.match(br'(.*?)\=\"(.*)(?=\")',file_contents[i])):
                # data fields within quotes
                field,value=re.findall(br'(.*?)\=\"(.*)(?=\")',file_contents[i]).pop()
                s_DSD_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()
            elif bool(re.match(br'(.*?)\=(.*)',file_contents[i])):
                # data fields without quotes
                field,value=re.findall(br'(.*?)\=(.*)',file_contents[i]).pop()
                s_DSD_fields[field.decode('utf-8')] = value.decode('utf-8').rstrip()

        # Return block name array to calling function
        return s_DSD_fields

    def cryosat_baseline_AB(self, fid, n_records):
        """
        Read L1b MDS variables for CryoSat Baselines A and B
        """
        n_SARIN_RW = 512
        n_SAR_RW = 128
        n_LRM_RW = 128
        n_blocks = 20
        n_BeamBehaviourParams = 50

        # Bind all the variables of the l1b_mds together into a single dictionary
        CS_l1b_mds = {}

        # CryoSat-2 Time and Orbit Group
        CS_l1b_mds['Location'] = {}
        # Time: day part
        CS_l1b_mds['Location']['Day'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32,fill_value=0)
        # Time: second part
        CS_l1b_mds['Location']['Second'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint32)
        # Time: microsecond part
        CS_l1b_mds['Location']['Micsec'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint32)
        # USO correction factor
        CS_l1b_mds['Location']['USO_Corr'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Mode ID
        CS_l1b_mds['Location']['Mode_ID'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint16)
        # Source sequence counter
        CS_l1b_mds['Location']['SSC'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint16)
        # Instrument configuration
        CS_l1b_mds['Location']['Inst_config'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint32)
        # Record Counter
        CS_l1b_mds['Location']['Rec_Count'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint32)
        # Lat: packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_mds['Location']['Lat'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Lon: packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_mds['Location']['Lon'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Alt: packed units (mm, 1e-3 m)
        # Altitude of COG above reference ellipsoid (interpolated value)
        CS_l1b_mds['Location']['Alt'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Instantaneous altitude rate derived from orbit: packed units (mm/s, 1e-3 m/s)
        CS_l1b_mds['Location']['Alt_rate'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Satellite velocity vector. In ITRF: packed units (mm/s, 1e-3 m/s)
        # ITRF= International Terrestrial Reference Frame
        CS_l1b_mds['Location']['Sat_velocity'] = np.ma.zeros((n_records,n_blocks,3),dtype=np.int32)
        # Real beam direction vector. In CRF: packed units (micro-m, 1e-6 m)
        # CRF= CryoSat Reference Frame.
        CS_l1b_mds['Location']['Real_beam'] = np.ma.zeros((n_records,n_blocks,3),dtype=np.int32)
        # Interferometric baseline vector. In CRF: packed units (micro-m, 1e-6 m)
        CS_l1b_mds['Location']['Baseline'] = np.ma.zeros((n_records,n_blocks,3),dtype=np.int32)
        # Measurement Confidence Data Flags
        # Generally the MCD flags indicate problems when set
        # If MCD is 0 then no problems or non-nominal conditions were detected
        # Serious errors are indicated by setting bit 31
        CS_l1b_mds['Location']['MCD'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint32)

        # CryoSat-2 Measurement Group
        # Derived from instrument measurement parameters
        CS_l1b_mds['Data'] = {}
        # Window Delay reference (two-way) corrected for instrument delays
        CS_l1b_mds['Data']['TD'] = np.ma.zeros((n_records,n_blocks),dtype=np.int64)
        # H0 Initial Height Word from telemetry
        CS_l1b_mds['Data']['H_0'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # COR2 Height Rate: on-board tracker height rate over the radar cycle
        CS_l1b_mds['Data']['COR2'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Coarse Range Word (LAI) derived from telemetry
        CS_l1b_mds['Data']['LAI'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Fine Range Word (FAI) derived from telemetry
        CS_l1b_mds['Data']['FAI'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Automatic Gain Control Channel 1: AGC gain applied on Rx channel 1.
        # Gain calibration corrections are applied (Sum of AGC stages 1 and 2
        # plus the corresponding corrections) (dB/100)
        CS_l1b_mds['Data']['AGC_CH1'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Automatic Gain Control Channel 2: AGC gain applied on Rx channel 2.
        # Gain calibration corrections are applied (dB/100)
        CS_l1b_mds['Data']['AGC_CH2'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Total Fixed Gain On Channel 1: gain applied by the RF unit. (dB/100)
        CS_l1b_mds['Data']['TR_gain_CH1'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Total Fixed Gain On Channel 2: gain applied by the RF unit. (dB/100)
        CS_l1b_mds['Data']['TR_gain_CH2'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Transmit Power in microWatts
        CS_l1b_mds['Data']['TX_Power'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Doppler range correction: Radial component (mm)
        # computed for the component of satellite velocity in the nadir direction
        CS_l1b_mds['Data']['Doppler_range'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Instrument Range Correction: transmit-receive antenna (mm)
        # Calibration correction to range on channel 1 computed from CAL1.
        CS_l1b_mds['Data']['TR_inst_range'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Instrument Range Correction: receive-only antenna (mm)
        # Calibration correction to range on channel 2 computed from CAL1.
        CS_l1b_mds['Data']['R_inst_range'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Instrument Gain Correction: transmit-receive antenna (dB/100)
        # Calibration correction to gain on channel 1 computed from CAL1
        CS_l1b_mds['Data']['TR_inst_gain'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Instrument Gain Correction: receive-only (dB/100)
        # Calibration correction to gain on channel 2 computed from CAL1
        CS_l1b_mds['Data']['R_inst_gain'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Internal Phase Correction (microradians)
        CS_l1b_mds['Data']['Internal_phase'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # External Phase Correction (microradians)
        CS_l1b_mds['Data']['External_phase'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Noise Power measurement (dB/100): converted from telemetry units to be
        # the noise floor of FBR measurement echoes.
        # Set to -9999.99 when the telemetry contains zero.
        CS_l1b_mds['Data']['Noise_power'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Phase slope correction (microradians)
        # Computed from the CAL-4 packets during the azimuth impulse response
        # amplitude (SARIN only). Set from the latest available CAL-4 packet.
        CS_l1b_mds['Data']['Phase_slope'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        CS_l1b_mds['Data']['Spares1'] = np.ma.zeros((n_records,n_blocks,4),dtype=np.int8)

        # CryoSat-2 External Corrections Group
        CS_l1b_mds['Geometry'] = {}
        # Dry Tropospheric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['dryTrop'] = np.ma.zeros((n_records),dtype=np.int32)
        # Wet Tropospheric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['wetTrop'] = np.ma.zeros((n_records),dtype=np.int32)
        # Inverse Barometric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['InvBar'] = np.ma.zeros((n_records),dtype=np.int32)
        # Delta Inverse Barometric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['DAC'] = np.ma.zeros((n_records),dtype=np.int32)
        # GIM Ionospheric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['Iono_GIM'] = np.ma.zeros((n_records),dtype=np.int32)
        # Model Ionospheric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['Iono_model'] = np.ma.zeros((n_records),dtype=np.int32)
        # Ocean tide Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['ocTideElv'] = np.ma.zeros((n_records),dtype=np.int32)
        # Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['lpeTideElv'] = np.ma.zeros((n_records),dtype=np.int32)
        # Ocean loading tide Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['olTideElv'] = np.ma.zeros((n_records),dtype=np.int32)
        # Solid Earth tide Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['seTideElv'] = np.ma.zeros((n_records),dtype=np.int32)
        # Geocentric Polar tide Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['gpTideElv'] = np.ma.zeros((n_records),dtype=np.int32)
        # Surface Type: enumerated key to classify surface at nadir
        # 0 = Open Ocean
        # 1 = Closed Sea
        # 2 = Continental Ice
        # 3 = Land
        CS_l1b_mds['Geometry']['Surf_type'] = np.ma.zeros((n_records),dtype=np.uint32)
        CS_l1b_mds['Geometry']['Spare1'] = np.ma.zeros((n_records,4),dtype=np.int8)
        # Corrections Status Flag
        CS_l1b_mds['Geometry']['Corr_status'] = np.ma.zeros((n_records),dtype=np.uint32)
        # Correction Error Flag
        CS_l1b_mds['Geometry']['Corr_error'] = np.ma.zeros((n_records),dtype=np.uint32)
        CS_l1b_mds['Geometry']['Spare2'] = np.ma.zeros((n_records,4),dtype=np.int8)

        # CryoSat-2 Average Waveforms Groups
        CS_l1b_mds['Waveform_1Hz'] = {}
        if (self.MODE == 'LRM'):
            # Low-Resolution Mode
            # Data Record Time (MDSR Time Stamp)
            CS_l1b_mds['Waveform_1Hz']['Day'] = np.zeros((n_records),dtype=np.int32)
            CS_l1b_mds['Waveform_1Hz']['Second'] = np.zeros((n_records),dtype=np.uint32)
            CS_l1b_mds['Waveform_1Hz']['Micsec'] = np.zeros((n_records),dtype=np.uint32)
            # Lat: packed units (0.1 micro-degree, 1e-7 degrees)
            CS_l1b_mds['Waveform_1Hz']['Lat'] = np.zeros((n_records),dtype=np.int32)
            # Lon: packed units (0.1 micro-degree, 1e-7 degrees)
            CS_l1b_mds['Waveform_1Hz']['Lon'] = np.zeros((n_records),dtype=np.int32)
            # Alt: packed units (mm, 1e-3 m)
            # Altitude of COG above reference ellipsoid (interpolated value)
            CS_l1b_mds['Waveform_1Hz']['Alt'] = np.zeros((n_records),dtype=np.int32)
            # Window Delay (two-way) corrected for instrument delays
            CS_l1b_mds['Waveform_1Hz']['TD'] = np.zeros((n_records),dtype=np.int64)
            # 1 Hz Averaged Power Echo Waveform
            CS_l1b_mds['Waveform_1Hz']['Waveform'] = np.zeros((n_records,n_LRM_RW),dtype=np.uint16)
            # Echo Scale Factor (to scale echo to watts)
            CS_l1b_mds['Waveform_1Hz']['Linear_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
            # Echo Scale Power (a power of 2)
            CS_l1b_mds['Waveform_1Hz']['Power2_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
            # Number of echoes averaged
            CS_l1b_mds['Waveform_1Hz']['N_avg_echoes'] = np.zeros((n_records),dtype=np.uint16)
            CS_l1b_mds['Waveform_1Hz']['Flags'] = np.zeros((n_records),dtype=np.uint16)
        elif (self.MODE == 'SAR'):
            # SAR Mode
            # Data Record Time (MDSR Time Stamp)
            CS_l1b_mds['Waveform_1Hz']['Day'] = np.zeros((n_records),dtype=np.int32)
            CS_l1b_mds['Waveform_1Hz']['Second'] = np.zeros((n_records),dtype=np.uint32)
            CS_l1b_mds['Waveform_1Hz']['Micsec'] = np.zeros((n_records),dtype=np.uint32)
            # Lat: packed units (0.1 micro-degree, 1e-7 degrees)
            CS_l1b_mds['Waveform_1Hz']['Lat'] = np.zeros((n_records),dtype=np.int32)
            # Lon: packed units (0.1 micro-degree, 1e-7 degrees)
            CS_l1b_mds['Waveform_1Hz']['Lon'] = np.zeros((n_records),dtype=np.int32)
            # Alt: packed units (mm, 1e-3 m)
            # Altitude of COG above reference ellipsoid (interpolated value)
            CS_l1b_mds['Waveform_1Hz']['Alt'] = np.zeros((n_records),dtype=np.int32)
            # Window Delay (two-way) corrected for instrument delays
            CS_l1b_mds['Waveform_1Hz']['TD'] = np.zeros((n_records),dtype=np.int64)
            # 1 Hz Averaged Power Echo Waveform
            CS_l1b_mds['Waveform_1Hz']['Waveform'] = np.zeros((n_records,n_SAR_RW),dtype=np.uint16)
            # Echo Scale Factor (to scale echo to watts)
            CS_l1b_mds['Waveform_1Hz']['Linear_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
            # Echo Scale Power (a power of 2)
            CS_l1b_mds['Waveform_1Hz']['Power2_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
            # Number of echoes averaged
            CS_l1b_mds['Waveform_1Hz']['N_avg_echoes'] = np.zeros((n_records),dtype=np.uint16)
            CS_l1b_mds['Waveform_1Hz']['Flags'] = np.zeros((n_records),dtype=np.uint16)
        elif (self.MODE == 'SIN'):
            # SARIN Mode
            # Same as the LRM/SAR groups but the waveform array is 512 bins instead of
            # 128 and the number of echoes averaged is different.
            # Data Record Time (MDSR Time Stamp)
            CS_l1b_mds['Waveform_1Hz']['Day'] = np.zeros((n_records),dtype=np.int32)
            CS_l1b_mds['Waveform_1Hz']['Second'] = np.zeros((n_records),dtype=np.uint32)
            CS_l1b_mds['Waveform_1Hz']['Micsec'] = np.zeros((n_records),dtype=np.uint32)
            # Lat: packed units (0.1 micro-degree, 1e-7 degrees)
            CS_l1b_mds['Waveform_1Hz']['Lat'] = np.zeros((n_records),dtype=np.int32)
            # Lon: packed units (0.1 micro-degree, 1e-7 degrees)
            CS_l1b_mds['Waveform_1Hz']['Lon'] = np.zeros((n_records),dtype=np.int32)
            # Alt: packed units (mm, 1e-3 m)
            # Altitude of COG above reference ellipsoid (interpolated value)
            CS_l1b_mds['Waveform_1Hz']['Alt'] = np.zeros((n_records),dtype=np.int32)
            # Window Delay (two-way) corrected for instrument delays
            CS_l1b_mds['Waveform_1Hz']['TD'] = np.zeros((n_records),dtype=np.int64)
            # 1 Hz Averaged Power Echo Waveform
            CS_l1b_mds['Waveform_1Hz']['Waveform'] = np.zeros((n_records,n_SARIN_RW),dtype=np.uint16)
            # Echo Scale Factor (to scale echo to watts)
            CS_l1b_mds['Waveform_1Hz']['Linear_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
            # Echo Scale Power (a power of 2)
            CS_l1b_mds['Waveform_1Hz']['Power2_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
            # Number of echoes averaged
            CS_l1b_mds['Waveform_1Hz']['N_avg_echoes'] = np.zeros((n_records),dtype=np.uint16)
            CS_l1b_mds['Waveform_1Hz']['Flags'] = np.zeros((n_records),dtype=np.uint16)

        # CryoSat-2 Waveforms Groups
        # Beam Behavior Parameters
        Beam_Behavior = {}
        # Standard Deviation of Gaussian fit to range integrated stack power.
        Beam_Behavior['SD'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
        # Stack Center: Mean of Gaussian fit to range integrated stack power.
        Beam_Behavior['Center'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
        # Stack amplitude parameter scaled in dB/100.
        Beam_Behavior['Amplitude'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
        # 3rd moment: providing the degree of asymmetry of the range integrated
        # stack power distribution.
        Beam_Behavior['Skewness'] = np.zeros((n_records,n_blocks),dtype=np.int16)
        # 4th moment: Measure of peakiness of range integrated stack power distribution.
        Beam_Behavior['Kurtosis'] = np.zeros((n_records,n_blocks),dtype=np.int16)
        Beam_Behavior['Spare'] = np.zeros((n_records,n_blocks,n_BeamBehaviourParams-5),dtype=np.int16)

        # CryoSat-2 mode specific waveforms
        CS_l1b_mds['Waveform_20Hz'] = {}
        if (self.MODE == 'LRM'):
            # Low-Resolution Mode
            # Averaged Power Echo Waveform [128]
            CS_l1b_mds['Waveform_20Hz']['Waveform'] = np.zeros((n_records,n_blocks,n_LRM_RW),dtype=np.uint16)
            # Echo Scale Factor (to scale echo to watts)
            CS_l1b_mds['Waveform_20Hz']['Linear_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
            # Echo Scale Power (a power of 2 to scale echo to Watts)
            CS_l1b_mds['Waveform_20Hz']['Power2_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
            # Number of echoes averaged
            CS_l1b_mds['Waveform_20Hz']['N_avg_echoes'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
            CS_l1b_mds['Waveform_20Hz']['Flags'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
        elif (self.MODE == 'SAR'):
            # SAR Mode
            # Averaged Power Echo Waveform [128]
            CS_l1b_mds['Waveform_20Hz']['Waveform'] = np.zeros((n_records,n_blocks,n_SAR_RW),dtype=np.uint16)
            # Echo Scale Factor (to scale echo to watts)
            CS_l1b_mds['Waveform_20Hz']['Linear_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
            # Echo Scale Power (a power of 2 to scale echo to Watts)
            CS_l1b_mds['Waveform_20Hz']['Power2_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
            # Number of echoes averaged
            CS_l1b_mds['Waveform_20Hz']['N_avg_echoes'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
            CS_l1b_mds['Waveform_20Hz']['Flags'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
            # Beam behaviour parameters
            CS_l1b_mds['Waveform_20Hz']['Beam'] = Beam_Behavior
        elif (self.MODE == 'SIN'):
            # SARIN Mode
            # Averaged Power Echo Waveform [512]
            CS_l1b_mds['Waveform_20Hz']['Waveform'] = np.zeros((n_records,n_blocks,n_SARIN_RW),dtype=np.uint16)
            # Echo Scale Factor (to scale echo to watts)
            CS_l1b_mds['Waveform_20Hz']['Linear_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
            # Echo Scale Power (a power of 2 to scale echo to Watts)
            CS_l1b_mds['Waveform_20Hz']['Power2_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
            # Number of echoes averaged
            CS_l1b_mds['Waveform_20Hz']['N_avg_echoes'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
            CS_l1b_mds['Waveform_20Hz']['Flags'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
            # Beam behaviour parameters
            CS_l1b_mds['Waveform_20Hz']['Beam'] = Beam_Behavior
            # Coherence [512]: packed units (1/1000)
            CS_l1b_mds['Waveform_20Hz']['Coherence'] = np.zeros((n_records,n_blocks,n_SARIN_RW),dtype=np.int16)
            # Phase Difference [512]: packed units (microradians)
            CS_l1b_mds['Waveform_20Hz']['Phase_diff'] = np.zeros((n_records,n_blocks,n_SARIN_RW),dtype=np.int32)

        # for each record in the CryoSat file
        for r in range(n_records):
            # CryoSat-2 Time and Orbit Group
            for b in range(n_blocks):
                CS_l1b_mds['Location']['Day'].data[r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Location']['Second'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Location']['Micsec'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Location']['USO_Corr'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Location']['Mode_ID'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                CS_l1b_mds['Location']['SSC'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                CS_l1b_mds['Location']['Inst_config'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Location']['Rec_Count'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Location']['Lat'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Location']['Lon'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Location']['Alt'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Location']['Alt_rate'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Location']['Sat_velocity'][r,b,:] = np.fromfile(fid,dtype='>i4',count=3)
                CS_l1b_mds['Location']['Real_beam'][r,b,:] = np.fromfile(fid,dtype='>i4',count=3)
                CS_l1b_mds['Location']['Baseline'][r,b,:] = np.fromfile(fid,dtype='>i4',count=3)
                CS_l1b_mds['Location']['MCD'][r,b] = np.fromfile(fid,dtype='>u4',count=1)

            # CryoSat-2 Measurement Group
            # Derived from instrument measurement parameters
            for b in range(n_blocks):
                CS_l1b_mds['Data']['TD'][r,b] = np.fromfile(fid,dtype='>i8',count=1)
                CS_l1b_mds['Data']['H_0'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['COR2'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['LAI'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['FAI'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['AGC_CH1'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['AGC_CH2'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['TR_gain_CH1'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['TR_gain_CH2'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['TX_Power'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['Doppler_range'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['TR_inst_range'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['R_inst_range'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['TR_inst_gain'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['R_inst_gain'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['Internal_phase'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['External_phase'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['Noise_power'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['Phase_slope'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['Spares1'][r,b,:] = np.fromfile(fid,dtype='>i1',count=4)

            # CryoSat-2 External Corrections Group
            CS_l1b_mds['Geometry']['dryTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['wetTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['InvBar'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['DAC'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['Iono_GIM'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['Iono_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['ocTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['lpeTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['olTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['seTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['gpTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['Surf_type'][r] = np.fromfile(fid,dtype='>u4',count=1)
            CS_l1b_mds['Geometry']['Spare1'][r,:] = np.fromfile(fid,dtype='>i1',count=4)
            CS_l1b_mds['Geometry']['Corr_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
            CS_l1b_mds['Geometry']['Corr_error'][r] = np.fromfile(fid,dtype='>u4',count=1)
            CS_l1b_mds['Geometry']['Spare2'][r,:] = np.fromfile(fid,dtype='>i1',count=4)

            # CryoSat-2 Average Waveforms Groups
            if (self.MODE == 'LRM'):
                # Low-Resolution Mode
                CS_l1b_mds['Waveform_1Hz']['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Second'][r] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Micsec'][r] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Lat'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Lon'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Alt'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['TD'][r] = np.fromfile(fid,dtype='>i8',count=1)
                CS_l1b_mds['Waveform_1Hz']['Waveform'][r,:] = np.fromfile(fid,dtype='>u2',count=n_LRM_RW)
                CS_l1b_mds['Waveform_1Hz']['Linear_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Power2_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['N_avg_echoes'][r] = np.fromfile(fid,dtype='>u2',count=1)
                CS_l1b_mds['Waveform_1Hz']['Flags'][r] = np.fromfile(fid,dtype='>u2',count=1)
            elif (self.MODE == 'SAR'):
                # SAR Mode
                CS_l1b_mds['Waveform_1Hz']['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Second'][r] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Micsec'][r] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Lat'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Lon'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Alt'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['TD'][r] = np.fromfile(fid,dtype='>i8',count=1)
                CS_l1b_mds['Waveform_1Hz']['Waveform'][r,:] = np.fromfile(fid,dtype='>u2',count=n_SAR_RW)
                CS_l1b_mds['Waveform_1Hz']['Linear_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Power2_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['N_avg_echoes'][r] = np.fromfile(fid,dtype='>u2',count=1)
                CS_l1b_mds['Waveform_1Hz']['Flags'][r] = np.fromfile(fid,dtype='>u2',count=1)
            elif (self.MODE == 'SIN'):
                # SARIN Mode
                CS_l1b_mds['Waveform_1Hz']['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Second'][r] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Micsec'][r] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Lat'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Lon'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Alt'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['TD'][r] = np.fromfile(fid,dtype='>i8',count=1)
                CS_l1b_mds['Waveform_1Hz']['Waveform'][r,:] = np.fromfile(fid,dtype='>u2',count=n_SARIN_RW)
                CS_l1b_mds['Waveform_1Hz']['Linear_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Power2_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['N_avg_echoes'][r] = np.fromfile(fid,dtype='>u2',count=1)
                CS_l1b_mds['Waveform_1Hz']['Flags'][r] = np.fromfile(fid,dtype='>u2',count=1)

            # CryoSat-2 Waveforms Groups
            if (self.MODE == 'LRM'):
                # Low-Resolution Mode
                for b in range(n_blocks):
                    CS_l1b_mds['Waveform_20Hz']['Waveform'][r,b,:] = np.fromfile(fid,dtype='>u2',count=n_LRM_RW)
                    CS_l1b_mds['Waveform_20Hz']['Linear_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Power2_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                    CS_l1b_mds['Waveform_20Hz']['N_avg_echoes'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Flags'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
            elif (self.MODE == 'SAR'):
                # SAR Mode
                for b in range(n_blocks):
                    CS_l1b_mds['Waveform_20Hz']['Waveform'][r,b,:] = np.fromfile(fid,dtype='>u2',count=n_SAR_RW)
                    CS_l1b_mds['Waveform_20Hz']['Linear_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Power2_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                    CS_l1b_mds['Waveform_20Hz']['N_avg_echoes'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Flags'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['SD'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Center'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Amplitude'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Skewness'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Kurtosis'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Spare'][r,b,:] = np.fromfile(fid,dtype='>i2',count=(n_BeamBehaviourParams-5))
            elif (self.MODE == 'SIN'):
                # SARIN Mode
                for b in range(n_blocks):
                    CS_l1b_mds['Waveform_20Hz']['Waveform'][r,b,:] = np.fromfile(fid,dtype='>u2',count=n_SARIN_RW)
                    CS_l1b_mds['Waveform_20Hz']['Linear_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Power2_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                    CS_l1b_mds['Waveform_20Hz']['N_avg_echoes'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Flags'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['SD'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Center'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Amplitude'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Skewness'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Kurtosis'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Spare'][r,b,:] = np.fromfile(fid,dtype='>i2',count=(n_BeamBehaviourParams-5))
                    CS_l1b_mds['Waveform_20Hz']['Coherence'][r,b,:] = np.fromfile(fid,dtype='>i2',count=n_SARIN_RW)
                    CS_l1b_mds['Waveform_20Hz']['Phase_diff'][r,b,:] = np.fromfile(fid,dtype='>i4',count=n_SARIN_RW)

        # set the mask from day variables
        mask_20Hz = CS_l1b_mds['Location']['Day'].data == CS_l1b_mds['Location']['Day'].fill_value
        Location_keys = [key for key in CS_l1b_mds['Location'].keys() if not re.search(r'Spare',key)]
        Data_keys = [key for key in CS_l1b_mds['Data'].keys() if not re.search(r'Spare',key)]
        Geometry_keys = [key for key in CS_l1b_mds['Geometry'].keys() if not re.search(r'Spare',key)]
        Wfm_1Hz_keys = [key for key in CS_l1b_mds['Waveform_1Hz'].keys() if not re.search(r'Spare',key)]
        Wfm_20Hz_keys = [key for key in CS_l1b_mds['Waveform_20Hz'].keys() if not re.search(r'Spare',key)]
        for key in Location_keys:
            CS_l1b_mds['Location'][key].mask = mask_20Hz.copy()
        for key in Data_keys:
            CS_l1b_mds['Data'][key].mask = mask_20Hz.copy()
        # return the output dictionary
        return CS_l1b_mds

    def cryosat_baseline_C(self, fid, n_records):
        """
        Read L1b MDS variables for CryoSat Baseline C
        """
        n_SARIN_BC_RW = 1024
        n_SARIN_RW = 512
        n_SAR_BC_RW = 256
        n_SAR_RW = 128
        n_LRM_RW = 128
        n_blocks = 20
        n_BeamBehaviourParams = 50

        # Bind all the variables of the l1b_mds together into a single dictionary
        CS_l1b_mds = {}

        # CryoSat-2 Time and Orbit Group
        CS_l1b_mds['Location'] = {}
        # Time: day part
        CS_l1b_mds['Location']['Day'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32,fill_value=0)
        # Time: second part
        CS_l1b_mds['Location']['Second'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint32)
        # Time: microsecond part
        CS_l1b_mds['Location']['Micsec'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint32)
        # USO correction factor
        CS_l1b_mds['Location']['USO_Corr'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Mode ID
        CS_l1b_mds['Location']['Mode_ID'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint16)
        # Source sequence counter
        CS_l1b_mds['Location']['SSC'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint16)
        # Instrument configuration
        CS_l1b_mds['Location']['Inst_config'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint32)
        # Record Counter
        CS_l1b_mds['Location']['Rec_Count'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint32)
        # Lat: packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_mds['Location']['Lat'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Lon: packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_mds['Location']['Lon'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Alt: packed units (mm, 1e-3 m)
        # Altitude of COG above reference ellipsoid (interpolated value)
        CS_l1b_mds['Location']['Alt'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Instantaneous altitude rate derived from orbit: packed units (mm/s, 1e-3 m/s)
        CS_l1b_mds['Location']['Alt_rate'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Satellite velocity vector. In ITRF: packed units (mm/s, 1e-3 m/s)
        # ITRF= International Terrestrial Reference Frame
        CS_l1b_mds['Location']['Sat_velocity'] = np.ma.zeros((n_records,n_blocks,3),dtype=np.int32)
        # Real beam direction vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
        # CRF= CryoSat Reference Frame.
        CS_l1b_mds['Location']['Real_beam'] = np.ma.zeros((n_records,n_blocks,3),dtype=np.int32)
        # Interferometric baseline vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
        CS_l1b_mds['Location']['Baseline'] = np.ma.zeros((n_records,n_blocks,3),dtype=np.int32)
        # Star Tracker ID
        CS_l1b_mds['Location']['ST_ID'] = np.ma.zeros((n_records,n_blocks),dtype=np.int16)
        # Antenna Bench Roll Angle (Derived from star trackers)
        # packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_mds['Location']['Roll'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Antenna Bench Pitch Angle (Derived from star trackers)
        # packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_mds['Location']['Pitch'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Antenna Bench Yaw Angle (Derived from star trackers)
        # packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_mds['Location']['Yaw'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Measurement Confidence Data Flags
        # Generally the MCD flags indicate problems when set
        # If MCD is 0 then no problems or non-nominal conditions were detected
        # Serious errors are indicated by setting bit 31
        CS_l1b_mds['Location']['MCD'] = np.ma.zeros((n_records,n_blocks),dtype=np.uint32)
        CS_l1b_mds['Location']['Spares'] = np.ma.zeros((n_records,n_blocks,2),dtype=np.int16)

        # CryoSat-2 Measurement Group
        # Derived from instrument measurement parameters
        CS_l1b_mds['Data'] = {}
        # Window Delay reference (two-way) corrected for instrument delays
        CS_l1b_mds['Data']['TD'] = np.ma.zeros((n_records,n_blocks),dtype=np.int64)
        # H0 Initial Height Word from telemetry
        CS_l1b_mds['Data']['H_0'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # COR2 Height Rate: on-board tracker height rate over the radar cycle
        CS_l1b_mds['Data']['COR2'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Coarse Range Word (LAI) derived from telemetry
        CS_l1b_mds['Data']['LAI'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Fine Range Word (FAI) derived from telemetry
        CS_l1b_mds['Data']['FAI'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Automatic Gain Control Channel 1: AGC gain applied on Rx channel 1.
        # Gain calibration corrections are applied (Sum of AGC stages 1 and 2
        # plus the corresponding corrections) (dB/100)
        CS_l1b_mds['Data']['AGC_CH1'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Automatic Gain Control Channel 2: AGC gain applied on Rx channel 2.
        # Gain calibration corrections are applied (dB/100)
        CS_l1b_mds['Data']['AGC_CH2'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Total Fixed Gain On Channel 1: gain applied by the RF unit. (dB/100)
        CS_l1b_mds['Data']['TR_gain_CH1'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Total Fixed Gain On Channel 2: gain applied by the RF unit. (dB/100)
        CS_l1b_mds['Data']['TR_gain_CH2'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Transmit Power in microWatts
        CS_l1b_mds['Data']['TX_Power'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Doppler range correction: Radial component (mm)
        # computed for the component of satellite velocity in the nadir direction
        CS_l1b_mds['Data']['Doppler_range'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Instrument Range Correction: transmit-receive antenna (mm)
        # Calibration correction to range on channel 1 computed from CAL1.
        CS_l1b_mds['Data']['TR_inst_range'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Instrument Range Correction: receive-only antenna (mm)
        # Calibration correction to range on channel 2 computed from CAL1.
        CS_l1b_mds['Data']['R_inst_range'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Instrument Gain Correction: transmit-receive antenna (dB/100)
        # Calibration correction to gain on channel 1 computed from CAL1
        CS_l1b_mds['Data']['TR_inst_gain'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Instrument Gain Correction: receive-only (dB/100)
        # Calibration correction to gain on channel 2 computed from CAL1
        CS_l1b_mds['Data']['R_inst_gain'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Internal Phase Correction (microradians)
        CS_l1b_mds['Data']['Internal_phase'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # External Phase Correction (microradians)
        CS_l1b_mds['Data']['External_phase'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Noise Power measurement (dB/100)
        CS_l1b_mds['Data']['Noise_power'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        # Phase slope correction (microradians)
        # Computed from the CAL-4 packets during the azimuth impulse response
        # amplitude (SARIN only). Set from the latest available CAL-4 packet.
        CS_l1b_mds['Data']['Phase_slope'] = np.ma.zeros((n_records,n_blocks),dtype=np.int32)
        CS_l1b_mds['Data']['Spares1'] = np.ma.zeros((n_records,n_blocks,4),dtype=np.int8)

        # CryoSat-2 External Corrections Group
        CS_l1b_mds['Geometry'] = {}
        # Dry Tropospheric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['dryTrop'] = np.ma.zeros((n_records),dtype=np.int32)
        # Wet Tropospheric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['wetTrop'] = np.ma.zeros((n_records),dtype=np.int32)
        # Inverse Barometric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['InvBar'] = np.ma.zeros((n_records),dtype=np.int32)
        # Delta Inverse Barometric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['DAC'] = np.ma.zeros((n_records),dtype=np.int32)
        # GIM Ionospheric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['Iono_GIM'] = np.ma.zeros((n_records),dtype=np.int32)
        # Model Ionospheric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['Iono_model'] = np.ma.zeros((n_records),dtype=np.int32)
        # Ocean tide Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['ocTideElv'] = np.ma.zeros((n_records),dtype=np.int32)
        # Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['lpeTideElv'] = np.ma.zeros((n_records),dtype=np.int32)
        # Ocean loading tide Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['olTideElv'] = np.ma.zeros((n_records),dtype=np.int32)
        # Solid Earth tide Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['seTideElv'] = np.ma.zeros((n_records),dtype=np.int32)
        # Geocentric Polar tide Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['gpTideElv'] = np.ma.zeros((n_records),dtype=np.int32)
        # Surface Type: enumerated key to classify surface at nadir
        # 0 = Open Ocean
        # 1 = Closed Sea
        # 2 = Continental Ice
        # 3 = Land
        CS_l1b_mds['Geometry']['Surf_type'] = np.ma.zeros((n_records),dtype=np.uint32)
        CS_l1b_mds['Geometry']['Spare1'] = np.ma.zeros((n_records,4),dtype=np.int8)
        # Corrections Status Flag
        CS_l1b_mds['Geometry']['Corr_status'] = np.ma.zeros((n_records),dtype=np.uint32)
        # Correction Error Flag
        CS_l1b_mds['Geometry']['Corr_error'] = np.ma.zeros((n_records),dtype=np.uint32)
        CS_l1b_mds['Geometry']['Spare2'] = np.ma.zeros((n_records,4),dtype=np.int8)

        # CryoSat-2 Average Waveforms Groups
        CS_l1b_mds['Waveform_1Hz'] = {}
        if (self.MODE == 'LRM'):
            # Low-Resolution Mode
            # Data Record Time (MDSR Time Stamp)
            CS_l1b_mds['Waveform_1Hz']['Day'] = np.zeros((n_records),dtype=np.int32)
            CS_l1b_mds['Waveform_1Hz']['Second'] = np.zeros((n_records),dtype=np.uint32)
            CS_l1b_mds['Waveform_1Hz']['Micsec'] = np.zeros((n_records),dtype=np.uint32)
            # Lat: packed units (0.1 micro-degree, 1e-7 degrees)
            CS_l1b_mds['Waveform_1Hz']['Lat'] = np.zeros((n_records),dtype=np.int32)
            # Lon: packed units (0.1 micro-degree, 1e-7 degrees)
            CS_l1b_mds['Waveform_1Hz']['Lon'] = np.zeros((n_records),dtype=np.int32)
            # Alt: packed units (mm, 1e-3 m)
            # Altitude of COG above reference ellipsoid (interpolated value)
            CS_l1b_mds['Waveform_1Hz']['Alt'] = np.zeros((n_records),dtype=np.int32)
            # Window Delay (two-way) corrected for instrument delays
            CS_l1b_mds['Waveform_1Hz']['TD'] = np.zeros((n_records),dtype=np.int64)
            # 1 Hz Averaged Power Echo Waveform
            CS_l1b_mds['Waveform_1Hz']['Waveform'] = np.zeros((n_records,n_LRM_RW),dtype=np.uint16)
            # Echo Scale Factor (to scale echo to watts)
            CS_l1b_mds['Waveform_1Hz']['Linear_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
            # Echo Scale Power (a power of 2 to scale echo to Watts)
            CS_l1b_mds['Waveform_1Hz']['Power2_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
            # Number of echoes averaged
            CS_l1b_mds['Waveform_1Hz']['N_avg_echoes'] = np.zeros((n_records),dtype=np.uint16)
            CS_l1b_mds['Waveform_1Hz']['Flags'] = np.zeros((n_records),dtype=np.uint16)
        elif (self.MODE == 'SAR'):
            # SAR Mode
            # Data Record Time (MDSR Time Stamp)
            CS_l1b_mds['Waveform_1Hz']['Day'] = np.zeros((n_records),dtype=np.int32)
            CS_l1b_mds['Waveform_1Hz']['Second'] = np.zeros((n_records),dtype=np.uint32)
            CS_l1b_mds['Waveform_1Hz']['Micsec'] = np.zeros((n_records),dtype=np.uint32)
            # Lat: packed units (0.1 micro-degree, 1e-7 degrees)
            CS_l1b_mds['Waveform_1Hz']['Lat'] = np.zeros((n_records),dtype=np.int32)
            # Lon: packed units (0.1 micro-degree, 1e-7 degrees)
            CS_l1b_mds['Waveform_1Hz']['Lon'] = np.zeros((n_records),dtype=np.int32)
            # Alt: packed units (mm, 1e-3 m)
            # Altitude of COG above reference ellipsoid (interpolated value)
            CS_l1b_mds['Waveform_1Hz']['Alt'] = np.zeros((n_records),dtype=np.int32)
            # Window Delay (two-way) corrected for instrument delays
            CS_l1b_mds['Waveform_1Hz']['TD'] = np.zeros((n_records),dtype=np.int64)
            # 1 Hz Averaged Power Echo Waveform
            CS_l1b_mds['Waveform_1Hz']['Waveform'] = np.zeros((n_records,n_SAR_RW),dtype=np.uint16)
            # Echo Scale Factor (to scale echo to watts)
            CS_l1b_mds['Waveform_1Hz']['Linear_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
            # Echo Scale Power (a power of 2 to scale echo to Watts)
            CS_l1b_mds['Waveform_1Hz']['Power2_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
            # Number of echoes averaged
            CS_l1b_mds['Waveform_1Hz']['N_avg_echoes'] = np.zeros((n_records),dtype=np.uint16)
            CS_l1b_mds['Waveform_1Hz']['Flags'] = np.zeros((n_records),dtype=np.uint16)
        elif (self.MODE == 'SIN'):
            # SARIN Mode
            # Same as the LRM/SAR groups but the waveform array is 512 bins instead of
            # 128 and the number of echoes averaged is different.
            # Data Record Time (MDSR Time Stamp)
            CS_l1b_mds['Waveform_1Hz']['Day'] = np.zeros((n_records),dtype=np.int32)
            CS_l1b_mds['Waveform_1Hz']['Second'] = np.zeros((n_records),dtype=np.uint32)
            CS_l1b_mds['Waveform_1Hz']['Micsec'] = np.zeros((n_records),dtype=np.uint32)
            # Lat: packed units (0.1 micro-degree, 1e-7 degrees)
            CS_l1b_mds['Waveform_1Hz']['Lat'] = np.zeros((n_records),dtype=np.int32)
            # Lon: packed units (0.1 micro-degree, 1e-7 degrees)
            CS_l1b_mds['Waveform_1Hz']['Lon'] = np.zeros((n_records),dtype=np.int32)
            # Alt: packed units (mm, 1e-3 m)
            # Altitude of COG above reference ellipsoid (interpolated value)
            CS_l1b_mds['Waveform_1Hz']['Alt'] = np.zeros((n_records),dtype=np.int32)
            # Window Delay (two-way) corrected for instrument delays
            CS_l1b_mds['Waveform_1Hz']['TD'] = np.zeros((n_records),dtype=np.int64)
            # 1 Hz Averaged Power Echo Waveform
            CS_l1b_mds['Waveform_1Hz']['Waveform'] = np.zeros((n_records,n_SARIN_RW),dtype=np.uint16)
            # Echo Scale Factor (to scale echo to watts)
            CS_l1b_mds['Waveform_1Hz']['Linear_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
            # Echo Scale Power (a power of 2 to scale echo to Watts)
            CS_l1b_mds['Waveform_1Hz']['Power2_Wfm_Multiplier'] = np.zeros((n_records),dtype=np.int32)
            # Number of echoes averaged
            CS_l1b_mds['Waveform_1Hz']['N_avg_echoes'] = np.zeros((n_records),dtype=np.uint16)
            CS_l1b_mds['Waveform_1Hz']['Flags'] = np.zeros((n_records),dtype=np.uint16)

        # CryoSat-2 Waveforms Groups
        # Beam Behavior Parameters
        Beam_Behavior = {}
        # Standard Deviation of Gaussian fit to range integrated stack power.
        Beam_Behavior['SD'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
        # Stack Center: Mean of Gaussian fit to range integrated stack power.
        Beam_Behavior['Center'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
        # Stack amplitude parameter scaled in dB/100.
        Beam_Behavior['Amplitude'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
        # 3rd moment: providing the degree of asymmetry of the range integrated
        # stack power distribution.
        Beam_Behavior['Skewness'] = np.zeros((n_records,n_blocks),dtype=np.int16)
        # 4th moment: Measure of peakiness of range integrated stack power distribution.
        Beam_Behavior['Kurtosis'] = np.zeros((n_records,n_blocks),dtype=np.int16)
        # Standard deviation as a function of boresight angle (microradians)
        Beam_Behavior['SD_boresight_angle'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
        # Stack Center angle as a function of boresight angle (microradians)
        Beam_Behavior['Center_boresight_angle'] = np.zeros((n_records,n_blocks),dtype=np.int16)
        Beam_Behavior['Spare'] = np.zeros((n_records,n_blocks,n_BeamBehaviourParams-7),dtype=np.int16)

        # CryoSat-2 mode specific waveform variables
        CS_l1b_mds['Waveform_20Hz'] = {}
        if (self.MODE == 'LRM'):
            # Low-Resolution Mode
            # Averaged Power Echo Waveform [128]
            CS_l1b_mds['Waveform_20Hz']['Waveform'] = np.zeros((n_records,n_blocks,n_LRM_RW),dtype=np.uint16)
            # Echo Scale Factor (to scale echo to watts)
            CS_l1b_mds['Waveform_20Hz']['Linear_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
            # Echo Scale Power (a power of 2)
            CS_l1b_mds['Waveform_20Hz']['Power2_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
            # Number of echoes averaged
            CS_l1b_mds['Waveform_20Hz']['N_avg_echoes'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
            CS_l1b_mds['Waveform_20Hz']['Flags'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
        elif (self.MODE == 'SAR'):
            # SAR Mode
            # Averaged Power Echo Waveform [256]
            CS_l1b_mds['Waveform_20Hz']['Waveform'] = np.zeros((n_records,n_blocks,n_SAR_BC_RW),dtype=np.uint16)
            # Echo Scale Factor (to scale echo to watts)
            CS_l1b_mds['Waveform_20Hz']['Linear_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
            # Echo Scale Power (a power of 2)
            CS_l1b_mds['Waveform_20Hz']['Power2_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
            # Number of echoes averaged
            CS_l1b_mds['Waveform_20Hz']['N_avg_echoes'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
            CS_l1b_mds['Waveform_20Hz']['Flags'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
            # Beam behaviour parameters
            CS_l1b_mds['Waveform_20Hz']['Beam'] = Beam_Behavior
        elif (self.MODE == 'SIN'):
            # SARIN Mode
            # Averaged Power Echo Waveform [1024]
            CS_l1b_mds['Waveform_20Hz']['Waveform'] = np.zeros((n_records,n_blocks,n_SARIN_BC_RW),dtype=np.uint16)
            # Echo Scale Factor (to scale echo to watts)
            CS_l1b_mds['Waveform_20Hz']['Linear_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
            # Echo Scale Power (a power of 2)
            CS_l1b_mds['Waveform_20Hz']['Power2_Wfm_Multiplier'] = np.zeros((n_records,n_blocks),dtype=np.int32)
            # Number of echoes averaged
            CS_l1b_mds['Waveform_20Hz']['N_avg_echoes'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
            CS_l1b_mds['Waveform_20Hz']['Flags'] = np.zeros((n_records,n_blocks),dtype=np.uint16)
            # Beam behaviour parameters
            CS_l1b_mds['Waveform_20Hz']['Beam'] = Beam_Behavior
            # Coherence [1024]: packed units (1/1000)
            CS_l1b_mds['Waveform_20Hz']['Coherence'] = np.zeros((n_records,n_blocks,n_SARIN_BC_RW),dtype=np.int16)
            # Phase Difference [1024]: packed units (microradians)
            CS_l1b_mds['Waveform_20Hz']['Phase_diff'] = np.zeros((n_records,n_blocks,n_SARIN_BC_RW),dtype=np.int32)

        # for each record in the CryoSat file
        for r in range(n_records):
            # CryoSat-2 Time and Orbit Group
            for b in range(n_blocks):
                CS_l1b_mds['Location']['Day'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Location']['Second'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Location']['Micsec'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Location']['USO_Corr'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Location']['Mode_ID'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                CS_l1b_mds['Location']['SSC'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                CS_l1b_mds['Location']['Inst_config'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Location']['Rec_Count'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Location']['Lat'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Location']['Lon'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Location']['Alt'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Location']['Alt_rate'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Location']['Sat_velocity'][r,b,:] = np.fromfile(fid,dtype='>i4',count=3)
                CS_l1b_mds['Location']['Real_beam'][r,b,:] = np.fromfile(fid,dtype='>i4',count=3)
                CS_l1b_mds['Location']['Baseline'][r,b,:] = np.fromfile(fid,dtype='>i4',count=3)
                CS_l1b_mds['Location']['ST_ID'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
                CS_l1b_mds['Location']['Roll'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Location']['Pitch'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Location']['Yaw'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Location']['MCD'][r,b] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Location']['Spares'][r,b,:] = np.fromfile(fid,dtype='>i2',count=2)

            # CryoSat-2 Measurement Group
            # Derived from instrument measurement parameters
            for b in range(n_blocks):
                CS_l1b_mds['Data']['TD'][r,b] = np.fromfile(fid,dtype='>i8',count=1)
                CS_l1b_mds['Data']['H_0'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['COR2'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['LAI'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['FAI'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['AGC_CH1'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['AGC_CH2'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['TR_gain_CH1'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['TR_gain_CH2'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['TX_Power'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['Doppler_range'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['TR_inst_range'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['R_inst_range'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['TR_inst_gain'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['R_inst_gain'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['Internal_phase'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['External_phase'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['Noise_power'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['Phase_slope'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Data']['Spares1'][r,b,:] = np.fromfile(fid,dtype='>i1',count=4)

            # CryoSat-2 External Corrections Group
            CS_l1b_mds['Geometry']['dryTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['wetTrop'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['InvBar'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['DAC'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['Iono_GIM'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['Iono_model'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['ocTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['lpeTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['olTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['seTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['gpTideElv'][r] = np.fromfile(fid,dtype='>i4',count=1)
            CS_l1b_mds['Geometry']['Surf_type'][r] = np.fromfile(fid,dtype='>u4',count=1)
            CS_l1b_mds['Geometry']['Spare1'][r,:] = np.fromfile(fid,dtype='>i1',count=4)
            CS_l1b_mds['Geometry']['Corr_status'][r] = np.fromfile(fid,dtype='>u4',count=1)
            CS_l1b_mds['Geometry']['Corr_error'][r] = np.fromfile(fid,dtype='>u4',count=1)
            CS_l1b_mds['Geometry']['Spare2'][r,:] = np.fromfile(fid,dtype='>i1',count=4)

            # CryoSat-2 Average Waveforms Groups
            if (self.MODE == 'LRM'):
                # Low-Resolution Mode
                CS_l1b_mds['Waveform_1Hz']['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Second'][r] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Micsec'][r] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Lat'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Lon'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Alt'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['TD'][r] = np.fromfile(fid,dtype='>i8',count=1)
                CS_l1b_mds['Waveform_1Hz']['Waveform'][r,:] = np.fromfile(fid,dtype='>u2',count=n_LRM_RW)
                CS_l1b_mds['Waveform_1Hz']['Linear_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Power2_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['N_avg_echoes'][r] = np.fromfile(fid,dtype='>u2',count=1)
                CS_l1b_mds['Waveform_1Hz']['Flags'][r] = np.fromfile(fid,dtype='>u2',count=1)
            elif (self.MODE == 'SAR'):
                # SAR Mode
                CS_l1b_mds['Waveform_1Hz']['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Second'][r] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Micsec'][r] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Lat'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Lon'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Alt'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['TD'][r] = np.fromfile(fid,dtype='>i8',count=1)
                CS_l1b_mds['Waveform_1Hz']['Waveform'][r,:] = np.fromfile(fid,dtype='>u2',count=n_SAR_RW)
                CS_l1b_mds['Waveform_1Hz']['Linear_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Power2_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['N_avg_echoes'][r] = np.fromfile(fid,dtype='>u2',count=1)
                CS_l1b_mds['Waveform_1Hz']['Flags'][r] = np.fromfile(fid,dtype='>u2',count=1)
            elif (self.MODE == 'SIN'):
                # SARIN Mode
                CS_l1b_mds['Waveform_1Hz']['Day'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Second'][r] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Micsec'][r] = np.fromfile(fid,dtype='>u4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Lat'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Lon'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Alt'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['TD'][r] = np.fromfile(fid,dtype='>i8',count=1)
                CS_l1b_mds['Waveform_1Hz']['Waveform'][r,:] = np.fromfile(fid,dtype='>u2',count=n_SARIN_RW)
                CS_l1b_mds['Waveform_1Hz']['Linear_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['Power2_Wfm_Multiplier'][r] = np.fromfile(fid,dtype='>i4',count=1)
                CS_l1b_mds['Waveform_1Hz']['N_avg_echoes'][r] = np.fromfile(fid,dtype='>u2',count=1)
                CS_l1b_mds['Waveform_1Hz']['Flags'][r] = np.fromfile(fid,dtype='>u2',count=1)

            # CryoSat-2 Waveforms Groups
            if (self.MODE == 'LRM'):
                # Low-Resolution Mode
                for b in range(n_blocks):
                    CS_l1b_mds['Waveform_20Hz']['Waveform'][r,b,:] = np.fromfile(fid,dtype='>u2',count=n_LRM_RW)
                    CS_l1b_mds['Waveform_20Hz']['Linear_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Power2_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                    CS_l1b_mds['Waveform_20Hz']['N_avg_echoes'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Flags'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
            elif (self.MODE == 'SAR'):
                # SAR Mode
                for b in range(n_blocks):
                    CS_l1b_mds['Waveform_20Hz']['Waveform'][r,b,:] = np.fromfile(fid,dtype='>u2',count=n_SAR_BC_RW)
                    CS_l1b_mds['Waveform_20Hz']['Linear_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Power2_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                    CS_l1b_mds['Waveform_20Hz']['N_avg_echoes'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Flags'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['SD'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Center'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Amplitude'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Skewness'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Kurtosis'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['SD_boresight_angle'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Center_boresight_angle'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Spare'][r,b,:] = np.fromfile(fid,dtype='>i2',count=(n_BeamBehaviourParams-7))
            elif (self.MODE == 'SIN'):
                # SARIN Mode
                for b in range(n_blocks):
                    CS_l1b_mds['Waveform_20Hz']['Waveform'][r,b,:] = np.fromfile(fid,dtype='>u2',count=n_SARIN_BC_RW)
                    CS_l1b_mds['Waveform_20Hz']['Linear_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Power2_Wfm_Multiplier'][r,b] = np.fromfile(fid,dtype='>i4',count=1)
                    CS_l1b_mds['Waveform_20Hz']['N_avg_echoes'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Flags'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['SD'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Center'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Amplitude'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Skewness'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Kurtosis'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['SD_boresight_angle'][r,b] = np.fromfile(fid,dtype='>u2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Center_boresight_angle'][r,b] = np.fromfile(fid,dtype='>i2',count=1)
                    CS_l1b_mds['Waveform_20Hz']['Beam']['Spare'][r,b,:] = np.fromfile(fid,dtype='>i2',count=(n_BeamBehaviourParams-7))
                    CS_l1b_mds['Waveform_20Hz']['Coherence'][r,b,:] = np.fromfile(fid,dtype='>i2',count=n_SARIN_BC_RW)
                    CS_l1b_mds['Waveform_20Hz']['Phase_diff'][r,b,:] = np.fromfile(fid,dtype='>i4',count=n_SARIN_BC_RW)

        # set the mask from day variables
        mask_20Hz = CS_l1b_mds['Location']['Day'].data == CS_l1b_mds['Location']['Day'].fill_value
        Location_keys = [key for key in CS_l1b_mds['Location'].keys() if not re.search(r'Spare',key)]
        Data_keys = [key for key in CS_l1b_mds['Data'].keys() if not re.search(r'Spare',key)]
        Geometry_keys = [key for key in CS_l1b_mds['Geometry'].keys() if not re.search(r'Spare',key)]
        Wfm_1Hz_keys = [key for key in CS_l1b_mds['Waveform_1Hz'].keys() if not re.search(r'Spare',key)]
        Wfm_20Hz_keys = [key for key in CS_l1b_mds['Waveform_20Hz'].keys() if not re.search(r'Spare',key)]
        for key in Location_keys:
            CS_l1b_mds['Location'][key].mask = mask_20Hz.copy()
        for key in Data_keys:
            CS_l1b_mds['Data'][key].mask = mask_20Hz.copy()
        # return the output dictionary
        return CS_l1b_mds

    def cryosat_baseline_D(self, full_filename, unpack=False):
        """
        Read L1b MDS variables for CryoSat Baseline D (netCDF4)
        """
        # open netCDF4 file for reading
        fid = netCDF4.Dataset(os.path.expanduser(full_filename),'r')
        # use original unscaled units unless unpack=True
        fid.set_auto_scale(unpack)
        # get dimensions
        ind_first_meas_20hz_01 = fid.variables['ind_first_meas_20hz_01'][:].copy()
        ind_meas_1hz_20_ku = fid.variables['ind_meas_1hz_20_ku'][:].copy()
        n_records = len(ind_first_meas_20hz_01)
        n_SARIN_D_RW = 1024
        n_SARIN_RW = 512
        n_SAR_D_RW = 256
        n_SAR_RW = 128
        n_LRM_RW = 128
        n_blocks = 20

        # Bind all the variables of the l1b_mds together into a single dictionary
        CS_l1b_mds = {}

        # CryoSat-2 Time and Orbit Group
        CS_l1b_mds['Location'] = {}
        # MDS Time
        CS_l1b_mds['Location']['Time'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Time'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        time_20_ku = fid.variables['time_20_ku'][:].copy()
        # Time: day part
        CS_l1b_mds['Location']['Day'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Day'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        # Time: second part
        CS_l1b_mds['Location']['Second'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Second'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        # Time: microsecond part
        CS_l1b_mds['Location']['Micsec'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Micsec'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        # USO correction factor
        CS_l1b_mds['Location']['USO_Corr'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['USO_Corr'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        uso_cor_20_ku = fid.variables['uso_cor_20_ku'][:].copy()
        # Mode ID
        CS_l1b_mds['Location']['Mode_ID'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Mode_ID'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        flag_instr_mode_op_20_ku =fid.variables['flag_instr_mode_op_20_ku'][:].copy()
        # Mode Flags
        CS_l1b_mds['Location']['Mode_flags'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Mode_flags'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        flag_instr_mode_flags_20_ku =fid.variables['flag_instr_mode_flags_20_ku'][:].copy()
        # Platform attitude control mode
        CS_l1b_mds['Location']['Att_control'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Att_control'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        flag_instr_mode_att_ctrl_20_ku =fid.variables['flag_instr_mode_att_ctrl_20_ku'][:].copy()
        # Instrument configuration
        CS_l1b_mds['Location']['Inst_config'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Inst_config'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        flag_instr_conf_rx_flags_20_ku = fid.variables['flag_instr_conf_rx_flags_20_ku'][:].copy()
        # acquisition band
        CS_l1b_mds['Location']['Inst_band'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Inst_band'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        flag_instr_conf_rx_bwdt_20_ku = fid.variables['flag_instr_conf_rx_bwdt_20_ku'][:].copy()
        # instrument channel
        CS_l1b_mds['Location']['Inst_channel'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Inst_channel'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        flag_instr_conf_rx_in_use_20_ku = fid.variables['flag_instr_conf_rx_in_use_20_ku'][:].copy()
        # tracking mode
        CS_l1b_mds['Location']['Tracking_mode'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Tracking_mode'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        flag_instr_conf_rx_trk_mode_20_ku = fid.variables['flag_instr_conf_rx_trk_mode_20_ku'][:].copy()
        # Source sequence counter
        CS_l1b_mds['Location']['SSC'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['SSC'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        seq_count_20_ku = fid.variables['seq_count_20_ku'][:].copy()
        # Record Counter
        CS_l1b_mds['Location']['Rec_Count'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Rec_Count'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        rec_count_20_ku = fid.variables['rec_count_20_ku'][:].copy()
        # Lat: packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_mds['Location']['Lat'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Lat'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        lat_20_ku = fid.variables['lat_20_ku'][:].copy()
        # Lon: packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_mds['Location']['Lon'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Lon'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        lon_20_ku = fid.variables['lon_20_ku'][:].copy()
        # Alt: packed units (mm, 1e-3 m)
        # Altitude of COG above reference ellipsoid (interpolated value)
        CS_l1b_mds['Location']['Alt'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Alt'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        alt_20_ku = fid.variables['alt_20_ku'][:].copy()
        # Instantaneous altitude rate derived from orbit: packed units (mm/s, 1e-3 m/s)
        CS_l1b_mds['Location']['Alt_rate'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Alt_rate'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        orb_alt_rate_20_ku = fid.variables['orb_alt_rate_20_ku'][:].copy()
        # Satellite velocity vector. In ITRF: packed units (mm/s, 1e-3 m/s)
        # ITRF= International Terrestrial Reference Frame
        CS_l1b_mds['Location']['Sat_velocity'] = np.ma.zeros((n_records,n_blocks,3))
        CS_l1b_mds['Location']['Sat_velocity'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        sat_vel_vec_20_ku = fid.variables['sat_vel_vec_20_ku'][:].copy()
        # Real beam direction vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
        # CRF= CryoSat Reference Frame.
        CS_l1b_mds['Location']['Real_beam'] = np.ma.zeros((n_records,n_blocks,3))
        CS_l1b_mds['Location']['Real_beam'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        beam_dir_vec_20_ku = fid.variables['beam_dir_vec_20_ku'][:].copy()
        # Interferometric baseline vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
        CS_l1b_mds['Location']['Baseline'] = np.ma.zeros((n_records,n_blocks,3))
        CS_l1b_mds['Location']['Baseline'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        inter_base_vec_20_ku = fid.variables['inter_base_vec_20_ku'][:].copy()
        # Star Tracker ID
        CS_l1b_mds['Location']['ST_ID'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['ST_ID'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        flag_instr_conf_rx_str_in_use_20_ku = fid.variables['flag_instr_conf_rx_str_in_use_20_ku'][:].copy()
        # Antenna Bench Roll Angle (Derived from star trackers)
        # packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_mds['Location']['Roll'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Roll'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        off_nadir_roll_angle_str_20_ku = fid.variables['off_nadir_roll_angle_str_20_ku'][:].copy()
        # Antenna Bench Pitch Angle (Derived from star trackers)
        # packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_mds['Location']['Pitch'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Pitch'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        off_nadir_pitch_angle_str_20_ku = fid.variables['off_nadir_pitch_angle_str_20_ku'][:].copy()
        # Antenna Bench Yaw Angle (Derived from star trackers)
        # packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_mds['Location']['Yaw'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['Yaw'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        off_nadir_yaw_angle_str_20_ku = fid.variables['off_nadir_yaw_angle_str_20_ku'][:].copy()
        # Measurement Confidence Data Flags
        # Generally the MCD flags indicate problems when set
        # If MCD is 0 then no problems or non-nominal conditions were detected
        # Serious errors are indicated by setting bit 31
        CS_l1b_mds['Location']['MCD'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Location']['MCD'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        flag_mcd_20_ku = fid.variables['flag_mcd_20_ku'][:].copy()

        # CryoSat-2 Measurement Group
        # Derived from instrument measurement parameters
        CS_l1b_mds['Data'] = {}
        # Window Delay reference (two-way) corrected for instrument delays
        CS_l1b_mds['Data']['TD'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['TD'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        window_del_20_ku = fid.variables['window_del_20_ku'][:].copy()
        # H0 Initial Height Word from telemetry
        CS_l1b_mds['Data']['H_0'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['H_0'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        h0_applied_20_ku = fid.variables['h0_applied_20_ku'][:].copy()
        # COR2 Height Rate: on-board tracker height rate over the radar cycle
        CS_l1b_mds['Data']['COR2'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['COR2'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        cor2_applied_20_ku = fid.variables['cor2_applied_20_ku'][:].copy()
        # Coarse Range Word (LAI) derived from telemetry
        CS_l1b_mds['Data']['LAI'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['LAI'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        h0_lai_word_20_ku = fid.variables['h0_lai_word_20_ku'][:].copy()
        # Fine Range Word (FAI) derived from telemetry
        CS_l1b_mds['Data']['FAI'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['FAI'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        h0_fai_word_20_ku = fid.variables['h0_fai_word_20_ku'][:].copy()
        # Automatic Gain Control Channel 1: AGC gain applied on Rx channel 1.
        # Gain calibration corrections are applied (Sum of AGC stages 1 and 2
        # plus the corresponding corrections) (dB/100)
        CS_l1b_mds['Data']['AGC_CH1'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['AGC_CH1'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        agc_ch1_20_ku = fid.variables['agc_ch1_20_ku'][:].copy()
        # Automatic Gain Control Channel 2: AGC gain applied on Rx channel 2.
        # Gain calibration corrections are applied (dB/100)
        CS_l1b_mds['Data']['AGC_CH2'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['AGC_CH2'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        agc_ch2_20_ku = fid.variables['agc_ch2_20_ku'][:].copy()
        # Total Fixed Gain On Channel 1: gain applied by the RF unit. (dB/100)
        CS_l1b_mds['Data']['TR_gain_CH1'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['TR_gain_CH1'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        tot_gain_ch1_20_ku = fid.variables['tot_gain_ch1_20_ku'][:].copy()
        # Total Fixed Gain On Channel 2: gain applied by the RF unit. (dB/100)
        CS_l1b_mds['Data']['TR_gain_CH2'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['TR_gain_CH2'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        tot_gain_ch2_20_ku = fid.variables['tot_gain_ch2_20_ku'][:].copy()
        # Transmit Power in microWatts
        CS_l1b_mds['Data']['TX_Power'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['TX_Power'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        transmit_pwr_20_ku = fid.variables['transmit_pwr_20_ku'][:].copy()
        # Doppler range correction: Radial component (mm)
        # computed for the component of satellite velocity in the nadir direction
        CS_l1b_mds['Data']['Doppler_range'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['Doppler_range'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        dop_cor_20_ku = fid.variables['dop_cor_20_ku'][:].copy()
        # Value of Doppler Angle for the first single look echo (1e-7 radians)
        CS_l1b_mds['Data']['Doppler_angle_start'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['Doppler_angle_start'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        dop_angle_start_20_ku = fid.variables['dop_angle_start_20_ku'][:].copy()
        # Value of Doppler Angle for the last single look echo (1e-7 radians)
        CS_l1b_mds['Data']['Doppler_angle_stop'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['Doppler_angle_stop'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        dop_angle_stop_20_ku = fid.variables['dop_angle_stop_20_ku'][:].copy()
        # Instrument Range Correction: transmit-receive antenna (mm)
        # Calibration correction to range on channel 1 computed from CAL1.
        CS_l1b_mds['Data']['TR_inst_range'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['TR_inst_range'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        instr_cor_range_tx_rx_20_ku = fid.variables['instr_cor_range_tx_rx_20_ku'][:].copy()
        # Instrument Range Correction: receive-only antenna (mm)
        # Calibration correction to range on channel 2 computed from CAL1.
        CS_l1b_mds['Data']['R_inst_range'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['R_inst_range'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        instr_cor_range_rx_20_ku = fid.variables['instr_cor_range_rx_20_ku'][:].copy()
        # Instrument Gain Correction: transmit-receive antenna (dB/100)
        # Calibration correction to gain on channel 1 computed from CAL1
        CS_l1b_mds['Data']['TR_inst_gain'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['TR_inst_gain'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        instr_cor_gain_tx_rx_20_ku = fid.variables['instr_cor_gain_tx_rx_20_ku'][:].copy()
        # Instrument Gain Correction: receive-only (dB/100)
        # Calibration correction to gain on channel 2 computed from CAL1
        CS_l1b_mds['Data']['R_inst_gain'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['R_inst_gain'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        instr_cor_gain_rx_20_ku = fid.variables['instr_cor_gain_rx_20_ku'][:].copy()
        # Internal Phase Correction (microradians)
        CS_l1b_mds['Data']['Internal_phase'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['Internal_phase'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        instr_int_ph_cor_20_ku = fid.variables['instr_int_ph_cor_20_ku'][:].copy()
        # External Phase Correction (microradians)
        CS_l1b_mds['Data']['External_phase'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['External_phase'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        instr_ext_ph_cor_20_ku = fid.variables['instr_ext_ph_cor_20_ku'][:].copy()
        # Noise Power measurement (dB/100)
        CS_l1b_mds['Data']['Noise_power'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['Noise_power'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        noise_power_20_ku = fid.variables['noise_power_20_ku'][:].copy()
        # Phase slope correction (microradians)
        # Computed from the CAL-4 packets during the azimuth impulse response
        # amplitude (SARIN only). Set from the latest available CAL-4 packet.
        CS_l1b_mds['Data']['Phase_slope'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Data']['Phase_slope'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        ph_slope_cor_20_ku = fid.variables['ph_slope_cor_20_ku'][:].copy()

        # CryoSat-2 External Corrections Group
        CS_l1b_mds['Geometry'] = {}
        # Data Record Time (MDSR Time Stamp)
        CS_l1b_mds['Geometry']['Time'] = fid.variables['time_cor_01'][:].copy()
        # Dry Tropospheric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['dryTrop'] = fid.variables['mod_dry_tropo_cor_01'][:].copy()
        # Wet Tropospheric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['wetTrop'] = fid.variables['mod_wet_tropo_cor_01'][:].copy()
        # Inverse Barometric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['InvBar'] = fid.variables['inv_bar_cor_01'][:].copy()
        # Delta Inverse Barometric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['DAC'] = fid.variables['hf_fluct_total_cor_01'][:].copy()
        # GIM Ionospheric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['Iono_GIM'] = fid.variables['iono_cor_gim_01'][:].copy()
        # Model Ionospheric Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['Iono_model'] = fid.variables['iono_cor_01'][:].copy()
        # Ocean tide Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['ocTideElv'] = fid.variables['ocean_tide_01'][:].copy()
        # Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['lpeTideElv'] = fid.variables['ocean_tide_eq_01'][:].copy()
        # Ocean loading tide Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['olTideElv'] = fid.variables['load_tide_01'][:].copy()
        # Solid Earth tide Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['seTideElv'] = fid.variables['solid_earth_tide_01'][:].copy()
        # Geocentric Polar tide Correction packed units (mm, 1e-3 m)
        CS_l1b_mds['Geometry']['gpTideElv'] = fid.variables['pole_tide_01'][:].copy()
        # Surface Type: enumerated key to classify surface at nadir
        # 0 = Open Ocean
        # 1 = Closed Sea
        # 2 = Continental Ice
        # 3 = Land
        CS_l1b_mds['Geometry']['Surf_type'] = fid.variables['surf_type_01'][:].copy()
        # Corrections Status Flag
        CS_l1b_mds['Geometry']['Corr_status'] = fid.variables['flag_cor_status_01'][:].copy()
        # Correction Error Flag
        CS_l1b_mds['Geometry']['Corr_error'] = fid.variables['flag_cor_err_01'][:].copy()

        # Same as the LRM/SAR groups but the waveform array is 512 bins instead of
        # 128 and the number of echoes averaged is different.
        CS_l1b_mds['Waveform_1Hz'] = {}
        # Data Record Time (MDSR Time Stamp)
        # Time (seconds since 2000-01-01)
        time_avg_01_ku = fid.variables['time_avg_01_ku'][:].copy()
        CS_l1b_mds['Waveform_1Hz']['Time'] = time_avg_01_ku.copy()
        # Time: day part
        CS_l1b_mds['Waveform_1Hz']['Day'] = np.array(time_avg_01_ku/86400.0, dtype=np.int32)
        # Time: second part
        CS_l1b_mds['Waveform_1Hz']['Second'] = np.array(time_avg_01_ku -
            CS_l1b_mds['Waveform_1Hz']['Day'][:]*86400.0, dtype=np.uint32)
        # Time: microsecond part
        CS_l1b_mds['Waveform_1Hz']['Micsec'] = np.array((time_avg_01_ku -
            CS_l1b_mds['Waveform_1Hz']['Day'][:]*86400.0 -
            CS_l1b_mds['Waveform_1Hz']['Second'][:])*1e6, dtype=np.uint32)
        # Lat: packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_mds['Waveform_1Hz']['Lat'] = fid.variables['lat_avg_01_ku'][:].copy()
        # Lon: packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_mds['Waveform_1Hz']['Lon'] = fid.variables['lon_avg_01_ku'][:].copy()
        # Alt: packed units (mm, 1e-3 m)
        # Altitude of COG above reference ellipsoid (interpolated value)
        CS_l1b_mds['Waveform_1Hz']['Alt'] = fid.variables['alt_avg_01_ku'][:].copy()
        # Window Delay (two-way) corrected for instrument delays
        CS_l1b_mds['Waveform_1Hz']['TD'] = fid.variables['window_del_avg_01_ku'][:].copy()
        # 1 Hz Averaged Power Echo Waveform
        CS_l1b_mds['Waveform_1Hz']['Waveform'] = fid.variables['pwr_waveform_avg_01_ku'][:].copy()
        # Echo Scale Factor (to scale echo to watts)
        CS_l1b_mds['Waveform_1Hz']['Linear_Wfm_Multiplier'] = fid.variables['echo_scale_factor_avg_01_ku'][:].copy()
        # Echo Scale Power (a power of 2 to scale echo to Watts)
        CS_l1b_mds['Waveform_1Hz']['Power2_Wfm_Multiplier'] = fid.variables['echo_scale_pwr_avg_01_ku'][:].copy()
        # Number of echoes averaged
        CS_l1b_mds['Waveform_1Hz']['N_avg_echoes'] = fid.variables['echo_numval_avg_01_ku'][:].copy()
        CS_l1b_mds['Waveform_1Hz']['Flags'] = fid.variables['flag_echo_avg_01_ku'][:].copy()

        # CryoSat-2 Waveforms Groups
        CS_l1b_mds['Waveform_20Hz'] = {}
        # Echo Scale Factor (to scale echo to watts)
        CS_l1b_mds['Waveform_20Hz']['Linear_Wfm_Multiplier'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Waveform_20Hz']['Linear_Wfm_Multiplier'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        echo_scale_factor_20_ku = fid.variables['echo_scale_factor_20_ku'][:].copy()
        # Echo Scale Power (a power of 2)
        CS_l1b_mds['Waveform_20Hz']['Power2_Wfm_Multiplier'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Waveform_20Hz']['Power2_Wfm_Multiplier'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        echo_scale_pwr_20_ku = fid.variables['echo_scale_pwr_20_ku'][:].copy()
        # Number of echoes averaged
        CS_l1b_mds['Waveform_20Hz']['N_avg_echoes'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Waveform_20Hz']['N_avg_echoes'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        echo_numval_20_ku = fid.variables['echo_numval_20_ku'][:].copy()
        # Flags for errors or information about 20Hz waveform
        CS_l1b_mds['Waveform_20Hz']['Flags'] = np.ma.zeros((n_records,n_blocks))
        CS_l1b_mds['Waveform_20Hz']['Flags'].mask = np.zeros((n_records,n_blocks),dtype=bool)
        flag_echo_20_ku = fid.variables['flag_echo_20_ku'][:].copy()
        # CryoSat-2 mode specific waveform variables
        if (self.MODE == 'LRM'):
            # Low-Resolution Mode
            # Averaged Power Echo Waveform [128]
            CS_l1b_mds['Waveform_20Hz']['Waveform'] = np.ma.zeros((n_records,n_blocks,n_LRM_RW))
            CS_l1b_mds['Waveform_20Hz']['Waveform'].mask = np.zeros((n_records,n_blocks,n_LRM_RW),dtype=bool)
            pwr_waveform_20_ku = fid.variables['pwr_waveform_20_ku'][:].copy()
        elif (self.MODE == 'SAR'):
            # SAR Mode
            # Averaged Power Echo Waveform [256]
            CS_l1b_mds['Waveform_20Hz']['Waveform'] = np.ma.zeros((n_records,n_blocks,n_SAR_D_RW))
            CS_l1b_mds['Waveform_20Hz']['Waveform'].mask = np.zeros((n_records,n_blocks,n_SAR_D_RW),dtype=bool)
            pwr_waveform_20_ku = fid.variables['pwr_waveform_20_ku'][:].copy()
        elif (self.MODE == 'SIN'):
            # SARIN Mode
            # Averaged Power Echo Waveform [1024]
            CS_l1b_mds['Waveform_20Hz']['Waveform'] = np.ma.zeros((n_records,n_blocks,n_SARIN_D_RW))
            CS_l1b_mds['Waveform_20Hz']['Waveform'].mask = np.zeros((n_records,n_blocks,n_SARIN_D_RW),dtype=bool)
            pwr_waveform_20_ku = fid.variables['pwr_waveform_20_ku'][:].copy()
            # Coherence [1024]: packed units (1/1000)
            CS_l1b_mds['Waveform_20Hz']['Coherence'] = np.ma.zeros((n_records,n_blocks,n_SARIN_D_RW))
            CS_l1b_mds['Waveform_20Hz']['Coherence'].mask = np.zeros((n_records,n_blocks,n_SARIN_D_RW),dtype=bool)
            coherence_waveform_20_ku = fid.variables['coherence_waveform_20_ku'][:].copy()
            # Phase Difference [1024]: packed units (microradians)
            CS_l1b_mds['Waveform_20Hz']['Phase_diff'] = np.ma.zeros((n_records,n_blocks,n_SARIN_D_RW))
            CS_l1b_mds['Waveform_20Hz']['Phase_diff'].mask = np.zeros((n_records,n_blocks,n_SARIN_D_RW),dtype=bool)
            ph_diff_waveform_20_ku = fid.variables['ph_diff_waveform_20_ku'][:].copy()

        # Beam Behavior Parameters
        if self.MODE in ('SAR','SIN'):
            CS_l1b_mds['Waveform_20Hz']['Beam'] = {}
            # Standard Deviation of Gaussian fit to range integrated stack power.
            CS_l1b_mds['Waveform_20Hz']['Beam']['SD'] = np.ma.zeros((n_records,n_blocks))
            CS_l1b_mds['Waveform_20Hz']['Beam']['SD'].mask = np.zeros((n_records,n_blocks),dtype=bool)
            stack_std_20_ku = fid.variables['stack_std_20_ku'][:].copy()
            # Stack Center: Mean of Gaussian fit to range integrated stack power.
            CS_l1b_mds['Waveform_20Hz']['Beam']['Center'] = np.ma.zeros((n_records,n_blocks))
            CS_l1b_mds['Waveform_20Hz']['Beam']['Center'].mask = np.zeros((n_records,n_blocks),dtype=bool)
            stack_centre_20_ku = fid.variables['stack_centre_20_ku'][:].copy()
            # Stack amplitude parameter scaled in dB/100.
            CS_l1b_mds['Waveform_20Hz']['Beam']['Amplitude'] = np.ma.zeros((n_records,n_blocks))
            CS_l1b_mds['Waveform_20Hz']['Beam']['Amplitude'].mask = np.zeros((n_records,n_blocks),dtype=bool)
            stack_scaled_amplitude_20_ku = fid.variables['stack_scaled_amplitude_20_ku'][:].copy()
            # 3rd moment: providing the degree of asymmetry of the range integrated
            # stack power distribution.
            CS_l1b_mds['Waveform_20Hz']['Beam']['Skewness'] = np.ma.zeros((n_records,n_blocks))
            CS_l1b_mds['Waveform_20Hz']['Beam']['Skewness'].mask = np.zeros((n_records,n_blocks),dtype=bool)
            stack_skewness_20_ku = fid.variables['stack_skewness_20_ku'][:].copy()
            # 4th moment: Measure of peakiness of range integrated stack power distribution.
            CS_l1b_mds['Waveform_20Hz']['Beam']['Kurtosis'] = np.ma.zeros((n_records,n_blocks))
            CS_l1b_mds['Waveform_20Hz']['Beam']['Kurtosis'].mask = np.zeros((n_records,n_blocks),dtype=bool)
            stack_kurtosis_20_ku = fid.variables['stack_kurtosis_20_ku'][:].copy()
            # Stack peakiness computed from the range integrated power of the single look echoes
            CS_l1b_mds['Waveform_20Hz']['Beam']['Peakiness'] = np.ma.zeros((n_records,n_blocks))
            CS_l1b_mds['Waveform_20Hz']['Beam']['Peakiness'].mask = np.zeros((n_records,n_blocks),dtype=bool)
            stack_peakiness_20_ku = fid.variables['stack_peakiness_20_ku'][:].copy()
            # Stack residuals of Gaussian that fits the range integrated power of the single look echoes
            CS_l1b_mds['Waveform_20Hz']['Beam']['RMS'] = np.ma.zeros((n_records,n_blocks))
            CS_l1b_mds['Waveform_20Hz']['Beam']['RMS'].mask = np.zeros((n_records,n_blocks),dtype=bool)
            stack_gaussian_fitting_residuals_20_ku = fid.variables['stack_gaussian_fitting_residuals_20_ku'][:].copy()
            # Standard deviation as a function of boresight angle (microradians)
            CS_l1b_mds['Waveform_20Hz']['Beam']['SD_boresight_angle'] = np.ma.zeros((n_records,n_blocks))
            CS_l1b_mds['Waveform_20Hz']['Beam']['SD_boresight_angle'].mask = np.zeros((n_records,n_blocks),dtype=bool)
            stack_std_angle_20_ku = fid.variables['stack_std_angle_20_ku'][:].copy()
            # Stack Center angle as a function of boresight angle (microradians)
            CS_l1b_mds['Waveform_20Hz']['Beam']['Center_boresight_angle'] = np.ma.zeros((n_records,n_blocks))
            CS_l1b_mds['Waveform_20Hz']['Beam']['Center_boresight_angle'].mask = np.zeros((n_records,n_blocks),dtype=bool)
            stack_centre_angle_20_ku = fid.variables['stack_centre_angle_20_ku'][:].copy()
            # Stack Center angle as a function of look angle (microradians)
            CS_l1b_mds['Waveform_20Hz']['Beam']['Center_look_angle'] = np.ma.zeros((n_records,n_blocks))
            CS_l1b_mds['Waveform_20Hz']['Beam']['Center_look_angle'].mask = np.zeros((n_records,n_blocks),dtype=bool)
            stack_centre_look_angle_20_ku = fid.variables['stack_centre_look_angle_20_ku'][:].copy()
            # Number of contributing beams in the stack before weighting
            CS_l1b_mds['Waveform_20Hz']['Beam']['Number'] = np.ma.zeros((n_records,n_blocks))
            CS_l1b_mds['Waveform_20Hz']['Beam']['Number'].mask = np.zeros((n_records,n_blocks),dtype=bool)
            stack_number_before_weighting_20_ku = fid.variables['stack_number_before_weighting_20_ku'][:].copy()
            # Number of contributing beams in the stack after weighting
            CS_l1b_mds['Waveform_20Hz']['Beam']['Weighted_Number'] = np.ma.zeros((n_records,n_blocks))
            CS_l1b_mds['Waveform_20Hz']['Beam']['Weighted_Number'].mask = np.zeros((n_records,n_blocks),dtype=bool)
            stack_number_after_weighting_20_ku = fid.variables['stack_number_after_weighting_20_ku'][:].copy()

        # for each record in the CryoSat file
        for r in range(n_records):
            # index for record r
            idx = ind_first_meas_20hz_01[r]
            # number of valid blocks in record r
            cnt = np.count_nonzero(ind_meas_1hz_20_ku == r)

            # CryoSat-2 Time and Orbit Group
            CS_l1b_mds['Location']['Time'].data[r,:cnt] = time_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Time'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Day'].data[r,:cnt] = np.array(time_20_ku[idx:idx+cnt]/86400.0, dtype=np.int64)
            CS_l1b_mds['Location']['Day'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Second'].data[r,:cnt] = np.array(time_20_ku[idx:idx+cnt] -
                CS_l1b_mds['Location']['Day'].data[r,:cnt]*86400.0, dtype=np.int64)
            CS_l1b_mds['Location']['Second'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Micsec'].data[r,:cnt] = np.array((time_20_ku[idx:idx+cnt] -
                CS_l1b_mds['Location']['Day'].data[r,:cnt]*86400.0 -
                CS_l1b_mds['Location']['Second'].data[r,:cnt])*1e6, dtype=np.uint32)
            CS_l1b_mds['Location']['Micsec'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['USO_Corr'].data[r,:cnt] = uso_cor_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['USO_Corr'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Mode_ID'].data[r,:cnt] = flag_instr_mode_op_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Mode_ID'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Mode_flags'].data[r,:cnt] = flag_instr_mode_flags_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Mode_flags'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Att_control'].data[r,:cnt] = flag_instr_mode_att_ctrl_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Att_control'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Inst_config'].data[r,:cnt] = flag_instr_conf_rx_flags_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Inst_config'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Inst_band'].data[r,:cnt] = flag_instr_conf_rx_bwdt_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Inst_band'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Inst_channel'].data[r,:cnt] = flag_instr_conf_rx_in_use_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Inst_channel'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Tracking_mode'].data[r,:cnt] = flag_instr_conf_rx_trk_mode_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Tracking_mode'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['SSC'].data[r,:cnt] = seq_count_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['SSC'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Rec_Count'].data[r,:cnt] = rec_count_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Rec_Count'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Lat'].data[r,:cnt] = lat_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Lat'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Lon'].data[r,:cnt] = lon_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Lon'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Alt'].data[r,:cnt] = alt_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Alt'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Alt_rate'].data[r,:cnt] = orb_alt_rate_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Alt_rate'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Sat_velocity'].data[r,:cnt,:] = sat_vel_vec_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Sat_velocity'].mask[r,:cnt,:] = False
            CS_l1b_mds['Location']['Real_beam'].data[r,:cnt,:] = beam_dir_vec_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Real_beam'].mask[r,:cnt,:] = False
            CS_l1b_mds['Location']['Baseline'].data[r,:cnt,:] = inter_base_vec_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Baseline'].mask[r,:cnt,:] = False
            CS_l1b_mds['Location']['ST_ID'].data[r,:cnt] = flag_instr_conf_rx_str_in_use_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['ST_ID'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Roll'].data[r,:cnt] = off_nadir_roll_angle_str_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Roll'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Pitch'].data[r,:cnt] = off_nadir_pitch_angle_str_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Pitch'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['Yaw'].data[r,:cnt] = off_nadir_yaw_angle_str_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['Yaw'].mask[r,:cnt] = False
            CS_l1b_mds['Location']['MCD'].data[r,:cnt] = flag_mcd_20_ku[idx:idx+cnt]
            CS_l1b_mds['Location']['MCD'].mask[r,:cnt] = False

            # CryoSat-2 Measurement Group
            # Derived from instrument measurement parameters
            CS_l1b_mds['Data']['TD'].data[r,:cnt] = window_del_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['TD'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['H_0'].data[r,:cnt] = h0_applied_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['H_0'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['COR2'].data[r,:cnt] = cor2_applied_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['COR2'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['LAI'].data[r,:cnt] = h0_lai_word_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['LAI'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['FAI'].data[r,:cnt] = h0_fai_word_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['FAI'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['AGC_CH1'].data[r,:cnt] = agc_ch1_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['AGC_CH1'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['AGC_CH2'].data[r,:cnt] = agc_ch2_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['AGC_CH2'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['TR_gain_CH1'].data[r,:cnt] = tot_gain_ch1_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['TR_gain_CH1'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['TR_gain_CH2'].data[r,:cnt] = tot_gain_ch2_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['TR_gain_CH2'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['TX_Power'].data[r,:cnt] = transmit_pwr_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['TX_Power'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['Doppler_range'].data[r,:cnt] = dop_cor_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['Doppler_range'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['Doppler_angle_start'].data[r,:cnt] = dop_angle_start_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['Doppler_angle_start'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['Doppler_angle_stop'].data[r,:cnt] = dop_angle_stop_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['Doppler_angle_stop'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['TR_inst_range'].data[r,:cnt] = instr_cor_range_tx_rx_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['TR_inst_range'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['R_inst_range'].data[r,:cnt] = instr_cor_range_rx_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['R_inst_range'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['TR_inst_gain'].data[r,:cnt] = instr_cor_gain_tx_rx_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['TR_inst_gain'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['R_inst_gain'].data[r,:cnt] = instr_cor_gain_rx_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['R_inst_gain'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['Internal_phase'].data[r,:cnt] = instr_int_ph_cor_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['Internal_phase'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['External_phase'].data[r,:cnt] = instr_ext_ph_cor_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['External_phase'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['Noise_power'].data[r,:cnt] = noise_power_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['Noise_power'].mask[r,:cnt] = False
            CS_l1b_mds['Data']['Phase_slope'].data[r,:cnt] = ph_slope_cor_20_ku[idx:idx+cnt]
            CS_l1b_mds['Data']['Phase_slope'].mask[r,:cnt] = False

            # CryoSat-2 Waveforms Groups
            CS_l1b_mds['Waveform_20Hz']['Linear_Wfm_Multiplier'].data[r,:cnt] = echo_scale_factor_20_ku[idx:idx+cnt]
            CS_l1b_mds['Waveform_20Hz']['Linear_Wfm_Multiplier'].mask[r,:cnt] = False
            CS_l1b_mds['Waveform_20Hz']['Power2_Wfm_Multiplier'].data[r,:cnt] = echo_scale_pwr_20_ku[idx:idx+cnt]
            CS_l1b_mds['Waveform_20Hz']['Power2_Wfm_Multiplier'].mask[r,:cnt] = False
            CS_l1b_mds['Waveform_20Hz']['N_avg_echoes'].data[r,:cnt] = echo_numval_20_ku[idx:idx+cnt]
            CS_l1b_mds['Waveform_20Hz']['N_avg_echoes'].mask[r,:cnt] = False
            CS_l1b_mds['Waveform_20Hz']['Flags'].data[r,:cnt] = flag_echo_20_ku[idx:idx+cnt]
            CS_l1b_mds['Waveform_20Hz']['Flags'].mask[r,:cnt] = False
            CS_l1b_mds['Waveform_20Hz']['Waveform'].data[r,:cnt,:] = pwr_waveform_20_ku[idx:idx+cnt,:]
            CS_l1b_mds['Waveform_20Hz']['Waveform'].mask[r,:cnt,:] = False
            # SARIN Mode parameters
            if (self.MODE == 'SIN'):
                CS_l1b_mds['Waveform_20Hz']['Coherence'].data[r,:cnt,:] = coherence_waveform_20_ku[idx:idx+cnt,:]
                CS_l1b_mds['Waveform_20Hz']['Coherence'].mask[r,:cnt,:] = False
                CS_l1b_mds['Waveform_20Hz']['Phase_diff'].data[r,:cnt,:] = ph_diff_waveform_20_ku[idx:idx+cnt,:]
                CS_l1b_mds['Waveform_20Hz']['Phase_diff'].mask[r,:cnt,:] = False
            # SAR/SARIN waveform beam parameters
            if self.MODE in ('SAR','SIN'):
                CS_l1b_mds['Waveform_20Hz']['Beam']['SD'].data[r,:cnt] = stack_std_20_ku[idx:idx+cnt]
                CS_l1b_mds['Waveform_20Hz']['Beam']['SD'].mask[r,:cnt] = False
                CS_l1b_mds['Waveform_20Hz']['Beam']['Center'].data[r,:cnt] = stack_centre_20_ku[idx:idx+cnt]
                CS_l1b_mds['Waveform_20Hz']['Beam']['Center'].mask[r,:cnt] = False
                CS_l1b_mds['Waveform_20Hz']['Beam']['Amplitude'].data[r,:cnt] = stack_scaled_amplitude_20_ku[idx:idx+cnt]
                CS_l1b_mds['Waveform_20Hz']['Beam']['Amplitude'].mask[r,:cnt] = False
                CS_l1b_mds['Waveform_20Hz']['Beam']['Skewness'].data[r,:cnt] = stack_skewness_20_ku[idx:idx+cnt]
                CS_l1b_mds['Waveform_20Hz']['Beam']['Skewness'].mask[r,:cnt] = False
                CS_l1b_mds['Waveform_20Hz']['Beam']['Kurtosis'].data[r,:cnt] = stack_kurtosis_20_ku[idx:idx+cnt]
                CS_l1b_mds['Waveform_20Hz']['Beam']['Kurtosis'].mask[r,:cnt] = False
                CS_l1b_mds['Waveform_20Hz']['Beam']['Peakiness'].data[r,:cnt] = stack_peakiness_20_ku[idx:idx+cnt]
                CS_l1b_mds['Waveform_20Hz']['Beam']['Peakiness'].mask[r,:cnt] = False
                CS_l1b_mds['Waveform_20Hz']['Beam']['RMS'].data[r,:cnt] = stack_gaussian_fitting_residuals_20_ku[idx:idx+cnt]
                CS_l1b_mds['Waveform_20Hz']['Beam']['RMS'].mask[r,:cnt] = False
                CS_l1b_mds['Waveform_20Hz']['Beam']['SD_boresight_angle'].data[r,:cnt] = stack_std_angle_20_ku[idx:idx+cnt]
                CS_l1b_mds['Waveform_20Hz']['Beam']['SD_boresight_angle'].mask[r,:cnt] = False
                CS_l1b_mds['Waveform_20Hz']['Beam']['Center_boresight_angle'].data[r,:cnt] = stack_centre_angle_20_ku[idx:idx+cnt]
                CS_l1b_mds['Waveform_20Hz']['Beam']['Center_boresight_angle'].mask[r,:cnt] = False
                CS_l1b_mds['Waveform_20Hz']['Beam']['Center_look_angle'].data[r,:cnt] = stack_centre_look_angle_20_ku[idx:idx+cnt]
                CS_l1b_mds['Waveform_20Hz']['Beam']['Center_look_angle'].mask[r,:cnt] = False
                CS_l1b_mds['Waveform_20Hz']['Beam']['Number'].data[r,:cnt] = stack_number_before_weighting_20_ku[idx:idx+cnt]
                CS_l1b_mds['Waveform_20Hz']['Beam']['Number'].mask[r,:cnt] = False
                CS_l1b_mds['Waveform_20Hz']['Beam']['Weighted_Number'].data[r,:cnt] = stack_number_after_weighting_20_ku[idx:idx+cnt]
                CS_l1b_mds['Waveform_20Hz']['Beam']['Weighted_Number'].mask[r,:cnt] = False

        # extract global attributes and assign as MPH and SPH metadata
        CS_l1b_mds['METADATA'] = dict(MPH={},SPH={},DSD={})
        # MPH attributes
        CS_l1b_mds['METADATA']['MPH']['PRODUCT'] = fid.product_name
        CS_l1b_mds['METADATA']['MPH']['DOI'] = fid.doi
        CS_l1b_mds['METADATA']['MPH']['PROC_STAGE'] =  fid.processing_stage
        CS_l1b_mds['METADATA']['MPH']['REF_DOC'] =  fid.reference_document
        CS_l1b_mds['METADATA']['MPH']['ACQUISITION_STATION'] = fid.acquisition_station
        CS_l1b_mds['METADATA']['MPH']['PROC_CENTER'] = fid.processing_centre
        CS_l1b_mds['METADATA']['MPH']['PROC_TIME'] = fid.creation_time
        CS_l1b_mds['METADATA']['MPH']['SOFTWARE_VER'] = fid.software_version
        CS_l1b_mds['METADATA']['MPH']['SENSING_START'] = fid.sensing_start
        CS_l1b_mds['METADATA']['MPH']['SENSING_STOP'] = fid.sensing_stop
        CS_l1b_mds['METADATA']['MPH']['PHASE'] = fid.phase
        CS_l1b_mds['METADATA']['MPH']['CYCLE'] = fid.cycle_number
        CS_l1b_mds['METADATA']['MPH']['REL_ORBIT'] = fid.rel_orbit_number
        CS_l1b_mds['METADATA']['MPH']['ABS_ORBIT'] = fid.abs_orbit_number
        CS_l1b_mds['METADATA']['MPH']['STATE_VECTOR_TIME'] = fid.state_vector_time
        CS_l1b_mds['METADATA']['MPH']['DELTA_UT1'] = fid.delta_ut1
        CS_l1b_mds['METADATA']['MPH']['X_POSITION'] = fid.x_position
        CS_l1b_mds['METADATA']['MPH']['Y_POSITION'] = fid.y_position
        CS_l1b_mds['METADATA']['MPH']['Z_POSITION'] = fid.z_position
        CS_l1b_mds['METADATA']['MPH']['X_VELOCITY'] = fid.x_velocity
        CS_l1b_mds['METADATA']['MPH']['Y_VELOCITY'] = fid.y_velocity
        CS_l1b_mds['METADATA']['MPH']['Z_VELOCITY'] = fid.z_velocity
        CS_l1b_mds['METADATA']['MPH']['VECTOR_SOURCE'] = fid.vector_source
        CS_l1b_mds['METADATA']['MPH']['LEAP_UTC'] = fid.leap_utc
        CS_l1b_mds['METADATA']['MPH']['LEAP_SIGN'] = fid.leap_sign
        CS_l1b_mds['METADATA']['MPH']['LEAP_ERR'] = fid.leap_err
        CS_l1b_mds['METADATA']['MPH']['PRODUCT_ERR'] = fid.product_err
        # SPH attributes
        CS_l1b_mds['METADATA']['SPH']['START_RECORD_TAI_TIME'] = fid.first_record_time
        CS_l1b_mds['METADATA']['SPH']['STOP_RECORD_TAI_TIME'] = fid.last_record_time
        CS_l1b_mds['METADATA']['SPH']['ABS_ORBIT_START'] = fid.abs_orbit_start
        CS_l1b_mds['METADATA']['SPH']['REL_TIME_ASC_NODE_START'] = fid.rel_time_acs_node_start
        CS_l1b_mds['METADATA']['SPH']['ABS_ORBIT_STOP'] = fid.abs_orbit_stop
        CS_l1b_mds['METADATA']['SPH']['REL_TIME_ASC_NODE_STOP'] = fid.rel_time_acs_node_stop
        CS_l1b_mds['METADATA']['SPH']['EQUATOR_CROSS_TIME_UTC'] = fid.equator_cross_time
        CS_l1b_mds['METADATA']['SPH']['EQUATOR_CROSS_LONG'] = fid.equator_cross_long
        CS_l1b_mds['METADATA']['SPH']['ASCENDING_FLAG'] = fid.ascending_flag
        CS_l1b_mds['METADATA']['SPH']['START_LAT'] = fid.first_record_lat
        CS_l1b_mds['METADATA']['SPH']['START_LONG'] = fid.first_record_lon
        CS_l1b_mds['METADATA']['SPH']['STOP_LAT'] = fid.last_record_lat
        CS_l1b_mds['METADATA']['SPH']['STOP_LONG'] = fid.last_record_lon
        CS_l1b_mds['METADATA']['SPH']['L0_PROC_FLAG'] = fid.l0_proc_flag
        CS_l1b_mds['METADATA']['SPH']['L0_PROCESSING_QUALITY'] = fid.l0_processing_quality
        CS_l1b_mds['METADATA']['SPH']['L0_PROC_THRESH'] = fid.l0_proc_thresh
        CS_l1b_mds['METADATA']['SPH']['L0_GAPS_FLAG'] = fid.l0_gaps_flag
        CS_l1b_mds['METADATA']['SPH']['L0_GAPS_NUM'] = fid.l0_gaps_num
        CS_l1b_mds['METADATA']['SPH']['INSTR_ID'] = fid.instr_id
        CS_l1b_mds['METADATA']['SPH']['OPEN_OCEAN_PERCENT'] = fid.open_ocean_percent
        CS_l1b_mds['METADATA']['SPH']['CLOSE_SEA_PERCENT'] = fid.close_sea_percent
        CS_l1b_mds['METADATA']['SPH']['CONTINENT_ICE_PERCENT'] = fid.continent_ice_percent
        CS_l1b_mds['METADATA']['SPH']['LAND_PERCENT'] = fid.land_percent
        CS_l1b_mds['METADATA']['SPH']['L1_PROD_STATUS'] = fid.l1b_prod_status
        CS_l1b_mds['METADATA']['SPH']['L1_PROC_FLAG'] = fid.l1b_proc_flag
        CS_l1b_mds['METADATA']['SPH']['L1_PROCESSING_QUALITY'] = fid.l1b_processing_quality
        CS_l1b_mds['METADATA']['SPH']['L1_PROC_THRESH'] = fid.l1b_proc_thresh
        CS_l1b_mds['METADATA']['SPH']['SIR_CONFIGURATION'] = fid.sir_configuration
        CS_l1b_mds['METADATA']['SPH']['SIR_OP_MODE'] = fid.sir_op_mode
        CS_l1b_mds['METADATA']['SPH']['ORBIT_FILE'] = fid.xref_orbit
        CS_l1b_mds['METADATA']['SPH']['PROC_CONFIG_PARAMS_FILE'] = fid.xref_pconf
        CS_l1b_mds['METADATA']['SPH']['CONSTANTS_FILE'] = fid.xref_constants
        CS_l1b_mds['METADATA']['SPH']['IPF_RA_DATABASE_FILE'] = fid.xref_siral_characterisation
        CS_l1b_mds['METADATA']['SPH']['DORIS_USO_DRIFT_FILE'] = fid.xref_uso
        CS_l1b_mds['METADATA']['SPH']['STAR_TRACKER_ATTREF_FILE'] = fid.xref_star_tracker_attref
        CS_l1b_mds['METADATA']['SPH']['SIRAL_LEVEL_0_FILE'] = fid.xref_siral_l0
        CS_l1b_mds['METADATA']['SPH']['CALIBRATION_TYPE_1_FILE'] = fid.xref_cal1
        CS_l1b_mds['METADATA']['SPH']['SIR_COMPLEX_CAL1_SARIN'] = fid.xref_cal1_sarin
        CS_l1b_mds['METADATA']['SPH']['CALIBRATION_TYPE_2_FILE'] = fid.xref_cal2
        CS_l1b_mds['METADATA']['SPH']['SURFACE_PRESSURE_FILE'] = fid.xref_surf_pressure
        CS_l1b_mds['METADATA']['SPH']['MEAN_PRESSURE_FILE'] = fid.xref_mean_pressure
        CS_l1b_mds['METADATA']['SPH']['WET_TROPOSPHERE_FILE'] = fid.xref_wet_trop
        CS_l1b_mds['METADATA']['SPH']['U_WIND_FILE'] = fid.xref_u_wind
        CS_l1b_mds['METADATA']['SPH']['V_WIND_FILE'] = fid.xref_v_wind
        CS_l1b_mds['METADATA']['SPH']['METEO_GRID_DEF_FILE'] = fid.xref_meteo
        CS_l1b_mds['METADATA']['SPH']['S1S2_PRESSURE_00H_MAP'] = fid.xref_s1s2_pressure_00h
        CS_l1b_mds['METADATA']['SPH']['S1S2_PRESSURE_06H_MAP'] = fid.xref_s1s2_pressure_06h
        CS_l1b_mds['METADATA']['SPH']['S1S2_PRESSURE_12H_MAP'] = fid.xref_s1s2_pressure_12h
        CS_l1b_mds['METADATA']['SPH']['S1S2_PRESSURE_18H_MAP'] = fid.xref_s1s2_pressure_18h
        CS_l1b_mds['METADATA']['SPH']['S1_TIDE_AMPLITUDE_MAP'] = fid.xref_s1_tide_amplitude
        CS_l1b_mds['METADATA']['SPH']['S1_TIDE_PHASE_MAP'] = fid.xref_s1_tide_phase
        CS_l1b_mds['METADATA']['SPH']['S2_TIDE_AMPLITUDE_MAP'] = fid.xref_s2_tide_amplitude
        CS_l1b_mds['METADATA']['SPH']['S2_TIDE_PHASE_MAP'] = fid.xref_s2_tide_phase
        CS_l1b_mds['METADATA']['SPH']['GPS_IONO_MAP'] = fid.xref_gim
        CS_l1b_mds['METADATA']['SPH']['IONO_COEFFICENTS_FILE'] = fid.xref_iono_cor
        CS_l1b_mds['METADATA']['SPH']['SAI_FILE'] = fid.xref_sai
        CS_l1b_mds['METADATA']['SPH']['OCEAN_TIDE_FILE'] = fid.xref_ocean_tide
        CS_l1b_mds['METADATA']['SPH']['TIDAL_LOADING_FILE'] = fid.xref_tidal_load
        CS_l1b_mds['METADATA']['SPH']['EARTH_TIDE_FILE'] = fid.xref_earth_tide
        CS_l1b_mds['METADATA']['SPH']['POLE_TIDE_FILE'] = fid.xref_pole_location
        CS_l1b_mds['METADATA']['SPH']['SURFACE_TYPE_FILE'] = fid.xref_surf_type

        # return the output dictionary
        return CS_l1b_mds

    def cryosat_scaling_factors(self):
        """
        Get scaling factors for converting original unpacked units in binary files
        """

        # dictionary of scale factors for CryoSat-2 variables
        CS_l1b_scale = {}

        # CryoSat-2 Time and Orbit Group
        CS_l1b_scale['Location'] = {}
        # Time: day part
        CS_l1b_scale['Location']['Day'] = 1.0
        # Time: second part
        CS_l1b_scale['Location']['Second'] = 1.0
        # Time: microsecond part
        CS_l1b_scale['Location']['Micsec'] = 1.0
        # USO correction factor
        CS_l1b_scale['Location']['USO_Corr'] = 1e-15
        # Mode ID
        CS_l1b_scale['Location']['Mode_ID'] = 1
        # Source sequence counter
        CS_l1b_scale['Location']['SSC'] = 1
        # Instrument configuration
        CS_l1b_scale['Location']['Inst_config'] = 1
        # Record Counter
        CS_l1b_scale['Location']['Rec_Count'] = 1
        # Lat: packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_scale['Location']['Lat'] = 1e-7
        # Lon: packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_scale['Location']['Lon'] = 1e-7
        # Alt: packed units (mm, 1e-3 m)
        # Altitude of COG above reference ellipsoid (interpolated value)
        CS_l1b_scale['Location']['Alt'] = 1e-3
        # Instantaneous altitude rate derived from orbit: packed units (mm/s, 1e-3 m/s)
        CS_l1b_scale['Location']['Alt_rate'] = 1e-3
        # Satellite velocity vector. In ITRF: packed units (mm/s, 1e-3 m/s)
        # ITRF= International Terrestrial Reference Frame
        CS_l1b_scale['Location']['Sat_velocity'] = 1e-3
        # Real beam direction vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
        # CRF= CryoSat Reference Frame.
        CS_l1b_scale['Location']['Real_beam'] = 1e-6
        # Interferometric baseline vector. In CRF: packed units (micro-m/s, 1e-6 m/s)
        CS_l1b_scale['Location']['Baseline'] = 1e-6
        # Star Tracker ID
        CS_l1b_scale['Location']['ST_ID'] = 1
        # Antenna Bench Roll Angle (Derived from star trackers)
        # packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_scale['Location']['Roll'] = 1e-7
        # Antenna Bench Pitch Angle (Derived from star trackers)
        # packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_scale['Location']['Pitch'] = 1e-7
        # Antenna Bench Yaw Angle (Derived from star trackers)
        # packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_scale['Location']['Yaw'] = 1e-7
        # Measurement Confidence Data Flags
        # Generally the MCD flags indicate problems when set
        # If MCD is 0 then no problems or non-nominal conditions were detected
        # Serious errors are indicated by setting bit 31
        CS_l1b_scale['Location']['MCD'] = 1
        CS_l1b_scale['Location']['Spares'] = 1

        # CryoSat-2 Measurement Group
        # Derived from instrument measurement parameters
        CS_l1b_scale['Data'] = {}
        # Window Delay reference (two-way) corrected for instrument delays
        CS_l1b_scale['Data']['TD'] = 1e-12
        # H0 Initial Height Word from telemetry
        CS_l1b_scale['Data']['H_0'] = 4.88e-11
        # COR2 Height Rate: on-board tracker height rate over the radar cycle
        CS_l1b_scale['Data']['COR2'] = 3.05e-12
        # Coarse Range Word (LAI) derived from telemetry
        CS_l1b_scale['Data']['LAI'] = 1.25e-8
        # Fine Range Word (FAI) derived from telemetry
        CS_l1b_scale['Data']['FAI'] = 12.5e-9/256.0
        # Automatic Gain Control Channel 1: AGC gain applied on Rx channel 1.
        # Gain calibration corrections are applied (Sum of AGC stages 1 and 2
        # plus the corresponding corrections) (dB/100)
        CS_l1b_scale['Data']['AGC_CH1'] = 1e-2
        # Automatic Gain Control Channel 2: AGC gain applied on Rx channel 2.
        # Gain calibration corrections are applied (dB/100)
        CS_l1b_scale['Data']['AGC_CH2'] = 1e-2
        # Total Fixed Gain On Channel 1: gain applied by the RF unit. (dB/100)
        CS_l1b_scale['Data']['TR_gain_CH1'] = 1e-2
        # Total Fixed Gain On Channel 2: gain applied by the RF unit. (dB/100)
        CS_l1b_scale['Data']['TR_gain_CH2'] = 1e-2
        # Transmit Power in microWatts
        CS_l1b_scale['Data']['TX_Power'] = 1e-6
        # Doppler range correction: Radial component (mm)
        # computed for the component of satellite velocity in the nadir direction
        CS_l1b_scale['Data']['Doppler_range'] = 1e-3
        # Instrument Range Correction: transmit-receive antenna (mm)
        # Calibration correction to range on channel 1 computed from CAL1.
        CS_l1b_scale['Data']['TR_inst_range'] = 1e-3
        # Instrument Range Correction: receive-only antenna (mm)
        # Calibration correction to range on channel 2 computed from CAL1.
        CS_l1b_scale['Data']['R_inst_range'] = 1e-3
        # Instrument Gain Correction: transmit-receive antenna (dB/100)
        # Calibration correction to gain on channel 1 computed from CAL1
        CS_l1b_scale['Data']['TR_inst_gain'] = 1e-2
        # Instrument Gain Correction: receive-only (dB/100)
        # Calibration correction to gain on channel 2 computed from CAL1
        CS_l1b_scale['Data']['R_inst_gain'] = 1e-2
        # Internal Phase Correction (microradians)
        CS_l1b_scale['Data']['Internal_phase'] = 1e-6
        # External Phase Correction (microradians)
        CS_l1b_scale['Data']['External_phase'] = 1e-6
        # Noise Power measurement (dB/100)
        CS_l1b_scale['Data']['Noise_power'] = 1e-2
        # Phase slope correction (microradians)
        # Computed from the CAL-4 packets during the azimuth impulse response
        # amplitude (SARIN only). Set from the latest available CAL-4 packet.
        CS_l1b_scale['Data']['Phase_slope'] = 1e-6
        CS_l1b_scale['Data']['Spares1'] = 1

        # CryoSat-2 External Corrections Group
        CS_l1b_scale['Geometry'] = {}
        # Dry Tropospheric Correction packed units (mm, 1e-3 m)
        CS_l1b_scale['Geometry']['dryTrop'] = 1e-3
        # Wet Tropospheric Correction packed units (mm, 1e-3 m)
        CS_l1b_scale['Geometry']['wetTrop'] = 1e-3
        # Inverse Barometric Correction packed units (mm, 1e-3 m)
        CS_l1b_scale['Geometry']['InvBar'] = 1e-3
        # Delta Inverse Barometric Correction packed units (mm, 1e-3 m)
        CS_l1b_scale['Geometry']['DAC'] = 1e-3
        # GIM Ionospheric Correction packed units (mm, 1e-3 m)
        CS_l1b_scale['Geometry']['Iono_GIM'] = 1e-3
        # Model Ionospheric Correction packed units (mm, 1e-3 m)
        CS_l1b_scale['Geometry']['Iono_model'] = 1e-3
        # Ocean tide Correction packed units (mm, 1e-3 m)
        CS_l1b_scale['Geometry']['ocTideElv'] = 1e-3
        # Long period equilibrium ocean tide Correction packed units (mm, 1e-3 m)
        CS_l1b_scale['Geometry']['lpeTideElv'] = 1e-3
        # Ocean loading tide Correction packed units (mm, 1e-3 m)
        CS_l1b_scale['Geometry']['olTideElv'] = 1e-3
        # Solid Earth tide Correction packed units (mm, 1e-3 m)
        CS_l1b_scale['Geometry']['seTideElv'] = 1e-3
        # Geocentric Polar tide Correction packed units (mm, 1e-3 m)
        CS_l1b_scale['Geometry']['gpTideElv'] = 1e-3
        # Surface Type: enumerated key to classify surface at nadir
        # 0 = Open Ocean
        # 1 = Closed Sea
        # 2 = Continental Ice
        # 3 = Land
        CS_l1b_scale['Geometry']['Surf_type'] = 1
        CS_l1b_scale['Geometry']['Spare1'] = 1
        # Corrections Status Flag
        CS_l1b_scale['Geometry']['Corr_status'] = 1
        # Correction Error Flag
        CS_l1b_scale['Geometry']['Corr_error'] = 1
        CS_l1b_scale['Geometry']['Spare2'] = 1

        # CryoSat-2 Average Waveforms Groups
        CS_l1b_scale['Waveform_1Hz'] = {}
        # Data Record Time (MDSR Time Stamp)
        CS_l1b_scale['Waveform_1Hz']['Day'] = 1.0
        CS_l1b_scale['Waveform_1Hz']['Second'] = 1.0
        CS_l1b_scale['Waveform_1Hz']['Micsec'] = 1.0
        # Lat: packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_scale['Waveform_1Hz']['Lat'] = 1e-7
        # Lon: packed units (0.1 micro-degree, 1e-7 degrees)
        CS_l1b_scale['Waveform_1Hz']['Lon'] = 1e-7
        # Alt: packed units (mm, 1e-3 m)
        # Altitude of COG above reference ellipsoid (interpolated value)
        CS_l1b_scale['Waveform_1Hz']['Alt'] = 1e-3
        # Window Delay (two-way) corrected for instrument delays
        CS_l1b_scale['Waveform_1Hz']['TD'] = 1e-12
        # 1 Hz Averaged Power Echo Waveform
        CS_l1b_scale['Waveform_1Hz']['Waveform'] = 1.0
        # Echo Scale Factor (to scale echo to watts)
        CS_l1b_scale['Waveform_1Hz']['Linear_Wfm_Multiplier'] = 1.0
        # Echo Scale Power (a power of 2 to scale echo to Watts)
        CS_l1b_scale['Waveform_1Hz']['Power2_Wfm_Multiplier'] = 1.0
        # Number of echoes averaged
        CS_l1b_scale['Waveform_1Hz']['N_avg_echoes'] = 1
        CS_l1b_scale['Waveform_1Hz']['Flags'] = 1

        # CryoSat-2 Waveforms Groups
        # Beam Behavior Parameters
        Beam_Behavior = {}
        # Standard Deviation of Gaussian fit to range integrated stack power.
        Beam_Behavior['SD'] = 1e-2
        # Stack Center: Mean of Gaussian fit to range integrated stack power.
        Beam_Behavior['Center'] = 1e-2
        # Stack amplitude parameter scaled in dB/100.
        Beam_Behavior['Amplitude'] = 1e-2
        # 3rd moment: providing the degree of asymmetry of the range integrated
        # stack power distribution.
        Beam_Behavior['Skewness'] = 1e-2
        # 4th moment: Measure of peakiness of range integrated stack power distribution.
        Beam_Behavior['Kurtosis'] = 1e-2
        # Standard deviation as a function of boresight angle (microradians)
        Beam_Behavior['SD_boresight_angle'] = 1e-6
        # Stack Center angle as a function of boresight angle (microradians)
        Beam_Behavior['Center_boresight_angle'] = 1e-6
        Beam_Behavior['Spare'] = 1

        # CryoSat-2 waveform variables
        CS_l1b_scale['Waveform_20Hz'] = {}
        # Averaged Power Echo Waveform
        CS_l1b_scale['Waveform_20Hz']['Waveform'] = 1.0
        # Echo Scale Factor (to scale echo to watts)
        CS_l1b_scale['Waveform_20Hz']['Linear_Wfm_Multiplier'] = 1.0
        # Echo Scale Power (a power of 2)
        CS_l1b_scale['Waveform_20Hz']['Power2_Wfm_Multiplier'] = 1.0
        # Number of echoes averaged
        CS_l1b_scale['Waveform_20Hz']['N_avg_echoes'] = 1
        CS_l1b_scale['Waveform_20Hz']['Flags'] = 1
        # Beam behaviour parameters
        CS_l1b_scale['Waveform_20Hz']['Beam'] = Beam_Behavior
        # Coherence [SARIN]: packed units (1/1000)
        CS_l1b_scale['Waveform_20Hz']['Coherence'] = 1e-3
        # Phase Difference [SARIN]: packed units (microradians)
        CS_l1b_scale['Waveform_20Hz']['Phase_diff'] = 1e-6

        # return the scaling factors
        return CS_l1b_scale
