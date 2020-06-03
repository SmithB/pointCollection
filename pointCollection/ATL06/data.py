#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 10:04:52 2019

@author: ben
"""

"""
ATL06 data structure
"""

import pointCollection as pc
import h5py
import numpy as np
import warnings

class data(pc.data):
    def __init__(self, fields=None, SRS_proj4=None, field_dict=None, beam_pair=2, columns=2):
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
        self.shape=(0,2)
        self.size=0
        self.beam_pair=beam_pair
        self.beam_type=[None, None]
        self.filename=None

    def __default_field_dict__(self):
        """
        Define the default fields that get read from the h5 file
        """
        field_dict={None:['delta_time','h_li','h_li_sigma','latitude','longitude','atl06_quality_summary','segment_id','sigma_geo_h'],
                    'ground_track':['x_atc', 'y_atc','seg_azimuth','sigma_geo_at','sigma_geo_xt'],
                    'fit_statistics':['dh_fit_dx','dh_fit_dx_sigma','h_mean', 'dh_fit_dy','h_rms_misfit','h_robust_sprd','n_fit_photons', 'signal_selection_source','snr_significance','w_surface_window_final'],
                    'geophysical':['bsnow_conf','bsnow_h','cloud_flg_asr','cloud_flg_atm','r_eff','tide_ocean'],
                    'derived':['valid', 'BP','LR', 'cycle_number', 'rgt']}
        return field_dict

    def from_h5(self, filename, group=None, field_dict=None, index_range=None):
        """
        Read data from a file.
        """
  
        if field_dict is not None:
            fields=[]
            for group in self.field_dict.keys():
                for field in self.field_dict[group]:
                        fields.append(field)
            self.fields=fields
            self.field_dict=field_dict
        self.file=filename
        h5_f=h5py.File(filename,'r')
        # generate the name for the hdf5 beam groups
        beam_names=['gt%d%s' %(self.beam_pair, b) for b in ['l','r']]
        if beam_names[0] not in h5_f or beam_names[1] not in h5_f:
            # return empty data structure
            for group in self.field_dict:
                for field in self.field_dict[group]:
                    setattr(self, field, np.zeros((0,2)))
                self.__update_size_and_shape__()
            return self
        if beam_names[0] not in h5_f.keys():
            return None
        # get the strong/weak beam info
        for count, beam_name in enumerate(beam_names):
            try:
                self.beam_type[count]=h5_f[beam_name].attrs['atlas_beam_type']
            except KeyError:
                pass  # leave the beam type as the default
        if index_range is None or index_range[1]==-1:
            index_range=[0, np.minimum(h5_f[beam_names[0]]['land_ice_segments']['h_li'].size, h5_f[beam_names[1]]['land_ice_segments']['h_li'].size)]
        n_vals=index_range[-1]-index_range[0]
        # read the orbit number
        try:
            self.rgt=int(h5_f['orbit_info']['rgt'][0])+np.zeros((n_vals,2))
            self.orbit=int(h5_f['orbit_info']['orbit_number'][0])
        except:
            pass

        for group in self.field_dict.keys():
            for field in self.field_dict[group]:
                if field not in self.fields:
                    self.fields.append(field)
                try:
                    #bad_val=-9999 # default value, works for most groups
                    if group is None:
                        #The None group corresponds to the top level of the land_ice_segments hirearchy
                        # most datasets have a __fillValue attribute marking good data, but some don't
                        try:
                            bad_val=h5_f[beam_names[0]]['land_ice_segments'][field].attrs['_FillValue']
                        except KeyError:
                            #latitude and longitude don't have fill values, but -9999 works OK.
                            bad_val=-9999
                        data=np.c_[
                            np.array(h5_f[beam_names[0]]['land_ice_segments'][field][index_range[0]:index_range[1]]).transpose(),
                            np.array(h5_f[beam_names[1]]['land_ice_segments'][field][index_range[0]:index_range[1]]).transpose()]
                    elif group == "orbit_info":
                        data=np.zeros((n_vals, 2))+h5_f['orbit_info'][field]
                        bad_val=-9999
                    elif group == "derived" and field == "valid":
                        data=np.ones((index_range[1]-index_range[0], 2), dtype='bool')
                        bad_val=0
                    elif group == "derived" and field == "BP":
                        data=np.zeros((index_range[1]-index_range[0], 2))+self.beam_pair
                        bad_val=0
                    elif group == "derived" and field == "LR":
                        data=np.ones((index_range[1]-index_range[0], 2))
                        data[:,0]=0
                        bad_val=-9999
                    elif group == "derived" and field == "n_pixels":
                        if self.beam_type[0]=='weak':
                            data=np.tile([4, 16], [n_vals, 1])
                        else:
                            data=np.tile([16, 4], [n_vals, 1])
                        bad_val=-9999
                    elif group == "derived" and field == "spot":
                        data=np.ones((index_range[1]-index_range[0], 2))
                        for bb in [0, 1]:
                            data[:,bb]=np.float64(h5_f[beam_names[bb]].attrs['atlas_spot_number'])
                    elif group == "derived" and field == "sigma_geo_r":
                        # Earlier versions of the data don't include a sigma_geo_r.  This is a fake value
                        data=np.zeros((index_range[1]-index_range[0], 2)) +0.03
                    elif group == "derived":
                        continue
                    else:
                        # All other groups are under the land_ice_segments/group hirearchy
                         try:
                            bad_val=h5_f[beam_names[0]]['land_ice_segments'][group][field].attrs['_FillValue']
                         except KeyError:
                            #flags don't have fill values, but -9999 works OK.
                            bad_val=-9999
                         try:
                             data=np.c_[
                                np.array(h5_f[beam_names[0]]['land_ice_segments'][group][field][index_range[0]:index_range[1]]).transpose(),
                                np.array(h5_f[beam_names[1]]['land_ice_segments'][group][field][index_range[0]:index_range[1]]).transpose()]
                         except KeyError:
                             print("missing hdf field for land_ice_segments/%s/%s in %s" % (group, field, filename))
                             data=np.zeros((index_range[1]+1-index_range[0], 2))+bad_val
                    # identify the bad data elements before converting the field to double
                    bad=data==bad_val
                    if field != 'valid':
                        data=data.astype(np.float64)
                    # mark data that are invalid (according to the h5 file) with NaNs
                    data[bad]=np.NaN
                    setattr(self, field, data)
                except KeyError:
                    print("could not read %s/%s" % (group, field))
                    setattr(self, field, np.zeros( [n_vals, 2])+np.NaN)
        self.__update_size_and_shape__()

        if "derived" in self.field_dict and "matlab_time" in self.field_dict['derived']:
            self.get_Matlab_time()
        self.__update_size_and_shape__()

        if "derived" in self.field_dict and "rss_along_track_dh" in self.field_dict['derived']:
            self.get_rss_along_track_dh()
        if "derived" in self.field_dict and "min_along_track_dh" in self.field_dict['derived']:            
            setattr(self, 'min_along_track_dh', self.calc_min_along_track_dh())

        # assign fields that must be copied from single-value attributes in the
        # h5 file
        if 'cycle_number' in self.fields:
            # get the cycle number
            try:
                cycle_number=h5_f['ancillary_data']['start_cycle']
            except KeyError:
                cycle_number=-9999
            setattr(self, 'cycle_number', cycle_number+np.zeros(self.shape, dtype=np.float64))
        h5_f.close()

        return self

    def get_Matlab_time(self):
        self.assign({'matlab_time': 737061. + self.delta_time/24./3600.})

    def calc_min_along_track_dh(self):
        min_along_track_dh=np.zeros(list(self.shape)+[2])
        n_pts=self.shape[0]
        if n_pts > 1:
            i0=slice(1, n_pts-1)
            for dim3, ii in enumerate([-1, 1]):
                i1=slice(1+ii, n_pts-1+ii)
                dx=self.x_atc[i0,:]-self.x_atc[i1,:]
                min_along_track_dh[i0,:, dim3] = np.abs(self.h_li[i0,:]-self.dh_fit_dx[i0,:]*dx-self.h_li[i1,:])
            # Note that np.nanmin will throw a RuntimeWarning in the (common)
            # case where both ends of the segment have an undefined difference.
            # Catch this error to avoid distractions:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                min_along_track_dh = np.nanmin(min_along_track_dh, axis=2)

            min_along_track_dh[0,:]=np.abs(self.h_li[1,:]-self.h_li[0,:] - (self.x_atc[1,:]-self.x_atc[0,:])*self.dh_fit_dx[0,:])
            min_along_track_dh[-1,:]=np.abs(self.h_li[-1,:]-self.h_li[-2,:] - (self.x_atc[-1,:]-self.x_atc[-2,:])*self.dh_fit_dx[-1,:])
        else:
            min_along_track_dh=np.zeros([1,2])+np.NaN
        return min_along_track_dh
        
    def get_rss_along_track_dh(self):
        self.rss_along_track_dh=np.zeros(self.shape)
        n_pts=self.shape[0]
        if n_pts > 1:
            i0=slice(1, n_pts-1)
            for ii in [-1, 1]:
                i1=slice(1+ii, n_pts-1+ii)
                dx=self.x_atc[i0,:]-self.x_atc[i1,:]
                self.rss_along_track_dh[i0,:] += (self.h_li[i0,:]-self.dh_fit_dx[i0,:]*dx-self.h_li[i1,:])**2
            self.rss_along_track_dh[0,:]=(self.h_li[1,:]-self.h_li[0,:] - (self.x_atc[1,:]-self.x_atc[0,:])*self.dh_fit_dx[0,:])**2
            self.rss_along_track_dh[-1,:]=(self.h_li[-1,:]-self.h_li[-2,:] - (self.x_atc[-1,:]-self.x_atc[-2,:])*self.dh_fit_dx[-1,:])**2
            self.rss_along_track_dh = np.sqrt(self.rss_along_track_dh)

def delta_t_to_Matlab(delta_t):
    return 730486 + delta_t/24./3600.
