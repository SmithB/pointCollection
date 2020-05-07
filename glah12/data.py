#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:44:26 2020

@author: ben
"""


import pointCollection as pc
import h5py
import numpy as np
import warnings

class data(pc.data):
    def __default_field_dict__(self):
        fd={
            'Elevation_Corrections':['d_GmC', 'd_satElevCorr'],
            'Elevation_Surfaces':['d_elev','d_IceSVar'],
            'Geolocation':['d_lat', 'd_lon'],
            'Geophysical':['d_DEM_elv', 'd_deltaEllip', 'd_ocElv'],
            'Quality':['elev_use_flg', 'sigma_att_flg'],
            'Reflectivity':['d_reflctUC'],
            'Time':['d_UTCTime_40'],
            'Waveform':['i_numPk']
            }
        temp= {'Data_40HZ/'+key:fd[key] for key in fd}
        temp.update({'derived':['campaign']})
        return temp
    def __internal_field_calc__(self, field_dict):
        if 'derived' in field_dict and 'campaign' in field_dict['derived']:
            with h5py.File(self.filename,'r') as h5f:
                campaign=h5f['ANCILLARY_DATA'].attrs['Campaign'].decode('ascii')
            campaign_list=['1A','2A','2B','2C','3A','3B','3C','3D','3E','3F',
                           '3G','3H','3I','3J','3K', '2D','2E','2F']
            self.assign({'campaign': np.ones_like(self.d_elev)+campaign_list.index(campaign)})

    def __apply_corrections__(self):
        #delta_ellip should be TP - WGS84, but it seems to be the opposite.  It needs to be subtracted
        z=self.elev-self.deltaEllip
        # remove the tide correction
        z[np.isfinite(self.ocElv)] -= self.ocElv[np.isfinite(self.ocElv)]
        # add the saturation elevation correction
        z[np.isfinite(self.satElevCorr)] += self.satElevCorr[np.isfinite(self.satElevCorr)]
        self.assign({'z':z})

    def __apply_filters__(self):
        good=(self.IceSVar < 0.035) & \
            (self.reflctUC > 0.05) & \
            (self.satElevCorr < 1) & \
            (self.i_numPk==1)
        self.index(good)

    def from_h5(self, filename, field_dict=None, lat_range=None, apply_corrections=True, apply_filters=True):

        if field_dict is None:
            field_dict=self.__default_field_dict__()
        # call pc.data() to read the file,
        D0=super().from_h5(filename, field_dict=field_dict)
        temp=D0.__dict__
        out={}
        for field in D0.fields:
            if field=='d_lat':
                out['latitude']=temp[field]
            elif field=='d_lon':
                out['longitude']=temp[field]
            elif field[0:2]=='d_':
                out[field[2:]]=temp[field]
            else:
                out[field]=temp[field]
        self.from_dict(out)
        if lat_range is not None:
            self.index((self.latitude>lat_range[0]) & (self.latitude < lat_range[1]))
        self.__internal_field_calc__(self.field_dict)
        if apply_corrections:
            self.__apply_corrections__()
        if apply_filters:
            self.__apply_filters__()
        self.assign({'time':self.UTCTime_40})
        return self