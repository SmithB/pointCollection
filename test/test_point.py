#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 08:32:55 2020

@author: ben
"""
import os
import inspect
import warnings
import pytest
import numpy as np
import pointCollection as pc

def test_point():
    x=np.arange(10)
    y=x**2
    z=y
    D1 = pc.data().from_dict({'x':x,'y':y,'z':z})
    D2 = pc.data().from_list([D1])
    assert all((getattr(D1,key) == getattr(D2,key)).all() for key in D1.fields)

def test_read_ATL06():
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    filepath = os.path.dirname(os.path.abspath(filename))
    ATL06_file = os.path.join(filepath,'..','test_data',
        'ATL06_20190205041106_05910210_209_01.h5')
    # general list of ATL06 fields with projected coordinate fields
    ATL06_fields = ['delta_time','h_li','h_li_sigma','latitude','longitude',
        'atl06_quality_summary','segment_id','sigma_geo_h','x_atc','y_atc',
        'seg_azimuth','sigma_geo_at','sigma_geo_xt','dh_fit_dx','dh_fit_dx_sigma',
        'h_mean','dh_fit_dy','h_rms_misfit','h_robust_sprd','n_fit_photons',
        'signal_selection_source','snr_significance','w_surface_window_final',
        'bsnow_conf','bsnow_h','cloud_flg_asr','cloud_flg_atm','r_eff',
        'tide_ocean','valid','BP','LR','cycle_number','rgt','x','y']
    # read ATL06 file and get projected coordinates
    D_06 = pc.ATL06.data().from_h5(os.path.relpath(ATL06_file))
    D_06.get_xy(EPSG=3031)
    assert (D_06.fields == ATL06_fields)
