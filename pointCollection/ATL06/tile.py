#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 18:45:18 2019

@author: ben
"""

import pointCollection as pc


class tile(pc.tile):

    def __default_field_dict__(self):
        return {None:['delta_time','h_li','h_li_sigma','latitude','longitude','atl06_quality_summary','segment_id','sigma_geo_h'],
            'fit_statistics':['dh_fit_dx'],
            'ground_track':['x_atc', 'sigma_geo_xt','sigma_geo_at'],
            'geophysical' : ['dac','tide_ocean'],
            'orbit_info':['rgt','cycle_number'],
            'derived':['valid','matlab_time', 'n_pixels','LR','BP','spot','rss_along_track_dh']}

    def __time_field__(self):
        return 'delta_time'

    def __z_field__(self):
        return 'h_li'
