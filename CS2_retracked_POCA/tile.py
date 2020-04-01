#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 18:45:18 2019

@author: ben
"""

import pointCollection as pc


class tile(pc.tile):

    def __default_field_dict__(self):
        return {None: ["DEM", "abs_orbit", "ambiguity", "block_h_spread", "burst", "coherence", "count", "dh_spline_dx", "geoid", "h", "phase", "power", "pt_num", "range_surf", "ret_count", "samp", "time", "x", "y"]}
    
    def __time_field__(self):
        return 'time'

    def __z_field__(self):
        return 'h'
