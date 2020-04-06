#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 18:45:18 2019

@author: ben
"""

import pointCollection as pc


class tile(pc.tile):

    def __default_field_dict__(self):
        return {None: ["AD", "DEM", "POCA_flag", "R_POCA", "R_offnadir", "abs_orbit", "ambiguity", "block_h_spread", "burst", "coherence", "count", "dRange_POCA", "dh_spline_dx", "error_composite", "geoid", "h", "phase", "power", "pt_num", "range_surf", "samp", "seg_ind", "time", "x", "y"]}
    
    def __time_field__(self):
        return 'time'

    def __z_field__(self):
        return 'h'
