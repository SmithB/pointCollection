#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 14:10:10 2020

@author: ben
"""

import pointCollection as pc

SRS_proj4='+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
gI=pc.geoIndex(SRS_proj4=SRS_proj4).for_file('/Volumes/ice2/ben/scf/AA_11/U07/ATL11_000110_0306_02_vU07.h5', file_type='ATL11')
