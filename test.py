#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 08:32:55 2020

@author: ben
"""

import pointCollection as pc
import numpy as np
xy0=(np.array(-360000.), np.array(-1140000.))
gi=pc.geoIndex().from_file('/Volumes/ice2/ben/CS2_data/GL/retrack_BC/2019/index_POCA.h5')
D=gi.query_xy(xy0, fields=['x','y', 'h'])
gi=None