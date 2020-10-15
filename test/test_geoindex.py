#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 08:32:55 2020

@author: ben
"""

import os
import warnings
import pytest
import pointCollection as pc
import numpy as np

def test_geoindex():
    xy0=(np.array(-360000.), np.array(-1140000.))
    geoIndex_file = '/Volumes/ice2/ben/CS2_data/GL/retrack_BC/2019/index_POCA.h5'
    # only run test if file exists
    try:
        gi=pc.geoIndex().from_file(geoIndex_file)
    except OSError:
        print('{0} not found on file system'.format(geoIndex_file))
    else:
        D=gi.query_xy(xy0, fields=['x','y', 'h'])
        gi=None
