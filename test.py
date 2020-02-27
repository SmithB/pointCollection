#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 16:25:38 2019

@author: ben
"""

import pointCollection as pc
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
ATL06_file='../pointCollection/test_data/ATL06_20190205041106_05910210_209_01.h5'
Dsub=pc.data().from_h5(ATL06_file, field_dict={'/gt1l/land_ice_segments':['h_li', 'segment_id']})
print(Dsub.h_li[0:100])
