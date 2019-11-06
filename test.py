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

test_dir='/Users/ben/git_repos/pointCollection/test_data/for_geoindex'

test_files=glob.glob(test_dir+'/*.h5')
gi_2=pc.geoIndex(delta=[10,10]).for_files(test_files,'h5')
