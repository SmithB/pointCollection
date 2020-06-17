# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 10:46:21 2017

Class to read and manipulate ATL06 data.  Currently set up for Ben-style fake data, should be modified to work with the official ATL06 prodct foramt

@author: ben
"""
import h5py
import numpy as np
from datetime import datetime, timedelta
from osgeo import osr
import matplotlib.pyplot as plt
from pointCollection import data
import re

class qFit(data):
    np.seterr(invalid='ignore')


    def __default_field_dict__(self):
        return {None:['latitude','longitude','elevation'], \
                'instrument_parameters':['azimuth','rel_time'],\
                    '__calc_internal__':['days_J2k']}

    def __internal_field_calc__(self, field_dict):
        # find the date and time number in filename
        if 'days_J2k' in field_dict['__calc_internal__']:
            m=re.search(r"ATM1B.*_(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2}).*.h5",self.filename)
            # create list of date variables
            this_time=[int(m.group(ind+1)) for ind in range(6)]

            t0=datetime(*this_time[0:3]) + timedelta(hours=this_time[3], minutes=this_time[4], seconds=this_time[5])-datetime(2000, 1, 1, 0, 0, 0)
            t0=t0.days+t0.seconds/24./3600
            self.days_J2K=(t0-datetime(2000,1, 1, 0, 0, 0)).days + self.seconds_of_day/24./3600.
