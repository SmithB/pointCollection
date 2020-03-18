#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  2 10:04:16 2020

pointCollection class for reading ATM Waveform geolocation information


@author: ben
"""
import numpy as np
import pointCollection as pc
from datetime import datetime, timedelta
import re

class data(pc.data):
    np.seterr(invalid='ignore')

    def __default_field_dict__(self):
        return {'footprint':['latitude','longitude','elevation'],
                'time':['seconds_of_day'],
                '__calc_internal__':['days_J2k']}

    def __internal_field_calc__(self, field_dict):
        # find the date and time number in filename
        if 'days_J2k' in field_dict['__calc_internal__']:
            m=re.search(r"\D*_(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2}).*.h5",self.filename)
			# create list of date variables
            this_time=[int(m.group(ind+1)) for ind in range(6)]

            t0=datetime(*this_time[0:3]) + timedelta(hours=this_time[3], minutes=this_time[4], seconds=this_time[5])-datetime(2000, 1, 1, 0, 0, 0)
            t0=np.float64(t0.days)

            self.days_J2k = t0 + self.seconds_of_day.astype(np.float64)/24./3600.
