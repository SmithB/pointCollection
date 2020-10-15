# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:51:44 2019

@author: ben
"""
import numpy as np
import re
from datetime import date

def WV_date(filename):
    date_re=re.compile(r'\d\d.*_(2\d\d\d)(\d\d)(\d\d)_')
    m=date_re.search(filename)
    if m is None:
        return np.NaN
    return date(int(m.group(1)), int(m.group(2)), int(m.group(3)))

def WV_year(filename):
    this_date=WV_date(filename)
    this_delta=this_date-date(2000, 1, 1)
    return 2000+this_delta.days/365.25

def WV_MatlabDate(filename):
    # matlab date gives days relative to date 0 0 0.  The earliest date that
    # a datetime can handle is 1, 1, 1, and the difference between the two is
    # 367 days.
    this_date=WV_date(filename)
    this_delta=this_date-date(1, 1, 1)
    this_delta=this_delta.days+this_delta.seconds/24./3600.+367.
    return this_delta
