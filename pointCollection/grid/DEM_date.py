# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:51:44 2019

@author: ben
"""
import numpy as np
import re
from datetime import datetime
import os
import glob
import json

def DEM_date(filename):
    """Parse a Worldview or tanDEM-x filename for a date."""

    temp=datetime_from_meta(filename)
    if temp is not None:
        return temp
    date_re=re.compile(r'\d.*_(2\d\d\d)(\d\d)(\d\d)_')
    m=date_re.search(filename)
    if m is not None:
        return datetime(int(m.group(1)), int(m.group(2)), int(m.group(3)), 0, 0, 0)
    date_re = re.compile(r'_(\d{4})(\d{2})(\d{2})\D(\d{2})(\d{2})(\d{2})_')
    m=date_re.search(filename)
    if m is not None:
        return datetime(*map(int, m.groups()))
    return np.nan

def DEM_year(filename):
    """Return the decimal year date for a Worldview or tamDEM-x filename."""
    this_date=DEM_date(filename)
    this_delta=this_date-datetime(2000, 1, 1, 0, 0, 0)
    return 2000 + 1/365.25*(this_delta.days + this_delta.seconds/24./3600.)

def DEM_MatlabDate(filename):
    """Return the Matlab-style date for a WorldView or tanDEM-x filename."""
    # matlab date gives days relative to date 0 0 0.  The earliest date that
    # a datetime can handle is 1, 1, 1, and the difference between the two is
    # 367 days.
    this_date=DEM_date(filename)
    this_delta=this_date-datetime(1, 1, 1, 0, 0, 0)
    this_delta=this_delta.days+this_delta.seconds/24./3600.+367.
    return this_delta

def datetime_from_meta(filename):

    """Look for a metadata file in PGC format"""
    thedir, thefile=os.path.split(filename)
    theparent = os.path.dirname(thedir)
    try:
        res=re.compile(r'_(\d+m)_').search(thefile).group(1)
        glob_str = thefile.replace(res,'*')
        glob_str=glob_str.replace('_dem_filt.tif','.json')
        if not glob_str.endswith('json'):
            return None
        temp=glob.glob(os.path.join(thedir, glob_str))
        if len(temp)==0:
            temp=glob.glob(os.path.join(thedir,'meta', glob_str))
            if len(temp)==0:
                temp=glob.glob(os.path.join(theparent, 'meta', glob_str))
                if len(temp)==0:
                    return None
        with open(temp[0]) as fh:
            thestr=json.load(fh)['properties']['datetime']
        return datetime.strptime(thestr,'%Y-%m-%dT%H:%M:%SZ')
    except Exception as e:
        print(f"DEM_date for {filename}: exception thrown:")
        print(e)
        return None
