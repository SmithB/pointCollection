#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 22:20:41 2019

@author: ben
"""

import datetime

def t_0():
    """
    return a datetime object corresponding to 2018 01 01 00:00:00
    """
    return datetime.datetime(2018, 1, 1, 0, 0, 0)

def to_delta_time(t_in_0):
    """
    convert a datetime or a time tuple to an icesat2 delta_time value
    """
    if isinstance(t_in_0, (tuple, list)):
        t_in = datetime.datetime(*tuple(t_in_0))
    else:
        t_in = t_in_0
    return (t_in-t_0()).total_seconds()

def from_delta_time(delta_time):
    """
    convert a delta_time to a datetime object
    """
    return t_0()+datetime.timedelta(seconds=delta_time)