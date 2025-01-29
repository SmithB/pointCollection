# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 10:46:21 2017

Class to read ATM data in the icessn format

@author: ben
"""
import numpy as np
import pointCollection as pc

class data(pc.data):
    np.seterr(invalid='ignore')

    from .from_text import from_text, get_ATM_time

    def __default_field_dict__(self):
        return {None:['x','y','latitude','longitude','z', 'time', 'track_identifier','sigma']}
