#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 09:48:40 2020

@author: ben
"""

import numpy as np

def campaign_bias_correction(campaign):
    '''
    calculate additive bias correction based on ICESat campaign number.

    the resulting correction is ADDED to the heights to correct for the bias.
    This is based on Helen's description of Adrian's correction:
        subtract 1.7 cm from Laser 2 and add 1.1 cm to laser 3

    Parameters
    ----------
    campaign : numpy array
        campaign number

    Returns
    -------
    correction value

    '''
    correction=np.zeros_like(campaign, dtype=float)+np.nan
    laser2=np.in1d(campaign.astype(int), 1+np.concatenate([np.arange(1,4), np.arange(15,18)]) )
    correction[laser2] -= 0.017
    laser3=np.in1d(campaign.astype(int), 1+np.arange(4, 15) )
    correction[laser3] += 0.011
