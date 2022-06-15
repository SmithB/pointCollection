#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 12:14:24 2022

@author: ben
"""

import matplotlib.pyplot as plt
import numpy as np
import pointCollection as pc

def in_axes(D, ax=None):
    '''
    identify the elements of D that are in a set of axes

    Parameters
    ----------
    D : pointCollection.data, iterable, or pointCollection.grid.data
        Object with coordinates.  If D is a pointCollection.data or an iterable
        containing two sets of coordinates, a single boolean will be returned 
        specifying the in-axes points.  If D is a pointCollection.grid.data, 
        one set of booleans will be returned for the y axis, one for the x axis
    ax : TYPE, optional
        The axes that will be queried.  If None, then gca() will be used.
        The default is None.

    Raises
    ------
    TypeError
        Error if the input data type is not understood.

    Returns
    -------
    boolean numpy array
        One (for point data) or two (for image data) boolean arrays indicating
        the in-axes elements.  For images, the y axis is listed first.

    '''
    if ax is None:
        ax=plt.gca()
    XR=ax.get_xlim()
    YR=ax.get_ylim()
    if isinstance(D, pc.grid.data):
        x=D.x
        y=D.y
        return  (y >= YR[0]) & (y< YR[1]), (x>=XR[0]) & (x <= XR[1])
    elif isinstance(D, pc.data):
        x=D.x
        y=D.y
    elif isinstance(D, (list, tuple)):
        x=D[0]
        y=D[1]
    else:
        raise TypeError('type of argument "D" unknown')
    return (x>=XR[0]) & (x <= XR[1]) & (y >= YR[0]) & (y< YR[1])