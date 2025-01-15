import numpy as np
#import matplotlib.pyplot as plt
import pyproj
import pointCollection as pc

"""
History
2/2024 : Written by Tyler Sutterley
3/2024 : Modified to include 2005 and 2008 epochs by Ben Smith
"""


# PURPOSE: convert heights between reference frames
def convert_ITRF(lon, lat, data, tdec, epochs=[2014, 2020], direction='inverse'):
    """Converts height data between ITRF2014 and ITRF2020

    Parameters
    ----------
    lon: np.ndarray
        Longitude (degrees)
    lat: np.ndarray
        Latitude (degrees)
    data: np.ndarray
        Height (m)
    tdec: np.ndarray
        Time (decimal years)
    direction: str
        Direction of transform

            - ``'forward'``: early to late
            - ``'inverse'``: late to early
    """
    # get the transform for converting from ITRF2014 to ITRF2020
    if tuple(epochs)==(2014, 2020):
        transform = wgs84_itrf2014_to_wgs84_itrf2020()
    elif tuple(epochs)==(2005, 2014):
        transform = wgs84_itrf2005_to_wgs84_itrf2014()
    elif tuple(epochs)==(2008, 2014):
        transform = wgs84_itrf2008_to_wgs84_itrf2014()
    # transform the data to ITRF
    return transform.transform(lon, lat, data, tdec, direction=direction)

# WGS84 Ellipsoid in ITRF2014 to WGS84 Ellipsoid in ITRF2020
def wgs84_itrf2014_to_wgs84_itrf2020():
    """``pyproj`` transform for WGS84 Ellipsoid in ITRF2014 to WGS84 Ellipsoid in ITRF2020"""
    pipeline = """+proj=pipeline
        +step +proj=unitconvert +xy_in=deg +z_in=m +xy_out=rad +z_out=m
        +step +proj=cart +ellps=WGS84
        +step +proj=helmert +x=0.0014 +y=0.0009 +z=-0.0014 +rx=0 +ry=0 +rz=0 +s=0.00042
            +dx=0 +dy=0.0001 +dz=-0.0002 +drx=0 +dry=0 +drz=0 +ds=0
            +t_epoch=2015 +convention=position_vector
        +step +inv +proj=cart +ellps=WGS84
        +step +proj=unitconvert +xy_in=rad +z_in=m +xy_out=deg +z_out=m"""
    return pyproj.Transformer.from_pipeline(pipeline)


# WGS84 Ellipsoid in ITRF2008 to WGS84 Ellipsoid in ITRF2014
def wgs84_itrf2008_to_wgs84_itrf2014():
    """``pyproj`` transform for WGS84 Ellipsoid in ITRF2008 to WGS84 Ellipsoid in ITRF2014"""
    pipeline = """+proj=pipeline
        +step +proj=unitconvert +xy_in=deg +z_in=m +xy_out=rad +z_out=m
        +step +proj=cart +ellps=WGS84
        +step +proj=helmert +x=-0.0016 +y=-0.0019 +z=-0.0024 +rx=0 +ry=0 +rz=0 +s=0.00002
            +dx=0 +dy=0.0000 +dz=0.0001 +drx=0 +dry=0 +drz=0 +ds=-0.00003
            +t_epoch=2015 +convention=position_vector
        +step +inv +proj=cart +ellps=WGS84
        +step +proj=unitconvert +xy_in=rad +z_in=m +xy_out=deg +z_out=m"""
    return pyproj.Transformer.from_pipeline(pipeline)


# WGS84 Ellipsoid in ITRF2005 to WGS84 Ellipsoid in ITRF2014
def wgs84_itrf2005_to_wgs84_itrf2014():
    """``pyproj`` transform for WGS84 Ellipsoid in ITRF2005 to WGS84 Ellipsoid in ITRF2014"""
    pipeline = """+proj=pipeline
        +step +proj=unitconvert +xy_in=deg +z_in=m +xy_out=rad +z_out=m
        +step +proj=cart +ellps=WGS84
        +step +proj=helmert +x=-0.0026 +y=-0.001 +z=0.0023 +rx=0 +ry=0 +rz=0 +s=-0.00092
            +dx=-0.0003 +dy=0.0 +dz=0.0001 +drx=0 +dry=0 +drz=0 +ds=-0.00003
            +t_epoch=2010 +convention=position_vector
        +step +inv +proj=cart +ellps=WGS84
        +step +proj=unitconvert +xy_in=rad +z_in=m +xy_out=deg +z_out=m"""
    return pyproj.Transformer.from_pipeline(pipeline)
