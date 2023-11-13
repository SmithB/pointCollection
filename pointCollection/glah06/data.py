#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:44:26 2020

@author: ben
"""


import pointCollection as pc
import h5py
import numpy as np
import warnings

class data(pc.data):
    def __default_field_dict__(self):
        fd={
            'Elevation_Corrections':['d_GmC', 'd_satElevCorr'],
            'Elevation_Surfaces':['d_elev'],
            'Geolocation':['d_lat', 'd_lon'],
            'Geophysical':['d_DEM_elv', 'd_deltaEllip', 'd_ocElv'],
            'Quality':['elev_use_flg', 'sigma_att_flg'],
            'Reflectivity':['d_reflctUC', 'd_sDevNsOb1'],
            'Time':['d_UTCTime_40'],
            'Waveform':['i_numPk']
            }
        temp= {'Data_40HZ/'+key:fd[key] for key in fd}
        temp.update({'derived':['campaign']})
        return temp
    def __internal_field_calc__(self, field_dict):
        if 'derived' in field_dict and 'campaign' in field_dict['derived']:
            with h5py.File(self.filename,'r') as h5f:
                campaign=h5f['ANCILLARY_DATA'].attrs['Campaign'].decode('ascii')
            campaign_list=['1A','2A','2B','2C','3A','3B','3C','3D','3E','3F',
                           '3G','3H','3I','3J','3K', '2D','2E','2F']
            self.assign({'campaign': np.ones_like(self.d_elev)+campaign_list.index(campaign)})

    def __apply_corrections__(self, convert=False):
        # calculate the heights and latitudes for different ellipsoids
        if convert:
            # semimajor axis (a) and flattening (f) for TP and WGS84 ellipsoids
            atop,ftop=(6378136.3,1.0/298.257)
            awgs,fwgs=(6378137.0,1.0/298.257223563)
            # calculate latitudes and heights for WGS84
            phi,z=self.convert_ellipsoid(atop, ftop, awgs, fwgs)
        else:
            # delta_ellip is TP (elev) - WGS84 (elev)
            # (not difference between ellipsoids)
            # delta_ellip values need to be subtracted to convert to WGS84
            z=self.elev-self.deltaEllip
            phi=np.copy(self.latitude)
        # remove the tide correction
        z[np.isfinite(self.ocElv)] -= self.ocElv[np.isfinite(self.ocElv)]
        # add the saturation elevation correction
        z[np.isfinite(self.satElevCorr)] += self.satElevCorr[np.isfinite(self.satElevCorr)]
        self.assign({'z':z,'latitude':phi})

    def __apply_filters__(self):
        good=(self.d_sDevNsOb1 < 0.035) & \
            (self.reflctUC > 0.05) & \
            (self.satElevCorr < 1) & \
            (self.i_numPk==1)
        self.index(good)

    def from_h5(self, filename, field_dict=None, lat_range=None, apply_corrections=True, apply_filters=True):

        if field_dict is None:
            field_dict=self.__default_field_dict__()
        # call pc.data() to read the file,
        D0=super().from_h5(filename, field_dict=field_dict)
        temp=D0.__dict__
        out={}
        for field in D0.fields:
            if field=='d_lat':
                out['latitude']=temp[field]
            elif field=='d_lon':
                out['longitude']=temp[field]
            elif field[0:2]=='d_':
                out[field[2:]]=temp[field]
            else:
                out[field]=temp[field]
        self.from_dict(out)
        if lat_range is not None:
            self.index((self.latitude>lat_range[0]) & (self.latitude < lat_range[1]))
        self.__internal_field_calc__(self.field_dict)
        if apply_corrections:
            self.__apply_corrections__()
        if apply_filters:
            self.__apply_filters__()
        self.assign({'time':self.UTCTime_40})
        return self

    def convert_ellipsoid(self, a1, f1, a2, f2, EPS=1e-12, ITMAX=10):
        """
        Convert latitudes and heights to a different ellipsoid using newton-raphson

        Inputs:
            a1: semi-major axis of input ellipsoid
            f1: flattening of input ellipsoid
            a2: semi-major axis of output ellipsoid
            f2: flattening of output ellipsoid

        Options:
            EPS: tolerance to prevent division by small numbers
                and to determine convergence
            ITMAX: maximum number of iterations to use in newton-raphson
        """
        # semiminor axis of input and output ellipsoid
        b1 = (1.0 - f1)*a1
        b2 = (1.0 - f2)*a2
        # initialize output arrays
        npts = len(self.d_lat)
        phi2 = np.zeros((npts))
        h2 = np.zeros((npts))
        # for each point
        for N in range(npts):
            # copy variables to prevent being modified
            phi1 = np.copy(self.latitude[N])
            h1 = np.copy(self.elev[N])
            # force phi1 into range -90 <= phi1 <= 90
            phi1 = -90.0 if (phi1 < -90.0) else phi1
            phi1 = 90.0 if (phi1 > 90.0) else phi1

            # handle special case near the equator
            # phi2 = phi1 (latitudes congruent)
            # h2 = h1 + a1 - a2
            if (np.abs(phi1) < EPS):
                phi2[N] = np.copy(phi1)
                h2[N] = h1 + a1 - a2
            # handle special case near the poles
            # phi2 = phi1 (latitudes congruent)
            # h2 = h1 + b1 - b2
            elif ((90.0 - np.abs(phi1)) < EPS):
                phi2[N] = np.copy(phi1)
                h2[N] = h1 + b1 - b2
            # handle case if latitude is within 45 degrees of equator
            elif (np.abs(phi1) <= 45):
                # convert phi1 to radians
                phi1r = phi1 * np.pi/180.0
                sinphi1 = np.sin(phi1r)
                cosphi1 = np.cos(phi1r)
                # prevent division by very small numbers
                cosphi1 = np.copy(EPS) if (cosphi1 < EPS) else cosphi1
                # calculate tangent
                tanphi1 = sinphi1 / cosphi1
                u1 = np.arctan(b1 / a1 * tanphi1)
                hpr1sin = b1 * np.sin(u1) + h1 * sinphi1
                hpr1cos = a1 * np.cos(u1) + h1 * cosphi1
                # set initial value for u2
                u2 = np.copy(u1)
                # setup constants
                k0 = b2 * b2 - a2 * a2
                k1 = a2 * hpr1cos
                k2 = b2 * hpr1sin
                # perform newton-raphson iteration to solve for u2
                # cos(u2) will not be close to zero since abs(phi1) <= 45
                for i in range(0, ITMAX+1):
                    cosu2 = np.cos(u2)
                    fu2 = k0 * np.sin(u2) + k1 * np.tan(u2) - k2
                    fu2p = k0 * cosu2 + k1 / (cosu2 * cosu2)
                    if (np.abs(fu2p) < EPS):
                        i = np.copy(ITMAX)
                    else:
                        delta = fu2 / fu2p
                        u2 -= delta
                        if (np.abs(delta) < EPS):
                            i = np.copy(ITMAX)
                # convert latitude to degrees and verify values between +/- 90
                phi2r = np.arctan(a2 / b2 * np.tan(u2))
                phi2[N] = phi2r*180.0/np.pi
                phi2[N] = -90.0 if (phi2[N] < -90.0) else phi2[N]
                phi2[N] = 90.0 if (phi2[N] > 90.0) else phi2[N]
                # calculate height
                h2[N] = (hpr1cos - a2 * np.cos(u2)) / np.cos(phi2r)
            # handle final case where latitudes are between 45 degrees and pole
            else:
                # convert phi1 to radians
                phi1r = phi1 * np.pi/180.0
                sinphi1 = np.sin(phi1r)
                cosphi1 = np.cos(phi1r)
                # prevent division by very small numbers
                cosphi1 = np.copy(EPS) if (cosphi1 < EPS) else cosphi1
                # calculate tangent
                tanphi1 = sinphi1 / cosphi1
                u1 = np.arctan(b1 / a1 * tanphi1)
                hpr1sin = b1 * np.sin(u1) + h1 * sinphi1
                hpr1cos = a1 * np.cos(u1) + h1 * cosphi1
                # set initial value for u2
                u2 = np.copy(u1)
                # setup constants
                k0 = a2 * a2 - b2 * b2
                k1 = b2 * hpr1sin
                k2 = a2 * hpr1cos
                # perform newton-raphson iteration to solve for u2
                # sin(u2) will not be close to zero since abs(phi1) > 45
                for i in range(0, ITMAX+1):
                    sinu2 = np.sin(u2)
                    fu2 = k0 * np.cos(u2) + k1 / np.tan(u2) - k2
                    fu2p =  -1 * (k0 * sinu2 + k1 / (sinu2 * sinu2))
                    if (np.abs(fu2p) < EPS):
                        i = np.copy(ITMAX)
                    else:
                        delta = fu2 / fu2p
                        u2 -= delta
                        if (np.abs(delta) < EPS):
                            i = np.copy(ITMAX)
                # convert latitude to degrees and verify values between +/- 90
                phi2r = np.arctan(a2 / b2 * np.tan(u2))
                phi2[N] = phi2r*180.0/np.pi
                phi2[N] = -90.0 if (phi2[N] < -90.0) else phi2[N]
                phi2[N] = 90.0 if (phi2[N] > 90.0) else phi2[N]
                # calculate height
                h2[N] = (hpr1sin - b2 * np.sin(u2)) / np.sin(phi2r)

        # return the latitude and height
        return (phi2, h2)
