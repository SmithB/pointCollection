#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 18:08:00 2019

@author: ben
"""
import pointCollection as pc
import ATL11
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import scipy.stats as sps
import argparse

def get_map(map_file, bounds=None):
    #if bounds is None:
    #    print("map file is " + map_file)
    #    backgrnd=pc.grid.data().from_geotif(map_file)
    #else:
    backgrnd=pc.grid.data().from_geotif(map_file, bounds=bounds)

    if backgrnd.z.ndim >2:
        backgrnd.z=backgrnd.z[:,:,0].astype(float)
    cr=sps.scoreatpercentile(backgrnd.z[np.isfinite(backgrnd.z) & (backgrnd.z > np.min(backgrnd.z)) ], [16, 84])
    backgrnd.show(cmap='gray', vmin=cr[0]-np.diff(cr), vmax=cr[1]+np.diff(cr))


class IndexPicker:
    def __init__(self, fig, index_file, map_file=None, W=4e4, datatype='ATL06'):
        self.xy=pc.geoIndex().from_file(index_file).bins_as_array()
        self.index_file=index_file
        self.map_file=map_file
        self.fig=fig
        self.datatype=datatype
        self.W=4.e4
        if fig is None:
            self.fig=plt.figure()
            if map_file is not None:
                get_map(map_file, None)
            print("about to plot %d points" % len(self.xy[0]))
            plt.plot(self.xy[0], self.xy[1],'.')
            print("starting picker")
        self.cid=self.fig.canvas.mpl_connect('button_press_event', self.click)
        plt.show(block=True)
        #self.fig.canvas.draw()

    def click(self, event):
        if not event.dblclick:
            return
        print('click', event)
        print(self)
        if self.datatype=='ATL06':
            field_dict={None:['delta_time','h_li','h_li_sigma','latitude','longitude','atl06_quality_summary','segment_id','sigma_geo_h'],
                        'fit_statistics':['dh_fit_dx'],
                        'ground_track':['x_atc', 'sigma_geo_xt','sigma_geo_at'],
                        'geophysical' : ['dac','tide_ocean','r_eff'],
                        'orbit_info':['rgt','cycle_number'],
                        'derived':['valid','matlab_time','LR','BP','spot','rss_along_track_dh']}
        elif self.datatype=='ATL11':
                field_dict={\
                            'corrected_h':['latitude','longitude','delta_time','h_corr','h_corr_sigma','ref_pt'],\
                            'cycle_stats':['ATL06_summary_zero_count'],\
                            'crossing_track_data':['delta_time','h_corr','h_corr_sigma','ref_pt','rgt',\
                               'spot_crossing', 'along_track_rss',\
                               'atl06_quality_summary','cycle_number'],
                            'ref_surf':['x_atc','y_atc']}
        elif self.datatype=='CS2':
            field_dict={'None':['x','y','time','h']}
        W=self.W
        self.fig.gca().plot(event.xdata, event.ydata,'m*')
        self.fig.canvas.draw()
        xy0=np.round(np.array([event.xdata, event.ydata])/1.e4)*1.e4
        print(xy0)
        D=pc.geoIndex().from_file(self.index_file).query_xy_box(\
                   event.xdata+np.array([-W/2, W/2]), event.ydata+np.array([-W/2, W/2]), fields=field_dict)

        #D=geo_index().from_file(self.index_file).query_xy_box((self.xy[0][best], self.xy[1][best]), fields=field_dict)
        print(f"click: index_file is {self.index_file}")

        TrackPicker(D, self.map_file, self.index_file, self.datatype)

class TrackPicker:
    def __init__(self, D, map_file, index_file, datatype):
        self.files=[]
        self.datatype=datatype
        self.map_file=map_file
        self.index_file=index_file
        print(f"TrackPicker: index_file is {index_file}")
        srs_proj4=pc.geoIndex().from_file(index_file).attrs['SRS_proj4']
        if datatype == 'ATL06':
            for ii, Di in enumerate(D):
                Di.get_xy(srs_proj4)
                Di.assign({'file_num':np.zeros_like(Di.x)+ii})
                self.files += [Di.filename]
            self.D_all=pc.data().from_list(D)
        elif datatype == 'ATL11':
           D_list=[]
           for ii, Di in enumerate(D):
               self.files += [Di.filename]
               Di.get_xy(srs_proj4)
               D_list.append(pc.data().from_dict({
                       'x':Di.x, 'y':Di.y, 'file_num':np.zeros_like(Di.x)+ii,
                       'h_corr':Di.corrected_h.h_corr[:,-1].ravel()}))
           self.D_all=pc.data().from_list(D_list)
        else:
            self.D_all=pc.data().from_list(D)
        self.D=D
        XR=[np.nanmin(self.D_all.x), np.nanmax(self.D_all.x)]
        YR=[np.nanmin(self.D_all.y), np.nanmax(self.D_all.y)]
        self.fig=plt.figure()
        if map_file is not None:
            get_map(map_file, bounds=[XR, YR])
        #for Di in D:
        #    plt.plot(Di.x, Di.y,'.')
        if self.datatype=='ATL06':
            plt.scatter(self.D_all.x, self.D_all.y, 6, c=self.D_all.r_eff, linewidth=0, vmax=1, vmin=0)
        elif self.datatype=='ATL11':
            plt.scatter(self.D_all.x, self.D_all.y, 6, c=self.D_all.h_corr); plt.colorbar()
        elif self.datatype=='ATM_waveform_fit':
            plt.scatter(self.D_all.x, self.D_all.y, 6, c=np.log10(self.D_all.K0))
        elif self.datatype=='CS2':
            plt.scatter(self.D_all.x, self.D_all.y, 6, c=self.D_all.h); plt.colorbar()
            hax=plt.axes([0.7, 0.05, 0.25, 0.25])
            hax.hist((self.D_all.time-730486)/365.25+2000, 100)
        self.cid=self.fig.canvas.mpl_connect('button_press_event', self)
        plt.show()

    def __call__(self, event):
        xyp=event.xdata+1j*event.ydata
        best = np.argmin(np.abs(self.D_all.x+1j*self.D_all.y - xyp))
        this = int(self.D_all.file_num[best])
        if self.datatype=='ATL06':
            self.ATL06_plot(self.D[this], this)
        else:
            self.ATL11_plot(self.D[this], this)

    def ATL06_plot(self, D, this):
        fig=plt.figure()
        plt.scatter(D.x_atc, D.h_li, c=np.log10(D.r_eff), cmap='rainbow', vmin=np.log10(0.25), vmax=np.log10(4))
        plt.colorbar()
        plt.title(self.files[this])
        print(self.files[this])
        fig.canvas.draw()
        plt.show()

    def ATL11_plot(self, D, this):
        print(self.files[this])

        fig=plt.figure()
        plt.subplot(211)
        plt.scatter(D.corrected_h.ref_pt*20, D.corrected_h.h_corr[:,2], 4, D.corrected_h.h_corr_sigma[:,2]); plt.colorbar()
        plt.subplot(212)
        plt.plot(D.corrected_h.ref_pt*20, D.corrected_h.h_corr[:,2], 'k')
        ref, xo, delta=D.get_xovers()
        print(delta)
        plt.scatter(ref.ref_pt*20, xo.h, 12, c=delta.h, linewidth=0); plt.colorbar()
        plt.plot()
        fig.canvas.draw()
        plt.show()



def main():
    #fig=plt.figure()
    fig=None
    parser=argparse.ArgumentParser("IndexPicker: class to make clickable maps of a geoIndex")
    parser.add_argument("index_file", type=str)
    parser.add_argument("data_type", type=str)
    parser.add_argument("--map_file",'-m', type=str)
    args=parser.parse_args()

    IndexPicker(fig, args.index_file, map_file=args.map_file, datatype=args.data_type)
    #plt.show(block=True)

if __name__=='__main__':
    main()
