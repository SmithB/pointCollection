# -*- coding: utf-8 -*-
"""
Created on Wed May  1 10:40:57 2019

@author: ben
"""

import pointCollection as pc

import numpy as np
import os
import h5py
#import matplotlib.pyplot as plt
import glob
#import re
#import sys


def ATL11_crossovers(index_file, xy0, W, EPSG=3031):
    field_dict_11={None:['latitude','longitude','delta_time',\
                        'h_corr','h_corr_sigma','h_corr_sigma_systematic', 'ref_pt'],\
                        '__calc_internal__' : ['rgt'],
                        'cycle_stats' : {'tide_ocean','dac'},
                        'ref_surf':['e_slope','n_slope', 'x_atc', 'fit_quality', 'dem_h']}
    D = pc.geoIndex().from_file(index_file).query_xy_box(
                xy0[0]+np.array([-1, 1])*W/2, xy0[1]+np.array([-1, 1])*W/2, \
                field_dict=field_dict_11)
    for Di in D:
        Di.get_xy(EPSG=EPSG)
        Di.assign({'pair':np.zeros_like(Di.rgt)+Di.pair})
    
    xover_list=list()
    #plt.clf()
    count=0
    for ii in np.arange(len(D)-1):
        for jj in np.arange(ii+1,len(D)):
            temp=[pc.data().from_dict({'x':Di.x[:,0], 'y':Di.y[:,0]}) for Di in [D[ii], D[jj]]]
            xyC, inds, L = pc.cross_tracks(temp, delta=60, delta_coarse=500)
            #print(xyC)
            if xyC is not None:
                xover_list.append({'xyC':xyC, 'data_0':D[ii][inds[0]], 'data_1':D[jj][inds[1]], 'L0':L[0], 'L1':L[1]})
    N_cols=xover_list[0]['data_0'].x.shape[1]
    v_fields=['latitude', 'longitude', 'h_corr', 'h_corr_sigma', 'h_corr_sigma_systematic', 'delta_time', 'quality_summary', 'ref_pt', 'rgt', 'pair', 'dem_h', 'x_atc', 'fit_quality', 'cycle_number', 'x', 'y']

    v_list=[{field:np.zeros([len(xover_list), N_cols])+np.NaN for field in v_fields} for ii in [0, 1]]
    for count, X in enumerate(xover_list):
        for data, li, vi in zip([X['data_0'],X['data_1']], [X['L0'],X['L1']], v_list):
            w=np.array([1-li, li])
            for field in v_fields:
                vi[field][count,:]=w.dot(getattr(data, field))
    v=[pc.data(columns=N_cols).from_dict(vi) for vi in v_list]
    
    ind=(np.abs(v[0].x -xy0[0])<=W/2) & (np.abs(v[0].y -xy0[1])<=W/2) 
    for ii in [0,1]:
        v[ii].index(ind)
    
    return v

def write_xovers(v, out_dir, xy0):
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    out_file=os.path.join(out_dir, f'E{int(np.round(xy0[0]/1000))}_N{int(np.round(xy0[0]/1000))}.h5')
    v[0].to_h5(out_file, group='data_0', replace=True)
    v[1].to_h5(out_file, group='data_1', replace=False)
    
def make_queue(index_file, W, hemisphere, out_dir):

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    
    xy = pc.geoIndex().from_file(index_file).bins_as_array()
    xyB = np.unique(np.round((xy[:,0]+1j*xy[:,1])/W))*W
    
    for Bi in xyB:
        print(f'cross_ATL11_tile.py --xy0 {int(np.real(Bi))} {int(np.imag(Bi))} -W {W} -H {hemisphere} --out_dir {out_dir}')
        
def main():
    import argparse
    parser=argparse.ArgumentParser(description='Find crossovers in an ATL11 tile')
    parser.add_argument('--index_file', type=str, help="glob which matches the tiles")
    parser.add_argument('--xy0', type=float, nargs=2, help="fit center location")
    parser.add_argument('--Width','-W',  type=float, help="fit width")
    parser.add_argument('--out_dir', type=str, help="output directory")
    parser.add_argument('--hemisphere', '-H', type=int, help="hemisphere, -1 for Antarctica, 1, for Arctic")
    parser.add_argument('--queue','-q', action="store_true")
    args=parser.parse_args()
    
    hemisphere=args.hemisphere
    out_dir=args.out_dir
    
    if hemisphere==1:
        EPSG=3413
    else:
        EPSG=3031
    
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    if args.queue:
        make_queue(args.index_file, args.W, hemisphere, out_dir)
        return

    v=ATL11_crossovers(args.index_file, args.xy0, args.W, EPSG=EPSG)
    if len(xover_list) > 0:
        write_xovers(v, args.out_dir, args.xy0)


if __name__=='__main__':
    main()
