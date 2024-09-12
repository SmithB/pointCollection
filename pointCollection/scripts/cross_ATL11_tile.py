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

def read_data_from_index(index_file, xy0, W, field_dict_11=None):

    D = pc.geoIndex().from_file(index_file).query_xy_box(
                xy0[0]+np.array([-1, 1])*W/2, xy0[1]+np.array([-1, 1])*W/2, \
                field_dict=field_dict_11)
    return D

def read_data_from_glob(glob_str, field_dict_11=None):
    D=[]
    for file in glob.glob(glob_str):
        for pair in [1, 2, 3]:
            try:

                D += [pc.ATL11.data(pair=pair).from_h5(file, field_dict = field_dict_11)]
            except Exception as e:
                print(e)
                print(f"no data for pair {pair}, file {file}")
    return D

def ATL11_crossovers(D, EPSG=3031):
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
    v_fields=['latitude', 'longitude', 'h_corr', 'h_corr_sigma', 'h_corr_sigma_systematic', 'delta_time', 'ref_pt', 'rgt', 'pair', 'dem_h', 'x_atc', 'fit_quality', 'cycle_number', 'x', 'y']

    if 'dh_geoloc' in xover_list[0]['data_0'].fields:
        v_fields += ['dh_geoloc']

    v_list=[{field:np.zeros([len(xover_list), N_cols])+np.nan for field in v_fields} for ii in [0, 1]]
    for count, X in enumerate(xover_list):
        for data, li, vi in zip([X['data_0'],X['data_1']], [X['L0'],X['L1']], v_list):
            w=np.array([1-li, li])
            for field in v_fields:
                vi[field][count,:]=w.dot(getattr(data, field))
    v=[pc.data(columns=N_cols).from_dict(vi) for vi in v_list]

    return v

def write_xovers(v, out_dir, xy0, out_file):
    if out_dir is not None and not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    if out_file is None and xy0 is not None:
        out_file=os.path.join(out_dir, f'E{int(np.round(xy0[0]/1000))}_N{int(np.round(xy0[0]/1000))}.h5')
    elif out_file is None:
        out_file=os.path.join(out_dir, 'all_xovers.h5')
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
    parser.add_argument('--glob_str','-g', type=str, help="glob string to find files to cross, to be specified if xy0 and W are not")
    parser.add_argument('--out_file', type=str, help="output file, to be specified if xy0 is not")
    parser.add_argument('--hemisphere', '-H', type=int, help="hemisphere, -1 for Antarctica, 1, for Arctic")
    parser.add_argument('--queue','-q', action="store_true")
    args=parser.parse_args()


    if args.hemisphere==1:
        EPSG=3413
    else:
        EPSG=3031

    if args.out_dir is not None and not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    if args.queue:
        make_queue(args.index_file, args.W, args.hemisphere, args.out_dir)
        return

    field_dict_11={None:['latitude','longitude','delta_time',\
                        'h_corr','h_corr_sigma','h_corr_sigma_systematic', 'ref_pt'],\
                        '__calc_internal__' : ['rgt'],
                        'cycle_stats' : {'tide_ocean','dac','dh_geoloc'},
                        'ref_surf':['e_slope','n_slope', 'x_atc', 'fit_quality', 'dem_h']}

    if args.glob_str is not None:
        D=read_data_from_glob(args.glob_str, field_dict_11=field_dict_11)
    else:
        D=read_data_from_index(args.index_file, args.xy0, args.W, field_dict_11=field_dict_11)

    v=ATL11_crossovers(D, EPSG=EPSG)

    if len(v) > 0:
        write_xovers(v, args.out_dir, args.xy0, args.out_file)

if __name__=='__main__':
    main()
