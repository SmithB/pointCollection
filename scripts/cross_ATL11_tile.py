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

def ATL11_crossovers(index_file, xy0, W):

    D = pc.geoIndex.from_file(index_file).query_xy_box(
        xy0[0]+np.array([-1, 1])*W/2, xy0[1]+np.array([-1, 1])*W/2)
 

    xover_list=list()
    #plt.clf()
    for ii in np.arange(len(D)):
        for jj in np.arange(len(D)):
            if ii == jj:
                continue 
            xyC, inds, L = pc.cross_tracks([D[ii], D[jj]], delta=60, delta_coarse=1000)
            if xyC is not None:
                try:
                    xover_list.append({'xyC':xyC, 'data_0':D[ii][inds[0]], 'data_1':D[jj][inds[1]], 'L0':L[0], 'L1':L[1]})
                except Exception as e:
                    print("HERE")
    return xover_list


def write_xovers(xover_list, out_file):
    if os.path.isfile(out_file):
        os.remove(out_file)
    with h5py.File(out_file,'w') as h5f:
        for key_D in ['data_0', 'data_1']:
            group='/'+key_D
            h5f.create_group(group)

            key_L=key_D.replace('data_','L')
            L=np.c_[[item[key_L] for item in xover_list]]

            Dtemp=[item[key_D] for item in xover_list]
            Dtemp=pc.data().from_list(Dtemp)
            shape=[np.int(Dtemp.size/2), 2]
            Dtemp.shape=shape
            for key in Dtemp.fields:
                temp=getattr(Dtemp, key)
                temp.shape=shape
                h5f.create_dataset(group+'/'+key, data=temp)
            h5f.create_dataset(group+'/W', data=np.c_[1-L, L])

        xy=np.c_[[item['xyC'] for item in xover_list]]
        h5f.create_dataset('/x', data=xy[:,0])
        h5f.create_dataset('/y', data=xy[:,1])
        try:
            h5f.create_dataset('/slope_x', data=np.array([item['slope_x'] for item in xover_list]))
            h5f.create_dataset('/slope_y', data=np.array([item['slope_y'] for item in xover_list]))
            h5f.create_dataset('/grounded', data=np.array([item['grounded'] for item in xover_list]))
        except Exception as e:
            pass
    return #xover_list

def read_xovers(xover_dir):

    tiles=glob.glob(xover_dir+'/*.h5')
    with h5py.File(tiles[0],'r') as h5f:
        fields=[key for key in h5f['D0'].keys()]
    D=[]
    X=[]
    for tile in glob.glob(xover_dir+'/*.h5'):
        D.append([pc.data().from_file(tile, field_dict={gr : fields}) for gr in ['D0','D1']])
        with h5py.open(tile,'r') as h5f:
            X.append(pc.data(fields=['x','y']).from_file(tile,field_dict={None:['x','y']}))
    return D, X

def make_queue(files, hemisphere, out_dir):

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    
    for file in files:
        if not os.path.isfile(out_dir+'/'+os.path.basename(file)):
            # N.B.  this requires that cross_ATL06_tile is on the unix path.
            if hemisphere is not None:
                print('cross_ATL06_tile.py %s %d %s' %   (file, hemisphere, out_dir))
            else:
                print('cross_ATL06_tile.py %s %s' %   (file, out_dir)) 

def main():
    import argparse
    parser=argparse.ArgumentParser(description='Find crossovers in an ATL11 tile')
    parser.add_argument('index_file', type=str, help="glob which matches the tiles")
    parser.add_argument('--xy0', type=float, nargs=2, help="fit center location")
    parser.add_argument('--Width','-W',  type=float, help="fit width")
    parser.add_argument('out_dir', type=str, help="output directory")
    parser.add_argument('--hemisphere', '-H', type=int, help="hemisphere, -1 for Antarctica, 1, for Arctic")
    parser.add_argument('--queue','-q', action="store_true")
    args=parser.parse_args()
    
    hemisphere=args.hemisphere
    out_dir=args.out_dir
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    if args.queue:
        make_queue(args.index_file, args.xy0, args.W, hemisphere, out_dir)
        return

    xover_list=ATL11_crossovers(args.index_file, args.xy0, args.W)
    if len(xover_list) > 0:
        write_xovers(xover_list, args.out_dir, args.xy0)


if __name__=='__main__':
    main()
