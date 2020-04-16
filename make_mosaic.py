#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_mosaic.py
Written by Tyler Sutterley (03/2020)

Create a weighted mosaic from a series of tiles

COMMAND LINE OPTIONS:
    --help: list the command line options
    -d X, --directory X: directory to run
    -g X, --glob_string X: quoted string to pass to glob to find the files
    -r X, --range X: valid range of tiles [xmin,xmax,ymin,ymax]
    -G X, --group X: input HDF5 group
    -F X, --field X: input HDF5 field map
    -p X, --pad X: pad width in meters for weights
    -f X, --feather X: feathering width in meters for weights
    -O X, --output X: output filename
    -v, --verbose: verbose output of run
    -m X, --mode X: Local permissions mode of the output mosaic

UPDATE HISTORY:
    Updated 03/2020: adding argparse bug fix for negative arguments
        made output filename a command line option
    Written 03/2020
"""
import os
import re
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pointCollection as pc

import sys
def main(argv):
    """
    Create a weighted mosaic from a series of tiles
    """

    # account for a bug in argparse that misinterprets negative arguments
    # this bug may be fixed in some versions of argparse
    # but this should keep compatibility
    for i, arg in enumerate(argv):
        if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

    parser=argparse.ArgumentParser()
    parser.add_argument('--directory','-d', type=str, default=os.getcwd(), help='directory to run')
    parser.add_argument('--glob_string','-g', type=str, default='/*/*.h5', help='quoted string to pass to glob to find the files')
    parser.add_argument('--range','-r', type=float, nargs='+', help='valid range of tiles [xmin,xmax,ymin,ymax]')
    parser.add_argument('--group','-G', type=str, default='dz/', help='input HDF5 group')
    parser.add_argument('--field','-F', type=str, nargs='+', default=['z','dz'], help='input HDF5 field map')
    parser.add_argument('--pad','-p', type=float, default=1e4, help='pad width in meters for weights')
    parser.add_argument('--feather','-f', type=float, default=2e4, help='feathering width in meters for weights')
    parser.add_argument('--output','-O', type=str, default='mosaic.h5',  help='output filename')
    parser.add_argument('--verbose','-v', default=False, action='store_true', help='verbose output of run')
    parser.add_argument('--show','-s', action="store_true")
    parser.add_argument('--mode','-m', default=0o775, help='permissions mode of output mosaic')
    args=parser.parse_args()

    # valid range of tiles
    if args.range:
        xmin,xmax,ymin,ymax = args.range
    else:
        xmin,xmax,ymin,ymax = [-np.inf,np.inf,-np.inf,np.inf]
    # convert field mapping from list to dict
    field_mapping = {args.field[0]:args.field[1]}

    # find list of valid files
    file_list = []
    for file in glob.glob(args.directory +'/'+args.glob_string):
        xc,yc=[int(item)*1.e3 for item in re.compile('E(.*)_N(.*).h5').search(file).groups()]
        if ((xc >=xmin) and (xc <= xmax) & (yc >= ymin) and (yc <= ymax)):
            file_list.append(file)

    # get bounds, grid spacing and dimensions of output mosaic
    mosaic=pc.grid.mosaic()
    for file in file_list:
        # read ATL14 grid from HDF5
        temp=pc.grid.mosaic().from_h5(file, group=args.group, field_mapping=field_mapping)
        # update grid spacing of output mosaic
        mosaic.update_spacing(temp)
        # update the extents of the output mosaic
        mosaic.update_bounds(temp)
        # update dimensions of output mosaic with new extents
        mosaic.update_dimensions(temp)

    # create output mosaic
    mosaic.data = np.zeros(mosaic.dimensions)
    mosaic.mask = np.ones(mosaic.dimensions,dtype=np.bool)
    mosaic.weight = np.zeros((mosaic.dimensions[0],mosaic.dimensions[1]))
    for file in file_list:
        # read ATL14 grid from HDF5
        temp=pc.grid.mosaic().from_h5(file, group=args.group, field_mapping=field_mapping)
        temp=temp.weights(pad=args.pad, feather=args.feather, apply=True)
        # get the image coordinates of the input file
        iy,ix = mosaic.image_coordinates(temp)
        for band in range(mosaic.dimensions[2]):
            mosaic.data[iy,ix,band] += temp.z[:,:,band]
            mosaic.mask[iy,ix,band] = False
        # add weights to total weight matrix
        mosaic.weight[iy,ix] += temp.weight[:,:]

    # find valid weights
    iy,ix = np.nonzero(mosaic.weight == 0)
    mosaic.mask[iy,ix,:] = True
    # normalize weights
    iy,ix = np.nonzero(mosaic.weight > 0)
    for band in range(mosaic.dimensions[2]):
        mosaic.data[iy,ix,band] /= mosaic.weight[iy,ix]
    # replace invalid points with fill_value
    mosaic.data[mosaic.mask] = mosaic.fill_value

    # output mosaic to HDF5
    fields=['x','y', 'data','weight']
    if mosaic.t.dtype in (int, float, np.int, np.float):
        fields += 't'
    mosaic.to_h5(os.path.join(args.directory,args.output),
        field_list=fields)
    os.chmod(os.path.join(args.directory,args.output), args.mode)

    if args.show:
        if len(mosaic.z.shape) > 2:
            plt.imshow(mosaic.z[:,:,-1]-mosaic.z[:,:,0], extent=mosaic.extent)
        else:
            mosaic.z.show()
        plt.colorbar()
        plt.show()
if __name__=='__main__':
    main(sys.argv)
