#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_mosaic.py
Written by Tyler Sutterley (07/2021)

Create a weighted mosaic from a series of tiles

COMMAND LINE OPTIONS:
    --help: list the command line options
    -d X, --directory X: directory to run
    -g X, --glob_string X: quoted string to pass to glob to find the files
    -r X, --range X: valid range of tiles to read [xmin,xmax,ymin,ymax]
    -G X, --group X: input HDF5 group
    -F X, --field X: input HDF5 field map
    -p X, --pad X: pad width in meters for weights
    -f X, --feather X: feathering width in meters for weights
    -S X, --spacing X: output grid spacing if creating from uniform tiles
    -c X, --crop X: crop mosaic to bounds [xmin,xmax,ymin,ymax]
    -O X, --output X: output filename
    -v, --verbose: verbose output of run
    -R, --replace: overwrite existing output files
    -s, --show: create plot of output mosaic
    -m X, --mode X: Local permissions mode of the output mosaic

UPDATE HISTORY:
    Updated 07/2021: added option replace for overwriting existing files
        added option for cropping output mosaic
    Updated 11/2020: added option spacing for setting output grid
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
    parser.add_argument('--directory','-d',
        type=os.path.expanduser, default=os.getcwd(),
        help='directory to run')
    parser.add_argument('--glob_string','-g', type=str,
        default='/*/*.h5',
        help='quoted string to pass to glob to find the files')
    parser.add_argument('--range','-r', type=float,
        nargs=4, default=[-np.inf,np.inf,-np.inf,np.inf],
        metavar=('xmin','xmax','ymin','ymax'),
        help='valid range of tiles to read')
    parser.add_argument('--in_group','-G',
        type=str, default='/',
        help='input HDF5 group')
    parser.add_argument('--out_group', type=str,
        help="""output hdf5 group (specify "/" for root;
        default is the same as the input group""")
    # parser.add_argument('--field','-F', type=str,
    #     nargs='+', default=['z','dz'],
    #     help='input HDF5 field map')
    parser.add_argument('--fields','-F', type=str,
        nargs='+', default=['z'],
        help='input HDF5 fields')
    parser.add_argument('--pad','-p', type=float,
        default=0,
        help='pad width in meters for weights')
    parser.add_argument('--feather','-f', type=float,
        default=0,
        help='feathering width in meters for weights')
    parser.add_argument('--spacing','-S', type=float,
        nargs=2, default=[None,None], metavar=('dx','dy'),
        help='output grid spacing if creating from uniform tiles')
    parser.add_argument('--crop','-c', type=float,
        nargs=4, metavar=('xmin','xmax','ymin','ymax'),
        help='crop mosaic to bounds')
    parser.add_argument('--output','-O',
        type=str, default='mosaic.h5',
        help='output filename')
    parser.add_argument('--verbose','-v',
        default=False, action='store_true',
        help='verbose output of run')
    parser.add_argument('--replace','-R',
        default=False, action='store_true',
        help='overwrite existing output files')
    parser.add_argument('--show','-s',
        action="store_true",
        help='create plot of output mosaic')
    parser.add_argument('--mode','-m',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output mosaic')
    args=parser.parse_args()

    # valid range of tiles
    xmin,xmax,ymin,ymax = args.range
    # convert field mapping from list to dict
    #field_mapping = {args.field[0]:args.field[1]}

    if args.out_group is None:
        args.out_group=args.in_group

    if args.verbose:
        print("searching glob string:"+"["+args.glob_string+"]")
    # find list of valid files
    file_list = []
    for file in glob.glob(args.directory +'/'+args.glob_string):
        try:
            xc,yc=[int(item)*1.e3 for item in re.compile(r'E(.*)_N(.*).h5').search(file).groups()]
        except Exception:
            continue
        if ((xc >= xmin) and (xc <= xmax) & (yc >= ymin) and (yc <= ymax)):
            file_list.append(file)

    # get bounds, grid spacing and dimensions of output mosaic
    mosaic=pc.grid.mosaic(spacing=args.spacing)
    for file in file_list.copy():
        # read tile grid from HDF5
        try:
            temp=pc.grid.data().from_h5(file, group=args.in_group, fields=[])
            # update grid spacing of output mosaic
            mosaic.update_spacing(temp)
            # update the extents of the output mosaic
            mosaic.update_bounds(temp)
            # update dimensions of output mosaic with new extents
            mosaic.update_dimensions(temp)
        except Exception:
            print(f"failed to read group {args.in_group} "+ file)
            file_list.remove(file)
    # create output mosaic
    mosaic.assign({field:np.zeros(mosaic.dimensions) for field in args.fields})
    #mosaic.data = np.zeros(mosaic.dimensions)
    mosaic.invalid = np.ones(mosaic.dimensions,dtype=np.bool)

    # read data grid from a single tile HDF5
    temp=pc.grid.mosaic().from_h5(file_list[0], group=args.in_group, fields=args.fields)
    these_fields=[field for field in args.fields if field in temp.fields]
    field_dims={field:getattr(temp, field).ndim for field in these_fields}
    # create the weights for a single file
    # as the tiles have the same dimensions, we only have to calculate
    # the first set.  After that we can just copy the first
    temp.weights(pad=args.pad, feather=args.feather, apply=False)

    # check if weight matrix exists
    weight_file = 'w0p{0:0.0f}f{1:0.0f}.h5'.format(args.pad,args.feather)
    if os.access(os.path.join(args.directory,weight_file),os.F_OK):
        # read input weight matrix
        mosaic.weight=pc.grid.data().from_h5(
            os.path.join(args.directory,weight_file),
            group=args.out_group, fields=['weight']).data
    else:
        # allocate for output weight matrix
        mosaic.weight=np.zeros((mosaic.dimensions[0],mosaic.dimensions[1]))
        tile_weight=temp.weight.copy()
        for count, file in enumerate(file_list):
            # read data grid from HDF5
            temp=pc.grid.mosaic().from_h5(file, group=args.in_group, fields=args.fields)
            these_fields=[field for field in args.fields if field in temp.fields]
            # copy weights for tile
            temp.weight=tile_weight.copy()
            # get the image coordinates of the input file
            iy,ix = mosaic.image_coordinates(temp)
            for field in these_fields:
                field_data=np.atleast_3d(getattr(temp, field))
                try:
                    for band in range(mosaic.dimensions[2]):
                        getattr(mosaic, field)[iy,ix,band] += field_data[:,:,band]*temp.weight
                        mosaic.invalid[iy,ix,band] = False
                except IndexError:
                    print(f"problem with field {field} in file {file}")
            # add weights to total weight matrix
            mosaic.weight[iy,ix] += temp.weight[:,:]
            # save weights to file
            pc.grid.data().from_dict({'x':mosaic.x,'y':mosaic.y,
                'weight':mosaic.weight}).to_h5(os.path.join(args.directory,weight_file), \
                group=args.out_group, replace=args.replace)

    # find valid weights
    iy,ix = np.nonzero(mosaic.weight == 0)
    mosaic.invalid[iy,ix,:] = True
    # normalize weights
    iy,ix = np.nonzero(mosaic.weight > 0)
    for band in range(mosaic.dimensions[2]):
        for field in mosaic.fields:
            getattr(mosaic, field)[iy,ix,band] /= mosaic.weight[iy,ix]
    # replace invalid points with fill_value
    for field in mosaic.fields:
        getattr(mosaic, field)[mosaic.invalid] = mosaic.fill_value
    # crop mosaic to bounds
    if np.any(args.crop):
        # x and y range (verify min and max order)
        XR = np.sort([args.range[0],args.range[1]])
        YR = np.sort([args.range[2],args.range[3]])
        mosaic.crop(XR, YR, fields=mosaic.fields)

    # output each field
    for field in mosaic.fields:
        if field_dims[field] == 2:
            pc.grid.data().from_dict({'x':mosaic.x,'y':mosaic.y,\
                               field:np.squeeze(getattr(mosaic,field)[:,:,0])})\
                .to_h5(os.path.join(args.directory,args.output), \
                       group=args.out_group, replace=args.replace)
        else:
            pc.grid.data().from_dict({'x':mosaic.x,'y':mosaic.y, 't': mosaic.t,\
                               field:getattr(mosaic,field)})\
                .to_h5(os.path.join(args.directory,args.output), \
                       group=args.out_group, replace=args.replace)

    if args.show:
        if len(mosaic.z.shape) > 2:
            plt.imshow(mosaic.z[:,:,-1]-mosaic.z[:,:,0], extent=mosaic.extent)
        else:
            mosaic.z.show()
        plt.colorbar()
        plt.show()
if __name__=='__main__':
    main(sys.argv)
