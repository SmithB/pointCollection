#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_mosaic.py
Written by Tyler Sutterley (05/2024)

Create a weighted mosaic from a series of tiles

COMMAND LINE OPTIONS:
    --help: list the command line options
    -d X, --directory X: directory to run
    -g X, --glob_string X: quoted string to pass to glob to find the files
    --proj4: projection string for the tiles and output mosaic
    -r X, --range X: valid range of tiles to read [xmin,xmax,ymin,ymax]
    -G X, --group X: input HDF5 group
    -F X, --field X: input HDF5 field map
    -p X, --pad X: pad width in meters for weights
    -f X, --feather X: feathering width in meters for weights
    -w, --weight: use a weighted summation scheme for calculating mosaic
    -S X, --spacing X: output grid spacing if creating from uniform tiles
    -c X, --crop X: crop mosaic to bounds
        [xmin,xmax,ymin,ymax] or [xmin,xmax,ymin,ymax,tmin,tmax]
    -O X, --output X: output filename
    -v, --verbose: verbose output of run
    -R, --replace: overwrite existing output files
    -s, --show: create plot of output mosaic
    -m X, --mode X: Local permissions mode of the output mosaic
    -N --ignore_Nodata: ignore nodata values (except Nan) in inputs

UPDATE HISTORY:
    Updated 05/2024: allow cropping in time for 3D fields
    Updated 01/2021: added option for setting projection attributes
    Updated 10/2021: added option for using a non-weighted summation
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
def main():
    """
    Create a weighted mosaic from a series of tiles
    """

    # account for a bug in argparse that misinterprets negative arguments
    # this bug may be fixed in some versions of argparse
    # but this should keep compatibility
    argv=sys.argv

    for i, arg in enumerate(argv):
        if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

    parser=argparse.ArgumentParser()
    parser.add_argument('--directory','-d',
        type=os.path.expanduser, default=os.getcwd(),
        help='directory to run')
    parser.add_argument('--glob_string','-g', type=str, nargs='+',
        default='/*/*.h5',
        help='quoted string to pass to glob to find the files')
    parser.add_argument('--proj4', type=str, default=None,
        help='projection string for the tiles and output mosaic')
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
    parser.add_argument('--fields','-F', type=str,
        nargs='+', default=None,
        help='input HDF5 fields')
    parser.add_argument('--pad','-p', type=float,
        default=0,
        help='pad width in meters for weights')
    parser.add_argument('--feather','-f', type=float,
        default=0,
        help='feathering width in meters for weights')
    parser.add_argument('--weight','-w',
        action="store_true",
        help='use a weighted summation scheme for calculating mosaic')
    parser.add_argument('--by_band',
        action="store_true",
        help="apply a weighted sum to each band individually, treating invalid pixels as missing data.  If False, invalid pixels in any input will overwrite all valid pixels for the same location")
    parser.add_argument('--spacing','-S', type=float,
        nargs=2, default=[None,None], metavar=('dx','dy'),
        help='output grid spacing if creating from uniform tiles')
    parser.add_argument('--crop','-c', type=float,
        nargs='+',
        help='crop mosaic to bounds')
    parser.add_argument('--t_range', type=float,
        nargs=2,
        help="crop time range to bounds")
    parser.add_argument('--output','-O',
        type=str, default='mosaic.h5',
        help='output filename')
    parser.add_argument('--verbose','-v',
        default=False, action='store_true',
        help='verbose output of run')
    parser.add_argument('--replace','-R',
        default=False, action='store_true',
        help='overwrite existing output files')
    parser.add_argument('--time', type=float,
        help='specify time for one-band DEM')
    parser.add_argument('--show','-s',
        action="store_true",
        help='create plot of output mosaic')
    parser.add_argument('--mode','-m',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output mosaic')
    try:
        assert(len(sys.argv)>1)
        args=parser.parse_args()
    except Exception:
        parser.print_usage()
        sys.exit(1)
    # valid range of tiles
    if np.any(np.isfinite(np.array(args.range))):
        xmin,xmax,ymin,ymax = args.range
    else:
        xmin,xmax,ymin,ymax = [None, None, None, None]
    # convert field mapping from list to dict
    #field_mapping = {args.fields[0]:args.field[1]}

    if args.out_group is None:
        args.out_group=args.in_group

    if args.verbose:
        print("searching in directory "+args.directory+" with  glob string:"+"["+str(args.glob_string)+"]")
    # find list of valid files

    if isinstance(args.glob_string, str):
        initial_file_list = glob.glob(args.directory +'/'+args.glob_string)
    else:
        initial_file_list = []
        for glob_string in args.glob_string:
            initial_file_list += glob.glob(args.directory +'/'+glob_string)

    if args.verbose:
        print(f"initial file list contains {len(initial_file_list)} files")
        #print(initial_file_list)

    if xmin is None:
        file_list=initial_file_list
    else:
        file_list = []
        for file in initial_file_list:
            try:
                xc,yc=[int(item)*1.e3 for item in re.compile(r'E([0-9-.]+)_N([0-9-.]+).h5').search(file).groups()]
            except Exception:
                continue
            if ((xc >= xmin) and (xc <= xmax) & (yc >= ymin) and (yc <= ymax)):
                file_list.append(file)

    if args.verbose:
        print(f"found {len(file_list)} files")

    if not args.weight:
        args.pad=None
        args.feather=None

    mosaic=pc.grid.mosaic(spacing=args.spacing).from_list(file_list,
                                      fields=args.fields,
                                      group=args.in_group,
                                      pad=args.pad,
                                      feather=args.feather,
                                      by_band=args.by_band)
    if isinstance(mosaic, str):
        if args.verbose:
            print(f"pc.grid.mosaic failed for group {args.in_group} and fields {args.fields} with message:")
            print(mosaic)
            return
    # crop mosaic to bounds
    if np.any(args.crop):
        # keyword arguments for cropping
        kwds = dict(fields=mosaic.fields)
        # x and y range (verify min and max order)
        XR = np.sort([args.crop[0],args.crop[1]])
        YR = np.sort([args.crop[2],args.crop[3]])
        # add time range if 3D
        if len(args.crop) == 6:
            kwds['TR'] = np.sort([args.crop[4],args.crop[5]])
        # crop the mosaic
        mosaic.crop(XR, YR, **kwds)

    if args.t_range is not None:
        kwds = dict(fields=mosaic.fields)
        kwds['TR']=args.t_range
        mosaic.crop(None, None, **kwds)

    if args.time is not None:
        mosaic.t = args.time
    # output each field
    for field in mosaic.fields:
        if mosaic.field_dims[field] == 2:
            pc.grid.data().from_dict({'x':mosaic.x,'y':mosaic.y,\
                                      field:getattr(mosaic,field)})\
                .to_h5(os.path.join(args.directory,args.output), \
                       group=args.out_group, replace=args.replace, \
                       srs_proj4=args.proj4)
        else:
            pc.grid.data().from_dict({'x':mosaic.x,'y':mosaic.y, 't': mosaic.t,\
                               field:getattr(mosaic,field)})\
                .to_h5(os.path.join(args.directory,args.output), \
                       group=args.out_group, replace=args.replace, \
                       srs_proj4=args.proj4)
        # only want the 'replace' argument on the first field
        args.replace=False
    if args.show:
        if len(mosaic.z.shape) > 2:
            plt.imshow(mosaic.z[:,:,-1]-mosaic.z[:,:,0], extent=mosaic.extent)
        else:
            mosaic.z.show()
        plt.colorbar()
        plt.show()

if __name__=='__main__':
    main()
