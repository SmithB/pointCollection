#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_mosaic.py
Written by Tyler Sutterley (01/2022)

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
    -c X, --crop X: crop mosaic to bounds [xmin,xmax,ymin,ymax]
    -O X, --output X: output filename
    -v, --verbose: verbose output of run
    -R, --replace: overwrite existing output files
    -s, --show: create plot of output mosaic
    -m X, --mode X: Local permissions mode of the output mosaic
    -N --ignore_Nodata: ignore nodata values (except Nan) in inputs

UPDATE HISTORY:
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
    parser.add_argument('--time', type=float,
        help='specify time for one-band DEM')
    parser.add_argument('--show','-s',
        action="store_true",
        help='create plot of output mosaic')
    parser.add_argument('--mode','-m',
        type=lambda x: int(x,base=8), default=0o775,
        help='permissions mode of output mosaic')
    args=parser.parse_args()

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

    #get bounds, grid spacing and dimensions of output mosaic
    mosaic=pc.grid.mosaic(spacing=args.spacing)
    for file in file_list.copy():
        # read tile grid from HDF5
        try:
            temp=pc.grid.data().from_h5(file, group=args.in_group, fields=[])
            #update the mosaic bounds to include this tile
            mosaic.update_bounds(temp)
        except Exception:
            print(f"failed to read group {args.in_group} "+ file)
            file_list.remove(file)

    mosaic.update_spacing(temp)
    mosaic.update_dimensions(temp)

    # create output mosaic, weights, and mask
    # read data grid from the first tile HDF5, use it to set the field dimensions
    temp=pc.grid.mosaic().from_h5(file_list[0], group=args.in_group, fields=args.fields)
    if len(temp.fields) ==0:
        print(f"make_mosaic.py: did not find fields {args.fields} in file {temp.filename}, exiting")
        return
    # if time is not specified, squeeze extra dimenstions out of inputs
    use_time=False
    for field in ['time','t']:
        if hasattr(temp, field) and getattr(temp, field) is not None:
            use_time=True
    if not use_time:
        for field in temp.fields:
            setattr(temp, field, np.squeeze(getattr(temp, field)))
    if args.fields is None:
        args.fields=temp.fields
    these_fields=[field for field in args.fields if field in temp.fields]
    field_dims={field:getattr(temp, field).ndim for field in these_fields}
    for field in these_fields:
        mosaic.assign({field:np.zeros(mosaic.dimensions[0:field_dims[field]])})
    mosaic.invalid = np.ones(mosaic.dimensions,dtype=bool)


    # check if using a weighted summation scheme for calculating mosaic
    if args.weight:
        mosaic.weight = np.zeros((mosaic.dimensions[0],mosaic.dimensions[1]))
        # create the weights for a single file
        # as the tiles have the same dimensions, we only have to calculate
        # the first set.  After that we can just copy the first
        temp.weights(pad=args.pad, feather=args.feather, apply=False)
        # allocate for output weight matrix
        mosaic.weight=np.zeros((mosaic.dimensions[0],mosaic.dimensions[1]))
        tile_weight=temp.weight.copy()
        # for each file in the list
        for file in file_list:
            # read data grid from HDF5
            temp=pc.grid.mosaic().from_h5(file, group=args.in_group, fields=args.fields)
            if not use_time:
                for field in temp.fields:
                    setattr(temp, field, np.squeeze(getattr(temp, field)))
            these_fields=[field for field in args.fields if field in temp.fields]
            # copy weights for tile
            temp.weight=tile_weight.copy()
            # get the image coordinates of the input file
            iy,ix = mosaic.image_coordinates(temp)
            for field in these_fields:
                try:
                    if field_dims[field]==3:
                        field_data=np.atleast_3d(getattr(temp, field))
                        bands=range(mosaic.dimensions[2])
                        for band in bands:
                            getattr(mosaic, field)[iy,ix,band] += field_data[:,:,band]*temp.weight
                            mosaic.invalid[iy,ix,band] = False
                    else:
                        field_data=getattr(temp, field)
                        getattr(mosaic, field)[iy,ix] += field_data*temp.weight
                        mosaic.invalid[iy,ix,0] = False
                except (IndexError, ValueError) as e:
                    print(f"problem with field {field} in group {args.in_group} in file {file}")
                    raise(e)
            # add weights to total weight matrix
            mosaic.weight[iy,ix] += temp.weight[:,:]
        # find valid weights
        iy,ix = np.nonzero(mosaic.weight == 0)
        mosaic.invalid[iy,ix,:] = True
        # normalize fields by weights
        iy,ix = np.nonzero(mosaic.weight > 0)
        for field in mosaic.fields:
            if field_dims[field]==3:
                for band in range(mosaic.dimensions[2]):
                    getattr(mosaic, field)[iy,ix,band] /= mosaic.weight[iy,ix]
            else:
                getattr(mosaic, field)[iy,ix] /= mosaic.weight[iy,ix]
    else:
        # overwrite the mosaic with each subsequent tile
        # for each file in the list
        for file in file_list:
            # read data grid from HDF5
            temp=pc.grid.mosaic().from_h5(file, group=args.in_group, fields=args.fields)
            these_fields=[field for field in args.fields if field in temp.fields]
            # get the image coordinates of the input file
            iy,ix = mosaic.image_coordinates(temp)
            for field in these_fields:

                try:
                    if field_dims[field]==3:
                        field_data=np.atleast_3d(getattr(temp, field))
                        for band in range(mosaic.dimensions[2]):
                            getattr(mosaic, field)[iy,ix,band] = field_data[:,:,band]
                            mosaic.invalid[iy,ix,band] = False
                    else:
                        field_data=getattr(temp, field)
                        getattr(mosaic, field)[iy,ix] = field_data
                        mosaic.invalid[iy, ix] = False
                except IndexError:
                    print(f"problem with field {field} in group {args.in_group} file {file}")

    # replace invalid points with fill_value
    for field in mosaic.fields:
        if field_dims[field]==3:
            getattr(mosaic, field)[mosaic.invalid] = mosaic.fill_value
        else:
            getattr(mosaic, field)[mosaic.invalid[:,:,0]] = mosaic.fill_value
    # crop mosaic to bounds
    if np.any(args.crop):
        # x and y range (verify min and max order)
        XR = np.sort([args.crop[0],args.crop[1]])
        YR = np.sort([args.crop[2],args.crop[3]])
        mosaic.crop(XR, YR, fields=mosaic.fields)

    if args.time is not None:
        mosaic.t = args.time
    # output each field
    for field in mosaic.fields:
        if field_dims[field] == 2:
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
    main(sys.argv)
