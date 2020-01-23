#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 10:24:00 2020

@author: ben
"""

import pointCollection as pc
import os
import argparse
import sys
import numpy as np
import json

def make_tile(index_file, xy0, tile_W=2.e5, srs_proj4=None, \
              file_type=None, verbose=False, out_dir='.', field_dict=None):
    '''
    Make a tile for a location based on a geoIndex

    Parameters
    ----------
    index_file : string
        Index file for which to make the tile
    xy0 : list of 2 floats, required
        Tile center. 
    tile_W : numeric, optional
        width of the tile. The default is 2.e5.
    srs_proj4 : TYPE, optional
        DESCRIPTION. The default is None.
    file_type : TYPE, optional
        DESCRIPTION. The default is None.
    verbose : TYPE, optional
        DESCRIPTION. The default is False.
    out_dir : TYPE, optional
        DESCRIPTION. The default is '.'.

    Returns
    -------
    None.

    '''
    
    tile=getattr(pc, file_type).tile
    
    fields=[]
    for key in field_dict:
        fields += field_dict[key]
    
    tile(tile_W=tile_W, SRS_proj4=srs_proj4, xy0=xy0)\
        .from_geoIndex(GI_file=index_file, field_dict=field_dict)\
        .write(out_dir, fields=fields)

def make_queue(index_file, queue_file, tile_W=2.e5, hemisphere=-1, \
               file_type='data', verbose=False, out_dir='.', field_dict_json=None):
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    xy0=pc.geoIndex().from_file(index_file).bins_as_array()
    b0=np.unique(np.round((xy0[0]+1j*xy0[1])/tile_W)*tile_W)
    with open(queue_file,'w') as qh:
        for xx in b0:
            this_str=f"{__file__} -i {index_file} -W {tile_W} -t {file_type} -H {hemisphere} --xy {np.real(xx)} {np.imag(xx)} -o {out_dir}"
            if field_dict_json is not None:
                this_str += f" -j {field_dict_json}"
            if verbose: 
                print(this_str)
            qh.write(this_str+'\n')

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--index_file','-i', type=str, help="index file, defaults to GeoIndex.h5 in the dirname of the glob string")
    parser.add_argument('--tile_W', '-W', type=int, default=2.e5,  help='tile width, default=200000')
    parser.add_argument('--hemisphere','-H', type=int, required=True, help='hemisphere, must be 1 or -1, required')
    parser.add_argument('--verbose','-v', action='store_true')
    parser.add_argument('--type','-t', type=str, required=True, help="file type.  See datatypes with 'index' methods in pointCollection")
    parser.add_argument('--xy', type=float, nargs=2, default=[None, None], help='tile center to be generated.')
    parser.add_argument('--out_dir','-o', type=str, default='.', help='destination directory for tile files')
    parser.add_argument('--queue_file', '-q', type=str, default=None, help='file in which to write queue of commands.  Either -xy or --queue must be specified')
    parser.add_argument('--field_dict_json','-j', type=str, default=None, help='json file from which to read the field dictionary')
    args=parser.parse_args()
    
    if args.hemisphere==1:
        srs_proj4 = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs '
    elif args.hemisphere==-1:
        srs_proj4 = '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs '
    
    if args.field_dict_json is not None:
        with open(args.field_dict_json,'r') as jh:
            field_dict=json.load(jh)
            for key in list(field_dict):
                if key == "null":
                    temp=field_dict.pop(key)
                    field_dict[None]=temp
    else:
        field_dict=None
    
    if args.xy[0] is not None:
        if args.queue_file is not None:
            print("must specify one of --queue_file or --xy, not both.")
            sys.exit(1)
        make_tile(args.index_file, xy0=args.xy, tile_W=args.tile_W, 
                  srs_proj4=srs_proj4, file_type=args.type, out_dir=args.out_dir,
                  field_dict=field_dict, verbose=args.verbose)
    if args.queue_file is not None:
        make_queue(args.index_file, args.queue_file, tile_W=args.tile_W, 
                   hemisphere=args.hemisphere, file_type=args.type, 
                   out_dir=args.out_dir, field_dict_json=args.field_dict_json,
                   verbose=args.verbose)

if __name__=='__main__':
    main()
