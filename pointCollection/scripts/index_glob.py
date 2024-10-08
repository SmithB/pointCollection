#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 09:42:00 2019

@author: ben
"""

import pointCollection as pc

import os
import glob
import numpy as np
import matplotlib.pyplot as plt

import argparse


def index_for_glob(glob_string, dir_root=None, index_file=None, file_type=None, group=None,\
                   verbose=False, delta=[1.e4, 1.e4], SRS_proj4=None, relative=False):
    """
    make a geoindex for the files in a glob string
    
    arguments:
        glob_string: a string that when passed to glob will return a list of files
        dir_root: directory name in which to root the index
        index_file: the name of the output index file (defaults to basename(glob_string)/GeoIndex.h5)
        verbose: should a lot be echoed?  defaults to None
    """

    if dir_root is not None and dir_root[-1] != '/':
        dir_root += '/'

    files=[]
    for thestr in glob_string.split(' '):
        files += glob.glob(thestr)

    if index_file is None:
        index_file=os.path.join(os.path.dirname(glob_string),'GeoIndex.h5')

    index_list=[];
    for file in files:
        if verbose:
            print(f"index_glob: indexing {file}")
        #if os.path.basename(file)[0:3]=='Geo':
        #    continue
        try:
            index_list += [pc.geoIndex(delta=delta, SRS_proj4=SRS_proj4).for_file(file, file_type, group=group) ]
        except Exception as e:
            print(f"index_glob: Exception thrown for file {file}:")
            print(e)
    if verbose:
        print(f"index_glob: making index file: {index_file}")
    gI=pc.geoIndex(delta=delta, SRS_proj4=SRS_proj4).from_list(index_list, \
               dir_root=dir_root)
    if relative:
        gI.attrs['dir_root']=None
    gI.to_file(index_file)
    
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--glob_string','-g', type=str, required=True, help="quoted string to pass to glob to find the files to index")
    parser.add_argument('--type','-t', type=str, required=True, help="file type.  See pointCollection.geoIndex for options")    
    parser.add_argument('--group', type=str, required=False, help="group within the file to read.  Defaults to None")
    parser.add_argument('--index_file','-i', type=str, help="index file, defaults to GeoIndex.h5 in the dirname of the glob string")
    parser.add_argument('--dir_root','-r', type=str, help="directory root for the index.  Defaults to None (absolute paths)")
    parser.add_argument('--bin_size', '-b', type=int, default=1.e4,  help='index bin size, default=1.e4')
    parser.add_argument('--hemisphere','-H', type=int, help='hemisphere, must be 1 or -1,')
    parser.add_argument('--proj4', type=str, help='proj4 to be used for the tiles')
    parser.add_argument('--verbose','-v', action='store_true')
    parser.add_argument('--Relative','-R', action='store_true')
    args=parser.parse_args()
    
    if args.hemisphere==1:
        srs_proj4 = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs '
    elif args.hemisphere==-1:
        srs_proj4 = '+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs '
    if args.proj4 is not None:
        srs_proj4 = args.proj4
    
    index_for_glob(args.glob_string, dir_root=args.dir_root, \
                   index_file=args.index_file, \
                   file_type=args.type,\
                   verbose=args.verbose,\
                   delta=[np.float64(args.bin_size), np.float64(args.bin_size)], \
                   SRS_proj4=srs_proj4,\
                   group=args.group,\
                   relative=args.Relative)

if __name__=='__main__':
    main()

# -H 1 -g "/Volumes/ice2/ben/scf/GL_11/U01/*.h5" --index_file /Volumes/ice2/ben/scf/GL_11/U01/GeoIndex.h5 --dir_root /Volumes/ice2/ben/scf/GL_11/U01 -t "h5_geoindex"
    
    
   #-H 1 -g "/Volumes/ice2/ben/scf/GL_11/U01*/ATL11*.h5" --index_file /Volumes/ice2/ben/scf/GL_11/GeoIndex_U01.h5 --dir_root /Volumes/ice2/ben/scf/GL_11/ -t "h5_geoindex" 
