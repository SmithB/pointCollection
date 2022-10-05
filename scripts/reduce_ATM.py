#! /usr/bin/env python


import numpy as np
import pointCollection as pc
import scipy.spatial as sps
import argparse
import os
import sys
def reduce_data(D, stats_radius=50, min_sigma=0.5):
    kdt = sps.KDTree(np.c_[D.x, D.y])

    D_out=pc.data().from_dict({
        'x':D.x,
        'y':D.y,
        'z':D.z,
        'time':D.time
        })
    D_out.assign({field:np.zeros_like(D_out.x)+np.NaN for field in ['bias_50m', 'noise_50m', 'N_50m', 'slope_x', 'slope_y']})

    for i0 in range(D_out.size):
        #if np.mod(i0, 100)==0:
        #    print([i0, D_out.size])
        i1=np.array(kdt.query_ball_point([D.x[i0], D.y[i0]], stats_radius))
        if len(i1) <= 5:
            continue
        b, n, sl, N = pc.calc_bias_and_R(D, i0, i1, stats_radius, min_sigma=min_sigma)
        D_out.bias_50m[i0]=b
        D_out.noise_50m[i0]=n
        D_out.N_50m[i0]=N
        D_out.slope_x[i0]=sl[0]
        D_out.slope_y[i0]=sl[1]
    return D_out


def main():
    parser=argparse.ArgumentParser(description='decimate an ATM file, report statistics of decimated points')
    parser.add_argument('Qfit_file', type=str, help="Qfit file to read")
    parser.add_argument('--decimate_to','-d', type=float, default=10, help="scale to which to decimate the data")
    parser.add_argument('--stats_radius','-s', type=float, default=50, help="scale at which to calculate the statistics")
    parser.add_argument('--min_sigma','-m', type=float, default=0.5, help="minimum sigma for the three-sigma edit step")
    parser.add_argument('--time_format', type=str, default='days_J2k', help='specify one of: days_J2k, matlab, years_CE')
    parser.add_argument('--out_dir', type=str, default='.')
    parser.add_argument('--out_name', type=str)
    args=parser.parse_args()
    
    if args.out_name is None:
        args.out_name=os.path.join(args.out_dir, os.path.basename(args.Qfit_file).replace('.h5','_filt.h5'))
    else:
        if args.out_dir is not None:
            args.out_name=os.path.join(args.out_dir, args.Qfit_file)

    D=pc.ATM_Qfit.data().from_h5(args.Qfit_file)
    
    if np.nanmax(D.latitude > 10):
        EPSG=3413
    else:
        EPSG=3031
        
    D.get_xy(EPSG=EPSG)
    D.assign({'z':D.elevation})
    if args.decimate_to > 0:
        im=pc.pt_blockmedian(D.x, D.y, D.z, args.decimate_to, break_ties=True, index_only=True)
        D.index(im)
        D.assign({'pt_num':im})
    
    if args.time_format=='matlab':
        # year = (ml-730486.)/365.25+2000.
        # ml-730486. = (year-2000)*365.25$
        # ml = days_j2k + 730486.
        D.assign({'time':D.days_J2k + 730486.})
    elif args.time_format == 'years_CE':
        D.assign({'time':D.days_J2k/365.25+2000})
    elif args.time_format == 'days_J2k':
        D.assign({'time':D.days_J2k})
    else:
        print(f'time_format: {args.time_format} not understood')
        raise(ValueError)
    D_out=reduce_data(D, stats_radius=args.stats_radius, min_sigma=args.min_sigma)
    
    if os.path.isfile(args.out_name):
        os.remove(args.out_name)
    
    pc.indexedH5.data().to_file(D_out,args.out_name)

if __name__=='__main__':
    main()

