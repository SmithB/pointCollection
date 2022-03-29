#! /usr/bin/env python

import h5py
import sys


with h5py.File(sys.argv[1],'r+') as h5f:
    print(sys.argv[1])
    if 'dir_root' in h5f['index'].attrs:
        print(f"\told dir_root attribute was {h5f['index'].attrs['dir_root']}")
        h5f['index'].attrs['dir_root']=''
    else:
        print("\tindex has no dir_root attribute")
