#! /usr/bin/env python

import h5py
import sys


with h5py.File(sys.argv[1],'r+') as h5f:
    print(f"old dir_root attribute was {h5f['index'].attrs['dir_root']}")
    h5f['index'].attrs['dir_root']=''

