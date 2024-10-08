#! /usr/bin/env python3

import h5py
import sys
import glob
import re

main():

    file_re=re.compile('file_\d+')

    for file in sys.argv[1:]:
        print(file)
        with h5py.File(file, 'r+') as h5f:
            if h5f['index'].attrs['dir_root'][-1] != '/':
    #            print('adding a slash')
                h5f['index'].attrs['dir_root'] += '/'
            for key, item in h5f['index'].attrs.items():
                if file_re.match(key):
    #                print([key, item])
                    if item[0] =='/':
                        temp=item
                        while temp[0] == '/':
                            temp=temp[1:]
                        h5f['index'].attrs[key]=temp

if __name__=='__main__':
    main()
