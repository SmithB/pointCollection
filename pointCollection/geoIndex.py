# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 20:58:19 2018

@author: ben

This is a class that lets us generate a coarse-resolution database of the
data-point locations in a point-data file, and to build hierarchical indices
to allow efficient searches for data.
Coordinates below are always provided as tuples or lists, the first member of which is an array of x coordinates, the second is an array of y coordinates
"""
import numpy as np
import re
#from osgeo import osr
#import matplotlib.pyplot as plt
import pointCollection as pc
import os
from warnings import warn

class geoIndex(dict):
    def __init__(self, delta=[1000,1000], SRS_proj4=None, data=None, DEBUG=False):
        dict.__init__(self)
        self.attrs={'delta':delta,'SRS_proj4':SRS_proj4, 'n_files':0, 'dir_root':''}
        self.data=data
        if self.data is not None:
            if hasattr(data,'x'):
                self.from_xy([data.x, data.y])
            elif hasattr(data,'latitude'):
                self.from_latlon(self.data.latitude, self.data.longitude)
        self.h5_file=None
        self.filename=None
        self.DEBUG=DEBUG

    def __repr__(self):
        out = f"{self.__class__} with {len(self.keys())} bins, referencing {self.attrs['n_files']} files"
        return out

    def __copy__(self):
        """
        copy method,
        """
        out=geoIndex()
        for attr in self.attrs:
            out.attrs[attr]=self.attrs[attr]
        for field in self.keys:
            out[field]=self[field].copy()
        out.data=self.data
        return out

    def copy_subset(self, xyBin=None, pad=None):
        """
        copy method, may specify which bins to copy
        """
        out=geoIndex()
        out.filename=self.filename
        for attr in self.attrs:
            out.attrs[attr]=self.attrs[attr]
        #out.attrs=self.attrs.copy()
        if xyBin is None:
            these_keys=self.keys()
        else:
            if pad is not None and pad >0:
                xyBin=pad_bins(xyBin, pad, self.attrs['delta'])
            these_keys=self.keys_from_xy(xyBin)
        for field in these_keys:
            if isinstance(self[field], dict):
                out[field]=self[field].copy()
            else:
                if field not in out:
                    out[field]=dict()
                # this copies the reference to a h5 group
                for key in self[field].keys():
                    out[field][key]=self[field][key][:]
        out.data=self.data
        return out

    def from_xy(self, xy,  filename=None, file_type=None, number=0, bin_function=np.round, fake_offset_val=None, first_last=None):
        """
        build a geoIndex from a list of x, y points, for a specified filename
        and file_type.  If the file_type is 'geoIndex', optionally specify a
        value for 'fake_offset_val'
        """
        delta=self.attrs['delta']
        self.filename=filename
        xy_bin = bin_function(np.c_[xy[0].ravel(), xy[1].ravel()]/delta).astype(int)
        if first_last is None:
            # If the inputs haven't provided the first and last index for each bin, need to calculate it:
            # sort by magnitude of (x,y), then by angle
            ordering=np.sqrt(np.sum(xy_bin**2, axis=1))+(np.arctan2(xy_bin[:,0], xy_bin[:,1])+np.pi)/2/np.pi
            uOrd, first=np.unique(ordering, return_index=True)
            uOrd, temp=np.unique(-ordering[::-1], return_index=True)
            last=len(ordering)-1-temp[::-1]
            keys=['%d_%d' % (delta[0]*xy_bin[first[ind],0], delta[1]*xy_bin[first[ind],1]) for ind  in range(len(first))]
        else:
            # assume that the first_last values match the bins
            first, last=first_last
            keys=['%d_%d' % (xy_bin[ind,0]*delta[0], xy_bin[ind,1]*delta[1]) for ind in range(first.size)]

        for ind, key in enumerate(keys):
            if fake_offset_val is None:
                self[key] = {'file_num':np.array(int(number), ndmin=1), 'offset_start':np.array(first[ind], ndmin=1), 'offset_end':np.array(last[ind], ndmin=1)}
            else:
                self[key] = {'file_num':np.array(int(number), ndmin=1), 'offset_start':np.array(fake_offset_val, ndmin=1), 'offset_end':np.array(fake_offset_val, ndmin=1)}
        #In some cases the files are predefined.  If this is not the case, use the current filename
        if 'file_0' not in self.attrs:
            self.attrs['file_0']=filename
            self.attrs['type_0']=file_type
            self.attrs['n_files']=1
        return self

    def from_latlon(self, lat, lon,  filename=None, file_type=None, number=0, fake_offset_val=None):

        temp=pc.data().from_dict({'latitude':lat,'longitude':lon}).get_xy(SRS_proj4=self.attrs['SRS_proj4'])
        x=temp.x
        y=temp.y
        #out_srs=osr.SpatialReference()
        #out_srs.ImportFromProj4(self.attrs['SRS_proj4'])
        #ll_srs=osr.SpatialReference()
        #ll_srs.ImportFromEPSG(4326)
        #if hasattr(osr,'OAMS_TRADITIONAL_GIS_ORDER'):
        #    ll_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        #ct=osr.CoordinateTransformation(ll_srs, out_srs).TransformPoint
        #xy=[ct(*xyz)[0:2] for xyz in zip(np.ravel(lon), np.ravel(lat), np.zeros_like(lat).ravel())]
        #x, y, z = list(zip(*[ct(*xy) for xy in zip(np.ravel(lon), np.ravel(lat), np.zeros_like(lat).ravel())]))
        return self.from_xy([np.array(x),np.array(y)], filename, file_type, number, fake_offset_val)

    def from_list(self, index_list, dir_root=''):
        """
        build a geoIndex from a list of geo_indices.
        Each bin in the resulting geoIndex contains information for reading
        the files indexed by the geo_indices in index_list
        """

        dir_root=strip_double_slashes(dir_root)
        if len(index_list)==0:
            return
        for key in ['dir_root', 'SRS_proj4']:
            if key in index_list[0].attrs:
                if index_list[0].attrs[key] is not None:
                    self.attrs[key]=index_list[0].attrs[key]
        if dir_root is not None and len(dir_root) > 0:
            self.attrs['dir_root']=dir_root
        # make a list of files in the destination index (self)
        fileListTo=list()

        for index in index_list:
            # check if a particular input file is alread in the output index, otherwise add it
            # keep track of how the filenames in the input index correspond to those in the output index
            num_out=dict()
            alreadyIn=list()
            for fileNum in range(index.attrs['n_files']):
                thisFileName=index.attrs['file_%d' % fileNum]
                if 'dir_root' in index.attrs and index.attrs['dir_root'] is not None:
                    thisFileName=os.path.join(index.attrs['dir_root'],thisFileName)
                if dir_root is not None:
                    thisFileName = os.path.relpath(thisFileName, dir_root)
                thisFileType=index.attrs['type_%d' % fileNum]
                if thisFileName not in fileListTo:
                    fileListTo.append(thisFileName)
                    self.attrs['file_%d' % (len(fileListTo)-1)] = thisFileName
                    self.attrs['type_%d' % (len(fileListTo)-1)] = thisFileType
                else:
                    alreadyIn.append(fileNum)
                num_out[fileNum]=fileListTo.index(thisFileName)
            # loop bver the bins in the current index
            for bin in index.keys():
                # If a particular filename is already in fileListTo, its corresponding
                # number is in alreadIn, and we'll skip it for this bin so we don't
                # end up with duplicate data
                newFileNums=index[bin]['file_num'].copy()
                keep=np.logical_not(np.isin(newFileNums, alreadyIn))
                if not np.any(keep):
                    continue
                newFileNums=newFileNums[keep]
                for row in range(newFileNums.shape[0]):
                    # translate the newFileNums to the file numbers for the output index
                    newFileNums[row]=num_out[newFileNums[row]]
                # if the bin is alreay in self, copy the infomation to it
                if bin in self:
                    append_data(self[bin],'file_num', newFileNums)
                    for field in ('offset_start','offset_end'):
                        append_data(self[bin], field, index[bin][field][keep])
                # Otherwise, make a new bin in self
                else:
                    self[bin]=dict()
                    self[bin]['file_num']=newFileNums
                    for field in ('offset_start','offset_end'):
                        self[bin][field]=index[bin][field][keep]
        self.attrs['n_files']=len(fileListTo)
        return self

    def from_file(self, index_file, read_file=False, group='index'):
        """
        read geoIndex info from file 'index_file.'
        If read_file is set to False, the file is not read, but the
        h5_file_index attribute of the resulting geoIndex is set to a
        reference to the hdf_file's 'index' attribute.  This seems to be
        faster than reading the whole file.
        """

        import h5py
        h5_f = h5py.File(os.path.expanduser(index_file),'r')
        h5_i = h5_f[group]
        if read_file:
            for bin in h5_i.keys():
                self[bin]=h5_i[bin]
        self.attrs=h5_i.attrs
        self.h5_file=h5_f
        self.h5_file_index=h5_f['index']
        self.filename=index_file
        return self

    def change_root(self, new_root, old_root=None):
        """
        changes the root path to a new path
        """
        if old_root is None:
            if self.attrs['dir_root'] is not None:
                old_root = os.path.normpath(self.attrs['dir_root'])
            else:
                old_root = ''
        new_root = os.path.normpath(new_root)
        file_re = re.compile(r'file_\d+')
        for key in self.attrs.keys():
            if file_re.match(key) is not None:
                temp = os.path.join(old_root, self.attrs[key])
                self.attrs[key] = os.path.relpath(temp, new_root)
        self.attrs['dir_root'] = new_root
        return self

    def to_file(self, filename):
        """
        write the current geoindex to h5 file 'filename'
        """
        import h5py
        indexF = h5py.File(os.path.expanduser(filename),'a', libver='latest')
        if 'index' in indexF:
            del indexF['index']
        indexGrp=indexF.create_group('index')
        if 'n_files' in self.attrs:
            indexGrp.attrs['n_files'] = self.attrs['n_files']
        else:
            indexGrp.attrs['n_files']=0
        if 'dir_root' in self.attrs and self.attrs['dir_root'] is not None:
            indexGrp.attrs['dir_root']=self.attrs['dir_root']
        indexGrp.attrs['delta'] = self.attrs['delta']
        if 'SRS_proj4' in self.attrs and self.attrs['SRS_proj4'] is not None:
            indexGrp.attrs['SRS_proj4'] = self.attrs['SRS_proj4']
        for key in self.keys():
            indexGrp.create_group(key)
            for field in ['file_num','offset_start','offset_end']:
                indexGrp[key].create_dataset(field,data=self[key][field])
        for ii in range(self.attrs['n_files']):
            this_key='file_%d' % ii
            indexGrp.attrs[this_key]=self.attrs[this_key]
            this_type='type_%d' % ii
            indexGrp.attrs[this_type]=self.attrs[this_type]
        indexF.close()
        return

    def for_file(self, filename, file_type, number=0, dir_root='', group=None,
                 self_contained=False):
        """
        make a geoIndex for file 'filename'

        Parameters
        ----------
        filename : string
            the file to index. For file_type='h5', reading goes through
            pc.data.from_h5(), which opens files with h5py -- this works for
            both '.h5' files and netCDF4-format '.nc' files (netCDF4 is an
            HDF5 container format), but not for classic/netCDF3 '.nc' files.
        file_type : string
            the type of file being indexed (e.g. 'h5', 'ATL06', 'ATL11', ...)
        number : int, optional
            the file number to assign this source within the index.
        dir_root : string, optional
            a directory prefix common to indexed files, stripped from the
            stored filename and re-applied at query time (see resolve_path()).
        group : string, optional
            for file_type='h5', the group within the file containing the
            'x' and 'y' fields to index.
        self_contained : bool, optional
            for file_type='h5' only. If True, don't record `filename` as a
            separate source; instead mark this source as living in group
            `group` inside whatever file this geoIndex is itself eventually
            written to via to_file(). This lets a single file hold both a
            pc.data object (in group `group`) and the geoIndex for it (in
            the 'index' group). Requires `group` to be set. Build order
            matters: write the data first (e.g.
            `D.to_h5(path, group=group)` -- its default `replace=True` would
            wipe a previously-written index if called second), then save
            the index into the *same* path with
            `pc.geoIndex(...).for_file(path, 'h5', group=group, self_contained=True).to_file(path)`.
            Sources built this way are meant to be read back directly via
            `from_file(path)`/`query_xy(...)`; avoid merging them with
            `for_files()`/`from_list()` into an index saved to a *different*
            path, since the embedded reference is tied to whatever
            `self.filename` is at query time.
        """
        if self_contained and file_type != 'h5':
            raise ValueError("for_file: self_contained=True is only supported for file_type='h5'.")
        dir_root=strip_double_slashes(dir_root)
        self.filename=filename
        if dir_root is not None:
            # eliminate the string in 'dir_root' from the filename
            filename_out=strip_double_slashes(filename).replace(dir_root,'')
        if file_type in ['ATL06']:
            temp=list()
            this_field_dict={None:('latitude','longitude','h_li','delta_time')}
            for beam_pair in (1, 2, 3):
                D=pc.ATL06.data(beam_pair=beam_pair, field_dict=this_field_dict).from_h5(filename)
                D.get_xy(self.attrs['SRS_proj4'])
                if D.latitude.shape[0] > 0:
                    temp.append(geoIndex(delta=self.attrs['delta'], SRS_proj4=\
                                          self.attrs['SRS_proj4']).\
                                from_xy([np.nanmean(D.x, axis=1), np.nanmean(D.y, axis=1)],
                                        '%s:pair%d' % (filename_out, beam_pair), 'ATL06', number=number))
            self.from_list(temp, dir_root=dir_root)
        if file_type in ['ATL11','ATL11_xo']:
            temp=list()
            for beam_pair in (1, 2, 3):
                if file_type=='ATL11':
                    field_dict={f'pt{beam_pair}':['latitude','longitude']}
                else:
                    field_dict={f'pt{beam_pair}/crossing_track_data':['latitude','longitude']}
                try:
                    D=pc.data().from_h5(filename, field_dict=field_dict)
                    D.get_xy(self.attrs['SRS_proj4'])
                    if D.x.shape[0] > 0:
                        temp.append(geoIndex(delta=self.attrs['delta'], \
                                          SRS_proj4=self.attrs['SRS_proj4']).from_xy([D.x, D.y], '%s:pair%d' % (filename_out, beam_pair), file_type, number=number))
                except Exception as e:
                    if self.DEBUG:
                        raise(e)
                    pass
            self.from_list(temp)
        if file_type in ['h5']:
            if self_contained and not group:
                raise ValueError("for_file: self_contained=True requires 'group' to be set.")
            D=pc.data().from_h5(filename, field_dict={group:['x','y']})
            if D.x.size > 0:
                index_filename = (':' + group) if self_contained else filename_out
                self.from_xy((D.x, D.y), filename=index_filename, file_type='h5', number=number)
        if file_type in ['ATM_Qfit']:
            D=pc.ATM_Qfit.data().from_h5(filename)
            if D.latitude.shape[0] > 0:
                self.from_latlon(D.latitude, D.longitude,  filename_out, 'ATM_Qfit', number=number)
        if file_type in ['ATM_waveform']:
            D=pc.ATMwaveform.data().from_h5(filename)
            if D.latitude.shape[0] > 0:
                self.from_latlon(D.latitude, D.longitude,  filename_out, 'ATM_waveform', number=number)
        if file_type in ['glah12']:
            if int(re.compile(r'lat_0=(\S+)').search(self.SRS_proj4).group(1))<0:
                D=pc.glah12.data().from_h5(filename, lat_range=[-90, -60])
            else:
                D=pc.glah12.data().from_h5(filename, lat_range=[60, 90])
            self.from_latlon(D.latitude, D.longitude, filename_out, 'glah12', number=number)
        if file_type in ['filtered_DEM', 'DEM', 'geotif'] :
            D=pc.grid.data().from_geotif(filename, bands=[1], min_res=self.attrs['delta'][0]/10).as_points()
            if D.size > 0:
                self.from_xy((D.x, D.y), filename=filename_out, file_type=file_type, number=number)
        if file_type in ['h5_geoindex']:
            # read the file as a collection of points
            temp_GI=geoIndex().from_file(filename)
            xy_bin=temp_GI.bins_as_array()
            # loop over a minimal set of attributes:
            for attr in ['delta','SRS_proj4','dir_root']:
                if attr in temp_GI.attrs:
                    self.attrs[attr]=temp_GI.attrs[attr]
            self.attrs['file_%d' % number] = filename_out
            self.attrs['type_%d' % number] = file_type
            if dir_root is not None:
                self.attrs['dir_root']=dir_root
            self.attrs['n_files']=1
            self.from_xy(xy_bin, filename=filename_out, file_type=file_type, number=number, fake_offset_val=-1)
        if file_type in ['indexed_h5']:
            import h5py
            h5f=h5py.File(filename,'r')
            if 'INDEX' in h5f:
                xy=[np.array(h5f['INDEX']['bin_x']), np.array(h5f['INDEX']['bin_y'])]
                if 'bin_index' in h5f['INDEX']:
                    # this is the type of indexed h5 that has all of the data in single datasets
                    i0_i1=h5f['INDEX']['bin_index']
                    first_last=[i0_i1[0,:].ravel(), i0_i1[1,:].ravel()]
                    fake_offset=None
                else:
                    first_last=None
                    fake_offset=-1
            else:
                # there is no index-- just a bunch of bins, maybe?
                first_last=None
                fake_offset=-1
                bin_re=re.compile(r"(.*)E_(.*)N");
                xy=[[], []]
                for key in h5f:
                    m=bin_re.match(key)
                    if m is None:
                        continue
                    xy[0].append(float(m.group(1)))
                    xy[1].append(float(m.group(2)))
                xy[0]=np.array(xy[0])
                xy[1]=np.array(xy[1])
            self.from_xy(xy, filename=filename_out, file_type=file_type, number=number, first_last=first_last, fake_offset_val=fake_offset)
            if dir_root is not None:
                self.attrs['dir_root']=dir_root
            h5f.close()
        if file_type in ['indexed_h5_from_matlab']:
            import h5py
            h5f=h5py.File(filename,'r')
            xy=[np.array(h5f['INDEX']['bin_x']), np.array(h5f['INDEX']['bin_y'])]
            first_last=None
            fake_offset=-1
            self.from_xy(xy, filename_out, file_type, number=number, first_last=first_last, fake_offset_val=fake_offset)
            h5f.close()
        return self

    def for_files(self, filename_list, file_type, SRS_proj4=None, dir_root=''):
        index_list=list()
        for filename in filename_list:
            index_list.append(geoIndex(SRS_proj4=SRS_proj4, delta=self.attrs['delta']).for_file(filename, file_type, dir_root=dir_root, number=0))
        self.SRS_proj4=SRS_proj4
        return self.from_list(index_list, dir_root=dir_root)

    def query_latlon(self, lat, lon, get_data=True, fields=None, error_action='warn'):
        """
        query the current geoIndex for all bins that match the bin locations
        provided in (lat, lon),  Optionally return data, with field query in 'fields'
        """
        #out_srs=osr.SpatialReference()
        #out_srs.ImportFromProj4(self.attribs['SRS_proj4'])
        #ll_srs=osr.SpatialReference()
        #ll_srs.ImportFromEPSG(4326)
        #if hasattr(osr,'OAMS_TRADITIONAL_GIS_ORDER'):
        #    ll_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        #ct=osr.CoordinateTransformation(ll_srs, out_srs)
        #x, y = list(zip(*[ct.TransformPoint(xy) for xy in zip(np.ravel(lon), np.ravel(lat))]))
        temp=pc.data().from_dict({'latitude':lat,'longitude':lon}).get_xy(SRS_proj4=self.attrs['SRS_proj4'])
        x=temp.x
        y=temp.y
        delta=self.attrs['delta']
        xb=np.round(x/delta[0])*delta[0]
        yb=np.round(y/delta[1])*delta[1]
        return self.query_xy([xb, yb], get_data=get_data, fields=fields, error_action=error_action)

    def query_xy_box(self, xr, yr, get_data=True, fields=None, dir_root='',
                     full_path=False, error_action='warn', remote_file=None, fs=None,
                     trim_last_point=False):
        """
        query the current geoIndex for all bins in the box specified by box [xr,yr]

        remote_file, fs, and trim_last_point are passed through to query_xy();
        see its docstring.
        """
        xy_bin=self.bins_as_array()
        these=(xy_bin[0] >= xr[0]) & (xy_bin[0] <= xr[1]) &\
            (xy_bin[1] >= yr[0]) & (xy_bin[1] <= yr[1])
        return self.query_xy([xy_bin[0][these], xy_bin[1][these]], get_data=get_data, \
                             fields=fields,
                             dir_root=dir_root,
                             bounds=[xr, yr],
                             full_path = full_path,
                             error_action = error_action,
                             remote_file = remote_file,
                             fs = fs,
                             trim_last_point = trim_last_point)

    def intersect(self, other, pad=[0, 0]):
        """
        given a pair of geo indexes return the subsets that are common between the two, optionally padding one or both
        """
        bins_both=set(self.keys()).intersection(other.keys())
        xyB=np.c_[[np.fromstring(key, sep='_') for key in bins_both]]
        if xyB.size==0:
            return None, None
        self_sub=self.copy_subset(xyBin=[xyB[:,0], xyB[:,1]], pad=pad[0])
        other_sub=other.copy_subset(xyBin=[xyB[:,0], xyB[:,1]], pad=pad[1])
        return self_sub, other_sub

    def query_xy(self, xyb, cleanup = True,
                 get_data = True,
                 full_path = True,
                 fields=None,
                 pad=None,
                 dir_root='',
                 strict=False,
                 bounds=None,
                 error_action='warn',
                 remote_file=None,
                 fs=None,
                 trim_last_point=False):
        """
        check if data exist within the current geo index for bins in lists/arrays
            xb and yb.
        If argument delta is provided, find the bins in the current geoIndex
            that round to (xb, yb)
        If 'delta' is provided, read the underlying data sources, possibly recursively
            otherwise return a query_result: a dict with one entry for each source file
            in the current geoIndex, giving the bin locations provided by that file,
            and the offsets in the file corresponding to each.
        If 'pad' is provided, include bins between xb-pad*delta and xp+pad*delta (inclusive)
            in the query (likewise for y)
        If 'remote_file' is provided, it overrides the resolved file identity (e.g. with
            a real s3:// URI) for every bin in this query, preserving any ':pairN' suffix.
            This only makes sense against an index that indexes a single physical file
            (e.g. a per-granule ATL11 index) -- it is not meant for indices spanning
            multiple distinct source files.
        trim_last_point is passed through to get_data(); see its docstring.
        """
        delta=self.attrs['delta']
        if isinstance(xyb[0], np.ndarray):
            xyb=[xyb[0].copy().ravel(), xyb[1].copy().ravel()]
        if pad is not None:
            xyb=pad_bins(xyb, pad, delta)
        if isinstance(xyb[0], float) or isinstance(xyb[0], int):
            # if scalars were provided, keep the 'zip' from choking by making them iterable
            xyb=[np.array(xyb[0].copy()).reshape([1]), np.array(xyb[1].copy()).reshape([1])]
        # round the input bins to the bin resolution
        for ii in [0, 1]:
            xyb[ii]=np.round(xyb[ii]/self.attrs['delta'][ii])*self.attrs['delta'][ii]
        # make a temporary geoIndex to hold the subset of the current geoindex
        # corresponding to xb and yb
        temp_gi=geoIndex(delta=self.attrs['delta'], SRS_proj4=self.attrs.get('SRS_proj4'))
        for bin in set(zip(xyb[0], xyb[1])):
           bin_name='%d_%d' % bin
           if bin_name in self:
               temp_gi[bin_name]=self[bin_name]
           elif hasattr(self, 'h5_file_index') and bin_name in self.h5_file_index:
               temp_gi[bin_name]=self.h5_file_index[bin_name]
        if len(temp_gi.keys())==0:
            return None
        temp_dict=dict()
        for field in ['file_num','offset_start','offset_end']:
           temp_dict[field]=np.concatenate([temp_gi[key][field] for key in sorted(temp_gi)])
        # build an array of x and y values for the bins in temp_gi
        xy0=np.concatenate([np.tile(np.fromstring(key, sep='_').astype(int),(temp_gi[key]['file_num'].size,1)) for key in sorted(temp_gi)], axis=0)
        out_file_nums=np.unique(temp_dict['file_num'])
        query_results=dict()
        for out_file_num in out_file_nums:
            these=temp_dict['file_num']==out_file_num
            i0=np.array(temp_dict['offset_start'][these], dtype=int)
            i1=np.array(temp_dict['offset_end'][these], dtype=int)
            xy=xy0[these,:]
            if cleanup:
                # clean up the output: when the start of the next segment is
                #before or adjacent to the end of the previous, stick them together
                ii=np.argsort(i0)
                i0=i0[ii]
                i1=i1[ii]
                xy=xy[ii,:]
                keep=np.zeros(len(i0), dtype=bool)
                this=0
                keep[this]=True
                for kk in np.arange(1,len(i0)):
                    if i0[kk]<=i1[this]+1 and i0[kk] > 0 and i0[kk] > 0:
                        keep[kk]=False
                        i1[this]=np.maximum(i1[this], i1[kk])
                    else:
                        this=kk
                        keep[kk]=True
                i0=i0[keep]
                i1=i1[keep]
                xy=xy[keep,:]
            # if the file_N attribute begins with ':', it's a group in the current
            # file (built by for_file(..., self_contained=True)); self.filename
            # already refers to it exactly as opened, so no directory-joining
            # is needed (see resolve_path()'s self-referential-filename guard)
            this_query_file = self.attrs['file_%d' % out_file_num]
            if this_query_file is not None and this_query_file[0] == ':':
                this_query_file = self.filename + this_query_file
            elif full_path:
                this_query_file = self.resolve_path(this_query_file, dir_root)
            if remote_file is not None:
                suffix = ''
                if this_query_file is not None and ':' in this_query_file:
                    suffix = ':' + this_query_file.split(':', 1)[1]
                this_query_file = remote_file + suffix
            query_results[this_query_file]={
            'type':self.attrs['type_%d' % out_file_num],
            'offset_start':i0,
            'offset_end':i1,
            'x':xy[:,0],
            'y':xy[:,1]}
        if get_data:
            query_results=self.get_data(query_results, fields=fields, dir_root=dir_root, bounds=bounds, error_action=error_action, already_resolved=full_path, fs=fs, trim_last_point=trim_last_point)
            if strict is True:
                # take the subset of data that rounds exactly to the query (OTW, may get data that extend outside)
                if not isinstance(query_results, list):
                    query_results=[query_results]
                for item in query_results:
                    if not hasattr(item,'x'):
                        item.get_xy(self.attrs['SRS_proj4'])
                    keep=np.zeros_like(item.x, dtype=bool)
                    xr=np.round(item.x/delta[0])*delta[0]
                    yr=np.round(item.y/delta[0])*delta[0]
                    for xybi in zip(xyb[0], xyb[1]):
                        ii=(xr==xybi[0]) & (yr==xybi[1])
                        keep[ii]=True
                    item.index(keep)
        return query_results

    def resolve_path(self, filename, dir_root=None):
        """
        figure out where to find a file based on a query result

        Parameters
        ----------
        filename : string
            a filename provided by an index
        dir_root : string or None
            a directory that can be prepended to subdirectories to help find files

        Returns
        -------
        string
            absolute path for the file to read
        """
        if pc.io_utils.is_remote_path(filename):
            return filename
        if dir_root is None:
            dir_root=''
        self_dir_root=''
        if 'dir_root' in self.attrs:
            self_dir_root=self.attrs['dir_root']
        # if the filename begins with '/', it is absolute
        if filename is not None and filename[0]==os.path.sep:
            return filename
        # if filename already refers to this index's own file (e.g. built from
        # self.filename for a ':group' self-contained-file entry), it's already
        # fully resolved -- resolving it again would double any relative prefix
        if self.filename is not None and filename is not None and filename.startswith(self.filename):
            return filename
        # if self.attrs['dir_root'] begins with '/', it is absolute, and overrides the dir_root argument
        if len(self_dir_root)>0 and self_dir_root==os.path.sep:
            return os.path.join(self.attrs['dir_root'], filename)
        # if self.attrs['dir_root'] does not begin with a '/', it is relative
        if len(self_dir_root) > 0:
            return os.path.join(dir_root, self.attrs['dir_root'], filename)
        # if dir_root is provided, prepend it to the filename
        if len(dir_root) >0 and dir_root[0]==os.path.sep:
            return os.path.join(dir_root,filename)
        # otherwise, if len(dir_root) is 0 and self.attrs['dir_root'] is None,
        # assume that files are relative to the index path
        if self.filename is not None:
            return os.path.join(os.path.dirname(self.filename), filename)
        # if nothing has happened yet, return the filename
        return filename

    def get_data(self, query_results, fields=None,  data=None, dir_root='',
                 bounds=None, function=None, error_action='warn', already_resolved=False, fs=None,
                 trim_last_point=False):
        """
        read the data from a set of query results
        Currently the function knows how to read:
        h5_geoindex
        indexed h5s
        Qfit data (waveform and plain)
        DEM data (filtered and not)
        ATL06 data.
        Append more cases as needed

        trim_last_point : bool, optional
            offset_start/offset_end (as built by from_xy()) are an *inclusive*
            (first, last) row-index pair, but readers (data.py, ATL06/data.py)
            slice with an *exclusive* stop -- so by default (False) this
            passes index_range=(offset_start, offset_end+1) to readers for
            the 'h5', 'ATL11', 'ATL06', and 'ATM_Qfit' types, so the last row
            of each segment is included. Set True to reproduce the old
            behavior (silently dropping that last row) for legacy
            comparisons. Does not affect 'indexed_h5'/'indexed_h5_from_matlab'
            (their offsets can be a -1 sentinel or come from an externally-
            built, unverified index) or a user-supplied `function` (which
            always receives the raw, unmodified offsets).
        """
        out_data=list()

        # some data types take a dictionary rather than a list of fields
        if isinstance(fields, dict):
            field_dict=fields
            field_list=None
        else:
            field_dict=None
            field_list=fields

        # if we are querying any DEM data, work out the bounds of the query so we don't have to read the whole DEMs
        all_types=[query_results[key]['type'] for key in query_results]
        if 'DEM' in all_types or 'filtered_DEM' in all_types and bounds is None:
            all_x=list()
            all_y=list()
            for key, result in query_results.items():
                all_x += result['x'].tolist()
                all_y += result['y'].tolist()
            delta=self.attrs['delta']
            bounds=[[np.min(all_x)-delta[0]/2, np.max(all_x)+delta[0]/2], \
                    [np.min(all_y)-delta[1]/2, np.max(all_y)+delta[1]/2]]

        # Types that route through pc.data.from_h5()-family readers, which
        # accept an externally-supplied, already-open h5_f handle. For these,
        # group all reads that target the same physical file -- e.g. ATL11's
        # up to 3 beam pairs, or multiple disjoint offset segments -- so a
        # single handle can be opened once, reused for every read against
        # that file, and then closed, before moving on to the next physical
        # file (at most one handle open at a time). Without this, each
        # pair/segment reopened the same remote file independently -- a full
        # extra round trip per read, for a file we'd already opened moments
        # before. Other types (rasters, indexed_h5's external/sentinel
        # offsets, a user-supplied `function`, etc.) are read exactly as
        # before, one query_results entry at a time.
        SHAREABLE_TYPES = ('h5', 'ATL11', 'ATM_Qfit')
        file_groups = {}
        other_items = []
        if function is None:
            for file_key, result in query_results.items():
                if result['type'] not in SHAREABLE_TYPES:
                    other_items.append((file_key, result))
                    continue
                this_file = file_key if already_resolved else self.resolve_path(file_key, dir_root)
                # offset_end is stored inclusive; readers expect an exclusive
                # stop, so add 1 here unless legacy (last-row-dropped)
                # behavior was requested -- see trim_last_point in the docstring.
                read_offset_end = result['offset_end'] if trim_last_point else result['offset_end'] + 1
                if result['type'] == 'h5':
                    # a ':group' suffix (from a self-contained-file entry) marks
                    # a group inside this_file rather than a separate file
                    if ':' in this_file:
                        physical_file, h5_group = this_file.split(':', 1)
                    else:
                        physical_file, h5_group = this_file, None
                    pair_num = 0
                elif result['type'] == 'ATL11':
                    physical_file, pair = this_file.split(':pair')
                    pair_num = int(pair)
                    h5_group = None
                else:  # 'ATM_Qfit'
                    physical_file, pair_num, h5_group = this_file, 0, None
                for i0, i1 in zip(result['offset_start'], read_offset_end):
                    file_groups.setdefault(physical_file, []).append({
                        'type': result['type'], 'pair_num': pair_num, 'group': h5_group,
                        'index_range': np.array([i0, i1]),
                    })
        else:
            other_items = list(query_results.items())

        for physical_file, tasks in file_groups.items():
            try:
                if not (pc.io_utils.path_exists(physical_file, fs=fs) or pc.io_utils.path_exists(physical_file.split(':')[0], fs=fs)):
                    print(f'geoIndex.get_data(): missing file {physical_file}')
                    continue
                # Read pt1 before pt2 before pt3 (etc.), and in ascending row
                # order within a pair, as a cheap proxy for the file's actual
                # on-disk layout -- ATL11 files are chunked and compressed, so
                # there's no single byte offset to sort by directly, but a
                # standard ATL11 file is written pt1, then pt2, then pt3, so
                # this keeps the shared handle's access pattern close to
                # monotonic rather than jumping around arbitrarily.
                tasks.sort(key=lambda t: (t['pair_num'], int(t['index_range'][0])))
                with pc.io_utils.open_h5(physical_file, fs=fs) as h5f:
                    for task in tasks:
                        try:
                            if task['type'] == 'h5':
                                Di = pc.data().from_h5(filename=physical_file, group=task['group'],
                                                        index_range=task['index_range'],
                                                        field_dict=field_dict, h5_f=h5f, fs=fs)
                            elif task['type'] == 'ATL11':
                                Di = pc.ATL11.data().from_h5(filename=physical_file,
                                                              index_range=task['index_range'],
                                                              pair=task['pair_num'], field_dict=field_dict,
                                                              h5_f=h5f, fs=fs)
                            else:  # 'ATM_Qfit'
                                Di = pc.ATM_Qfit.data().from_h5(physical_file, index_range=task['index_range'],
                                                                 h5_f=h5f, fs=fs)
                            if Di.filename is None:
                                Di.filename = physical_file
                            out_data.append(Di)
                        except Exception as e:
                            if error_action == 'warn':
                                warn(f'geoIndex.py: caught exception reading {physical_file} '
                                     f'(type={task["type"]}, pair={task["pair_num"]}): {e}')
                            else:
                                print(f'geoindex.py: exception for file:{physical_file}')
                                raise(e)
            except Exception as e:
                if error_action == 'warn':
                    warn(f'geoIndex.py: caught exception attempting to read file: {physical_file}')
                    print(e)
                else:
                    print(f'geoindex.py: exception for file:{physical_file}')
                    raise(e)

        for file_key, result in other_items:
            if already_resolved:
                this_file = file_key
            else:
                this_file = self.resolve_path(file_key, dir_root)
            try:
                if not (pc.io_utils.path_exists(this_file, fs=fs) or pc.io_utils.path_exists(this_file.split(':')[0], fs=fs)):
                    print(f'geoIndex.get_data(): missing file {this_file}')
                    continue
                # offset_end is stored inclusive; readers expect an exclusive
                # stop, so add 1 here unless legacy (last-row-dropped) behavior
                # was requested. Only used by the 'ATL06' branch below -- see
                # trim_last_point in the docstring.
                read_offset_end = result['offset_end'] if trim_last_point else result['offset_end'] + 1
                if function is not None:
                    # user has provided a function to read the data
                    D=[function(filename=this_file, index_range=temp, field_dict=field_dict, bounds=bounds) for temp in zip(result['offset_start'], result['offset_end'])]
                elif result['type'] == 'h5_geoindex':
                    D=geoIndex().from_file(this_file).query_xy((result['x'], result['y']), fields=fields, get_data=True, dir_root=dir_root, error_action=error_action)
                elif result['type'] == 'ATL06':
                    this_file, pair = this_file.split(':pair')
                    if fields is None:
                        fields={None:(u'latitude',u'longitude',u'h_li',u'delta_time')}
                    D=[pc.ATL06.data(beam_pair=int(pair), fields=field_list, field_dict=field_dict).from_h5(\
                        filename=this_file, index_range=np.array(temp)) \
                        for temp in zip(result['offset_start'], read_offset_end)]
                elif result['type'] == 'ATM_waveform':
                    D=[pc.atmWaveform(filename=this_file, index_range=np.array(temp), waveform_format=True) for temp in zip(result['offset_start'], result['offset_end'])]
                elif result['type'] == 'geotif':
                    # assume single band
                    D=pc.grid.data().from_geotif(this_file, bounds=bounds, bands=[1], date_format='year').as_points(keep_all=True)
                    D.index(D, np.isfinite(D.z))
                elif result['type'] == 'DEM':
                    D=pc.grid.data().from_geotif(this_file, bounds=bounds, bands=[1], date_format='year')
                    if D.shape is None:
                        continue
                    D=D.as_points(keep_all=True)
                    D.index(np.isfinite(D.z))
                elif result['type'] == 'filtered_DEM':
                    try:
                        D=pc.grid.data().from_geotif(this_file, bounds=bounds, bands=[1], date_format='year')
                        if D.shape is None:
                            continue
                        D=D.as_points(keep_all=True)
                        try:
                            D1=pc.grid.data().from_geotif(this_file, bounds=bounds, bands=[2], date_format='year').as_points(keep_all=True)
                            D.assign({'sigma':D1.z})
                            D.index(np.isfinite(D.z) & np.isfinite(D.sigma))
                        except AttributeError:
                            D.index(np.isfinite(D.z))
                        except TypeError:
                            # this catches missing band 2
                            D.assign(sigma=np.ones_like(D.z))
                            D.index(np.isfinite(D.z))
                        D.filename=this_file
                    except IndexError as e:
                        warn(f"pointCollection.geoIndex: failed to read {this_file}:"+str(e))
                        continue
                elif result['type'] == 'indexed_h5':
                    D = [pc.indexedH5.data(filename=this_file).read([result['x'], result['y']],  fields=fields, index_range=[result['offset_start'], result['offset_end']])]
                elif result['type'] == 'indexed_h5_from_matlab':
                    D = [ pc.indexedH5.data(filename=this_file).read([result['x']/1000, result['y']/1000],  fields=fields, index_range=[result['offset_start'], result['offset_end']])]
                if result['type'] is None:
                    D = [data[np.arange(temp[0], temp[1])] for temp in zip(result['offset_start'], result['offset_end'])]
                # add data to list of results.  May be a list or a single result
                if isinstance(D,list):
                    for Di in D:
                        if Di.filename is None:
                            Di.filename=this_file
                    out_data += D
                else:
                    if D is None:
                        continue
                    if D.filename is None:
                        D.filename=this_file
                    out_data.append(D)
            except Exception as e:
                if error_action=='warn':
                    warn(f'geoIndex.py: caught exception attempting to read file: {this_file}')
                    print(e)
                else:
                    print(f'geoindex.py: exception for file:{this_file}')
                    raise(e)
        return out_data

    def bins_as_array(self):
        """
        return an array containing the locations for all the bins in an index
        """
        if len(self)>0:
            xy_bin=np.c_[[np.fromstring(key, sep='_') for key in self.keys()]]
        else:
            try:
                xy_bin=np.c_[[np.fromstring(key, sep='_') for key in self.h5_file_index.keys()]]
            except AttributeError as e:
                print("AttributeError in bins_as_array, continuing:")
                print(e)
                xy_bin=np.zeros(0)
        if xy_bin.size > 0:
            return (xy_bin[:,0].ravel(), xy_bin[:,1].ravel())
        else:
            return (np.zeros(0), np.zeros(0))

    def bin_latlon(self):
        xy_bin=self.bins_as_array()
        temp=pc.data().from_dict({'x':xy_bin[0],'y':xy_bin[1]}).get_latlon(SRS_proj4=self.attrs['SRS_proj4'])

        #internal_srs=osr.SpatialReference()
        #internal_srs.ImportFromProj4(self.attrs['SRS_proj4'])
        #ll_srs=osr.SpatialReference()
        #ll_srs.ImportFromEPSG(4326)
        #if hasattr(osr,'OAMS_TRADITIONAL_GIS_ORDER'):
        #    ll_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        #ct=osr.CoordinateTransformation( internal_srs, ll_srs)
        #lon, lat, z0 = list(zip(*[ct.TransformPoint(*xy) for xy in zip(np.ravel(xy_bin[:,0]), np.ravel(xy_bin[:,1]), np.ravel(np.zeros_like(xy_bin[:,1])))]))
        #return (lat, lon)
        return (temp.latitude, temp.longitude)

    def keys_from_xy(self, xy):
        pts=unique_points((xy[0].ravel(), xy[1].ravel()), delta=self.attrs['delta'])
        result=[p1 for p1 in ['%d_%d' % p0 for p0 in zip(pts[0], pts[1])] if p1 in self]
        return result

def unique_points(xy, delta=[1, 1]):
    xr=(np.round(np.array(xy[0])/delta[0])*delta[0]).ravel().tolist()
    yr=(np.round(np.array(xy[1])/delta[1])*delta[1]).ravel().tolist()
    xyb=np.concatenate([np.array(xybi).reshape([1,2]) for xybi in set(zip(xr, yr))], axis=0)
    return [xyb[:,0], xyb[:,1]]

def pad_bins(xyb, pad, delta):
    [xp,yp]=np.meshgrid(np.arange(-pad, pad+1)*delta[0], np.arange(-pad, pad+1)*delta[1])
    xp=xp.ravel(); yp=yp.ravel();
    if isinstance(xyb[0],int) or isinstance(xyb[0],float):
        xyb[0]=np.array([xpi+xyb[0] for xpi in xp])
        xyb[1]=np.array([ypi+xyb[1] for ypi in yp])
    else:
        xyb[0]=np.concatenate([xpi+xyb[0] for xpi in xp]).ravel()
        xyb[1]=np.concatenate([ypi+xyb[1] for ypi in yp]).ravel()

    # keep only the unique members of xb and yb
    xyb = unique_points(xyb, delta)
    return xyb


def append_data(group, field, newdata):
    """
    utility function that can append data either to an hdf5 field or a dict of numpy arrays
    """
    try:
        old_shape=np.array(group[field].shape)
        new_shape=old_shape.copy()
        new_shape[0]+=newdata.shape[0]
        group[field].reshape((new_shape))
        group[field][old_shape[0]:new_shape[0],:]=newdata
    except:
        group[field]=np.concatenate((group[field], newdata), axis=0)
    return

def index_list_for_files(filename_list, file_type, delta, SRS_proj4, dir_root=''):
    index_list=list()
    for filename in filename_list:
        index_list.append(geoIndex(SRS_proj4=SRS_proj4, delta=delta).for_file(filename, file_type, dir_root=dir_root, number=0))
    return index_list

def strip_double_slashes(thestr):
    if thestr is None:
        return thestr
    while '//' in thestr:
        thestr=thestr.replace('//','/')
    return thestr
