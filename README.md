# pointCollection
Utilities for organizing and manipulating point data

Please make sure that gdal is installed before you install this package. 

## introduction

This package is designed to allow efficient reading and writing of geospatial data from Python.  It covers some of the same ground as do community tools such as geopandas and pdal, but offers an abbreviated syntax and a modular design that may be helpful to some projects.  Benefits include:

- Flexible input and output from hdf5 and netCDF files: Data can be read from file hierarchies, and field names can be reassigned within the function calls (from_h5 and to_h5 functions in _data_ and _grid.data_ objects, _from_nc_ and _to_nc_ functions in _grid.data_ objects)

- Efficient spatial subsetting of objects: This is can be carried out during reading or in memory (.crop, .cropped methods of _data_ and _grid.data_ objects)

- Fast transforms between coordinate systems:  Transormation routines are inherited from the proj4 library (to_xy and to_latlon functions in _data_ objects)

- Data exchange between point data and gridded data: _grid.data_ objects can be converted to collections of points, and collections of points can be mapped to grids (

- Subclassible data objects: The _data_ and _grid.data_ objects are designed to be subclassible to create objects with custom readers and internal field names (see examples for altimetry point formats: ATL06, ATL11, ATM_Qfit, CS2)

- GeoIndex class:  The GeoIndex class can summarize the geolocation information in large numbers of files on disk, and allows efficient access to geospatial data stored within these files

- IndexedH5 data format: In this format point data are stored so that spatailly adjacent data stored in contiguous blocks within the file.  This can greatly improve the efficiency of data reads from the file.

## Installation

An efficient way to add pointCollection to your python environment is to clone the repository, change into the repository directory, and install the package with pip:
```bash
git clone https://github.com/SmithB/pointCollection.git
cd pointCollection
pip install -e .
```
(I like to install packages in developer mode.)

# What this package provides
The package provides two main base classes, and a set of helper classes designed for efficient data access of files on disk.  In general use, it is usually imported with the abbreviation:

```python
import pointCollection as pc
```

# Data classes

### pointCollection.data() class
This is the main class for this project.  It is intended to contain data, and has methods that read from and write to specific formats (primarily from hdf5).  It is constructed either by reading from a file, or from a dictionary interface:

```python
import pointCollection as pc
import numpy as np
D=pc.data().from_dict({'x':np.arange(5),'y':np.arange(5),'z':np.arange(5)**2})
```

This creates a pointCollection.data instance with fields 'x','y', and 'z', and we can access the data using dot syntax:

```python
print(D.z)
--[0 1 2 3 4]
```

All the fields in the object can also be indexed by applying standard Python slicing/indexing syntax to the object; for example:

```python
D.z[0:3]
```
returns a new pointCollection.data object containing first three elements of all fields of D.

pointCollection.data objects can also read data from hdf5 files, using the .from_h5() method, which accepts arguments to specify the fields and groups to be read.  

#### Methods:

pointCollection is designed to allow intelligent handling of subsets of data.  pc.data objects can be cropped, and indexed (affecting all fields simultaneously).  They can be read from and written to hdf5 files, and their coordinates can be transformed from gegraphic (latitude/longitude) to projected coordinates, and vice versa.  The _bounds_ method provides the spatial extent of the object.

#### Properties:

pointCollection.data objects have _size_ and _shape_ attributes that encode the dimensions of the fields in the object.  They have a _fields_ attribute that lists the fields that have been assigned to the object.  Objects that have been read from files retain a _filename_ attribute that specifies their source.

### pointCollection.grid.data() subclass

The pointCollection grid subclass reads and writes data from gridded files.  It is intended to contain 'x' and 'y' fields (and sometimes 't' or 'time' for 3-D data) specifying the corrdinates of the grid, and data fields, the default for which is 'z'.  pc.grid.data() objects can be indexed using 2-d or 3-d slicing operations.  Methods are provided to read and write data from hdf-5 files and geotifs.

### pointCollection.GeoIndex() class

pointCollection also provides the geoIndex class, which is indended to organize data from a variety of different datasets that is contained in individual files.  Once a geoIndex has been created for a set of files, it allows data to be read transparently from all the files, returning a list of pc.data objects for a specified input area.

### Subclassing pointCollection:

One strength of pointCollection is that the pointCollection.data and pointCollection.grid.data classes can be adapted to contain a wide variety of data.  Subclasses of pointCollection data objects can be adapted to handle reading and translation of different data formats, usually by creating a custom from_file() or from_h5() method for the subclass.  The subclasses then inherit the subsetting and projection methods from the parent class.  The subclass definition helps users read a standard set of field from a file, optionally translating the field names used in the file into a more standard name (i.e. translating a "time_of_day" field in the file into a "t" field in the pointCollection object).

Subclasses of pointCollection.data are provided for a few different ice-sheet altimetry data formats:
* ICESat-2 ATL06 and ATL11
* ICESat GLAH12
* IceBridge ATM Qfit, LVIS
* Cryosat-2 data
The definitions for these formats should serve as a model for how to generate a new subclass.

[![Language](https://img.shields.io/badge/python-v3.6-green.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/SmithB/pointCollection/blob/master/LICENSE)

