# pointCollection
Utilities for organizing and manipulating point data

Please make sure that gdal is installed before you install this package. 

## introduction

This package is designed to allow efficient reading and writing of geospatial data from Python.  It covers some of the same ground as do community tools such as geopandas and pdal, but offers an abbreviated syntax and a modular design that may be helpful to some projects.


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

pointCollection.data objects can also read data from hdf5 files, using the .from_h5() method, which accepts arguments to specify the fields and groups to be read.  Subclasses of pointCollection.data are provided for a few different ice-sheet altimetry data formats:
* ICESat-2 ATL06 and ATL11
* ICESat GLAH12
* IceBridge ATM Qfit, LVIS
* Cryosat-2 data

### pointCollection.grid.data() class

The pointCollection grid format reads and writes data from gridded files.  It is intended to contain 'x' and 'y' fields (and sometimes 't' or 'time' for 3-D data) specifying the corrdinates of the grid, and data fields, the default for which is 'z'.  pc.grid.data() objects can be indexed using 2-d or 3-d slicing operations.  Methods are provided to read and write data from hdf-5 files and geotifs.


### pointCollection.GeoIndex() class

pointCollection also provides the geoIndex class, which is indended to organize data from a variety of different datasets that is contained in individual files.  Once a geoIndex has been created for a set of files, it allows data to be read transparently from all the files, returning a list of pc.data objects for a specified input area.


[![Language](https://img.shields.io/badge/python-v3.6-green.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/SmithB/pointCollection/blob/master/LICENSE)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/SmithB/pointCollection/master)
[![Binder](https://binder.pangeo.io/badge.svg)](https://binder.pangeo.io/v2/gh/SmithB/pointCollection/master)

