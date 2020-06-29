Setup and Installation
======================

Presently pointCollection is only available for use as a [GitHub repository](https://github.com/SmithB/pointCollection).
The contents of the repository can be download as a [zipped file](https://github.com/SmithB/pointCollection/archive/master.zip)  or cloned.
To use this repository, please fork into your own account and then clone onto your system.
```bash
git clone https://github.com/SmithB/pointCollection.git
```
#### Dependencies
pointCollection is dependent on some open source programs:
- [gdal](https://gdal.org/index.html)
- [pyproj](https://download.osgeo.org/proj)
- [HDF5](https://www.hdfgroup.org)
- [netCDF](https://www.unidata.ucar.edu/software/netcdf)

The version of GDAL used within pointCollection will match the version of the installed C program.  The path to the C program that will be used with pointCollection is given by:
```bash
gdal-config --datadir
```
The pointCollection installation uses the `gdal-config` routines to set the GDAL package version.  

#### Installation
Can then install using `setuptools`:
```bash
python setup.py install
```
or `pip`:
```bash
python3 -m pip install --user .
```
Alternatively can install the pointCollection utilities directly from GitHub with `pip`:
```
python3 -m pip install --user git+https://github.com/SmithB/pointCollection.git
```
Executable versions of this repository can also be tested using [Binder](https://mybinder.org/v2/gh/SmithB/pointCollection/master) and [Pangeo](https://binder.pangeo.io/v2/gh/SmithB/pointCollection/master).
