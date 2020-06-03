from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='pointCollection',
    version='1.0.0.0',
    description='Utilities for organizing and manipulating point data',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/SmithB/pointCollection',
    author='Ben Smith',
    author_email='besmith@uw.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
    packages=find_packages(),
    install_requires=['numpy','scipy','pyproj','h5py','netCDF4','matplotlib','gdal'],
    dependency_links=['https://github.com/suzanne64/ATL11/tarball/master',
        'https://github.com/tsutterley/pyTMD/tarball/master'],
)
