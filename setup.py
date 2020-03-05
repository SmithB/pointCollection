from setuptools import setup, find_packages
setup(
    name='Utilities for organizing and manipulating point data',
    version='1.0.0.0',
    description='Python',
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
)
