from distutils.core import setup
import setuptools

setup(name='goes_ortho',
      version='0.2',
      description='Functions for downloading GOES-R ABI imagery, orthorectifying with a DEM, creating timeseries for a single point from a stack of ABI images',
      packages=setuptools.find_packages(),
     )