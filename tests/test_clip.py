'''Functions for testing goes_clip module.'''

from os import path
import xarray as xr
import goes_clip

### test functions for CONUS images ###

def goes_clip_conus():
    filepath = './tests/resources/C/'
    bounds = [30, 40, -110, -100]
    newfilepath = path.abspath('./tmp.nc')
    goes_clip.subsetNetCDF(filepath,bounds,newfilepath)

def test_goes_clip_conus():
    # run goes clip function
    goes_clip_conus()
    # do we have a dataset?
    assert type(xr.open_dataset('tmp.nc')) == xr.core.dataset.Dataset, "unable to open file"
    ds = xr.open_dataset('tmp.nc')
    # does it have non-zero dimensions?
    assert ds.dims['y'] != 0, "y dimension should be non-zero, check the latitude bounds you've provided"
    assert ds.dims['x'] != 0, "x dimension should be non-zero, check the longitude bounds you've provided"
    ds.close()

### test functions for Full Disk images ###

def goes_clip_fulldisk():   
    filepath = './tests/resources/F/'
    bounds = [-40, 40, -170, 20]
    newfilepath = path.abspath('./tmp.nc')
    goes_clip.subsetNetCDF(filepath,bounds,newfilepath)

def test_goes_clip_fulldisk():
    # run goes clip function
    goes_clip_fulldisk()
    # do we have a dataset?
    assert type(xr.open_dataset('tmp.nc')) == xr.core.dataset.Dataset, "unable to open file"
    ds = xr.open_dataset('tmp.nc')
    # does it have non-zero dimensions?
    assert ds.dims['y'] != 0, "y dimension should be non-zero, check the latitude bounds you've provided"
    assert ds.dims['x'] != 0, "x dimension should be non-zero, check the longitude bounds you've provided"
    ds.close()

### run tests ###

test_goes_clip_conus()

test_goes_clip_fulldisk()