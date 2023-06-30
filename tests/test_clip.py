'''Functions for testing goes_clip module.'''

from os import path
import xarray as xr
import goes_ortho as go

### test function for CONUS images ###

def test_goes_clip_conus(setup_session):
    # run goes clip function
    filepath = './tests/resources/spestana-goes-ortho-data-bc4e02e/data/C/'
    bounds = [30, 40, -110, -100]
    newfilepath = path.abspath('./tmp.nc')
    go.clip.subsetNetCDF(filepath,bounds,newfilepath)
    # do we have a dataset?
    assert type(xr.open_dataset('tmp.nc')) == xr.core.dataset.Dataset, "unable to open file"
    ds = xr.open_dataset('tmp.nc')
    # does it have non-zero dimensions?
    assert ds.dims['y'] != 0, "y dimension should be non-zero, check the latitude bounds you've provided"
    assert ds.dims['x'] != 0, "x dimension should be non-zero, check the longitude bounds you've provided"
    ds.close()

### test function for Full Disk images ###

def test_goes_clip_fulldisk(setup_session):
    # run goes clip function
    filepath = './tests/resources/spestana-goes-ortho-data-bc4e02e/data/F/'
    bounds = [-40, 40, -170, 20]
    newfilepath = path.abspath('./tmp.nc')
    go.clip.subsetNetCDF(filepath,bounds,newfilepath)
    # do we have a dataset?
    assert type(xr.open_dataset('tmp.nc')) == xr.core.dataset.Dataset, "unable to open file"
    ds = xr.open_dataset('tmp.nc')
    # does it have non-zero dimensions?
    assert ds.dims['y'] != 0, "y dimension should be non-zero, check the latitude bounds you've provided"
    assert ds.dims['x'] != 0, "x dimension should be non-zero, check the longitude bounds you've provided"
    ds.close()
