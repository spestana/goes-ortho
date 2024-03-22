'''Functions for testing clip module.'''

from os import path
import glob
import xarray as xr
import goes_ortho as go

### test function for CONUS images ###

def test_goes_clip_conus(setup_session):
    # setup for clip function
    file_list = glob.glob('./tests/resources/spestana-goes-ortho-data-*/data/C/*.nc')
    filepath = file_list[0]
    bounds = [30, 40, -110, -100]
    newfilepath = path.abspath('./tmp.nc')
    # run goes clip function
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
    # setup for clip function
    file_list = glob.glob('./tests/resources/spestana-goes-ortho-data-*/data/C/*.nc')
    filepath = file_list[0]
    bounds = [-170, -40, 40,  20] # [min_lon, min_lat, max_lon, max_lat]
    newfilepath = path.abspath('./tmp.nc')
    # run goes clip function
    go.clip.subsetNetCDF(filepath,bounds,newfilepath)
    # do we have a dataset?
    assert type(xr.open_dataset('tmp.nc')) == xr.core.dataset.Dataset, "unable to open file"
    ds = xr.open_dataset('tmp.nc')
    # does it have non-zero dimensions?
    assert ds.dims['y'] != 0, "y dimension should be non-zero, check the latitude bounds you've provided"
    assert ds.dims['x'] != 0, "x dimension should be non-zero, check the longitude bounds you've provided"
    ds.close()
