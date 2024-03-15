'''Functions for testing goes_geometry module.'''

from os import path
import sys
from pathlib import Path
import xarray as xr
from goes_ortho import goes_geometry
from goes_ortho import get_data

### test functions for CONUS images ###

def goes_clip_conus():
    filepath = './resources/spestana-goes-ortho-data-bc4e02e/data/C/'
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

def read_test_data() -> list:   
    filepath = './resources/spestana-goes-ortho-data-bc4e02e/data/F/'
    bounds = [-40, 40, -170, 20]
    newfilepath = path.abspath('./tmp.nc')
    goes_clip.subsetNetCDF(filepath,bounds,newfilepath)
    
    return [H, req, rpol, e, lon_0_deg]

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

get_data.download_example_data()

test_goes_clip_conus()

test_goes_clip_fulldisk()
=======

def test_LonLat2ABIangle(setup_session):
    # get geometry metadata from one of our CONUS test/example datasets
    filepaths = glob.glob(f'./tests/resources/spestana-goes-ortho-data-*/data/C/OR_ABI-L2-ACHAC-M6_G17_s20221210056177_e20221210058550_c20221210101260.nc')
    with xr.open_dataset(filepaths[0]) as ds:
        req = ds.goes_imager_projection.semi_major_axis
        rpol = ds.goes_imager_projection.semi_minor_axis
        H = ds.goes_imager_projection.perspective_point_height + ds.goes_imager_projection.semi_major_axis
        lon_0 = ds.goes_imager_projection.longitude_of_projection_origin
        e = 0.0818191910435 # GRS-80 eccentricity
    # run LonLat2ABIangle
    (x,y) = go.geometry.LonLat2ABIangle(-100, 45, 1000, H, req, rpol, e, lon_0)
    assert x == 0.06993880298923384, "ABI Fixed Grid x coordinate (radians) should be a numpy.float64"
    assert y == 0.11588291952114481, "ABI Fixed Grid y coordinate (radians) should be a numpy.float64"


def test_both(setup_session):
    # get geometry metadata from one of our Full Disk test/example datasets
    filepaths = glob.glob(f'./tests/resources/spestana-goes-ortho-data-*/data/F/OR_ABI-L1b-RadF-M6C16_G18_s20230370050213_e20230370059532_c20230370059568.nc')
    print(filepaths[0])
    with xr.open_dataset(filepaths[0]) as ds:
        req = ds.goes_imager_projection.semi_major_axis
        rpol = ds.goes_imager_projection.semi_minor_axis
        H = ds.goes_imager_projection.perspective_point_height + ds.goes_imager_projection.semi_major_axis
        lon_0 = ds.goes_imager_projection.longitude_of_projection_origin
        e = 0.0818191910435 # GRS-80 eccentricity
    # test both LonLat2ABIangle and ABIangle2LonLat
    original_lat = -45
    original_lon = 170
    (x,y) = go.geometry.LonLat2ABIangle(original_lon, original_lat, 0, H, req, rpol, e, lon_0)
    (lon, lat) = go.geometry.ABIangle2LonLat(x, y, H, req, rpol, lon_0)
    # allowing a very small difference (1e-9) due to rounding errors/precision limits in the functions
    assert abs(lon - original_lon) < 1e-9, "Longitude does not match original longitude after conversions"
    assert abs(lat - original_lat) < 1e-9, "Latitude does not match original latitude after conversions"

