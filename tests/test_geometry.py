'''Functions for testing goes_geometry module.'''

import glob
import xarray as xr
import numpy as np
import goes_ortho as go

def test_LonLat2ABIangle(setup_session):
    # get geometry metadata from one of our CONUS test/example datasets
    filepaths = glob.glob(f'./tests/resources/spestana-goes-ortho-data-bc4e02e/data/C/OR_ABI-L2-ACHAC-M6_G17_s20221210056177_e20221210058550_c20221210101260.nc')
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
    filepaths = glob.glob(f'./tests/resources/spestana-goes-ortho-data-bc4e02e/data/F/OR_ABI-L1b-RadF-M6C16_G18_s20230370050213_e20230370059532_c20230370059568.nc')
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
