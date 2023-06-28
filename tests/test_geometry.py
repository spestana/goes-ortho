'''Functions for testing goes_geometry module.'''

import glob
import xarray as xr
import numpy as np
from goes_ortho import goes_geometry
from goes_ortho import get_data

def test_LonLat2ABIangle():
    # get geometry metadata from one of our CONUS test/example datasets
    filepaths = glob.glob(f'./resources/spestana-goes-ortho-data-bc4e02e/data/C/*.nc')
    with xr.open_dataset(filepaths[0]) as ds:
        req = ds.goes_imager_projection.semi_major_axis
        rpol = ds.goes_imager_projection.semi_minor_axis
        H = ds.goes_imager_projection.perspective_point_height + ds.goes_imager_projection.semi_major_axis
        lon_0 = ds.goes_imager_projection.longitude_of_projection_origin
        e = 0.0818191910435 # GRS-80 eccentricity
    # run LonLat2ABIangle
    (x,y) = goes_geometry.LonLat2ABIangle(-100, 45, 1000, H, req, rpol, e, lon_0)
    assert type(x) == np.float64, "ABI Fixed Grid x coordinate (radians) should be a numpy.float64"
    assert type(y) == np.float64, "ABI Fixed Grid y coordinate (radians) should be a numpy.float64"


get_data.download_example_data()

# test for CONUS
test_LonLat2ABIangle()