"""
Functions for clipping GOES ABI imagery to smaller areas
"""

import xarray as xr

from goes_ortho.geometry import LonLat2ABIangle


def subsetNetCDF(filepath, bounds, newfilepath=None):
    """
    Crop a GOES-R ABI NetCDF file to latitude/longitude bounds.

    Parameters
    ------------
    filepath : str
        path to a NetCDF file
    bounds : list
        list or array containing latitude/longitude bounds like [min_lon, min_lat, max_lon, max_lat]
    newfilepath : str
        path and filename for a new NetCDF file to write out to, defaults to None where it will overwrite the input NetCDF file

    Returns
    ------------
    None

    Examples
    ------------
    Subset a GOES ABI CONUS image so that we only have the western half of CONUS within latitudes 30 and 50, and longitudes -125 and -105.

    >>> bounds = [-125, 30, -105, 50]
    >>> subsetNetCDF('CONUS.nc',bounds,newfilepath='westernCONUS.nc')

    """

    # get bounds: [min_lon, min_lat, max_lon, max_lat]
    lon_west = bounds[0]
    lat_south = bounds[1]
    lon_east = bounds[2]
    lat_north = bounds[3]

    with xr.open_dataset(filepath) as file:
        print(file)
        f = file.load()
        # Values needed for geometry calculations
        req = f.goes_imager_projection.semi_major_axis
        rpol = f.goes_imager_projection.semi_minor_axis
        H = (
            f.goes_imager_projection.perspective_point_height
            + f.goes_imager_projection.semi_major_axis
        )
        lon_0 = f.goes_imager_projection.longitude_of_projection_origin
        e = 0.0818191910435  # GRS-80 eccentricity

        # find corresponding look angles for the four corners
        x_rad_sw, y_rad_sw = LonLat2ABIangle(
            lon_west, lat_south, 0, H, req, rpol, e, lon_0
        )
        print("SW Corner: {}, {}".format(x_rad_sw, y_rad_sw))
        x_rad_se, y_rad_se = LonLat2ABIangle(
            lon_east, lat_south, 0, H, req, rpol, e, lon_0
        )
        print("SE Corner: {}, {}".format(x_rad_se, y_rad_se))
        x_rad_nw, y_rad_nw = LonLat2ABIangle(
            lon_west, lat_north, 0, H, req, rpol, e, lon_0
        )
        print("NW Corner: {}, {}".format(x_rad_nw, y_rad_nw))
        x_rad_ne, y_rad_ne = LonLat2ABIangle(
            lon_east, lat_north, 0, H, req, rpol, e, lon_0
        )
        print("NE Corner: {}, {}".format(x_rad_ne, y_rad_ne))
        # choose the bounds that cover the largest extent
        y_rad_s = min(
            y_rad_sw, y_rad_se
        )  # choose southern-most coordinate (scan angle in radians)
        y_rad_n = max(y_rad_nw, y_rad_ne)  # northern-most
        x_rad_e = max(x_rad_se, x_rad_ne)  # eastern-most (scan angle in radians)
        x_rad_w = min(x_rad_sw, x_rad_nw)  # western-most
        print(
            "Corner coords chosen: N: {}, S: {}; E: {}, W: {}".format(
                y_rad_n, y_rad_s, x_rad_e, x_rad_w
            )
        )
        # Use these coordinates to subset the whole dataset
        y_rad_bnds, x_rad_bnds = [y_rad_n, y_rad_s], [x_rad_w, x_rad_e]
        ds = f.sel(x=slice(*x_rad_bnds), y=slice(*y_rad_bnds))
    # Close the original file
    f.close()
    if newfilepath is None:
        # Overwrite the original file
        ds.to_netcdf(
            filepath, "w", encoding={"x": {"dtype": "float"}, "y": {"dtype": "float"}}
        )
    else:
        # save to new file
        ds.to_netcdf(
            newfilepath,
            "w",
            encoding={"x": {"dtype": "float"}, "y": {"dtype": "float"}},
        )

    return None
