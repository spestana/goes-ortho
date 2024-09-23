"""
Geometry functions for GOES-R ABI imagery
"""

import numpy as np


def ABIangle2LonLat(x, y, H, req, rpol, lon_0_deg):
    """
    Computes the latitude and longitude (degrees) of a ground point given GOES-R ABI Fixed Grid image coordinates

    Parameters
    ------------
    x : float
        coordinate in the ABI Fixed Grid, scan angle [radians]
    y : float
        coordinate in the ABI Fixed Grid, elevation angle [radians]
    H : float
        satellite distance to Earth center [km]
    req : float
        Earth semi-major axis of GRS80 ellipsoid (equatorial radius) [km]
    rpol : float
        Earth semi-minor axis of GRS80 ellipsoid (polar radius) [km]
    lon_0_deg : float
        longitude of projection origin (longitude of sub-satellite point) [degrees]

    Returns
    -----------
    lon : float
        longitude of ground point [degrees]
    lat : float
        latitude of ground point [degrees]

    Examples
    -----------

    """

    # intermediate calculations
    a = np.sin(x) ** 2 + (
        np.cos(x) ** 2 * (np.cos(y) ** 2 + (req**2 / rpol**2) * np.sin(y) ** 2)
    )
    b = -2 * H * np.cos(x) * np.cos(y)
    c = H**2 - req**2

    rs = (-b - np.sqrt(b**2 - 4 * a * c)) / (
        2 * a
    )  # distance from satellite point (S) to P

    # solve for rc on the ellipsoid
    # _rc = c*cos(A) ± √[ a2 - c2 sin2 (A) ]
    # add elevation z to rc
    # compute new rs value

    Sx = rs * np.cos(x) * np.cos(y)
    Sy = -rs * np.sin(x)
    Sz = rs * np.cos(x) * np.sin(y)

    # calculate lat and lon
    lat = np.arctan((req**2 / rpol**2) * (Sz / np.sqrt((H - Sx) ** 2 + Sy**2)))
    lat = np.degrees(lat)  # *
    lon = lon_0_deg - np.degrees(np.arctan(Sy / (H - Sx)))

    # handle when longidue is further west of -180 degrees
    if lon < -180:
        lon = lon + 360

    return (lon, lat)


def LonLat2ABIangle(lon_deg, lat_deg, z, H, req, rpol, e, lon_0_deg):
    """
    Computes the GOES-R ABI Fixed Grid image coordinates given latitude and longitude (degrees) of a ground point.

    Parameters
    ------------
    lon_deg : float
        longitude of ground point [degrees]
    lat_deg : float
        latitude of ground point [degrees]
    z : float
        elevation of ground point above GRS80 ellipsoid [meters]
    H : float
        satellite distance to Earth center [km]
    req : float
        Earth semi-major axis of GRS80 ellipsoid (equatorial radius) [km]
    rpol : float
        Earth semi-minor axis of GRS80 ellipsoid (polar radius) [km]
    e : float
        eccentricity of ellipsoid (e=0.0818191910435 for GRS80) [unitless]
    lon_0_deg : float
        longitude of projection origin (longitude of sub-satellite point) [degrees]

    Returns
    ------------
    x : float
        ABI Fixed Grid x coordinate (scan angle) [radians]
    y : float
        ABI Fixed Grid y coordinate (elevation angle) [radians]

    Examples
    ------------

    """

    # convert lat and lon from degrees to radians
    lon = np.radians(lon_deg)
    lat = np.radians(lat_deg)
    lon_0 = np.radians(lon_0_deg)

    # geocentric latitude
    lat_geo = np.arctan((rpol**2 / req**2) * np.tan(lat))

    # geocentric distance to point on the ellipsoid
    rc = rpol / np.sqrt(
        1 - (e**2) * (np.cos(lat_geo) ** 2)
    )  # this is rc if point is on the ellipsoid
    if ~isinstance(z, int):
        rc = (
            rc + z
        )  # this is rc if the point is offset from the ellipsoid by z (meters)

    # intermediate calculations
    Sx = H - rc * np.cos(lat_geo) * np.cos(lon - lon_0)
    Sy = -rc * np.cos(lat_geo) * np.sin(lon - lon_0)
    Sz = rc * np.sin(lat_geo)

    # calculate x and y scan angles
    y = np.arctan(Sz / Sx)
    x = np.arcsin(-Sy / np.sqrt(Sx**2 + Sy**2 + Sz**2))

    ## determine if this point is visible to the satellite
    # condition = ( H * (H-Sx) ) < ( Sy**2 + (req**2 / rpol**2)*Sz**2 )
    # if condition == True:
    #    print('Point at {},{} not visible to satellite.'.format(lon_deg,lat_deg))
    #    return (np.nan, np.nan)
    # else:
    #    return (x,y)
    return (x, y)


def calcLookAngles(lon_deg, lat_deg, lon_0_deg):
    """
    Calculate azimuth and elevation angles for a geostationary satellite viewed from Earth's surface.

    Parameters
    ------------
    lon_deg : float
        longitude of ground point [degrees]
    lat_deg : float
        latitude of ground point [degrees]
    lon_0_deg : float
        longitude of projection origin (longitude of sub-satellite point) [degrees]

    Returns
    ------------
    az : float
        azimuth angle [degrees]
    el : float
        elevation angle [degrees]

    Examples
    ------------

    """

    # convert lat and lon from degrees to radians
    lon = np.radians(lon_deg)
    lat = np.radians(lat_deg)
    lon_0 = np.radians(lon_0_deg)

    s = lon_0 - lon

    el = np.arctan(
        ((np.cos(s) * np.cos(lat)) - 0.1512)
        / (np.sqrt(1 - ((np.cos(s) ** 2) * (np.cos(lat) ** 2))))
    )

    az = np.arctan(np.tan(s) / np.sin(lat))

    return (np.degrees(az) + 180, np.degrees(el))


def goes_lza(lat_ssp, lon_ssp, lat, lon, H=42164.16, r_eq=6378.137):
    """
    Compute the Locan Zenith Angle for a point on Earth surface to a GOES-R geostationary satellite. See more details from NOAA here: https://www.ncdc.noaa.gov/sites/default/files/attachments/GOES-R_ABI_local_zenith_angle_description.docx

    Parameters
    ------------
    lat_ssp : float
       sub-satellite point latitude [degrees]
    lon_ssp : float
       sub-satellite point longitude [degrees]
    lat : float
       view point latitude on Earth's surfaace [degrees]
    lon : float
       view point longitude on Earth's surface [degrees]
    elev : float
       view point elevation (height above GRS80 ellispoid) [km]
    H : float
       satellite distance to Earth center [km] (defaults to 42164.16 km)
    r_eq : float
        Earth semi-major axis (GRS80 ellipsoid) [km] (defaults to 6378.137 km)

    Returns
    ------------
    LZA : float
        local zenith angle [degrees]
    is_point_visible : bool
        True/False flag indicating if the ground point is actually visible to the satellite

    Examples
    ------------

    """

    # intermediate calculation
    B = np.arccos(
        np.cos(np.radians(lat) - np.radians(lat_ssp))
        * np.cos(np.radians(lon) - np.radians(lon_ssp))
    )

    # determine if point is visible to the satellite
    is_point_visible = B < np.arccos(r_eq / (H + r_eq))

    # compute LZA
    LZA_radians = np.arcsin(
        (H * np.sin(B)) / (np.sqrt(H**2 + r_eq**2 - 2 * H * r_eq * np.cos(B)))
    )

    # convert LZA from radians to degrees
    LZA = LZA_radians * 180 / np.pi

    return LZA, is_point_visible


def goes_azi(lat_ssp, lon_ssp, lat, lon):
    """
    Compute azimuth for geostationary satellite, not GOES specific, spherical Earth assumption. See also: http://tiij.org/issues/issues/3_2/3_2e.html

    Parameters
    ------------
    lat_ssp : float
        sub-satellite point latitude [degrees]
    lon_ssp : float
        sub-satellite point longitude [degrees]
    lat : float
        view point latitude on Earth's surfaace [degrees]
    lon : float
        view point longitude on Earth's surface [degrees]

    Returns
    ------------
    azi : float
        azimuth angle [degrees]

    Examples
    ------------

    """

    azi = 180 + np.degrees(
        np.arctan(np.tan(np.radians(lon_ssp - lon)) / np.sin(np.radians(lat)))
    )

    return azi.T


def get_nested_coords(ds, x_rad, y_rad):
    """
    Given the coordinates of a single point in the ABI Fixed Grid coordinates (x_rad and y_rad, in radians) find within a GOES ABI-L1b-Rad dataset, (any of the 2km bands) the coordinates of the nearest "2 km" (56 urad) pixel center, the coordinates of each of the pixel centers of the four "1 km" (28 urad) pixels, and the sixteen "500 m" (14 urad) pixels that are nested within the "2 km" pixel.

    Parameters
    ------------
    ds : xarray.Dataset
        xarray dataset read from a GOES ABI-L1b-Rad NetCDF file of any of the "2 km" bands
    x_rad : float
        x coordinate in the ABI Fixed Grid, scan angle [radians]
    y_rad : float
        y coordinate in the ABI Fixed Grid, elevation angle [radians]

    Returns
    ------------
    nearest_xs_2km : float
        pixel-centered x coordinate of 2km pixel
    nearest_ys_2km : float
        pixel-centered y coordinate of 2km pixel
    nearest_xs_1km : float
        pixel-centered x coordinates of nested 1km pixels
    nearest_ys_1km : float
        pixel-centered y coordinates of nested 1km pixels
    nearest_xs_500m : float
        pixel-centered x coordinates of nested 500 m pixels
    nearest_ys_500m : float
        pixel-centered y coordinates of nested 500 m pixels

    Examples
    ------------

    """

    # "2 km" pixel coordinate
    nearest_xs = ds.sel(x=x_rad, y=y_rad, method="nearest").x
    nearest_ys = ds.sel(x=x_rad, y=y_rad, method="nearest").y
    nearest_xs_2km, nearest_ys_2km = np.meshgrid(nearest_xs, nearest_ys)

    # "1 km" pixel coordinates
    nearest_xs_1km, nearest_ys_1km = np.meshgrid(
        np.linspace(
            nearest_xs_2km[0][0] - (28e-6) * 0.5,
            nearest_xs_2km[0][0] + (28e-6) * 0.5,
            num=2,
        ),
        np.linspace(
            nearest_ys_2km[0][0] - (28e-6) * 0.5,
            nearest_ys_2km[0][0] + (28e-6) * 0.5,
            num=2,
        ),
    )

    # "500 m" pixel coordinates
    nearest_xs_500m, nearest_ys_500m = np.meshgrid(
        np.linspace(
            nearest_xs_2km[0][0] - (14e-6) * 1.5,
            nearest_xs_2km[0][0] + (14e-6) * 1.5,
            num=4,
        ),
        np.linspace(
            nearest_ys_2km[0][0] - (14e-6) * 1.5,
            nearest_ys_2km[0][0] + (14e-6) * 1.5,
            num=4,
        ),
    )

    return (
        nearest_xs_2km,
        nearest_ys_2km,
        nearest_xs_1km,
        nearest_ys_1km,
        nearest_xs_500m,
        nearest_ys_500m,
    )
