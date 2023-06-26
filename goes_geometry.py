"""
Geometry functions for GOES-R ABI imagery
"""

import numpy as np
import pandas as pd
import xarray as xr

def ABIangle2LonLat(x, y, H, req, rpol, lon_0_deg):
    '''This function finds the latitude and longitude (degrees) of point P 
    given x and y, the ABI elevation and scanning angle (radians)'''
    
    # intermediate calculations
    a = np.sin(x)**2 + ( np.cos(x)**2 * ( np.cos(y)**2 + ( req**2 / rpol**2 ) * np.sin(y)**2 ) )
    b = -2 * H * np.cos(x) * np.cos(y)
    c = H**2 - req**2

    rs = ( -b - np.sqrt( b**2 - 4*a*c ) ) / ( 2 * a ) # distance from satellite point (S) to P
    
    # solve for rc on the ellipsoid
    #_rc = c*cos(A) ± √[ a2 - c2 sin2 (A) ]
    # add elevation z to rc
    # compute new rs value
    
    Sx = rs * np.cos(x) * np.cos(y)
    Sy = -rs * np.sin(x)
    Sz = rs * np.cos(x) * np.sin(y)
    
    # calculate lat and lon
    lat = np.arctan( ( req**2 / rpol**2 ) * ( Sz / np.sqrt( ( H - Sx )**2 + Sy**2 ) ) )
    lat = np.degrees(lat) #*
    lon = lon_0_deg - np.degrees( np.arctan( Sy / ( H - Sx )) )
    
    return (lon,lat)

def LonLat2ABIangle(lon_deg, lat_deg, z, H, req, rpol, e, lon_0_deg):
    '''This function finds the ABI elevation (y) and scanning (x) angles (radians) of point P, 
    given a latitude and longitude (degrees)'''
    
    # convert lat and lon from degrees to radians
    lon = np.radians(lon_deg)
    lat = np.radians(lat_deg)
    lon_0 = np.radians(lon_0_deg)
      
    # geocentric latitude
    lat_geo = np.arctan( (rpol**2 / req**2) * np.tan(lat) )

    # geocentric distance to point on the ellipsoid
    rc = rpol / np.sqrt(1 - (e**2)*(np.cos(lat_geo)**2)) # this is rc if point is on the ellipsoid
    if z != 0:
        rc = rc + z # this is rc if the point is offset from the ellipsoid by z (meters)

    # intermediate calculations
    Sx = H - rc * np.cos(lat_geo) * np.cos(lon - lon_0)
    Sy = -rc * np.cos(lat_geo) * np.sin(lon - lon_0)
    Sz = rc * np.sin(lat_geo)
    
    # calculate x and y scan angles
    y = np.arctan( Sz / Sx )
    x = np.arcsin( -Sy / np.sqrt( Sx**2 + Sy**2 + Sz**2 ) )
    
    ## determine if this point is visible to the satellite
    #condition = ( H * (H-Sx) ) < ( Sy**2 + (req**2 / rpol**2)*Sz**2 )
    #if condition == True:
    #    print('Point at {},{} not visible to satellite.'.format(lon_deg,lat_deg))
    #    return (np.nan, np.nan)
    #else:
    #    return (x,y)
    return (x,y)

def calcLookAngles(lon_deg, lat_deg, lon_0_deg):
    '''Calculate azimuth and elevation angles (view from Earth's surface to satellite position)'''
    # convert lat and lon from degrees to radians
    lon = np.radians(lon_deg)
    lat = np.radians(lat_deg)
    lon_0 = np.radians(lon_0_deg)
    
    s = lon_0 - lon
    
    el = np.arctan( ((np.cos(s)*np.cos(lon)) - 0.1512) / (np.sqrt(1 - ((np.cos(s)**2)*(np.cos(lon)**2)))) )
    
    az = np.arctan( np.tan(s)/np.sin(lon) )
    
    return(np.degrees(az) + 180, np.degrees(el))

def goes_lza(lat_ssp, lon_ssp, lat, lon, H=42164.16, r_eq=6378.137):
    
    '''
    Compute the Locan Zenith Angle for a point on Earth surface to a GOES-R geostationary satellite.
        See more details from NOAA here: 
        https://www.ncdc.noaa.gov/sites/default/files/attachments/GOES-R_ABI_local_zenith_angle_description.docx
    
    Inputs:
        GOES-R satellite position
            lat_ssp: sub-satellite point latitude [degrees]
            lon_ssp: sub-satellite point longitude [degrees]
    
        View point (on Earth's surface) position
            lat: view point latitude on Earth's surfaace [degrees]
            lon: view point longitude on Earth's surface [degrees]
            elev: view point elevation (heigh above GRS80 ellispoid) [km]
            
        Earth model parameters (optional)
            H: satellite distance to Earth center [km] (defaults to 42164.16 km)
            r_eq: Earth semi-major axis (GRS80 ellipsoid) [km] (defaults to 6378.137 km)
            
    Returns:
        LZA: local zenith angle [degrees]
        is_point_visible: True/False flag indicating if the ground point is actually visible to the satellite
    
    '''

    # intermediate calculation
    B = np.arccos( np.cos(np.radians(lat)-np.radians(lat_ssp)) * np.cos(np.radians(lon)-np.radians(lon_ssp)) )

    # determine if point is visible to the satellite
    is_point_visible = (B < np.arccos(r_eq / (H+r_eq)))

    # compute LZA
    LZA_radians = np.arcsin( (H * np.sin(B) ) / ( np.sqrt( H**2 + r_eq**2 - 2*H*r_eq*np.cos(B) ) ) )
    
    # convert LZA from radians to degrees
    LZA = LZA_radians * 180/np.pi
    
    return LZA, is_point_visible

def goes_azi(lat_ssp, lon_ssp, lat, lon):
    
    '''quick calculation of azimuth for geostationary satellite, not GOES specific, spherical Earth assumption
    http://tiij.org/issues/issues/3_2/3_2e.html'''
    
    azi = 180 + np.degrees( np.arctan(np.tan(np.radians(lon_ssp - lon))/np.sin(np.radians(lat))) )
    
    return azi.T

def get_nested_coords(ds, x_rad, y_rad):
    
    ''' given the coordinates of a single point in radians: x_rad, y_rad (ABI Fixed Grid coordinates)
        find within ds (GOES ABI-L1b-Rad, any of the 2km bands) the coordinates of the nearest "2 km" (56 urad) pixel center,
        the coordinates of each of the pixel centers of the four "1 km" (28 urad) pixels,
        and the sixteen "500 m" (14 urad) pixels that are nested within the "2 km" pixel. '''
    
    # "2 km" pixel coordinate
    nearest_xs = ds.sel(x=x_rad, y=y_rad, method="nearest").x
    nearest_ys = ds.sel(x=x_rad, y=y_rad, method="nearest").y
    nearest_xs_2km, nearest_ys_2km = np.meshgrid(nearest_xs, nearest_ys)
    
    # "1 km" pixel coordinates
    nearest_xs_1km, nearest_ys_1km = np.meshgrid(np.linspace(nearest_xs_2km[0][0]-(28e-6)*0.5, nearest_xs_2km[0][0]+(28e-6)*0.5, num=2), \
                                                 np.linspace(nearest_ys_2km[0][0]-(28e-6)*0.5, nearest_ys_2km[0][0]+(28e-6)*0.5, num=2))

    # "500 m" pixel coordinates
    nearest_xs_500m, nearest_ys_500m = np.meshgrid(np.linspace(nearest_xs_2km[0][0]-(14e-6)*1.5, nearest_xs_2km[0][0]+(14e-6)*1.5, num=4), \
                                                   np.linspace(nearest_ys_2km[0][0]-(14e-6)*1.5, nearest_ys_2km[0][0]+(14e-6)*1.5, num=4))


    return nearest_xs_2km, nearest_ys_2km, nearest_xs_1km, nearest_ys_1km, nearest_xs_500m, nearest_ys_500m