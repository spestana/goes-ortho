"""
Functions to orthorectify GOES-R ABI images using a DEM
"""

# To do:
# - instead of specifying a DEM .tif file, use elevation
# - what about parts of the land surface we can't see?

#-------------------------------------------------------#

import numpy as np
import pandas as pd
import xarray as xr


#-------------------------------------------------------#

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
    _rc = rpol / np.sqrt(1 - (e**2)*(np.cos(lat_geo)**2)) # this is rc if point is on the ellipsoid
    rc = _rc + z # this is rc if the point is offset from the ellipsoid by z (meters)

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
    
	

def ABIpixelMap(abi_grid_x, abi_grid_y):
    '''
    Converts an array of continuous ABI scan angles into discrete pixel center locations 
    (in scan angle coordinates, incrimenting by the pixel IFOV)
	# NOTE: This function isn't needed for the applying the mapping to a GOES ABI image, 
    # but we can still use this to make some visualizations of what we're doing.
    '''
	
    # IFOV values for GOES ABI bands ("500 m" 14 urad; "1 km" 28 urad; "2 km" 56 urad)
    ifov=np.array([14e-6, 28e-6, 56e-6])
    
    # Convert from scan angle to pixel row/column coordinates 
    x_px = np.array([np.divide(abi_grid_x,i) for i in ifov])
    y_px = np.array([np.divide(abi_grid_y,i) for i in ifov])

    # Get the center coordinate of the pixel each grid cell lies within
    center_x = ((np.floor(np.abs(x_px))+0.5)*np.sign(x_px)) * ifov[:,None,None]
    center_y = ((np.floor(np.abs(y_px))+0.5)*np.sign(y_px)) * ifov[:,None,None]

    # Get the pixel coordinate (row/column) that each grid cell lies within
    #center_col = (np.floor(x_px))
    #center_row = (np.floor(y_px))

    return center_x, center_y#, center_col, center_row



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
	


#-------------------------------------------------------#
# Values needed for geometry calculations
req = 6378137 #g16_b02.goes_imager_projection.semi_major_axis
rpol = 6356752.31414 #g16_b02.goes_imager_projection.semi_minor_axis
H = 42164160 #g16_b02.goes_imager_projection.perspective_point_height + g16_b02.goes_imager_projection.semi_major_axis
lon_0 = -75 #g16_b02.goes_imager_projection.longitude_of_projection_origin
e = 0.0818191910435 # GRS-80 eccentricity

#-------------------------------------------------------#
# Load DEM
dem_filename = 'dem.tif'
dem = xr.open_rasterio(dem_filename) #or grandmesa_dem.tif
dem = dem.where(dem!=dem.nodatavals[0])[0,:,:] # replace nodata with nans
dem = dem.where(dem!=0) # replace zeros with nans


#-------------------------------------------------------#
# For the entire DEM, determine what ABI pixel every DEM grid cell belongs to, taking into account terrain parallax
# map every DEM grid cell to a corresponding ABI pixel (using center of pixel scan angle x and y)

# for each grid cell in the DEM, compute the corresponding ABI scan angle (x and y, radians)
X, Y = np.meshgrid(dem.x,dem.y) # Lon and Lat of each DEM grid cell
Z = dem.values # elevation of each DEM grid cell

# with parallax
abi_grid_x, abi_grid_y = LonLat2ABIangle(X,Y,Z,H,req,rpol,e,lon_0)
center_x, center_y = ABIpixelMap(abi_grid_x, abi_grid_y) # NOTE: this step is no longer required for our final product

# without parallax
_abi_grid_x, _abi_grid_y = LonLat2ABIangle(X,Y,np.zeros_like(Z),H,req,rpol,e,lon_0)
_center_x, _center_y = ABIpixelMap(_abi_grid_x, _abi_grid_y)  # NOTE: this step is no longer required for our final product

#-------------------------------------------------------#
# Calculate look angles
az, el = calcLookAngles(X,Y,lon_0)

#-------------------------------------------------------#

#Create metadata about this map (should probably clean up metadata, adhere to some set of standards)
#
#For every DEM grid cell (our finest resolution unit of information) we have a:
#
#    Longitude
#    Latitude
#    ABI scan angle to this grid cell
#    Elevation

# Create metadata dictionary:
metadata = {
            # Information about the projection geometry:
            'longitude_of_projection_origin': lon_0,
            'semi_major_axis': req,
            'semi_minor_axis': rpol,
            'satellite_height': H,
            'grs80_eccentricity': e,
    
            'longitude_of_projection_origin_info': 'longitude of geostationary satellite orbit',
            'semi_major_axis_info': 'semi-major axis of GRS 80 reference ellipsoid',
            'semi_minor_axis_info': 'semi-minor axis of GRS 80 reference ellipsoid',
            'satellite_height_info': 'distance from center of ellipsoid to satellite (perspective_point_height + semi_major_axis_info)',
            'grs80_eccentricity_info': 'eccentricity of GRS 80 reference ellipsoid',

            # Information about the DEM source file
            'dem_file': dem_filename,
            'dem_crs' : dem.crs,
            'dem_transform' : dem.transform,
            'dem_res' : dem.res,
            'dem_ifov': -9999, # TO DO
    
            'dem_file_info': 'filename of dem file used to create this mapping',
            'dem_crs_info' : 'coordinate reference system from DEM geotiff',
            'dem_transform_info' : 'transform matrix from DEM geotiff', 
            'dem_res_info' : 'resolution of DEM geotiff',
            'dem_ifov_info': 'instantaneous field of view (angular size of DEM grid cell)',
    
            # For each DEM grid cell, we have...
            'dem_px_angle_x_info': 'DEM grid cell X coordinate (east/west) scan angle in the ABI Fixed Grid',
            'dem_px_angle_y_info': 'DEM grid cell Y coordinate (north/south) scan angle in the ABI Fixed Grid',
            'longitude_info': 'longitude from DEM file',
            'latitude_info': 'latitude from DEM file',
            'elevation_info': 'elevation from DEM file'
}


#-------------------------------------------------------#

# Export this pixel map along with the metadata (NetCDF with xarray)
ds = xr.Dataset({    
                'elevation':          (['y', 'x'], dem.values)
                },
    
                coords={'longitude':  (['x'], dem.x),
                        'latitude':   (['y'], dem.y),
                        'dem_px_angle_x':     (['y', 'x'],  abi_grid_x),
                        'dem_px_angle_y':     (['y', 'x'],  abi_grid_y)},
                
                attrs=metadata)
                 
ds.to_netcdf('test_map_tuol_2019.nc',mode='w')

