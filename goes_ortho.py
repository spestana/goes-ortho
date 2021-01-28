"""
Functions to orthorectify GOES-R ABI images using a DEM
"""

# To do:
# - instead of specifying a DEM .tif file, use elevation module
# - what about parts of the land surface we can't see?

#-------------------------------------------------------#

import numpy as np
import pandas as pd
import xarray as xr
import os
import glob
from asp_binder_utils import get_dem, run_bash_command


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


def subset_abi_netcdf(filepath,bounds,new_filepath=None):
    '''Function to crop a GOES ABI netcdf file to lat/lon bounds.
        Inputs:
            - filepath: path to a netcdf file
            - bounds: list or array containing lat/lon bounds like [min_lat, max_lat, min_lon, max_lon]
                      where latitude is between -90 and 90, longitude between -180 and 180
            - new_filepath: path and filename of new file to save the cropped image to. If not provided or set to None, this will overwrite the original file at "filepath"
            '''

    # Show us the bounds we'll crop images to
    print('Subsetting \n{filepath}\n to these bounds:'.format(filepath=filepath))
    print('\t({w},{n}).\t.({e},{n})\n\n\n\n\t({w},{s}).\t.({e},{s})\n'.format(n=bounds[1],w=bounds[2],e=bounds[3],s=bounds[0]))

    
    # get bounds: Lat_min Lat_max Lon_min Lon_max
    lat_south = bounds[0]
    lat_north = bounds[1]
    lon_west = bounds[2]
    lon_east = bounds[3]
    

    with xr.open_dataset(filepath, decode_times=False) as file:
        # NOTE: for some reason (?) I sometimes get an error "ValueError: unable to decode time units 'seconds since 2000-01-01 12:00:00' with the default calendar. Try opening your dataset with decode_times=False." so I've added decode_times=False here.
        f = file.load()
        # Values needed for geometry calculations
        req = f.goes_imager_projection.semi_major_axis
        rpol = f.goes_imager_projection.semi_minor_axis
        H = f.goes_imager_projection.perspective_point_height + f.goes_imager_projection.semi_major_axis
        lon_0 = f.goes_imager_projection.longitude_of_projection_origin
        e = 0.0818191910435 # GRS-80 eccentricity
        # find corresponding look angles
        x_rad_w, y_rad_s = LonLat2ABIangle(lon_west,lat_south,0,H,req,rpol,e,lon_0)
        #print('SW Corner: {}, {}'.format(x_rad_w, y_rad_s))
        x_rad_e, y_rad_n = LonLat2ABIangle(lon_east,lat_north,0,H,req,rpol,e,lon_0)
        #print('NE Corner: {}, {}'.format(x_rad_e, y_rad_n))
        # Use these coordinates to subset the whole dataset
        y_rad_bnds, x_rad_bnds = [y_rad_n, y_rad_s], [x_rad_w, x_rad_e]
        ds = f.sel(x=slice(*x_rad_bnds), y=slice(*y_rad_bnds))
        
        # Close the original file
        f.close()
        
        # If no new_filepath is provided, overwrite the existing file
        if new_filepath == None:
            new_filepath = filepath
        
        # Write the new file
        ds.to_netcdf(new_filepath,'w',encoding={'x': {'dtype': 'float'},'y': {'dtype': 'float'}}) #
    
    return None



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
    

def make_ortho_map(goes_filepath, dem_filepath, out_filepath=None):
    '''For the entire DEM, determine the ABI scan angle coordinates for every DEM grid cell, 
    taking into account the underlying terrain and satellite's viewing geometry.
    
    Create the mapping between GOES-R ABI pixels (netCDF input file) and a DEM grid (geotiff input file)
    '''
    
    print('\nRUNNING: make_ortho_map()')
    
    # Open the GOES ABI image
    print('\nOpening GOES ABI image...')
    abi_image = xr.open_dataset(goes_filepath)
    # Get inputs: projection information from the ABI radiance product (values needed for geometry calculations)
    print('\nGet inputs: projection information from the ABI radiance product')
    req = abi_image.goes_imager_projection.semi_major_axis
    rpol = abi_image.goes_imager_projection.semi_minor_axis
    H = abi_image.goes_imager_projection.perspective_point_height + abi_image.goes_imager_projection.semi_major_axis
    lon_0 = abi_image.goes_imager_projection.longitude_of_projection_origin
    e = 0.0818191910435 # GRS-80 eccentricity
    print('...done')
    
    
    # Load DEM
    print('\nOpening DEM file...')
    dem = xr.open_rasterio(dem_filepath)
    dem = dem.where(dem!=dem.nodatavals[0])[0,:,:] # replace nodata with nans
    dem = dem.where(dem!=0) # replace zeros with nans
    # TO DO: use elevation library or something similar to grab an SRTM3 DEM from online for this step
    # Create 2D arrays of longitude and latitude from the DEM
    print('\nCreate 2D arrays of longitude and latitude from the DEM')
    X, Y = np.meshgrid(dem.x,dem.y) # Lon and Lat of each DEM grid cell
    Z = dem.values # elevation of each DEM grid cell
    print('...done')
    
    # For each grid cell in the DEM, compute the corresponding ABI scan angle (x and y, radians)
    print('\nFor each grid cell in the DEM, compute the corresponding ABI scan angle (x and y, radians)')
    abi_grid_x, abi_grid_y = LonLat2ABIangle(X,Y,Z,H,req,rpol,e,lon_0)
    print('...done')
    
    # Create metadata dictionary about this map (should probably clean up metadata, adhere to some set of standards)
    print('\nCreate metadata dictionary about this map')
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
                'dem_file': dem_filepath,
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
    print('...done')
    
    # Create pixel map dataset
    print('\nCreate pixel map dataset')
    ds = xr.Dataset({    
                    'elevation':          (['latitude', 'longitude'], dem.values)
                    },
        
                    coords={'longitude':  (['longitude'], dem.x),
                            'latitude':   (['latitude'], dem.y),
                            'dem_px_angle_x':     (['latitude', 'longitude'],  abi_grid_x),
                            'dem_px_angle_y':     (['latitude', 'longitude'],  abi_grid_y)},
                    
                    attrs=metadata)
    print(ds)
    print('...done')
                     
    if out_filepath != None:
        print('\nExport this pixel map along with the metadata (NetCDF with xarray)')
        # Export this pixel map along with the metadata (NetCDF with xarray)
        ds.to_netcdf(out_filepath,mode='w')
        print('...done')
    
    # Return the pixel map dataset
    print('\nReturn the pixel map dataset.')
    
    return ds

def orthorectify_abi_rad(goes_filepath, pixel_map, out_filename=None):
    '''Using the pixel mapping for a specific ABI viewing geometry over a particular location,
    orthorectify the ABI radiance values and return an xarray dataarray with those values.'''
    print('\nRUNNING: orthorectify_abi_rad()')
    
    # First check, Does the projection info in the image match our mapping?
    print('\nDoes the projection info in the image match our mapping?')
    # Open the GOES ABI image
    print('\nOpening GOES ABI image...\t\t\tABI image value\tPixel map value')
    abi_image = xr.open_dataset(goes_filepath)
    print('perspective_point_height + semi_major_axis:\t{}\t{}'.format(abi_image.goes_imager_projection.perspective_point_height 
                                                                       + abi_image.goes_imager_projection.semi_major_axis,
                                                          pixel_map.satellite_height))
    print('semi_major_axis:\t\t\t\t{}\t{}'.format(abi_image.goes_imager_projection.semi_major_axis,
                                                          pixel_map.semi_major_axis))
    print('semi_minor_axis:\t\t\t\t{}\t{}'.format(abi_image.goes_imager_projection.semi_minor_axis,
                                                          pixel_map.semi_minor_axis))
    print('longitude_of_projection_origin:\t\t\t{}\t\t{}'.format(abi_image.goes_imager_projection.longitude_of_projection_origin,
                                                          pixel_map.longitude_of_projection_origin))
    print('...done')
    
    # Map (orthorectify) and clip the image to the pixel map
    print('\nMap (orthorectify) and clip the image to the pixel map')
    abi_rad_values = abi_image.sel(x=pixel_map.dem_px_angle_x, y=pixel_map.dem_px_angle_y, method='nearest').Rad.values
    print('...done')
	
    # Output this result to a new NetCDF file
	#Create a new xarray dataset with the orthorectified ABI radiance values, 
    #Lat, Lon, Elevation, and metadata from the pixel map. 
    print('\nOutput this result to a new NetCDF file')
    if out_filename == None:
        out_filename=abi_image.dataset_name+'_ortho.nc'
    print('Saving file as: {}'.format(out_filename))
    pixel_map['rad'] = (['latitude', 'longitude'], abi_rad_values)
    pixel_map['tb'] = goesBrightnessTemp(pixel_map.rad, abi_image.planck_fk1.values, abi_image.planck_fk2.values, abi_image.planck_bc1.values, abi_image.planck_bc2.values)
    pixel_map.to_netcdf(out_filename)
    print('...done')
    
    return pixel_map


def orthorectify_abi(goes_filepath, pixel_map, data_vars, out_filename=None):
    '''Using the pixel mapping for a specific ABI viewing geometry over a particular location,
    orthorectify the ABI radiance values and return an xarray dataarray with those values.'''
    print('\nRUNNING: orthorectify_abi_rad()')
    
    # First check, Does the projection info in the image match our mapping?
    print('\nDoes the projection info in the image match our mapping?')
    # Open the GOES ABI image
    print('\nOpening GOES ABI image...\t\t\tABI image value\tPixel map value')
    abi_image = xr.open_dataset(goes_filepath)
    print('perspective_point_height + semi_major_axis:\t{}\t{}'.format(abi_image.goes_imager_projection.perspective_point_height 
                                                                       + abi_image.goes_imager_projection.semi_major_axis,
                                                          pixel_map.satellite_height))
    print('semi_major_axis:\t\t\t\t{}\t{}'.format(abi_image.goes_imager_projection.semi_major_axis,
                                                          pixel_map.semi_major_axis))
    print('semi_minor_axis:\t\t\t\t{}\t{}'.format(abi_image.goes_imager_projection.semi_minor_axis,
                                                          pixel_map.semi_minor_axis))
    print('longitude_of_projection_origin:\t\t\t{}\t\t{}'.format(abi_image.goes_imager_projection.longitude_of_projection_origin,
                                                          pixel_map.longitude_of_projection_origin))
    print('...done')
    
    # Map (orthorectify) and clip the image to the pixel map for each data variable we want
    for var in data_vars:
        print('\nMap (orthorectify) and clip the image to the pixel map for {}'.format(var))
        abi_var_values = abi_image.sel(x=pixel_map.dem_px_angle_x, y=pixel_map.dem_px_angle_y, method='nearest')[var].values
        print('...done')
        
        #Create a new xarray dataset with the orthorectified ABI radiance values, 
        #Lat, Lon, Elevation, and metadata from the pixel map. 
        pixel_map[var] = (['latitude', 'longitude'], abi_var_values)
        # If we are looking at an ABI-L1b-Rad product, create either a reflectance (bands 1-6) or brightness temperautre (bands 7-16) dataset
        if var == 'Rad':
            # if we are looking at bands 1-6, compute reflectance
            if abi_image.band_id.values[0] <= 6:
                pixel_map['ref'] = ref_or_tb = goesReflectance(pixel_map[var], abi_image.kappa0.values)
            # else, compute brightness temperature for bands 7-16
            else:
                pixel_map['tb'] = goesBrightnessTemp(pixel_map[var], abi_image.planck_fk1.values, abi_image.planck_fk2.values, abi_image.planck_bc1.values, abi_image.planck_bc2.values)
        
    # Output this result to a new NetCDF file
    print('\nOutput this result to a new NetCDF file')
    if out_filename == None:
        out_filename=abi_image.dataset_name+'_ortho.nc'
    print('Saving file as: {}'.format(out_filename))
    
    pixel_map.to_netcdf(out_filename)
    print('...done')
    
    return pixel_map


def ortho(goes_image_path, data_vars, bounds, new_goes_filename, dem_filepath=None, demtype='SRTMGL3', keep_dem=True):
    '''Wraps around get_dem(), make_ortho_map(), orthorectify_abi()'''
    
    if dem_filepath == None:
        dem_filepath = 'temp_{demtype}_DEM.tif'.format(demtype=demtype)
    get_dem(demtype=demtype, 
            bounds=bounds, 
            out_fn=dem_filepath, 
            proj='+proj=lonlat +datum=GRS80') # make sure to convert to GRS80 ellipsoid model GOES ABI fixed grid uses
    
    # create the mapping between scan angle coordinates and lat/lon given the GOES satellite position and our DEM
    goes_ortho_map = make_ortho_map(goes_image_path, 
                                    dem_filepath)
    
    # Apply the "ortho map" and save a new NetCDF file with data variables from the original file
    goes_ds = orthorectify_abi(goes_image_path, 
                               goes_ortho_map,
                               data_vars,
                               out_filename=new_goes_filename)
    
    # If keep_dem is False, delete the temporary DEM file we downloaded
    if keep_dem == False:
        os.remove(dem_filepath)
    
    return None

def output_ortho_netcdf(abi_rad_values, pixel_map, out_filename):

    print('\nRUNNING: output_ortho_netcdf()')
    
    # some metadata for this
    metadata = {'rad' : 'units'}
    
    # make the data array
    #rad_da = xr.DataArray(abi_rad_values, 
    #                      dims=('y','x'),
    #                      coords={'latitude': (['y'], pixel_map.latitude),
    #                            'longitude': (['x'], pixel_map.longitude)},
    #                     attrs=metadata)

    
    return None


def goesBrightnessTemp(rad, fk1, fk2, bc1, bc2): 
    # Convert Radiance to Brightness Temperature for GOES-R ABI emissive bands
    Tb = ( fk2 / (np.log((fk1 / rad) + 1)) - bc1 ) / bc2
    return Tb
	
def goesReflectance(rad, kappa): 
    # Convert Radiance to Reflectance for GOES-R ABI reflective bands
    ref = kappa * rad
    return ref
	
def make_abi_timeseries(directory, product, data_vars, lon, lat, z, outfilepath=None):
    '''Given a directory of GOES ABI products, create a timeseries of data variables (specified in data_vars) for a single point (at lon, lat, elevation).
	   Returns a pandas dataframe, optional output to a csv file.'''

    path = '{directory}/**/*{product}*.nc'.format(directory=directory, product=product)
    file_list = glob.glob(path, recursive=True) 
    
    # create empty dataframe to hold the data variables we want plus a timestamp
    df_columns = [var for var in data_vars]
    df_columns.append('time')
    # if Radiance is one of the data variables we are interested in
    if 'Rad' in data_vars:
        # create a new column for reflectance (for bands 1-6) or brightness temperature (for band 7-16)
        df_columns.append('ref_or_tb')
    # create the data frame we will populate with values
    df = pd.DataFrame(columns=df_columns)
    
    print('Creating a timeseries of {data_vars} from {product} at ({lat}, {lon}, {z})'.format(data_vars=data_vars,
                                                                                                 product=product,
                                                                                                 lat=lat,
                                                                                                 lon=lon,
                                                                                                 z=z))
    print('Reading:')
    for filename in file_list:
        print('{}'.format(filename), end='\r')
    
        with xr.open_dataset(filename, decode_times=False) as f:
            # I've included "decode_times=False" to this xr.open_dataset because I've encountered some ABI-L2-ACMC files where the timestamp couldn't be read
            # and xarray gave a "ValueError: unable to decode time units 'seconds since 2000-01-01 12:00:00' with the default calendar. Try opening your dataset with decode_times=False."
            # I've also switched which timestamp from the ABI files I'm reading (was f.time_bounds.values.min(), now f.time_coverage_start)
                       
            # Read goes_imager_projection values needed for geometry calculations           
            # and compute the corresponding look angles (in radiance) for the lat, lon, elevation we are interested in
            x_rad, y_rad = LonLat2ABIangle(lon,
                                           lat,
                                           z,
                                           f.goes_imager_projection.perspective_point_height + f.goes_imager_projection.semi_major_axis,
                                           f.goes_imager_projection.semi_major_axis,
                                           f.goes_imager_projection.semi_minor_axis,
                                           0.0818191910435, # GRS-80 eccentricity
                                           f.goes_imager_projection.longitude_of_projection_origin)
        
            # get the timestamp for this observation (these should all be UTC, but I am removing timezone info because not all timestamps are converting the same way, and I was getting a "Cannot compare tz-naive and tz-aware timestamps" error)
            timestamp = pd.Timestamp(f.time_coverage_start).replace(tzinfo=None)
            
            # create an empty dictionary we will populate with values from file f
            this_row_dict = {}
            
            # create an empty list of the same length as data_vars to hold each variable's value
            values = ['' for v in data_vars]
            
            # For each variable we are interested, specified in the list "data_vars"
            for i, var in enumerate(data_vars):
                            
                # find corresponding pixel data_var value nearest to these scan angles y_rad and x_rad
                values[i] = f[var].sel(y=y_rad, x=x_rad, method='nearest').values.mean()
                
                # For all other products set ref_or_tb to None
                ref_or_tb = None
                # For ABI-L1b-Rad products only:
                if var == 'Rad':
                    # If we are looking at a reflective band (bands 1-6), convert Radiance to Reflectance
                    if f.band_id.values <= 6:
                        ref_or_tb = goesReflectance(values[i], f.kappa0.values)
                    # If we are looking at an emissive band (bands 7-16), convert Radiance to Brightness Temperature (K)
                    else:
                        ref_or_tb = goesBrightnessTemp(values[i], f.planck_fk1.values, f.planck_fk2.values, f.planck_bc1.values, f.planck_bc2.values)
            
            
            # create a dictionary for this row of values (where each row is a GOES-R observation time)
            this_row_dict = dict( zip(data_vars, values ))
            # add our time stamp to this dict before we update the dataframe
            this_row_dict['time'] = timestamp
            
            # If we have reflectance or brightness temperature to add to our dataframe
            if ref_or_tb is not None:
                # add reflectance or brightness temperature to this row's update dict
                this_row_dict['ref_or_tb'] = ref_or_tb
            
            # Finally, append this_row_dict to our dataframe for this one GOES-R observation time
            df = df.append(this_row_dict, ignore_index=True)
            
    
    # set the dataframe intext to the timestamp column
    df.set_index('time', inplace = True, verify_integrity = True)
         
    # if an output filepath was provided, save the dataframe as a csv
    if outfilepath is not None:
        print('Saving csv file to: {}'.format(outfilepath))
        df.to_csv(outfilepath)
    
    return df



