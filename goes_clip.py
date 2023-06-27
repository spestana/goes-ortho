'''
Functions for clipping GOES ABI imagery to smaller areas
'''

import xarray as xr
import os
from goes_geometry import LonLat2ABIangle

def subsetNetCDF(filepath,bounds,newfilepath=None):
    '''Function to crop a GOES ABI netcdf file to lat/lon bounds.
        Inputs:
            - filepath: path to a netcdf file
            - bounds: list or array containing lat/lon bounds like [min_lat, max_lat, min_lon, max_lon]'''
    # get all the files we just downloaded to "filepath"
    file_list = getListOfFiles(filepath)
    
    # get bounds: Lat_min Lat_max Lon_min Lon_max
    lat_south = bounds[0]
    lat_north = bounds[1]
    lon_west = bounds[2]
    lon_east = bounds[3]
    
    for file_name in file_list:
        with xr.open_dataset(file_name) as file:
            f = file.load()
            # Values needed for geometry calculations
            req = f.goes_imager_projection.semi_major_axis
            rpol = f.goes_imager_projection.semi_minor_axis
            H = f.goes_imager_projection.perspective_point_height + f.goes_imager_projection.semi_major_axis
            lon_0 = f.goes_imager_projection.longitude_of_projection_origin
            e = 0.0818191910435 # GRS-80 eccentricity

            # find corresponding look angles for the four corners
            x_rad_sw, y_rad_sw = LonLat2ABIangle(lon_west,lat_south,0,H,req,rpol,e,lon_0)
            print('SW Corner: {}, {}'.format(x_rad_sw, y_rad_sw))
            x_rad_se, y_rad_se = LonLat2ABIangle(lon_east,lat_south,0,H,req,rpol,e,lon_0)
            print('SE Corner: {}, {}'.format(x_rad_se, y_rad_se))
            x_rad_nw, y_rad_nw = LonLat2ABIangle(lon_west,lat_north,0,H,req,rpol,e,lon_0)
            print('NW Corner: {}, {}'.format(x_rad_nw, y_rad_nw))
            x_rad_ne, y_rad_ne = LonLat2ABIangle(lon_east,lat_north,0,H,req,rpol,e,lon_0)
            print('NE Corner: {}, {}'.format(x_rad_ne, y_rad_ne))
            # choose the bounds that cover the largest extent
            y_rad_s = min(y_rad_sw, y_rad_se) # choose southern-most coordinate (scan angle in radians)
            y_rad_n = max(y_rad_nw, y_rad_ne) # northern-most
            x_rad_e = max(x_rad_se, x_rad_ne) # eastern-most (scan angle in radians)
            x_rad_w = min(x_rad_sw, x_rad_nw) # western-most
            print('Corner coords chosen: N: {}, S: {}; E: {}, W: {}'.format(y_rad_n, y_rad_s, x_rad_e, x_rad_w))
            # Use these coordinates to subset the whole dataset
            y_rad_bnds, x_rad_bnds = [y_rad_n, y_rad_s], [x_rad_w, x_rad_e]
            ds = f.sel(x=slice(*x_rad_bnds), y=slice(*y_rad_bnds))
        # Close the original file
        f.close()
        if newfilepath == None:
            # Overwrite the original file
            ds.to_netcdf(file_name,'w',encoding={'x': {'dtype': 'float'},'y': {'dtype': 'float'}})
        else:
            # save to new file
            ds.to_netcdf(newfilepath,'w',encoding={'x': {'dtype': 'float'},'y': {'dtype': 'float'}})
    
    return None

def getListOfFiles(dirName):
    # create a list of file and sub directories 
    # names in the given directory 
    # https://thispointer.com/python-how-to-get-list-of-files-in-directory-and-sub-directories/
    listOfFile = os.listdir(dirName)
    allFiles = list()
    # Iterate over all the entries
    for entry in listOfFile:
        # Create full path
        fullPath = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory 
        if os.path.isdir(fullPath):
            allFiles = allFiles + getListOfFiles(fullPath)
        else:
            allFiles.append(fullPath)
                
    return allFiles 