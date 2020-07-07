import numpy as np
import pandas as pd
import xarray as xr
import os
from goespy.Downloader import ABI_Downloader # https://github.com/palexandremello/goes-py
import sys, argparse

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
    
def LonLat2ABIangle_ellipsoid(lon_deg, lat_deg, H, req, rpol, e, lon_0_deg):
    '''Find the ABI elevation (y) and scanning (x) angles (radians) of point P , given a latitude and longitude (degrees)
    NOTE: this assumes the ideal ellipsoidal shape of the Earth defined by GRS80'''

    # convert lat and lon from degrees to radians
    lon = np.radians(lon_deg)
    lat = np.radians(lat_deg)
    lon_0 = np.radians(lon_0_deg)
      
    # geocentric latitude
    lat_geo = np.arctan( (rpol**2 / req**2) * np.tan(lat) )

    # geocentric distance to point on the ellipsoid
    rc = rpol / np.sqrt(1 - (e**2)*(np.cos(lat_geo)**2)) # this is rc if point is on the ellipsoid

    # intermediate calculations
    Sx = H - rc * np.cos(lat_geo) * np.cos(lon - lon_0)
    Sy = -rc * np.cos(lat_geo) * np.sin(lon - lon_0)
    Sz = rc * np.sin(lat_geo)
    
    # calculate x and y scan angles
    y = np.arctan( Sz / Sx )
    x = np.arcsin( -Sy / np.sqrt( Sx**2 + Sy**2 + Sz**2 ) )
    
    # determine if this point is visible to the satellite
    condition = ( H * (H-Sx) ) < ( Sy**2 + (req**2 / rpol**2)*Sz**2 )
    if condition == True:
        print('Point at {},{} not visible to satellite.'.format(lon_deg,lat_deg))
        return (np.nan, np.nan)
    else:
        return (x,y)
        


def subsetNetCDF(filepath,bounds):
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

            # find corresponding look angles
            x_rad_w, y_rad_s = LonLat2ABIangle_ellipsoid(lon_west,lat_south,H,req,rpol,e,lon_0)
            #print('SW Corner: {}, {}'.format(x_rad_w, y_rad_s))
            x_rad_e, y_rad_n = LonLat2ABIangle_ellipsoid(lon_east,lat_north,H,req,rpol,e,lon_0)
            #print('NE Corner: {}, {}'.format(x_rad_e, y_rad_n))

            # Use these coordinates to subset the whole dataset
            y_rad_bnds, x_rad_bnds = [y_rad_n, y_rad_s], [x_rad_w, x_rad_e]
            ds = f.sel(x=slice(*x_rad_bnds), y=slice(*y_rad_bnds))
        # Close the original file
        f.close()
        # Overwrite the original file
        ds.to_netcdf(file_name,'w',encoding={'x': {'dtype': 'float'},'y': {'dtype': 'float'}}) #
    
    return None

##############################################################

#---------------------------- COMMAND LINE ARGUMENTS ----------------------------#
# Set up argument and error handling
parser = argparse.ArgumentParser(description='Produces a timeseries of GOES ABI Radiance observations for a single location given a directory of GOES ABI files')
parser.add_argument('-B','--bucket', required=True, type=str, help='AWS S3 Bucket for GOES (e.g. noaa-goes16')
parser.add_argument('-Y','--year', required=True, type=int, help='Specify time range to search for GOES ABI imagery (year)')
parser.add_argument('-M','--month', required=True, type=int, help='Specify time range to search for GOES ABI imagery (month)')
parser.add_argument('-D','--days', required=True, type=int, nargs=2, help='Specify time range to search for GOES ABI imagery (start day, stop day)')
parser.add_argument('-p','--product', required=True, type=str, help='GOES ABI Product (e.g. ABI-L1b-RadC)')
parser.add_argument('-c','--channel', required=True, type=str, help='GOES ABI channel/band (e.g. C14)')
parser.add_argument('-b','--bounds', required=True, type=int, nargs=4, help='Bounds to crop GOES ABI image to (min_lat max_lat min_lon max_lon')
parser.add_argument('-d','--dir', required=True, help='Directory to save GOES ABI files (.nc)')
args = parser.parse_args()

#-----------------------SET ARGUMENTS TO VARIABLES----------------------------#
indir = args.dir


bucket = args.bucket # AWS S3 Bucket for GOES
satellite = bucket[5:] # get the last part of the bucket name which contains satellite name (goes16 or goes17)
# Specify time range to search for GOES ABI imagery (year, month, days)
year = args.year
month = args.month
start_day = args.days[0]
stop_day = args.days[1]
days=np.linspace(start_day,stop_day,stop_day-start_day+1,dtype=np.int16) # create list of days from start to stop date
hours=['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23'] # all hours (can we just use linspace of ints here instead of list of strings?)
# Specify GOES ABI product, channel, lat/lon bounds, directory path for storing files
product = args.product
channel = args.channel # e.g. 'C14' is the 11.2 micron channel, "Longwave window"
bounds = args.bounds # Define a bounding box to crop to: [Lat_min Lat_max Lon_min Lon_max] (e.g. 30, 50, -125, -105 for western CONUS)
storage_path = args.dir # Local path where data will be stored (to do: check if files already exist in this directory)

# Show us the bounds we'll crop images to
print('\nFiles will be downloaded and then cropped to these bounds:')
print('\t({w},{n}).\t.({e},{n})\n\n\n\n\t({w},{s}).\t.({e},{s})\n'.format(n=bounds[1],w=bounds[2],e=bounds[3],s=bounds[0]))


##############################################################
# For each S3 bucket, download the corresponding observations if we don't have them already
filepath = []; # store filepaths of the files we download
print('For each S3 bucket, download the corresponding observations')
i = 0
for d in range(len(days)):
    for h in range(len(hours)):
        filepath.append('{}/{}/{}/{}/{}/{}/{}/{}/'.format(storage_path,satellite,year,month,days[d],product,hours[h],channel))
        if not os.path.exists(filepath[i]):
            ABI = ABI_Downloader(storage_path,bucket,year,month,days[d],hours[h],product,channel)
        
            # now try and crop these so they don't take up so much space - this is very inefficient but oh well it's what I have right now
            if os.path.exists(filepath[i]): # we have to make sure the path exists (meaning we downloaded something) before running the subsetNetCDF function
                print('\nSubsetting files in...{}'.format(filepath[i]))
                subsetNetCDF(filepath[i],bounds)
        i+=1



print("Done")
