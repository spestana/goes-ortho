import numpy as np
import pandas as pd
import xarray as xr
import os
from goespy.Downloader import ABI_Downloader # https://github.com/palexandremello/goes-py
import sys, argparse


##############################################################
def getListOfFiles(dirName, channel):
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
            allFiles = allFiles + getListOfFiles(fullPath,channel)
        else:
            if channel in fullPath:
                allFiles.append(fullPath)
                
    return allFiles 


##############################################################
def LonLat2ABIangle(lon_deg, lat_deg, z, H, req, rpol, e, lon_0_deg):
    '''Find the ABI elevation (y) and scanning (x) angles (radians) of point P , given a latitude and longitude (degrees)'''
    
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
    
    # determine if this point is visible to the satellite
    condition = ( H * (H-Sx) ) < ( Sy**2 + (req**2 / rpol**2)*Sz**2 )
    if condition == True:
        print('Point at {},{} not visible to satellite.'.format(lon_deg,lat_deg))
        return (np.nan, np.nan)
    else:
        return (x,y)

##############################################################    
def goesBrightnessTemp(rad, fk1, fk2, bc1, bc2): 
    # for GOES emissive bands
    # TODO: Update this to work with the other emissive bands
    #fk1 = 8.51022e+03
    #fk2 = 1.28627e+03
    #bc1 = 0.22516
    #bc2 = 0.99920
    T = ( fk2 / (np.log((fk1 / rad) + 1)) - bc1 ) / bc2
    return T


#---------------------------- COMMAND LINE ARGUMENTS ----------------------------#
# Set up argument and error handling
parser = argparse.ArgumentParser(description='Produces a timeseries of GOES ABI Radiance observations for a single location given a directory of GOES ABI files')
parser.add_argument('-d','--dir', required=True, help='directory to search for GOES ABI files (.nc)')
parser.add_argument('-c','--channel', required=True, type=str, help='GOES ABI channel/band (e.g. C14)')
parser.add_argument('-l','--loc', required=True, nargs=3, type=float, help='Latitude, Longitude, and elevation of location to extract timeseries for')
parser.add_argument('-o','--ofile', required=True, nargs=1, type=str, help='Output filepath and name')
args = parser.parse_args()

#-----------------------SET ARGUMENTS TO VARIABLES----------------------------#
# Input directory
indir = args.dir

# GOES ABI Channel
channel = args.channel

# Lat, Lon, Elev.
lat_obs = args.loc[0]
lon_obs = args.loc[1]
z_obs = args.loc[2]

# Output filepath/filename
out_filepath_filename = args.ofile[0]

print(out_filepath_filename)



##############################################################
# for each path (load all the observations)
print('Load all observations from the directory provided')
path = indir
file_list = []
try:
    file_list.append(getListOfFiles(path,channel))
except FileNotFoundError:
    print('Could not find files at {}'.format(path))


##############################################################
# flatten into a single list
print('Flatten into single list')
file_list = [item for sublist in file_list for item in sublist]


##############################################################
# open this entire dataset as a "multi-file dataset"
print('Open as multi-file dataset')
g = xr.open_mfdataset(file_list, concat_dim='t', combine='nested')

##############################################################
# Get location from user input
lat = lat_obs
lon= lon_obs
z = z_obs
print('Extracting radiance values for: {}, {}, {}'.format(lon, lat, z))
##############################################################
print('Get values needed for geometry calculations')
with xr.open_dataset(file_list[0]) as f:
    # Values needed for geometry calculations
    req = f.goes_imager_projection.semi_major_axis
    rpol = f.goes_imager_projection.semi_minor_axis
    H = f.goes_imager_projection.perspective_point_height + f.goes_imager_projection.semi_major_axis
    lon_0 = f.goes_imager_projection.longitude_of_projection_origin
    e = 0.0818191910435 # GRS-80 eccentricity
       
    # find corresponding look angles
    x_rad, y_rad = LonLat2ABIangle(lon,lat,z,H,req,rpol,e,lon_0)
    print('Extracting radiance values for: {}, {}'.format(x_rad, y_rad))

##############################################################
# Take a look at how radiance changes at this single point
print('Take a look at how radiance changes at this single point')
values = []
Tb = []
time = []
for filename in file_list:
    with xr.open_dataset(filename) as f:
        # find corresponding pixel Radiance value nearest to these angles
        value = f.Rad.sel(y=y_rad, x=x_rad, method='nearest').values.mean()
        values.append(value)
        # convert Radiance to Tb
        # Values needed for radiometric conversion
        fk1 = f.planck_fk1.values
        fk2 = f.planck_fk2.values
        bc1 = f.planck_bc1.values
        bc2 = f.planck_bc2.values
        # calculate brightness temperature for this band (K)
        Tb.append(goesBrightnessTemp(value, fk1, fk2, bc1, bc2))
        start_time = f.time_bounds.values.min()
        time.append(start_time)


##############################################################

print('Make dataframe')
# values are in mW for the emissive bands
radiance = np.array(values)
# convert from K to C for brightness temperature
Tb = np.array(Tb) - 273.15

# make a pandas dataframe out of this
d = {'time': time, 'rad': radiance, 'tb': Tb}
df = pd.DataFrame(d)
print(df)
print(type(df))

##############################################################
# save this informaiton out to a pickle file
print('save this informaiton out to a pickle file')
print(out_filepath_filename)
df.to_pickle(out_filepath_filename, protocol=3)