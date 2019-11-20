import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import os
from goespy.Downloader import ABI_Downloader # https://github.com/palexandremello/goes-py

###############################
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
def goesBrightnessTemp(rad): 
    # for GOES emissive bands (only for band 14 right now)
    # TODO: Update this to work with the other emissive bands
    fk1 = 8.51022e+03
    fk2 = 1.28627e+03
    bc1 = 0.22516
    bc2 = 0.99920
    T = ( fk2 / (np.log((fk1 / rad) + 1)) - bc1 ) / bc2
    return T
	
##############################################################
# AWS S3 Bucket for GOES-17
bucket = 'noaa-goes16'
# Specify date, time, product, band (channel)
year='2017'
month='04'
start_day = 19
stop_day = 22
days=np.linspace(start_day,stop_day,stop_day-start_day+1,dtype=np.int16)
hours=['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23']
product='ABI-L1b-RadC'
channel='C14' # 11.2 micron channel
# Local paths where data will be stored
path = r'C:\Users\steve\goes16'
filepath = []; # store filepaths of the files we download


##############################################################

# For each S3 bucket, download the corresponding observations
i = 0
for d in range(len(days)):
    for h in range(len(hours)):
        filepath.append('{}/{}/{}/{}/{}/{}/{}/'.format(path,year,month,days[d],product,hours[h],channel))
        if not os.path.exists(filepath[i]):
            ABI = ABI_Downloader(bucket,year,month,days[d],hours[h],product,channel)
        i+=1



##############################################################
# for each path (load all the observations)
file_list = []
for path in filepath:
    file_list.append(getListOfFiles(path))



# flatten into a single list
file_list = [item for sublist in file_list for item in sublist]


##############################################################
# open this entire dataset as a "multi-file dataset"
g = xr.open_mfdataset(file_list,concat_dim='t')

##############################################################
# Approx location of Gaylor Pit
lat = 37.88 
lon= -119.31
z = 0
print(lon, lat, z)
##############################################################
with xr.open_dataset(file_list[0]) as f:
    # Values needed for geometry calculations
    req = f.goes_imager_projection.semi_major_axis
    rpol = f.goes_imager_projection.semi_minor_axis
    H = f.goes_imager_projection.perspective_point_height + f.goes_imager_projection.semi_major_axis
    lon_0 = f.goes_imager_projection.longitude_of_projection_origin
    e = 0.0818191910435 # GRS-80 eccentricity
    
    # find corresponding look angles
    x_rad, y_rad = LonLat2ABIangle(lon,lat,z,H,req,rpol,e,lon_0)
    print(x_rad, y_rad)

##############################################################
# Take a look at how radiance changes at this single point
values = []
time = []
for filename in file_list:
    with xr.open_dataset(filename) as f:
        # find corresponding pixel Radiance value nearest to these angles
        value = f.Rad.sel(y=y_rad, x=x_rad, method='nearest').values.mean()
        values.append(value)
        start_time = f.time_bounds.values.min()
        time.append(start_time)


##############################################################


# values are in mW for the emissive bands
radiance = np.array(values)
# calculate brightness temperature for this band, convert from K to C
Tb = goesBrightnessTemp(radiance) - 273.15

# make a pandas dataframe out of this
d = {'time': time, 'tb': Tb}
df = pd.DataFrame(d)


##############################################################
# save this informaiton out to a pickle file
df1.to_pickle('GOES_B14_BrightnessTemperature_19-22_Feb_2017_GaylorPit_z.pkl')