"""
Script to orthorectify GOES-R ABI images that correspond with specific ASTER observation times using a DEM
"""

#-------------------------------------------------------#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import os
import goes_ortho

#-------------------------------------------------------#
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
    
def nearest(items, pivot):
    # https://stackoverflow.com/questions/32237862/find-the-closest-date-to-a-given-date
    return min(items, key=lambda x: abs(x - pivot))



#-------------------------------------------------------#

'''Specify directories for inputs and outputs'''

# directory containing ASTER L1T geotiffs, we'll grab the datetime of each from their filenames
aster_directory = "/storage/spestana/ASTER/AST_L1T/geotiff/T/T_band14_Tuolumne-and-CUES"

# the top level directory that contains all the GOES-16 ABI files we want to search
goes_directory = "/storage/GOES/goes16"

# the top level directory that we will output all orthorectified GOES images to
# folder structure within will mirror our source goes_directory
output_directory = "/storage/GOES/orthorectified"

# pick a DEM to orthorectify with
dem_filename = 'dem\dem2.tif'
# The DEM I'm using here is merged/cropped around the upper Tuolumne River basin 
# and includes Mammoth Mountain, from four tiles of SRTM (1 Arc-Second Global) 
# retrieved from USGS EarthExplorer.

#-------------------------------------------------------#

'''Get the times of all our ASTER observations (in UTC)'''

# Find all our ASTER files
aster_files = getListOfFiles(aster_directory)

## Parse the date and time from ASTER filenames
aster_datetimes = []
aster_datetimes_UTC = []
for fpath in aster_files:
    fn = fpath.split('\\')[-1] # non-re method
    MM = fn.split('_')[2][3:5]
    DD = fn.split('_')[2][5:7]
    YYYY = fn.split('_')[2][7:11]
    hh = fn.split('_')[2][11:13]
    mm = fn.split('_')[2][13:15]
    ss = fn.split('_')[2][15:17]
    aster_datetimes_UTC.append(pd.Timestamp('{}-{}-{} {}:{}:{}'.format(YYYY, MM, DD, hh, mm, ss),tz='UTC'))
    aster_datetimes.append(pd.Timestamp('{}-{}-{} {}:{}:{}'.format(YYYY, MM, DD, hh, mm, ss),tz='UTC')- pd.Timedelta(8, unit='hours'))

aster = pd.DataFrame({'datetime': aster_datetimes, 'datetimeUTC': aster_datetimes_UTC, 'filepath': aster_files})
aster.sort_values('datetime',inplace=True)
aster.reset_index(inplace=True, drop=True)

#-------------------------------------------------------#

'''Find GOES ABI images that are closest to these ASTER observations and orthorectify'''

aster_counter = 0
# for every ASTER datetime (in UTC)
for aster_datetime_UTC in aster.datetimeUTC:
    # count
    aster_counter += 1
    print('\n\n File {} of {}'.format(aster_counter, aster.shape[0])
    print('\nFor ASTER observation at {}'.format(aster_datetime_UTC))
    # find the GOES subdirectory for the corresponding year-month-day and hour
    goes_subdir = r"/{year}/{month}/{day}/{product}/{hour}/{channel}/".format(
                            goes_directory=goes_directory, 
                            year=aster_datetime_UTC.strftime('%Y'), 
                            month=aster_datetime_UTC.strftime('%m'), 
                            day=aster_datetime_UTC.day, 
                            product='ABI-L1b-RadC', 
                            hour=aster_datetime_UTC.strftime('%H'), 
                            channel='C14')
    # now within this subdirectory, the same hour of this ASTER observation
    print('\nSearching for GOES ABI imagery within:\n{}{}'.format(goes_directory,goes_subdir))
    # get the filenames of each GOES ABI image in this subdirectory
    goes_files = getListOfFiles(os.path.normpath(goes_directory+goes_subdir))
    goes_datetimes_UTC_list = []
    goes_files_list = []
    
    # create an empty dictionary we'll fill with filenames and timestamps for each subdirectory we search
    goes_dict = {}
    
    for this_goes_file in goes_files:
        this_goes_filename = this_goes_file.split('\\')[-1]
        #print('\t{}'.format(this_goes_filename))
        
        # parse the timstamp in the filename 
        this_goes_datetime_UTC = this_goes_filename.split('_')[-1].split('.')[0][1:-1]
        this_goes_datetime_UTC = pd.to_datetime(this_goes_datetime_UTC, format="%Y%j%H%M%S")
        this_goes_datetime_UTC = pd.Timestamp(this_goes_datetime_UTC, tz='UTC')
        #print('\t{}'.format(this_goes_datetime_UTC))
        
        # add these to our dictionary, use the date as the key
        goes_dict[this_goes_datetime_UTC] = {}
        goes_dict[this_goes_datetime_UTC]['filepath'] = this_goes_file
    
    # now find the one closest to our ASTER observation
    nearest_goes_datetime_UTC = nearest(list(goes_dict.keys()), aster_datetime_UTC)
    #print('\t{} -- {}'.format(aster_datetime_UTC, nearest_goes_datetime_UTC))
    nearest_goes_filepath = goes_dict[nearest_goes_datetime_UTC]['filepath']
    #print(goes_dict[nearest_goes_datetime_UTC]['filepath'])
    print('\n\tFound nearest GOES ABI image:\n\t\tASTER datetime:\t{}\n\t\tGOES datetime:\t{}\n\t\tGOES filepath:\t{}'.format(
            aster_datetime_UTC, nearest_goes_datetime_UTC, nearest_goes_filepath))
    
    # create the output directory if it does not already exist
    output_subdir = r"{}{}".format(output_directory,goes_subdir)
    print('\n\tPreparing to output files to:\n\t{}'.format(output_subdir))
    if not os.path.exists(output_subdir):
        os.makedirs(output_subdir)
        
    # create a new filename for the orthorectified image
    new_file_name = nearest_goes_filepath.split('\\')[-1].split('.')[0] + '_orthorectified'
    print('\n\tNew files will be called:\n\t{}.*'.format(new_file_name))
    
    # Generate the pixel mapping that relates GOES ABI pixels to points on the DEM surface
    pixel_map = goes_ortho.make_ortho_map(nearest_goes_filepath, dem_filename)
    
    # Apply the pixel mapping to orthorectify the GOES ABI image
    ds = goes_ortho.orthorectify_abi_rad(nearest_goes_filepath, 
                                     pixel_map, 
                                     out_filename='{}{}.nc'.format(output_subdir,new_file_name))
    
    # Save a copy of the new orthorectified GOES ABI image (brightness temperature) as a GeoTIFF
    print('\nSave a copy as a brightness temperature GeoTIFF')
    new_nc_file_path =  '{}{}.nc'.format(output_subdir,new_file_name)
    new_gtiff_file_path = '{}{}.tif'.format(output_subdir,new_file_name)
    !gdal_translate -a_srs EPSG:4326 -of GTiff NETCDF:$new_nc_file_path:tb $new_gtiff_file_path



print('\n\n\nALL DONE!')