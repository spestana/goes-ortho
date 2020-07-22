"""
Script to convert a bunch of netcdfs to geotiffs wrapping around gdal_translate
Use the separate python script goes-orthorectify-aster.py to run the orthorectification routine for all the dates/times we want.
This script I'll just use to convert the resulting netCDF files to geotiffs if I need them.
"""

#-------------------------------------------------------#

import os

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
	
#-------------------------------------------------------#

directory = r"/storage/GOES/orthorectified"

# get all the files here
filepaths = getListOfFiles(directory)

for i, filepath in enumerate(filepaths):
    
	# only get netcdf files
	if os.path.splitext(os.path.normpath(filepath))[-1] == '.nc':
	
		print('Converting file {} of {}'.format(i+1, len(filepaths)))
		
		# get the orthorectified GOES NetCDF filepath
		netcdf_filepath = os.path.normpath(filepath)
		
		# make a filepath for the new tif file
		geotiff_filepath = '{}.tif'.format(os.path.splitext(netcdf_filepath)[0])
		
		# using gdal in the command line here to convert
		os.system('gdal_translate -a_srs {a_srs} -of GTiff NETCDF:{netcdf_filepath}:{band} {geotiff_filepath}'.format(
			a_srs = 'EPSG:4326',
			netcdf_filepath = netcdf_filepath,
			geotiff_filepath = geotiff_filepath,
			band = 'tb'))