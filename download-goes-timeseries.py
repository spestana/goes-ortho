import numpy as np
import pandas as pd
import xarray as xr
import os
import sys, argparse
from goespy.Downloader import ABI_Downloader # https://github.com/palexandremello/goes-py



#---------------------------- FUNCTIONS ----------------------------#
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


#---------------------------- COMMAND LINE ARGUMENTS ----------------------------#
# Set up argument and error handling
parser = argparse.ArgumentParser(description='Downloads a timeseries of GOES ABI radiance observations form AWS and returns a list of the filepaths.')
parser.add_argument('-s','--satellite', required=True, help='goes-16 or goes-17')
parser.add_argument('-b','--band', required=True, help='GOES ABI band number')
parser.add_argument('-t','--times', required=True, nargs=2, help='Start and stop times')
args = parser.parse_args()

#-----------------------SET ARGUMENTS TO VARIABLES----------------------------#
# Make sure that satellite selected is valid
satellite = args.satellite
print(satellite)

# Make sure that band selected is valid
band = args.band
print(band)

# Make sure the start and stop times are properly formatted
starttime = pd.to_datetime(args.times[0])
stoptime = pd.to_datetime(args.times[1])

print(starttime.month)
print(stoptime)


##############################################################
# TODO: How to best let the user input this criteria?
# -s satellite (specifies bucket we want to look at)
# -b band (specified which radiance product we want)
# -t time (start and stop dates, and then assume we want all observations within that interval)
###
# AWS S3 Buckets for GOES-16 and -17
buckets = ['noaa-goes16']
# Specify date, time, product, band (channel)
years=['2017']
months=['04']
start_day = starttime.day
stop_day = stoptime.day
days=np.linspace(start_day,stop_day,stop_day-start_day+1,dtype=np.int16)
hours=['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23']
products=['ABI-L1b-RadC']
channels=['C14'] # 11.2 micron channel
# Local paths where data will be stored
homepath = str(os.path.expanduser("~")) # home directory
filepath = []; # store filepaths of the files we download


##############################################################
# For each S3 bucket, download the corresponding observations if we don't have them already
print('For each S3 bucket, download the corresponding observations')
i = 0
for b in range(len(buckets)):
	satellite = buckets[b][5:] # get the last part of the bucket name
	for p in range(len(products)):
		for c in range(len(channels)):
			for y in range(len(years)):
				for m in range(len(months)):
					for d in range(len(days)):
						for h in range(len(hours)):
							print('Downloading: ',satellite,years[y],months[m],days[d],products[p],hours[h],channels[c])
							filepath.append('{}/{}/{}/{}/{}/{}/{}/{}/'.format(homepath,satellite,years[y],months[m],days[d],products[p],hours[h],channels[c]))
							if not os.path.exists(filepath[i]):
								ABI = ABI_Downloader(buckets[b],years[y],months[m],days[d],hours[h],products[p],channels[c])
							i+=1



##############################################################
# for each path (load all the observations)
print('Load all observations from each path')
file_list = []
for path in filepath:
    file_list.append(getListOfFiles(path))


##############################################################
# flatten into a single list
print('Flatten into single list')
file_list = [item for sublist in file_list for item in sublist]

# TODO:
# if this has been called from another python script, return 'file_list'
# otherwise, if this has been called from the command line, don't return anything