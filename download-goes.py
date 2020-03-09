import numpy as np
import pandas as pd
import xarray as xr
import os
from goespy.Downloader import ABI_Downloader # https://github.com/palexandremello/goes-py


##############################################################
###############  download-goes-timeseries.py ################# <-- the below should sit on its own
##############################################################

##############################################################
# AWS S3 Bucket for GOES-17
bucket = 'noaa-goes16'
satellite = bucket[5:] # get the last part of the bucket name
# Specify date, time, product, band (channel)
year='2017'
month='01'
start_day = 15
stop_day = 31
days=np.linspace(start_day,stop_day,stop_day-start_day+1,dtype=np.int16)
hours=['00','01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23']
product='ABI-L1b-RadC'
channel='C14' # 11.2 micron channel, "Longwave window"
# Local paths where data will be stored
homepath = '/storage/GOES' #str(os.path.expanduser("~")) # home directory
filepath = []; # store filepaths of the files we download


##############################################################
# For each S3 bucket, download the corresponding observations if we don't have them already
print('For each S3 bucket, download the corresponding observations')
i = 0
for d in range(len(days)):
    for h in range(len(hours)):
        filepath.append('{}/{}/{}/{}/{}/{}/{}/{}/'.format(homepath,satellite,year,month,days[d],product,hours[h],channel))
        if not os.path.exists(filepath[i]):
            ABI = ABI_Downloader(bucket,year,month,days[d],hours[h],product,channel)
        i+=1



print("Done")