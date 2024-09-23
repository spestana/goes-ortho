import argparse
import os
from glob import glob

import numpy as np
from goespy.Downloader import (
    ABI_Downloader,  # https://github.com/palexandremello/goes-py
)

from goes_ortho.clip import subsetNetCDF

##############################################################

# ---------------------------- COMMAND LINE ARGUMENTS ----------------------------#
# Set up argument and error handling
parser = argparse.ArgumentParser(
    description="Produces a timeseries of GOES ABI Radiance observations for a single location given a directory of GOES ABI files"
)
parser.add_argument(
    "-B",
    "--bucket",
    required=True,
    type=str,
    help="AWS S3 Bucket for GOES (e.g. noaa-goes16",
)
parser.add_argument(
    "-Y",
    "--year",
    required=True,
    type=int,
    help="Specify time range to search for GOES ABI imagery (year)",
)
parser.add_argument(
    "-M",
    "--month",
    required=True,
    type=int,
    help="Specify time range to search for GOES ABI imagery (month)",
)
parser.add_argument(
    "-D",
    "--days",
    required=True,
    type=int,
    nargs=2,
    help="Specify time range to search for GOES ABI imagery (start day, stop day)",
)
parser.add_argument(
    "-p",
    "--product",
    required=True,
    type=str,
    help="GOES ABI Product (e.g. ABI-L1b-RadC)",
)
parser.add_argument(
    "-c", "--channel", required=True, type=str, help="GOES ABI channel/band (e.g. C14)"
)
parser.add_argument(
    "-b",
    "--bounds",
    required=True,
    type=float,
    nargs=4,
    help="Bounds to crop GOES ABI image to (min_lat max_lat min_lon max_lon",
)
parser.add_argument(
    "-d", "--dir", required=True, help="Directory to save GOES ABI files (.nc)"
)
args = parser.parse_args()

# -----------------------SET ARGUMENTS TO VARIABLES----------------------------#
indir = args.dir


bucket = args.bucket  # AWS S3 Bucket for GOES
satellite = bucket[
    5:
]  # get the last part of the bucket name which contains satellite name (goes16 or goes17)
# Specify time range to search for GOES ABI imagery (year, month, days)
year = args.year
month = args.month
start_day = args.days[0]
stop_day = args.days[1]
days = np.linspace(
    start_day, stop_day, stop_day - start_day + 1, dtype=np.int16
)  # create list of days from start to stop date
hours = [
    "00",
    "01",
    "02",
    "03",
    "04",
    "05",
    "06",
    "07",
    "08",
    "09",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "17",
    "18",
    "19",
    "20",
    "21",
    "22",
    "23",
]  # all hours (can we just use linspace of ints here instead of list of strings?)
# Specify GOES ABI product, channel, lat/lon bounds, directory path for storing files
product = args.product
channel = args.channel  # e.g. 'C14' is the 11.2 micron channel, "Longwave window"
bounds = args.bounds  # Define a bounding box to crop to: [Lat_min Lat_max Lon_min Lon_max] (e.g. 30, 50, -125, -105 for western CONUS)
storage_path = args.dir  # Local path where data will be stored (to do: check if files already exist in this directory)

# Show us the bounds we'll crop images to
print("\nFiles will be downloaded and then cropped to these bounds:")
print(
    "\t({w},{n}).\t.({e},{n})\n\n\n\n\t({w},{s}).\t.({e},{s})\n".format(
        n=bounds[1], w=bounds[2], e=bounds[3], s=bounds[0]
    )
)


##############################################################
# For each S3 bucket, download the corresponding observations if we don't have them already
filepath = []  # store filepaths of the files we download
print("For each S3 bucket, download the corresponding observations")
i = 0
for d in range(len(days)):
    for h in range(len(hours)):
        filepath.append(
            "{}/{}/{}/{}/{}/{}/{}/{}/".format(
                storage_path,
                satellite,
                year,
                month,
                days[d],
                product,
                hours[h],
                channel,
            )
        )
        if not os.path.exists(filepath[i]):
            ABI = ABI_Downloader(
                storage_path, bucket, year, month, days[d], hours[h], product, channel
            )

            # now try and crop these so they don't take up so much space - this is very inefficient but oh well it's what I have right now
            if os.path.exists(
                filepath[i]
            ):  # we have to make sure the path exists (meaning we downloaded something) before running the subsetNetCDF function
                print("\nSubsetting files in...{}".format(filepath[i]))
                for file in glob(filepath[i] + "*.nc"):
                    subsetNetCDF(file, bounds)
        i += 1


print("Done")
