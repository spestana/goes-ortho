'''Get test data for tests and/or examples'''

# based on https://github.com/GlacioHack/xdem/blob/d91bf1cc9b3f36d77f3729649bc8e9edc6b42f9f/xdem/examples.py#L33

import os
import tarfile
import urllib.request
import shutil
import subprocess
import sys
import goes_ortho as go
from goespy.Downloader import ABI_Downloader
import json
from glob import glob
from dateutil import rrule, parser


def download_abi(downloadRequest_filepath):
    '''Download GOES ABI imagery as specified by an input JSON file. (this function wraps around goespy.ABIDownloader())'''

    # load json file that specifies what we'd like to download and parse its contents
    with open(downloadRequest_filepath, "r") as f:
        downloadRequest = json.load(f)

    startDatetime = parser.parse(downloadRequest['dateRange']['startDatetime'])
    endDatetime = parser.parse(downloadRequest['dateRange']['endDatetime'])
    bounds = [  downloadRequest['bounds']["min_lat"],
                downloadRequest['bounds']["max_lat"],
                downloadRequest['bounds']["min_lon"],
                downloadRequest['bounds']["max_lon"] 
            ]
    satellite = downloadRequest['satellite']
    bucket = 'noaa-' + satellite
    product = downloadRequest['product']
    outDir = downloadRequest['outputDirectory']

    # parse channels/bands and start a download for each
    for channel in downloadRequest['bands']:
        #print(channel)
        if type(channel) == int:
            channel = 'C{:02}'.format(channel) # correct channel to a string formatted like "C02" if we were provided with an integer channel/band number
            #print(channel)  
        # Show us the bounds we'll crop images to
        print('\nFiles will be downloaded and then cropped to these bounds:')
        print('\t({w},{n}).\t.({e},{n})\n\n\n\n\t({w},{s}).\t.({e},{s})\n'.format(n=bounds[1],w=bounds[2],e=bounds[3],s=bounds[0]))
        # For each S3 bucket, download the corresponding observations if we don't have them already
        filepath = []; # store filepaths of the files we download
        print('For each S3 bucket, download the corresponding observations')
        i = 0
        for dt in rrule.rrule(rrule.HOURLY, dtstart=startDatetime, until=endDatetime):
            filepath.append('{}/{}/{}/{}/{}/{}/{}/{}/'.format(outDir,satellite,dt.year,dt.month,dt.day,product,f'{dt.hour:02}',channel))
            if not os.path.exists(filepath[i]):
                ABI = ABI_Downloader(outDir,bucket,dt.year,dt.month,dt.day,f'{dt.hour:02}',product,channel)
                # now try and crop these so they don't take up so much space - this is very inefficient but oh well it's what I have right now
                if os.path.exists(filepath[i]): # we have to make sure the path exists (meaning we downloaded something) before running the subsetNetCDF function
                    print('\nSubsetting files in...{}'.format(filepath[i]))
                    for file in glob(filepath[i]+'*.nc'):
                        go.clip.subsetNetCDF(file,bounds)
                i+=1
    print("Done")
    return None

def get_dem(demtype, bounds, api_key, out_fn=None, proj='EPSG:4326'):
    """
    download a DEM of choice from OpenTopography World DEM (modified by Shashank Bhushan, first written by David Shean)
    
    Parameters
    ------------
    demtype : str
        type of DEM to fetch (e.g., SRTMGL1, SRTMGL1_E, SRTMGL3 etc)
    bounds : list
        geographic aoi extent in format (minx,miny,maxx,maxy)
    out_fn : str
        path to output filename
    t_srs : str
        output DEM projection
    
    Returns
    -----------
    out_DEM : str
        path to output DEM (useful if the downloaded DEM is reprojected to custom proj)

    Examples
    ------------
    
    """
    import requests
    from distutils.spawn import find_executable
    ### From David Shean
    base_url="https://portal.opentopography.org/API/globaldem?demtype={}&west={}&south={}&east={}&north={}&outputFormat=GTiff&API_Key={}"
    if out_fn is None:
        out_fn = '{}.tif'.format(demtype)
    if not os.path.exists(out_fn):
        #Prepare API request url
        #Bounds should be [minlon, minlat, maxlon, maxlat]
        url = base_url.format(demtype, *bounds, api_key)
        print(url)
        #Get
        response = requests.get(url)
        #Check for 200
        #Write to disk
        open(out_fn, 'wb').write(response.content)
    if proj != 'EPSG:4326':
        #Could avoid writing to disk and direclty reproject with rasterio, using gdalwarp for simplicity
        proj_fn = os.path.splitext(out_fn)[0]+'_proj.tif'
        if not os.path.exists(proj_fn):
            output_res = 30
            gdalwarp = find_executable('gdalwarp')
            gdalwarp_call = f"{gdalwarp} -r cubic -co COMPRESS=LZW -co TILED=YES -co BIGTIFF=IF_SAFER -tr {output_res} {output_res} -t_srs '{proj}' {out_fn} {proj_fn}"
            print(gdalwarp_call)
            run_bash_command(gdalwarp_call)
        out_DEM = proj_fn
    else:
        out_DEM = out_fn
    return out_DEM

def run_bash_command(cmd):
    #written by Scott Henderson
    # move to asp_binder_utils
    """Call a system command through the subprocess python module."""
    print(cmd)
    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode < 0:
            print("Child was terminated by signal", -retcode, file=sys.stderr)
        else:
            print("Child returned", retcode, file=sys.stderr)
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)


def download_example_data() -> None:
    """
    Fetch the GOES ABI example files.
    """

    # Static commit hash (to be updated as needed)
    #commit = "16756d3aff6ca41ebb0be999a82d2f66930e7851"
    # The URL from which to download the tarball
    url = f"https://github.com/spestana/goes-ortho-data/tarball/main" ##commit={commit}"

    # Make resources directory
    tmp_dir = "./tests/resources/"
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
		
    # Path and filename for tarball
    tar_path = tmp_dir + "data.tar.gz"
    if not os.path.exists(tar_path):
        response = urllib.request.urlopen(url)
        # If the response was right, download the tarball
        if response.getcode() == 200:
            with open(tar_path, "wb") as outfile:
                outfile.write(response.read())
        else:
            raise ValueError(f"Example GOES ABI data fetch gave non-200 response: {response.status_code}")
    
        # Extract the tarball
        with tarfile.open(tar_path) as tar:
            tar.extractall(tmp_dir)

def remove_example_data() -> None:
    tmp_dir = "./tests/resources/"
    shutil.rmtree(tmp_dir)
