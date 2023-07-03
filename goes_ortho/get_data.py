'''Get test data for tests and/or examples'''

# based on https://github.com/GlacioHack/xdem/blob/d91bf1cc9b3f36d77f3729649bc8e9edc6b42f9f/xdem/examples.py#L33

import os
import tarfile
import tempfile
import urllib.request
import shutil

def get_dem(demtype, bounds, api_key, out_fn=None, proj='EPSG:4326'):
    """
    download a DEM of choice from OpenTopography World DEM
    (modified by Shashank Bhushan, first written by David Shean)
    Parameters
    ------------
    demtype: str
        type of DEM to fetch (e.g., SRTMGL1, SRTMGL1_E, SRTMGL3 etc)
    bounds: list
        geographic aoi extent in format (minx,miny,maxx,maxy)
    out_fn: str
        path to output filename
    t_srs: str
        output DEM projection
    Returns
    -----------
    out_DEM: str
        path to output DEM (useful if the downloaded DEM is reprojected to custom proj)
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
    commit = "bc4e02ee93a9a8ca10738df3bcec2b829c838d69"
    # The URL from which to download the tarball
    url = f"https://github.com/spestana/goes-ortho-data/tarball/main#commit={commit}"

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