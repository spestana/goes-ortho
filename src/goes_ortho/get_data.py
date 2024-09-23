"""Get test data for tests and/or examples"""

# based on https://github.com/GlacioHack/xdem/blob/d91bf1cc9b3f36d77f3729649bc8e9edc6b42f9f/xdem/examples.py#L33

import datetime as dt
import json
import os
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

import xarray as xr
import zarr
from dateutil import parser, rrule

import goes_ortho as go


def build_zarr(downloadRequest_filepath):
    # download requested imagery
    print("download requested imagery")
    image_path_list = download_abi(downloadRequest_filepath)

    # parse json request file
    print("parse json request file")
    _, _, bounds, _, _, _, _, variables, apiKey, _, outputFilepath = parse_json(
        downloadRequest_filepath
    )

    # orthorectify all images
    print("orthorectify all images")
    new_image_path_list = []
    print(image_path_list)
    for goes_image_path in image_path_list:
        print("filename: ", goes_image_path)
        new_goes_filename = goes_image_path.with_name(
            goes_image_path.stem + "_o"
        ).with_suffix(".nc")
        # new_goes_filename = goes_image_path.split('.')[:-1][0] + '_o.nc'
        print("renamed to: ", new_goes_filename)
        new_image_path_list.append(new_goes_filename)
        go.orthorectify.ortho(
            goes_image_path, variables, bounds, apiKey, new_goes_filename, keep_dem=True
        )
    print(new_image_path_list)
    # add time dimension, fix CRS, build zarr file
    print("add time dimension, fix CRS, build zarr file")
    for variable in variables:
        print("add_datetime_crs")
        new_image_path_list, datetimes_list = add_datetime_crs(
            new_image_path_list, variable
        )

        # start Dask cluster
        print("start Dask cluster")
        _ = go.io.dask_start_cluster(
            workers=6,
            threads=2,
            open_browser=False,
            verbose=True,
        )

        print(new_image_path_list)
        # nc_files = sorted(
        #    new_image_path_list,
        #    key=datetimes_list
        # )
        nc_files = [
            img_path
            for img_dt, img_path in sorted(zip(datetimes_list, new_image_path_list))
        ]
        print(nc_files)
        # Open all the raster files as a single dataset (combining them together)
        # Why did we choose chunks = 500? 100MB?
        # https://docs.xarray.dev/en/stable/user-guide/dask.html#optimization-tips
        print("open all rasters")
        ds = xr.open_mfdataset(nc_files, chunks={"time": 500})
        #
        ## if 'Rad' is our variable, check if we should add reflectance 'ref', or brightness temperature 'tb' to the list too
        # if variable == 'Rad':
        #    if ds.band_id.values[0] <= 6:
        #        print('adding ref to variables list')
        #        variables.append('ref')
        #    else:
        #        print('adding tb to variables list')
        #        variables.append('tb')
        #
        # rechunk along time dimension
        # Dask's rechunk documentation: https://docs.dask.org/en/stable/generated/dask.array.rechunk.html
        # 0:-1 specifies that we want the dataset to be chunked along the 0th dimension -- the time dimension, which means that each chunk will have all 40 thousand values in time dimension
        # 1:'auto', 2:'auto' and balance=True specifies that dask can freely rechunk along the latitude and longitude dimensions to attain blocks that have a uniform size
        print("rechunk")
        ds[variable].data.rechunk(
            {0: -1, 1: "auto", 2: "auto"}, block_size_limit=1e8, balance=True
        )
        # Assign the dimensions of a chunk to variables to use for encoding afterwards
        t, y, x = (
            ds[variable].data.chunks[0][0],
            ds[variable].data.chunks[1][0],
            ds[variable].data.chunks[2][0],
        )
        # Create an output zarr file and write these chunks to disk
        # remove if file already exists
        shutil.rmtree(outputFilepath, ignore_errors=True)
        # chunk
        ds[variable].encoding = {"chunks": (t, y, x)}
        # output zarr file
        print("saving zarr file")
        ds.to_zarr(outputFilepath)
        # Display
        source_group = zarr.open(outputFilepath)
        source_array = source_group[variable]
        print(source_group.tree())
        print(source_array.info)
        del source_group
        del source_array
    print("Done.")
    return None


def make_request_json(
    workflowName,
    startDatetime,
    endDatetime,
    bounds,
    satellite,
    product,
    band,
    variable,
    apiKey,
):
    """For running through github actions, make a request json file from github user input to be read by the build_zarr function"""
    request_dict = {
        "dateRange": {"startDatetime": startDatetime, "endDatetime": endDatetime},
        "bounds": {
            "min_lon": bounds[0],
            "min_lat": bounds[1],
            "max_lon": bounds[2],
            "max_lat": bounds[3],
        },
        "satellite": satellite,
        "product": product,
        "bands": [band],
        "variables": [variable],
        "downloadDirectory": "./",
        "outputFilepath": "./{}.zarr".format(workflowName),
        "apiKey": apiKey,
    }
    filename = workflowName + ".json"
    with open(filename, "w") as f:
        json.dump(request_dict, f)


def get_start_date_from_abi_filename(s):
    return s.split("_s")[1].split("_")[0]


def add_datetime_crs(files, variable, crs="EPSG:4326"):
    print(files)
    print(variable)
    print(crs)
    new_files = []
    datetimes = [
        dt.datetime.strptime(
            f.stem.split("_")[3][
                1:-1
            ],  # parse the start time (the part "s2022__________" in the file name)
            "%Y%j%H%M%S",
        )
        for f in files
    ]
    print(datetimes)
    for i, file in enumerate(files):
        print(f"Processing {i} of {len(files)}...")
        print(datetimes[i])
        try:
            ds = xr.open_dataset(file)
            ds = ds.assign_coords({"time": datetimes[i]})
            ds = ds.expand_dims("time")
            ds = ds.reset_coords(drop=True)
            da = ds[variable]
            new_file_name = file.with_name(
                file.stem + "_{}".format(variable)
            ).with_suffix(".nc")
            # new_file_name = file.replace(
            #    ".nc",
            #    "_{}.nc".format(variable),
            # )
            da = da.rio.write_crs(crs)
            da.to_netcdf(new_file_name)
            new_files.append(new_file_name)
        except Exception as err:
            print(f"Failed on {file}")
            print(f"Error: {err}")
    return new_files, datetimes


def parse_json(downloadRequest_filepath):
    # load json file that specifies what we'd like to download and parse its contents
    with open(downloadRequest_filepath, "r") as f:
        downloadRequest = json.load(f)

    startDatetime = parser.parse(downloadRequest["dateRange"]["startDatetime"])
    endDatetime = parser.parse(downloadRequest["dateRange"]["endDatetime"])
    bounds = [
        downloadRequest["bounds"]["min_lon"],
        downloadRequest["bounds"]["min_lat"],
        downloadRequest["bounds"]["max_lon"],
        downloadRequest["bounds"]["max_lat"],
    ]  # bounds = [min_lon, min_lat, max_lon, max_lat]
    satellite = downloadRequest["satellite"]
    bucket = "noaa-" + satellite
    product = downloadRequest["product"]
    channels = downloadRequest["bands"]
    variables = downloadRequest["variables"]
    apiKey = downloadRequest["apiKey"]
    outDir = downloadRequest["downloadDirectory"]
    outputFilepath = downloadRequest["outputFilepath"]

    return (
        startDatetime,
        endDatetime,
        bounds,
        satellite,
        bucket,
        product,
        channels,
        variables,
        apiKey,
        outDir,
        outputFilepath,
    )


def download_abi(downloadRequest_filepath):
    """Download GOES ABI imagery as specified by an input JSON file. (this function wraps around goespy.ABIDownloader())"""

    (
        startDatetime,
        endDatetime,
        bounds,
        satellite,
        bucket,
        product,
        channels,
        _,
        _,
        outDir,
        _,
    ) = parse_json(downloadRequest_filepath)

    output_filepaths = []

    # parse channels/bands and start a download for each
    for channel in channels:
        # print(channel)
        if isinstance(channel, int):
            channel = "C{:02}".format(
                channel
            )  # correct channel to a string formatted like "C02" if we were provided with an integer channel/band number
            # print(channel)
        # Show us the bounds we'll crop images to
        print("\nFiles will be downloaded and then cropped to these bounds:")
        print(
            "\t({w},{n}).\t.({e},{n})\n\n\n\n\t({w},{s}).\t.({e},{s})\n".format(
                n=bounds[3], w=bounds[0], e=bounds[2], s=bounds[1]
            )
        )
        # For each S3 bucket, download the corresponding observations if we don't have them already
        download_filepaths = []  # store filepaths where our downloaded files are
        print("For each S3 bucket, download the corresponding observations")
        i = 0
        for this_datetime in rrule.rrule(
            rrule.HOURLY, dtstart=startDatetime, until=endDatetime
        ):
            if ("L1b-Rad" in product) or ("L2-CMIP" in product):
                this_filepath = (
                    Path(outDir)
                    / satellite
                    / str(this_datetime.year)
                    / str(this_datetime.month)
                    / str(this_datetime.day)
                    / product
                    / "{:02}".format(this_datetime.hour)
                    / channel
                )
            else:
                this_filepath = (
                    Path(outDir)
                    / satellite
                    / str(this_datetime.year)
                    / str(this_datetime.month)
                    / str(this_datetime.day)
                    / product
                    / "{:02}".format(this_datetime.hour)
                )
            print(this_filepath)
            download_filepaths.append(
                this_filepath
            )  #'{}/{}/{}/{}/{}/{}/{}/{}/'.format(outDir,satellite,dt.year,dt.month,dt.day,product,'{:02}'.format(dt.hour),channel)
            if not Path.is_dir(
                download_filepaths[i]
            ):  # os.path.exists(download_filepaths[i]):
                _ = go.Downloader.ABI_Downloader(
                    outDir,
                    bucket,
                    this_datetime.year,
                    this_datetime.month,
                    this_datetime.day,
                    f"{this_datetime.hour:02}",
                    product,
                    channel,
                )
                # now try and crop these so they don't take up so much space - this is very inefficient but oh well it's what I have right now
                if Path.is_dir(
                    download_filepaths[i]
                ):  # we have to make sure the path exists (meaning we downloaded something) before running the subsetNetCDF function
                    print("\nSubsetting files in...{}".format(download_filepaths[i]))
                    for file in download_filepaths[i].glob("*.nc"):
                        print(file)
                        output_filepaths.append(file)
                        go.clip.subsetNetCDF(file, bounds)
                i += 1
    print("Done")
    return output_filepaths


def get_dem(demtype, bounds, api_key, out_fn=None, proj="EPSG:4326"):
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
    from distutils.spawn import find_executable

    import requests

    ### From David Shean
    base_url = "https://portal.opentopography.org/API/globaldem?demtype={}&west={}&south={}&east={}&north={}&outputFormat=GTiff&API_Key={}"
    if out_fn is None:
        out_fn = "{}.tif".format(demtype)
    if not os.path.exists(out_fn):
        # Prepare API request url
        # Bounds should be [minlon, minlat, maxlon, maxlat]
        url = base_url.format(demtype, *bounds, api_key)
        print(url)
        # Get
        response = requests.get(url)
        # Check for 200
        # Write to disk
        open(out_fn, "wb").write(response.content)
    if proj != "EPSG:4326":
        # Could avoid writing to disk and directly reproject with rasterio, using gdalwarp for simplicity
        proj_fn = os.path.splitext(out_fn)[0] + "_proj.tif"
        if not os.path.exists(proj_fn):
            output_res = 30
            gdalwarp = find_executable("gdalwarp")
            gdalwarp_call = f"{gdalwarp} -r cubic -co COMPRESS=LZW -co TILED=YES -co BIGTIFF=IF_SAFER -tr {output_res} {output_res} -t_srs '{proj}' {out_fn} {proj_fn}"
            print(gdalwarp_call)
            run_bash_command(gdalwarp_call)
        out_DEM = proj_fn
    else:
        out_DEM = out_fn
    return out_DEM


def run_bash_command(cmd):
    # written by Scott Henderson
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
    # commit = "16756d3aff6ca41ebb0be999a82d2f66930e7851"
    # The URL from which to download the tarball
    url = "https://github.com/spestana/goes-ortho-data/tarball/main"  ##commit={commit}"

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
            raise ValueError(
                f"Example GOES ABI data fetch gave non-200 response: {response.status_code}"
            )

        # Extract the tarball
        with tarfile.open(tar_path) as tar:
            tar.extractall(tmp_dir)


def remove_example_data() -> None:
    tmp_dir = "./tests/resources/"
    shutil.rmtree(tmp_dir)
