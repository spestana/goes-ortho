"""Functions for downloading data"""

import datetime as dt
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path
from typing import List, Tuple, Union

import geojson
import goes2go as g2g
import numpy as np
import xarray as xr
import zarr
from dateutil import parser, rrule
from shapely.geometry import shape
from tqdm import tqdm

import goes_ortho as go


def build_zarr(download_request_filepath, downloader="goes2go"):
    """
    For running through github actions. Download and build a zarr file from a download request json file
    """
    # download requested imagery
    logging.info("download requested imagery")
    if downloader == "goespy":
        image_path_list = download_abi_goespy(download_request_filepath)
    if downloader == "goes2go":
        image_path_list = download_abi_goes2go(download_request_filepath)

    # parse json request file
    logging.info("parse json request file")
    _, _, bounds, _, _, _, _, variables, apiKey, _, outputFilepath = parse_json(
        download_request_filepath
    )

    # orthorectify all images
    logging.info("orthorectify all images")
    new_image_path_list = []
    logging.info(image_path_list)
    for goes_image_path in image_path_list:
        logging.info("filename: ", goes_image_path)
        new_goes_filename = goes_image_path.with_name(
            goes_image_path.stem + "_o"
        ).with_suffix(".nc")
        # new_goes_filename = goes_image_path.split('.')[:-1][0] + '_o.nc'
        logging.info("renamed to: ", new_goes_filename)
        new_image_path_list.append(new_goes_filename)
        go.orthorectify.ortho(
            goes_image_path, variables, bounds, apiKey, new_goes_filename, keep_dem=True
        )
    logging.info(new_image_path_list)
    # add time dimension, fix CRS, build zarr file
    logging.info("add time dimension, fix CRS, build zarr file")
    for variable in variables:
        logging.info("add_datetime_crs")
        new_image_path_list, datetimes_list = add_datetime_crs(
            new_image_path_list, variable
        )

        # start Dask cluster
        logging.info("start Dask cluster")
        _ = go.io.dask_start_cluster(
            workers=6,
            threads=2,
            open_browser=False,
            verbose=True,
        )

        logging.info(new_image_path_list)
        # nc_files = sorted(
        #    new_image_path_list,
        #    key=datetimes_list
        # )
        nc_files = [
            img_path
            for img_dt, img_path in sorted(zip(datetimes_list, new_image_path_list))
        ]
        logging.info(nc_files)
        # Open all the raster files as a single dataset (combining them together)
        # Why did we choose chunks = 500? 100MB?
        # https://docs.xarray.dev/en/stable/user-guide/dask.html#optimization-tips
        logging.info("open all rasters")
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
        logging.info("rechunk")
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
        logging.info("saving zarr file")
        ds.to_zarr(outputFilepath)
        # Display
        source_group = zarr.open(outputFilepath)
        source_array = source_group[variable]
        logging.info(source_group.tree())
        logging.info(source_array.info)
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
) -> str:
    """
    For running through github actions, make a request json file from github user input to be read by the build_zarr function
    """
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
    return filename


def get_start_date_from_abi_filename(s: str) -> str:
    """
    Read the start time from an ABI filename string.

    Parameters
    ------------
    s : str
        filename for an ABI image product

    Returns
    -----------
    time_str : str
        start time from the filename
    """
    time_str = s.split("_s")[1].split("_")[0]
    return time_str


def add_datetime_crs(
    files: List[str], variable: str, crs: str = "EPSG:4326"
) -> Tuple[List[str], List[dt.datetime]]:
    """
    Add datetime and CRS information to a GOES ABI image product file
    """
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


def parse_json(download_request_filepath: str) -> Tuple:
    """
    Parses a json file that specifies what should be downloaded, created with make_request_json()

    Parameters
    ------------
    download_request_filepath : str
        file path to a download request json file, created with make_request_json()

    Returns
    -----------
    parsed : Tuple
        Tuple of parsed information (startDatetime, endDatetime, bounds, satellite, bucket, product, channels, variables, apiKey, outDir, outputFilepath)

    Examples
    ------------
    (startDatetime, endDatetime, bounds, satellite, bucket, product, channels, variables, apiKey, outDir, outputFilepath) = parse_json("my_download_request.json")
    """
    # load json file that specifies what we'd like to download and parse its contents
    with open(download_request_filepath, "r") as f:
        download_request = json.load(f)

    startDatetime = parser.parse(download_request["dateRange"]["startDatetime"])
    endDatetime = parser.parse(download_request["dateRange"]["endDatetime"])
    bounds = [
        download_request["bounds"]["min_lon"],
        download_request["bounds"]["min_lat"],
        download_request["bounds"]["max_lon"],
        download_request["bounds"]["max_lat"],
    ]  # bounds = [min_lon, min_lat, max_lon, max_lat]
    satellite = download_request["satellite"]
    bucket = "noaa-" + satellite
    product = download_request["product"]
    channels = download_request["bands"]
    variables = download_request["variables"]
    apiKey = download_request["apiKey"]
    outDir = download_request["downloadDirectory"]
    outputFilepath = download_request["outputFilepath"]

    parsed = (
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

    return parsed


def download_abi_goespy(download_request_filepath: str) -> List[Path]:
    """
    Download GOES ABI imagery as specified by an input JSON file using the goespy package for downloading.

    Parameters
    ------------
    download_request_filepath : str
        file path to a download request json file, created with make_request_json()

    Returns
    -----------
    output_filepaths : List[Path]
        path to output DEM (useful if the downloaded DEM is reprojected to custom proj)

    Examples
    ------------
    filepaths = download_abi_goespy("my_download_request.json")
    """

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
    ) = parse_json(download_request_filepath)

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
            # print(this_filepath)
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


def download_abi_goes2go(download_request_filepath: str) -> List[Path]:
    """
    Download GOES ABI imagery as specified by an input JSON file using the goes2go package for downloading.

    Parameters
    ------------
    download_request_filepath : str
        file path to a download request json file, created with make_request_json()

    Returns
    -----------
    output_filepaths : List[Path]
        path to output DEM (useful if the downloaded DEM is reprojected to custom proj)

    Examples
    ------------
    filepaths = download_abi_goes2go("my_download_request.json")
    """

    (
        start_datetime,
        end_datetime,
        bounds,
        satellite,
        bucket,
        product_and_domain,
        channels,
        _,
        _,
        outDir,
        _,
    ) = parse_json(download_request_filepath)

    output_filepaths = []

    # get satellite integer (e.g. the 16, 17, or 18 part of goes16, goes17, goes18
    satellite_int = int(re.search(r"\d+", satellite).group())
    # separate product and domain (C, F, M)
    product = product_and_domain[:-1]
    domain = product_and_domain[-1]
    # check domain
    valid_domains = ["M", "C", "F"]
    assert any(
        domain in d for d in valid_domains
    ), "Invalid GOES-R ABI product name. Product name should end with a domain identifier: M, C, or F. (e.g. ABI-L2-LSTC)"

    this_start_datetime = start_datetime

    # step through and download n hours at a time
    n = 3  # try 3 hours at a time
    n_steps = int(
        np.ceil((end_datetime - start_datetime).total_seconds() / (3600 * n))
    )  # how many steps we'll have to take
    print(f"Estimated {n_steps} batches to download")
    batch_number = 1
    while this_start_datetime < end_datetime:
        print(f"Batch number {batch_number}")
        this_end_datetime = this_start_datetime + dt.timedelta(hours=n)
        print(
            f"Download batch of imagery from {this_start_datetime} to {this_end_datetime}"
        )

        # set up goes2go object for this satellite, product, and domain
        G = g2g.GOES(satellite=satellite_int, product=product, domain=domain)

        # see what is available to download
        # df = G.df(start=startDatetime.strftime("%Y-%m-%d %H:%M"), end=endDatetime.strftime("%Y-%m-%d %H:%M"))

        try:
            # download this batch
            df = G.timerange(
                start=this_start_datetime.strftime("%Y-%m-%d %H:%M"),
                end=this_end_datetime.strftime("%Y-%m-%d %H:%M"),
            )
        except FileNotFoundError as e:
            print(
                f"FileNotFoundError encountered. The requested image may not exist. Because this searched a time window of {n} hours, there may be some valid imagery within the time window. Try a smaller time window to search for valid imagery.\n{e}"
            )
            # set the next start datetime to the end of this batch
            this_start_datetime = this_end_datetime
            batch_number += 1
            continue

        # get the filepaths of this batch
        batch_filepaths = [
            Path(g2g.config["default"]["save_dir"], filepath) for filepath in df["file"]
        ]

        # now try and crop these so they don't take up so much space - this is very inefficient but oh well it's what I have right now
        print(f"Cropping image batch to {bounds}")
        pbar = tqdm(total=len(batch_filepaths))
        for _, this_filepath in enumerate(batch_filepaths):
            if Path.exists(
                this_filepath
            ):  # we have to make sure the path exists (meaning we downloaded something) before running the subsetNetCDF function
                logging.info("\nSubsetting file {}".format(this_filepath))
                go.clip.subsetNetCDF(this_filepath, bounds)
                pbar.update(1)

        pbar.close()

        # append the batch filepaths to our main filepaths list
        output_filepaths.append(batch_filepaths)

        # set the next start datetime to the end of this batch
        this_start_datetime = this_end_datetime

        batch_number += 1

    # flatten the list of lists that we've compiled
    output_filepaths = [
        filepath for batch_filepaths in output_filepaths for filepath in batch_filepaths
    ]

    print("Done")
    return output_filepaths


def get_dem(demtype, bounds, api_key, out_fn=None, proj="EPSG:4326"):
    """
    download a DEM of choice from OpenTopography World DEM

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

    # written by Shashank Bhushan, David Shean

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
    """
    Call a system command through the subprocess python module.
    """
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
    # based on https://github.com/GlacioHack/xdem/blob/d91bf1cc9b3f36d77f3729649bc8e9edc6b42f9f/xdem/examples.py#L33
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
    """
    Remove example files downloaded by download_example_data().
    """
    tmp_dir = "./tests/resources/"
    shutil.rmtree(tmp_dir)


def bounds_from_geojson(geojson_filepath: str) -> Union[List[List[float]], List[float]]:
    """
    Get a bounding box around a polygon from a geojson file
    """

    with open(geojson_filepath) as f:
        geojson_data = geojson.load(f)
        bounds = []
        for feature in geojson_data["features"]:
            bounds.append(shape(feature["geometry"]).bounds)  # [minx, miny, maxx, maxy]

        if len(bounds) > 1:
            print(
                "geojson file contains more than one feature, returning a list of bounds, one set of bounds per feature"
            )
        if len(bounds) == 1:
            bounds = bounds[
                0
            ]  # if there is only one set of bounds, remove it from the list

    return bounds
