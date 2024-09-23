import logging
import re
import shutil
import webbrowser
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from pathlib import Path
from subprocess import PIPE, STDOUT, Popen

import fsspec
import rioxarray
import xarray as xr
import zarr
from dask.distributed import Client, LocalCluster
from rasterio.enums import Resampling

"""
Basic io functions.
"""


@contextmanager
def redirect_stdout_stderr(stdout_fn=None, stderr_fn=None):
    """
    Writes stdout and/or stderr to file.
    Use os.devnull as the file name to silence entirely.
    """
    Path(stdout_fn).parent.mkdir(parents=True, exist_ok=True)
    Path(stderr_fn).parent.mkdir(parents=True, exist_ok=True)

    if stdout_fn and stderr_fn:
        with open(stdout_fn, "w") as stdout, open(stderr_fn, "w") as stderr:
            with redirect_stdout(stdout) as out, redirect_stderr(stderr) as err:
                yield (err, out)
    elif stdout_fn:
        with open(stdout_fn, "w") as stdout:
            with redirect_stdout(stdout) as out:
                yield out
    elif stderr_fn:
        with open(stderr_fn, "w") as stderr:
            with redirect_stdout(stderr) as err:
                yield err


def parse_urls_from_S3_bucket(
    s3_bucket_name,
    aws_server_url="s3.amazonaws.com",
    folder="",
    extension="tif",
):
    fs = fsspec.filesystem("s3", anon=True)
    bucket = "s3://" + Path(s3_bucket_name, folder).as_posix()
    base_url = Path(s3_bucket_name + "." + aws_server_url, folder).as_posix()
    file_names = [x.split("/")[-1] for x in fs.ls(bucket) if extension in x]
    urls = ["http://" + Path(base_url, x).as_posix() for x in file_names]

    return urls


def _get_test_sites(s3_bucket_name):
    fs = fsspec.filesystem("s3", anon=True)
    bucket = "s3://" + s3_bucket_name
    sites = [x.split("/")[-1] for x in fs.ls(bucket)]
    return sites


def run_command(command, verbose=True):
    """
    Run something from the command line.

    Example 1
    call = ['command', 'input', output']
    run_command(call)

    Example 2
    call = 'command "input" output'
    run_command(call)

    Use a space separated string if your command contains nested strings.
    """

    if isinstance(command, type(str())):
        if verbose:
            print(command)
        shell = True
    else:
        if verbose:
            print(*command)
        shell = False

    p = Popen(command, stdout=PIPE, stderr=STDOUT, shell=shell)

    while p.poll() is None:
        if verbose:
            try:
                line = (p.stdout.readline()).decode("ASCII").rstrip("\n")
            except Exception:
                line = p.stdout.read()
                pass
            print(line)


def parse_timestamps(
    file_list,
    date_string_pattern="....-..-..",
):
    tmp = re.compile(date_string_pattern)
    results = []
    for x in file_list:
        try:
            results.append(tmp.search(x).group(0))
        except AttributeError:
            print("pattern not found in", x)
    return results


def dask_start_cluster(
    workers,
    threads=1,
    ip_address=None,
    port=":8786",
    open_browser=False,
    verbose=True,
):
    """
    Starts a dask cluster. Can provide a custom IP or URL to view the progress dashboard.
    This may be necessary if working on a remote machine.
    """
    cluster = LocalCluster(
        n_workers=workers,
        threads_per_worker=threads,
        silence_logs=logging.ERROR,
        dashboard_address=port,
    )

    client = Client(cluster)

    if ip_address:
        if ip_address[-1] == "/":
            ip_address = ip_address[:-1]  # remove trailing '/' in case it exists
        port = str(cluster.dashboard_link.split(":")[-1])
        url = ":".join([ip_address, port])
        if verbose:
            print("\n" + "Dask dashboard at:", url)
    else:
        if verbose:
            print("\n" + "Dask dashboard at:", cluster.dashboard_link)
        url = cluster.dashboard_link

    if port not in url:
        if verbose:
            print("Port", port, "already occupied")

    if verbose:
        print("Workers:", workers)
        print("Threads per worker:", threads, "\n")

    if open_browser:
        webbrowser.open(url, new=0, autoraise=True)

    return client


def dask_get_mapped_tasks(dask_array):
    """
    Finds tasks associated with chunked dask array.
    """
    # TODO There has to be a better way to do this...
    txt = dask_array._repr_html_()
    idx = txt.find("Tasks")
    strings = txt[idx - 20 : idx].split(" ")
    tasks_count = max([int(i) for i in strings if i.isdigit()])
    return tasks_count


def xr_read_geotif(geotif_file_path, chunks="auto", masked=True):
    """
    Reads in single or multi-band GeoTIFF as dask array.
    Inputs
    ----------
    GeoTIFF_file_path : GeoTIFF file path
    Returns
    -------
    ds : xarray.Dataset
        Includes rioxarray extension to xarray.Dataset
    """

    da = rioxarray.open_rasterio(geotif_file_path, chunks=chunks, masked=True)

    # Extract bands and assign as variables in xr.Dataset()
    ds = xr.Dataset()
    for i, v in enumerate(da.band):
        da_tmp = da.sel(band=v)
        da_tmp.name = "band" + str(i + 1)

        ds[da_tmp.name] = da_tmp

    # Delete empty band coordinates.
    # Need to preserve spatial_ref coordinate, even though it appears empty.
    # See spatial_ref attributes under ds.coords.variables used by rioxarray extension.
    del ds.coords["band"]

    # Preserve top-level attributes and extract single value from value iterables e.g. (1,) --> 1
    ds.attrs = da.attrs
    for key, value in ds.attrs.items():
        try:
            if len(value) == 1:
                ds.attrs[key] = value[0]
        except TypeError:
            pass

    return ds


def xr_stack_geotifs(  # noqa: C901
    geotif_files_list,
    datetimes_list,
    reference_geotif_file,
    resampling="bilinear",
    save_to_nc=False,
    nc_out_dir=None,
    overwrite=True,
    cleanup=False,
    verbose=True,
):
    """
    Stack single or multi-band GeoTiFFs to reference_geotiff.
    Returns out-of-memory dask array, unless resampling occurs.

    Optionally, set save_to_nc true when resmapling is required to
    return an out-of-memory dask array.
    Inputs
    ----------
    geotif_files_list     : list of GeoTIFF file paths
    datetimes_list        : list of datetime objects for each GeoTIFF
    reference_geotif_file : GeoTIFF file path
    Returns
    -------
    ds : xr.Dataset()
    """
    ## TODO: Parameterize crs, res, bounds, transform
    ## TODO: rewrite with dask delayed https://tutorial.dask.org/03_dask.delayed.html

    if save_to_nc and nc_out_dir:
        nc_out_dir = Path(nc_out_dir)
        nc_out_dir.mkdir(parents=True, exist_ok=True)

    ## Check each geotiff has a datetime associated with it.
    if len(datetimes_list) == len(geotif_files_list):
        pass
    else:
        print("length of datetimes does not match length of GeoTIFF list")
        print("datetimes:", len(datetimes_list))
        print("geotifs:", len(geotif_files_list))
        return None

    ## Choose resampling method. Defaults to bilinear.
    if isinstance(resampling, type(Resampling.bilinear)):
        resampling = resampling
    elif resampling == "bilinear":
        resampling = Resampling.bilinear
    elif resampling == "nearest":
        resampling = Resampling.nearest
    elif resampling == "cubic":
        resampling = Resampling.cubic
    else:
        resampling = Resampling.bilinear

    ## Get target object with desired crs, res, bounds, transform
    ref = xr_read_geotif(reference_geotif_file)

    ## Stack geotifs and dimension in time
    datasets = []
    nc_files = []
    out_dirs = []

    c = 0
    for index, file_name in enumerate(geotif_files_list):
        if not nc_out_dir:
            out_fn = str(Path(file_name).with_suffix("")) + ".nc"
        else:
            out_fn = str(Path(nc_out_dir, Path(file_name).with_suffix("").name + ".nc"))

        if Path(out_fn).exists() and not overwrite:
            nc_files.append(out_fn)
            out_dir = str(Path(out_fn).parents[0])
            out_dirs.append(out_dir)
            src = xr.open_dataset(out_fn, chunks="auto")
            datasets.append(src)

        else:
            Path(out_fn).unlink(missing_ok=True)
            src = xr_read_geotif(file_name)
            #             if not check_xr_rio_ds_match(src, ref):
            src = src.rio.reproject_match(ref, resampling=resampling)
            c += 1
            src = src.assign_coords({"time": datetimes_list[index]})
            src = src.expand_dims("time")
            if save_to_nc:
                src.to_netcdf(out_fn)
                nc_files.append(out_fn)
                out_dir = str(Path(out_fn).parents[0])
                out_dirs.append(out_dir)
            datasets.append(src)

    # check if anything was resampled
    if c != 0:
        if verbose:
            print(
                "Resampled",
                c,
                "of",
                len(geotif_files_list),
                "dems to match reference DEM spatial_ref, crs, transform, bounds, and resolution.",
            )

    # Optionally ensure data are returned as dask array.
    if save_to_nc:
        if verbose:
            print("Reading files from", ",".join([str(i) for i in list(set(out_dirs))]))
        ds = xr.open_mfdataset(nc_files, chunks="auto")
        ds = ds.sortby("time")
        ds.rio.write_crs(ref.rio.crs, inplace=True)
        return ds.chunk("auto", balance=True)

    ds = xr.concat(datasets, dim="time", combine_attrs="no_conflicts")
    ds = ds.sortby("time")
    ds.rio.write_crs(ref.rio.crs, inplace=True)
    return ds.chunk("auto", balance=True)


def check_xr_rio_ds_match(ds1, ds2):
    """
    Checks if spatial attributes, crs, bounds, and transform match.
    Inputs
    ----------
    ds1 : xarray.Dataset with rioxarray extension
    ds2 : xarray.Dataset with rioxarray extension
    Returns
    -------
    bool
    """

    if (
        (ds1["spatial_ref"].attrs == ds2["spatial_ref"].attrs)
        & (ds1.rio.crs == ds2.rio.crs)
        & (ds1.rio.transform() == ds2.rio.transform())
        & (ds1.rio.bounds() == ds2.rio.bounds())
        & (ds1.rio.resolution() == ds2.rio.resolution())
    ):
        return True
    else:
        return False


def create_zarr_stack(  # noqa: C901
    xarray_dataset,
    output_directory="./",
    variable_name="band1",
    zarr_stack_file_name="stack.zarr",
    overwrite=False,
    verbose=True,
    cleanup=False,
):
    ds = xarray_dataset
    crs = ds.rio.crs
    print(crs)

    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)

    zarr_stack_fn = Path(output_directory, zarr_stack_file_name)
    zarr_stack_tmp = Path(output_directory, "stack_tmp.zarr")

    if overwrite:
        shutil.rmtree(zarr_stack_fn, ignore_errors=True)
        shutil.rmtree(zarr_stack_tmp, ignore_errors=True)
    if zarr_stack_fn.exists():
        if cleanup:
            if verbose:
                print("Removing temporary zarr stack")
            shutil.rmtree(zarr_stack_tmp, ignore_errors=True)

        ds = xr.open_dataset(zarr_stack_fn, chunks="auto", engine="zarr")
        if verbose:
            print("Zarr file already exists")
            print("Zarr file info")
            source_group = zarr.open(zarr_stack_fn)
            source_array = source_group[variable_name]
            print(source_group.tree())
            print(source_array.info)
            del source_group
            del source_array

        tc, yc, xc = determine_optimal_chuck_size(ds, verbose=verbose)
        ds = xr.open_dataset(
            zarr_stack_fn, chunks={"time": tc, "y": yc, "x": xc}, engine="zarr"
        )

    else:
        if zarr_stack_tmp.exists():
            shutil.rmtree(zarr_stack_tmp, ignore_errors=True)
        # remove attributes that zarr doesn't like
        try:
            ds = ds.drop(["spatial_ref"])
            for i in ds.data_vars:
                try:
                    del ds[i].attrs["grid_mapping"]
                except Exception:
                    pass
        except Exception:
            pass

        if verbose:
            print("Creating temporary zarr stack")

        ds[variable_name].data = ds[variable_name].data.rechunk(
            {0: "auto", 1: "auto", 2: "auto"}, block_size_limit=1e8, balance=True
        )
        arr = ds[variable_name].data
        t, y, x = arr.chunks[0][0], arr.chunks[1][0], arr.chunks[2][0]
        ds[variable_name].encoding = {"chunks": (t, y, x)}
        ds.to_zarr(zarr_stack_tmp)

        if verbose:
            source_group = zarr.open(zarr_stack_tmp)
            source_array = source_group[variable_name]
            print(source_group.tree())
            print(source_array.info)
            del source_group
            del source_array
            print("Rechunking temporary zarr stack and saving as")
            print(str(zarr_stack_fn))

        arr = ds[variable_name].data.rechunk(
            {0: -1, 1: "auto", 2: "auto"}, block_size_limit=1e8, balance=True
        )
        t, y, x = arr.chunks[0][0], arr.chunks[1][0], arr.chunks[2][0]
        ds = xr.open_dataset(
            zarr_stack_tmp, chunks={"time": t, "y": y, "x": x}, engine="zarr"
        )
        ds[variable_name].encoding = {"chunks": (t, y, x)}
        ds.rio.write_crs(crs, inplace=True)
        ds.attrs["crs"] = crs.to_wkt()
        ds.to_zarr(zarr_stack_fn)

        if verbose:
            print("Rechunked zarr file info")
            source_group = zarr.open(zarr_stack_fn)
            source_array = source_group[variable_name]
            print(source_group.tree())
            print(source_array.info)
            del source_group
            del source_array
        if cleanup:
            if verbose:
                print("Removing temporary zarr stack")
            shutil.rmtree(zarr_stack_tmp, ignore_errors=True)

        tc, yc, xc = determine_optimal_chuck_size(ds, verbose=verbose)
        ds = xr.open_dataset(
            zarr_stack_fn, chunks={"time": tc, "y": yc, "x": xc}, engine="zarr"
        )

    if verbose:
        print("Zarr file at", zarr_stack_fn)

    ds.rio.write_crs(crs, inplace=True)
    ds.attrs["crs"] = crs.to_wkt()
    return ds


def determine_optimal_chuck_size(
    ds, variable_name="band1", x_dim="x", y_dim="y", verbose=True
):
    if verbose:
        print("Dask chunk size:")
    ## set chunk size to 5 MB if single time series array < 1 MB in size
    ## else increase to max of 1 GB chunk sizes.

    time_series_array_size = (
        ds[variable_name]
        .sel(
            {
                x_dim: ds[variable_name][x_dim].values[0],
                y_dim: ds[variable_name][y_dim].values[0],
            }
        )
        .nbytes
    )
    MB = 1048576
    if time_series_array_size < 1e6:
        chunk_size_limit = 2 * MB
    elif time_series_array_size < 1e7:
        chunk_size_limit = 20 * MB
    elif time_series_array_size < 1e8:
        chunk_size_limit = 200 * MB
    else:
        chunk_size_limit = 1000 * MB
    arr = ds[variable_name].data.rechunk(
        {0: -1, 1: "auto", 2: "auto"}, block_size_limit=chunk_size_limit, balance=True
    )
    tc, yc, xc = arr.chunks[0][0], arr.chunks[1][0], arr.chunks[2][0]
    chunksize = ds[variable_name][:tc, :yc, :xc].nbytes / 1e6
    if verbose:
        print("Chunk shape:", "(" + ",".join([str(x) for x in [tc, yc, xc]]) + ")")
        print(
            "Chunk size:",
            ds[variable_name][:tc, :yc, :xc].nbytes,
            "(" + str(round(chunksize, 1)) + "MB)",
        )

    return tc, yc, xc
