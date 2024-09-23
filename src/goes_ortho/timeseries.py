"""
Functions for extracting timeseries from directories of GOES ABI imagery
"""

import glob

import pandas as pd
import xarray as xr

import goes_ortho as go

# def df_from_zarr(zarrFilepath, variable, point_lat_lon, outFilepath=None):

#     ds = xr.open_dataset(
#         zarrFilepath,
#         chunks={'time': 40785, 'latitude': 50, 'longitude': 50},
#         engine='zarr'
#     )
#     # When we pass in a chunks argument, the dataset opened will be filled with Dask arrays

#     point_timeseries = ds[variable].sel(latitude = point_lat_lon[0], longitude = point_lat_lon[1], method='nearest')

#     # Convert the timeseries into a pandas dataframe and save in a .csv file
#     df = point_timeseries.to_dataframe().drop(columns=['latitude', 'longitude'])

#     if outFilepath != None:
#         df.to_csv(outFilepath)

#     return df


def make_abi_timeseries(directory, product, data_vars, lon, lat, z, outfilepath=None):
    """Given a directory of GOES ABI products, create a timeseries of data variables (specified in data_vars) for a single point (at lon, lat, elevation).
    Returns a pandas dataframe, optional output to a csv file."""

    path = "{directory}/**/*{product}*.nc".format(directory=directory, product=product)
    file_list = glob.glob(path, recursive=True)

    # create empty dataframe to hold the data variables we want plus a timestamp
    df_columns = list(data_vars)
    df_columns.append("time")
    # if Radiance is one of the data variables we are interested in
    if "Rad" in data_vars:
        # create a new column for reflectance (for bands 1-6) or brightness temperature (for band 7-16)
        df_columns.append("ref_or_tb")
    # create the data frame we will populate with values
    df = pd.DataFrame(columns=df_columns)

    print(
        "Creating a timeseries of {data_vars} from {product} at ({lat}, {lon}, {z})".format(
            data_vars=data_vars, product=product, lat=lat, lon=lon, z=z
        )
    )
    print("Reading:")
    for filename in file_list:
        try:
            print("{}".format(filename), end="\r")

            with xr.open_dataset(filename, decode_times=False) as f:
                # I've included "decode_times=False" to this xr.open_dataset because I've encountered some ABI-L2-ACMC files where the timestamp couldn't be read
                # and xarray gave a "ValueError: unable to decode time units 'seconds since 2000-01-01 12:00:00' with the default calendar. Try opening your dataset with decode_times=False."
                # I've also switched which timestamp from the ABI files I'm reading (was f.time_bounds.values.min(), now f.time_coverage_start)

                # Read goes_imager_projection values needed for geometry calculations
                # and compute the corresponding look angles (in radiance) for the lat, lon, elevation we are interested in
                x_rad, y_rad = go.geometry.LonLat2ABIangle(
                    lon,
                    lat,
                    z,
                    f.goes_imager_projection.perspective_point_height
                    + f.goes_imager_projection.semi_major_axis,
                    f.goes_imager_projection.semi_major_axis,
                    f.goes_imager_projection.semi_minor_axis,
                    0.0818191910435,  # GRS-80 eccentricity
                    f.goes_imager_projection.longitude_of_projection_origin,
                )

                # get the timestamp for this observation (these should all be UTC, but I am removing timezone info because not all timestamps are converting the same way, and I was getting a "Cannot compare tz-naive and tz-aware timestamps" error)
                timestamp = pd.Timestamp(f.time_coverage_start).replace(tzinfo=None)

                # create an empty dictionary we will populate with values from file f
                this_row_dict = {}

                # create an empty list of the same length as data_vars to hold each variable's value
                values = ["" for v in data_vars]

                # For each variable we are interested, specified in the list "data_vars"
                for i, var in enumerate(data_vars):
                    # find corresponding pixel data_var value nearest to these scan angles y_rad and x_rad
                    values[i] = (
                        f[var].sel(y=y_rad, x=x_rad, method="nearest").values.mean()
                    )

                    # For all other products set ref_or_tb to None
                    ref_or_tb = None
                    # For ABI-L1b-Rad products only:
                    if var == "Rad":
                        # If we are looking at a reflective band (bands 1-6), convert Radiance to Reflectance
                        if f.band_id.values <= 6:
                            ref_or_tb = go.rad.goesReflectance(
                                values[i], f.kappa0.values
                            )
                        # If we are looking at an emissive band (bands 7-16), convert Radiance to Brightness Temperature (K)
                        else:
                            ref_or_tb = go.rad.goesBrightnessTemp(
                                values[i],
                                f.planck_fk1.values,
                                f.planck_fk2.values,
                                f.planck_bc1.values,
                                f.planck_bc2.values,
                            )

                # create a dictionary for this row of values (where each row is a GOES-R observation time)
                this_row_dict = dict(zip(data_vars, values))
                # add our time stamp to this dict before we update the dataframe
                this_row_dict["time"] = timestamp

                # If we have reflectance or brightness temperature to add to our dataframe
                if ref_or_tb is not None:
                    # add reflectance or brightness temperature to this row's update dict
                    this_row_dict["ref_or_tb"] = ref_or_tb

                # Finally, append this_row_dict to our dataframe for this one GOES-R observation time
                this_row_df = pd.DataFrame(this_row_dict, index=[0])
                df = pd.concat([df, this_row_df], ignore_index=True)

        except AttributeError as e:
            print(e)
            pass

    # drop duplicates if there are any, keep the first one
    df.drop_duplicates(["time"], keep="first", inplace=True)

    # set the dataframe intext to the timestamp column
    df.set_index("time", inplace=True, verify_integrity=True)

    # if an output filepath was provided, save the dataframe as a csv
    if outfilepath is not None:
        print("Saving csv file to: {}".format(outfilepath))
        df.to_csv(outfilepath)

    return df


def make_nested_abi_timeseries(
    directory, product, data_vars, lon, lat, z, outfilepath=None
):
    """Given a directory of GOES ABI products, create a timeseries of data variables (specified in data_vars) for a single point (at lon, lat, elevation).
    Retrieves all pixels nested within larger "2 km" ABI Fixed Grid cell.
    Returns a pandas dataframe, optional output to a csv file."""

    path = "{directory}/**/*{product}*.nc".format(directory=directory, product=product)
    file_list = glob.glob(path, recursive=True)

    path = "{directory}/**/*{product}*.nc".format(directory=directory, product=product)
    file_list = glob.glob(path, recursive=True)
    print(f"Found {len(file_list)} files in {path}")

    print(
        "Creating a timeseries of {data_vars} from {product} at ({lat}, {lon}, {z})\n".format(
            data_vars=data_vars, product=product, lat=lat, lon=lon, z=z
        )
    )
    # row_dicts = {}
    data_list = []
    # eSun_list = []
    print(f"Reading {len(file_list)} files from {path}\n")
    counter = 1
    for filename in file_list:
        try:
            print(
                "file {} of {}: {}".format(counter, len(file_list), filename), end="\r"
            )
            counter += 1
            with xr.open_dataset(filename, decode_times=False) as f:
                # I've included "decode_times=False" to this xr.open_dataset because I've encountered some ABI-L2-ACMC files where the timestamp couldn't be read
                # and xarray gave a "ValueError: unable to decode time units 'seconds since 2000-01-01 12:00:00' with the default calendar. Try opening your dataset with decode_times=False."
                # I've also switched which timestamp from the ABI files I'm reading (was f.time_bounds.values.min(), now f.time_coverage_start)
                # print(filename)
                # Read goes_imager_projection values needed for geometry calculations
                # and compute the corresponding look angles (in radiance) for the lat, lon, elevation we are interested in
                x_rad, y_rad = go.geometry.LonLat2ABIangle(
                    lon,
                    lat,
                    z,
                    f.goes_imager_projection.perspective_point_height
                    + f.goes_imager_projection.semi_major_axis,
                    f.goes_imager_projection.semi_major_axis,
                    f.goes_imager_projection.semi_minor_axis,
                    0.0818191910435,  # GRS-80 eccentricity
                    f.goes_imager_projection.longitude_of_projection_origin,
                )
                (
                    nearest_xs_2km,
                    nearest_ys_2km,
                    nearest_xs_1km,
                    nearest_ys_1km,
                    nearest_xs_500m,
                    nearest_ys_500m,
                ) = go.geometry.get_nested_coords(f, x_rad, y_rad)

                # get the timestamp for this observation (these should all be UTC, but I am removing timezone info because not all timestamps are converting the same way, and I was getting a "Cannot compare tz-naive and tz-aware timestamps" error)
                timestamp = (
                    pd.Timestamp(f.time_coverage_start)
                    .replace(tzinfo=None)
                    .round("min")
                )

                band = f.band_id.values[0]
                # band_formatted = "{:02.0f}".format(band)
                if band in [2]:
                    # print(f'Found band {f.band_id.values[0]} file...')
                    # print(f'Using pixel coordinates for 500m pixels: {nearest_xs_500m}, {nearest_ys_500m}')
                    # find corresponding pixel 'Rad' value nearest to these scan angles y_rad and x_rad
                    rad_values = (
                        f["Rad"]
                        .sel(
                            y=nearest_ys_500m[:, 0],
                            x=nearest_xs_500m[0, :],
                            method="nearest",
                        )
                        .rename("rad")
                    )  # .rename({'x': 'x05','y': 'y05'})
                    # If we are looking at a reflective band (bands 1-6), convert Radiance to Reflectance
                    ref_or_tb = go.rad.goesReflectance(
                        rad_values, f.kappa0.values
                    ).rename("ref")
                if band in [1, 3, 5]:
                    # print(f'Found band {f.band_id.values[0]} file...')
                    # print(f'Using pixel coordinates for 1km pixels: {nearest_xs_1km}, {nearest_ys_1km}')
                    # find corresponding pixel 'Rad' value nearest to these scan angles y_rad and x_rad
                    rad_values = (
                        f["Rad"]
                        .sel(
                            y=nearest_ys_1km[:, 0],
                            x=nearest_xs_1km[0, :],
                            method="nearest",
                        )
                        .rename("rad")
                    )  # .rename({'x': 'x1','y': 'y1'})
                    # If we are looking at a reflective band (bands 1-6), convert Radiance to Reflectance
                    ref_or_tb = go.rad.goesReflectance(
                        rad_values, f.kappa0.values
                    ).rename("ref")
                if band in [4, 6]:
                    # print(f'Found band {f.band_id.values[0]} file...')
                    # print(f'Using pixel coordinates for 1km pixels: {nearest_xs_2km}, {nearest_ys_2km}')
                    # find corresponding pixel 'Rad' value nearest to these scan angles y_rad and x_rad
                    rad_values = (
                        f["Rad"]
                        .sel(
                            y=nearest_ys_2km[:, 0],
                            x=nearest_xs_2km[0, :],
                            method="nearest",
                        )
                        .rename("rad")
                    )  #
                    # If we are looking at a reflective band (bands 1-6), convert Radiance to Reflectance
                    ref_or_tb = go.rad.goesReflectance(
                        rad_values, f.kappa0.values
                    ).rename("ref")
                if band in [7, 8, 9, 10, 11, 12, 13, 14, 15, 16]:
                    # print(f'Found band {f.band_id.values[0]} file...')
                    # print(f'Using pixel coordinates for 2km pixels: {nearest_xs_2km}, {nearest_ys_2km}')
                    # find corresponding pixel 'Rad' value nearest to these scan angles y_rad and x_rad
                    rad_values = (
                        f["Rad"]
                        .sel(
                            y=nearest_ys_2km[:, 0],
                            x=nearest_xs_2km[0, :],
                            method="nearest",
                        )
                        .rename("rad")
                    )  # .rename({'x': 'x2','y': 'y2'})
                    # If we are looking at an emissive band (bands 7-16), convert Radiance to Brightness Temperature (K)
                    ref_or_tb = go.rad.goesBrightnessTemp(
                        rad_values,
                        f.planck_fk1.values,
                        f.planck_fk2.values,
                        f.planck_bc1.values,
                        f.planck_bc2.values,
                    ).rename("tb")
                # append to list
                rad_values["t"] = timestamp.round("min")
                ref_or_tb["t"] = timestamp.round("min")
                data_list.append(
                    rad_values.expand_dims(dim={"t": 1})
                    .expand_dims(dim={"band": 1})
                    .assign_coords(band=("band", [band]))
                )
                data_list.append(
                    ref_or_tb.expand_dims(dim={"t": 1})
                    .expand_dims(dim={"band": 1})
                    .assign_coords(band=("band", [band]))
                )
        except (AttributeError, OSError) as e:
            print(e)
            pass

    df = data_list_to_df(data_list)

    # if an output filepath was provided, save the dataframe as a csv
    if outfilepath is not None:
        print("Saving csv file to: {}".format(outfilepath))
        df.to_csv(outfilepath)

    return df


def data_list_to_df(data_list):
    this_dict = {}
    counter = 1
    for i in range(len(data_list)):
        print("dataset {} of {}".format(counter, len(data_list)), end="\r")
        counter += 1
        if data_list[i].t.values[0] not in this_dict.keys():
            this_dict[
                data_list[i].t.values[0]
            ] = {}  # create new dict entry if it does not exist

        # now update that dict entry
        this_dict[data_list[i].t.values[0]]["t"] = data_list[i].t.values[0]
        this_dict[data_list[i].t.values[0]]["x_2km"] = data_list[i].x_image.values
        this_dict[data_list[i].t.values[0]]["y_2km"] = data_list[i].y_image.values
        if data_list[i].band.values == 2:  # 500m band
            this_dict[data_list[i].t.values[0]]["x_500m_WW"] = data_list[i].x.values[0]
            this_dict[data_list[i].t.values[0]]["x_500m_W"] = data_list[i].x.values[1]
            this_dict[data_list[i].t.values[0]]["x_500m_E"] = data_list[i].x.values[2]
            this_dict[data_list[i].t.values[0]]["x_500m_EE"] = data_list[i].x.values[3]
            this_dict[data_list[i].t.values[0]]["y_500m_SS"] = data_list[i].y.values[0]
            this_dict[data_list[i].t.values[0]]["y_500m_S"] = data_list[i].y.values[1]
            this_dict[data_list[i].t.values[0]]["y_500m_N"] = data_list[i].y.values[2]
            this_dict[data_list[i].t.values[0]]["y_500m_NN"] = data_list[i].y.values[3]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_500m_NW_NW"
            ] = data_list[i].values.ravel()[12]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_500m_NW_NE"
            ] = data_list[i].values.ravel()[13]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_500m_NW_SW"
            ] = data_list[i].values.ravel()[8]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_500m_NW_SE"
            ] = data_list[i].values.ravel()[9]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_500m_NE_NW"
            ] = data_list[i].values.ravel()[14]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_500m_NE_NE"
            ] = data_list[i].values.ravel()[15]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_500m_NE_SW"
            ] = data_list[i].values.ravel()[10]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_500m_NE_SE"
            ] = data_list[i].values.ravel()[11]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_500m_SW_NW"
            ] = data_list[i].values.ravel()[4]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_500m_SW_NE"
            ] = data_list[i].values.ravel()[5]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_500m_SW_SW"
            ] = data_list[i].values.ravel()[0]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_500m_SW_SE"
            ] = data_list[i].values.ravel()[1]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_500m_SE_NW"
            ] = data_list[i].values.ravel()[6]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_500m_SE_NE"
            ] = data_list[i].values.ravel()[7]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_500m_SE_SW"
            ] = data_list[i].values.ravel()[2]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_500m_SE_SE"
            ] = data_list[i].values.ravel()[3]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_2km"
            ] = data_list[i].values.mean()
        elif data_list[i].band.values in [1, 3, 5]:  # 1km bands
            this_dict[data_list[i].t.values[0]]["x_1km_W"] = data_list[i].x.values[0]
            this_dict[data_list[i].t.values[0]]["x_1km_E"] = data_list[i].x.values[1]
            this_dict[data_list[i].t.values[0]]["y_1km_N"] = data_list[i].y.values[1]
            this_dict[data_list[i].t.values[0]]["y_1km_S"] = data_list[i].y.values[0]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_1km_NW"
            ] = data_list[i].values.ravel()[0]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_1km_NE"
            ] = data_list[i].values.ravel()[1]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_1km_SW"
            ] = data_list[i].values.ravel()[2]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_1km_SE"
            ] = data_list[i].values.ravel()[3]
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_2km"
            ] = data_list[i].values.mean()
        else:  # 2km bands
            this_dict[data_list[i].t.values[0]][
                f"b{data_list[i].band.values[0]}_{data_list[i].name}_2km"
            ] = data_list[i].values.ravel()[0]

    # drop duplicates if there are any, keep the first one
    # df.drop_duplicates(['time'], keep='first', inplace=True)

    df = pd.DataFrame.from_dict(this_dict, orient="index")

    # set the dataframe intext to the timestamp column
    # df.set_index('time', inplace = True, verify_integrity = True)

    return df
