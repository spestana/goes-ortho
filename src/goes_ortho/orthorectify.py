"""
Functions to orthorectify GOES-R ABI images using a DEM
"""

import os

import numpy as np
import rioxarray
import xarray as xr

from goes_ortho.geometry import LonLat2ABIangle
from goes_ortho.get_data import get_dem
from goes_ortho.rad import goesBrightnessTemp, goesReflectance


def ABIpixelMap(abi_grid_x, abi_grid_y):
    """
    Converts an array of continuous ABI scan angles into discrete pixel center locations (in scan angle coordinates, incrimenting by the pixel IFOV) NOTE: This function isn't needed for the applying the mapping to a GOES ABI image, but we can still use this to make some visualizations of what we're doing.

    Parameters
    ------------
    abi_grid_x : np.array
        2-dimensional array of x coordinates (scan angle) in ABI Fixed Grid [radians]
    abi_grid_y : np.array
        2-dimensional array of y coordinates (elevation angle) in ABI Fixed Grid [radians]

    Returns
    ------------
    center_x : np.array
        pixel center x coordinates (scan angle) in ABI Fixed Grid [radians]
    center_y : np.array
        pixel center y coordinates (elevation angle) in ABI Fixed Grid [radians]

    Examples
    ------------

    """

    # IFOV values for GOES ABI bands ("500 m" 14 urad; "1 km" 28 urad; "2 km" 56 urad)
    ifov = np.array([14e-6, 28e-6, 56e-6])

    # Convert from scan angle to pixel row/column coordinates
    x_px = np.array([np.divide(abi_grid_x, i) for i in ifov])
    y_px = np.array([np.divide(abi_grid_y, i) for i in ifov])

    # Get the center coordinate of the pixel each grid cell lies within
    center_x = ((np.floor(np.abs(x_px)) + 0.5) * np.sign(x_px)) * ifov[:, None, None]
    center_y = ((np.floor(np.abs(y_px)) + 0.5) * np.sign(y_px)) * ifov[:, None, None]

    # Get the pixel coordinate (row/column) that each grid cell lies within
    # center_col = (np.floor(x_px))
    # center_row = (np.floor(y_px))

    return center_x, center_y  # , center_col, center_row


def make_ortho_map(goes_filepath, dem_filepath, out_filepath=None):
    """
    For the entire DEM, determine the ABI scan angle coordinates for every DEM grid cell, taking into account the underlying terrain and satellite's viewing geometry. Create the mapping between GOES-R ABI pixels (netCDF input file) and a DEM grid (geotiff input file)

    Parameters
    ------------
    goes_filepath : str
        filepath to GOES ABI NetCDF file
    dem_filepath : str
        filepath to digital elevation model (DEM), GeoTiff file
    out_filepath : str
        optional filepath and filename to save this map to, defaults to None

    Returns
    ------------
    ds : xarray.Dataset
        dataset of the map relating ABI Fixed Grid coordinates to latitude and longitude

    Examples
    ------------

    """

    print("\nRUNNING: make_ortho_map()")

    # Open the GOES ABI image
    print("\nOpening GOES ABI image...")
    abi_image = xr.open_dataset(goes_filepath, decode_times=False)
    # NOTE: for some reason (?) I sometimes get an error "ValueError: unable to decode time units 'seconds since 2000-01-01 12:00:00' with the default calendar. Try opening your dataset with decode_times=False." so I've added decode_times=False here.
    # Get inputs: projection information from the ABI radiance product (values needed for geometry calculations)
    print("\nGet inputs: projection information from the ABI radiance product")
    req = abi_image.goes_imager_projection.semi_major_axis
    rpol = abi_image.goes_imager_projection.semi_minor_axis
    H = (
        abi_image.goes_imager_projection.perspective_point_height
        + abi_image.goes_imager_projection.semi_major_axis
    )
    lon_0 = abi_image.goes_imager_projection.longitude_of_projection_origin
    e = 0.0818191910435  # GRS-80 eccentricity
    print("...done")

    # Load DEM
    print("\nOpening DEM file...")
    dem = rioxarray.open_rasterio(dem_filepath)
    dem = dem.where(dem != dem.attrs["_FillValue"])[0, :, :]  # replace nodata with nans
    dem = dem.fillna(
        0
    )  # fill nans with zeros for the ocean (temporary fix for fog project)
    # dem = dem.where(dem!=0) # replace zeros with nans
    # Create 2D arrays of longitude and latitude from the DEM
    print("\nCreate 2D arrays of longitude and latitude from the DEM")
    X, Y = np.meshgrid(dem.x, dem.y)  # Lon and Lat of each DEM grid cell
    Z = dem.values  # elevation of each DEM grid cell
    print("...done")

    # For each grid cell in the DEM, compute the corresponding ABI scan angle (x and y, radians)
    print(
        "\nFor each grid cell in the DEM, compute the corresponding ABI scan angle (x and y, radians)"
    )
    abi_grid_x, abi_grid_y = LonLat2ABIangle(X, Y, Z, H, req, rpol, e, lon_0)
    print("...done")

    # Create metadata dictionary about this map (should probably clean up metadata, adhere to some set of standards)
    print("\nCreate metadata dictionary about this map")
    metadata = {
        # Information about the projection geometry:
        "longitude_of_projection_origin": lon_0,
        "semi_major_axis": req,
        "semi_minor_axis": rpol,
        "satellite_height": H,
        "grs80_eccentricity": e,
        "longitude_of_projection_origin_info": "longitude of geostationary satellite orbit",
        "semi_major_axis_info": "semi-major axis of GRS 80 reference ellipsoid",
        "semi_minor_axis_info": "semi-minor axis of GRS 80 reference ellipsoid",
        "satellite_height_info": "distance from center of ellipsoid to satellite (perspective_point_height + semi_major_axis_info)",
        "grs80_eccentricity_info": "eccentricity of GRS 80 reference ellipsoid",
        # Information about the DEM source file
        "dem_file": dem_filepath,
        #'dem_crs' : dem.crs,
        #'dem_transform' : dem.transform,
        #'dem_res' : dem.res,
        #'dem_ifov': -9999, # TO DO
        "dem_file_info": "filename of dem file used to create this mapping",
        "dem_crs_info": "coordinate reference system from DEM geotiff",
        "dem_transform_info": "transform matrix from DEM geotiff",
        "dem_res_info": "resolution of DEM geotiff",
        "dem_ifov_info": "instantaneous field of view (angular size of DEM grid cell)",
        # For each DEM grid cell, we have...
        "dem_px_angle_x_info": "DEM grid cell X coordinate (east/west) scan angle in the ABI Fixed Grid",
        "dem_px_angle_y_info": "DEM grid cell Y coordinate (north/south) scan angle in the ABI Fixed Grid",
        "longitude_info": "longitude from DEM file",
        "latitude_info": "latitude from DEM file",
        "elevation_info": "elevation from DEM file",
    }
    print("...done")

    # Create pixel map dataset
    print("\nCreate pixel map dataset")
    ds = xr.Dataset(
        {"elevation": (["latitude", "longitude"], dem.values)},
        coords={
            "longitude": (["longitude"], dem.x.data),
            "latitude": (["latitude"], dem.y.data),
            "dem_px_angle_x": (["latitude", "longitude"], abi_grid_x),
            "dem_px_angle_y": (["latitude", "longitude"], abi_grid_y),
        },
        attrs=metadata,
    )
    print(ds)
    print("...done")

    if out_filepath is not None:
        print("\nExport this pixel map along with the metadata (NetCDF with xarray)")
        # Export this pixel map along with the metadata (NetCDF with xarray)
        ds.to_netcdf(out_filepath, mode="w")
        print("...done")

    # Return the pixel map dataset
    print("\nReturn the pixel map dataset.")

    return ds


def orthorectify_abi(goes_filepath, pixel_map, data_vars, out_filename=None):
    """
    Using the pixel mapping for a specific ABI viewing geometry over a particular location,
    orthorectify the ABI radiance values and return an xarray dataarray with those values.

    Parameters
    ------------
    goes_filepath : str
        filepath to GOES ABI NetCDF file
    pixel_map : xarray.Dataset
        dataset of the map relating ABI Fixed Grid coordinates to latitude and longitude
    data_vars : list
        list of variable names from the GOES ABI NetCDF file we wish to extract
    out_filename : str
        optional filepath and filename to save the orthorectified image to, defaults to None

    Returns
    ------------
    pixel_map : xarray.Dataset
        dataset of the orthorectified GOES ABI image

    Examples
    ------------

    """
    print("\nRUNNING: orthorectify_abi_rad()")

    # First check, Does the projection info in the image match our mapping?
    print("\nDoes the projection info in the image match our mapping?")
    # Open the GOES ABI image
    print("\nOpening GOES ABI image...\t\t\tABI image value\tPixel map value")
    abi_image = xr.open_dataset(goes_filepath, decode_times=False)
    print(
        "perspective_point_height + semi_major_axis:\t{}\t{}".format(
            abi_image.goes_imager_projection.perspective_point_height
            + abi_image.goes_imager_projection.semi_major_axis,
            pixel_map.satellite_height,
        )
    )
    print(
        "semi_major_axis:\t\t\t\t{}\t{}".format(
            abi_image.goes_imager_projection.semi_major_axis, pixel_map.semi_major_axis
        )
    )
    print(
        "semi_minor_axis:\t\t\t\t{}\t{}".format(
            abi_image.goes_imager_projection.semi_minor_axis, pixel_map.semi_minor_axis
        )
    )
    print(
        "longitude_of_projection_origin:\t\t\t{}\t\t{}".format(
            abi_image.goes_imager_projection.longitude_of_projection_origin,
            pixel_map.longitude_of_projection_origin,
        )
    )
    print("...done")

    # Map (orthorectify) and clip the image to the pixel map for each data variable we want
    for var in data_vars:
        print(
            "\nMap (orthorectify) and clip the image to the pixel map for {}".format(
                var
            )
        )
        abi_var_values = abi_image.sel(
            x=pixel_map.dem_px_angle_x, y=pixel_map.dem_px_angle_y, method="nearest"
        )[var].values
        print("...done")

        # Create a new xarray dataset with the orthorectified ABI radiance values,
        # Lat, Lon, Elevation, and metadata from the pixel map.
        pixel_map[var] = (["latitude", "longitude"], abi_var_values)
        # If we are looking at an ABI-L1b-Rad product, create either a reflectance (bands 1-6) or brightness temperautre (bands 7-16) dataset
        if var == "Rad":
            # if we are looking at bands 1-6, compute reflectance
            if abi_image.band_id.values[0] <= 6:
                pixel_map["ref"] = goesReflectance(
                    pixel_map[var], abi_image.kappa0.values
                )
            # else, compute brightness temperature for bands 7-16
            else:
                pixel_map["tb"] = goesBrightnessTemp(
                    pixel_map[var],
                    abi_image.planck_fk1.values,
                    abi_image.planck_fk2.values,
                    abi_image.planck_bc1.values,
                    abi_image.planck_bc2.values,
                )

    # Map (orthorectify) the original ABI Fixed Grid coordinate values to the new pixels for reference
    print(
        "\nMap (orthorectify) and clip the image to the pixel map for ABI Fixed Grid coordinates"
    )
    abi_fixed_grid_x_values = abi_image.sel(
        x=pixel_map.dem_px_angle_x.values.ravel(), method="nearest"
    ).x.values
    abi_fixed_grid_y_values = abi_image.sel(
        y=pixel_map.dem_px_angle_y.values.ravel(), method="nearest"
    ).y.values
    abi_fixed_grid_x_values_reshaped = np.reshape(
        abi_fixed_grid_x_values, pixel_map.dem_px_angle_x.shape
    )
    abi_fixed_grid_y_values_reshaped = np.reshape(
        abi_fixed_grid_y_values, pixel_map.dem_px_angle_y.shape
    )
    pixel_map["abi_fixed_grid_x"] = (
        ("latitude", "longitude"),
        abi_fixed_grid_x_values_reshaped,
    )
    pixel_map["abi_fixed_grid_y"] = (
        ("latitude", "longitude"),
        abi_fixed_grid_y_values_reshaped,
    )
    print("...done")

    # drop DEM from dataset
    # pixel_map = pixel_map.drop(['elevation'])

    print(
        "\nCreate zone labels for each unique pair of ABI Fixed Grid coordinates (for each orthorectified pixel footprint)"
    )
    # Found this clever solution here: https://stackoverflow.com/a/32326297/11699349
    # Create unique values for every "zone" (the GOES ABI pixel footprints) with the same ABI Fixed Grid X and Y values
    unique_values = (
        pixel_map.abi_fixed_grid_x.values
        * (pixel_map.abi_fixed_grid_y.values.max() + 1)
        + pixel_map.abi_fixed_grid_y.values
    )
    # Find the index of all unique values we just created
    _, idx = np.unique(unique_values, return_inverse=True)
    # Use these indices, reshaped to the original shape, as our zone labels
    zone_labels = idx.reshape(pixel_map.abi_fixed_grid_y.values.shape)
    # Add the zone_labels to the dataset
    pixel_map["zone_labels"] = (("latitude", "longitude"), zone_labels)
    print("...done")

    # Output this result to a new NetCDF file
    print("\nOutput this result to a new NetCDF file")
    if out_filename is None:
        out_filename = abi_image.dataset_name + "_ortho.nc"
    print("Saving file as: {}".format(out_filename))

    pixel_map.to_netcdf(out_filename)
    print("...done")

    return pixel_map


def ortho(
    goes_image_path,
    data_vars,
    bounds,
    api_key,
    new_goes_filename,
    dem_filepath=None,
    demtype="SRTMGL3",
    keep_dem=True,
):
    """
    Wraps around get_dem(), make_ortho_map(), orthorectify_abi()

    Parameters
    ------------
    goes_image_path : str
        filepath to GOES ABI NetCDF file
    data_vars : list
        list of variable names from the GOES ABI NetCDF file we wish to extract
    bounds : list
        longitude and latitude bounds to clip and orthorectify GOES ABI image, like [min_lon, min_lat, max_lon, max_lat]
    api_key : str
        Opentopography.org API key, can be created at https://portal.opentopography.org/requestService?service=api
    new_goes_filename : str
        new filepath and filename to save the orthorectified image to
    dem_filepath : str
        filepath to save DEM to, defaults to None
    demtype : str
        DEM from Opentopography.org, see documentation in get_data.get_dem()
    keep_dem : bool
        option to save DEM file or delete after use

    Returns
    ------------
    None

    Examples
    ------------

    """

    if dem_filepath is None:
        dem_filepath = "temp_{demtype}_DEM.tif".format(
            demtype=demtype,
        )
    get_dem(
        demtype=demtype,
        bounds=bounds,
        api_key=api_key,
        out_fn=dem_filepath,
        proj="+proj=lonlat +datum=GRS80",
    )  # make sure to convert to GRS80 ellipsoid model GOES ABI fixed grid uses

    # create the mapping between scan angle coordinates and lat/lon given the GOES satellite position and our DEM
    goes_ortho_map = make_ortho_map(goes_image_path, dem_filepath)

    # Apply the "ortho map" and save a new NetCDF file with data variables from the original file
    _ = orthorectify_abi(
        goes_image_path, goes_ortho_map, data_vars, out_filename=new_goes_filename
    )

    # If keep_dem is False, delete the temporary DEM file we downloaded
    if not keep_dem:
        os.remove(dem_filepath)

    return None
