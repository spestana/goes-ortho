import numpy as np
import pandas as pd

"""
Functions for converting GOES ABI Radiance values
"""

"""Steven Pestana, March 2021 (spestana@uw.edu)

For background information, see:
 - http://cimss.ssec.wisc.edu/goes/calibration/#eqw
 - http://cimss.ssec.wisc.edu/goes/calibration/Converting_AHI_RadianceUnits_24Feb2015.pdf
"""


def goesBrightnessTemp(rad, fk1, fk2, bc1, bc2):
    """
    Convert Radiance to Brightness Temperature for GOES-R ABI emissive bands (7-16)

    Parameters
    ------------
    rad : float, np.array, or xarray.DataArray
        radiance [mW / m^2 sr cm^-1]
    fk1 : float
        Planck function coefficient 1, from GOES ABI product metadata
    fk2 : float
        Planck function coefficient 2, from GOES ABI product metadata
    bc1 : float
        spectral response function offset correction term, from GOES ABI product metadata
    bc2 : float
        spectral response function scale correction term, from GOES ABI product metadata

    Returns
    ------------
    Tb : float, np.array, or xarray.DataArray
        brightness temperature [K]

    Examples
    ------------

    """
    Tb = (fk2 / (np.log((fk1 / rad) + 1)) - bc1) / bc2
    return Tb


def goesReflectance(rad, kappa):
    """
    Convert Radiance to Reflectance for GOES-R ABI reflective bands (1-6)

    Parameters
    ------------
    rad : float, np.array, or xarray.DataArray
        radiance [mW / m^2 sr cm^-1]
    kappa : float
        incident Lambertian equivalent radiance, from GOES ABI product metadata

    Returns
    ------------
    ref : float, np.array, or xarray.DataArray
        reflectance factor

    Examples
    ------------

    """
    ref = kappa * rad
    return ref


def abi_radiance_wavenumber_to_wavelength(goes, channel, rad_wn):
    """
    Convert GOES ABI Radiance units from [mW / m^2 sr cm^-1] to [W / m^2 sr um]

    Parameters
    ------------
    goes : int
        16, 17, or 18 to select GOES-16, GOES-17, or GOES-18
    channel : int
        1-16 to select GOES ABI channel/band
    rad_wn : float, np.array, or xarray.DataArray
        GOES ABI Radiance in "wavenumber" units [mW / m^2 sr cm^-1]

    Returns
    ------------
    rad_wl :  float, np.array, or xarray.DataArray
        GOES ABI Radiance in "wavelength" units [W / m^2 sr um]

    Examples
    ------------

    """

    # Read in Band Equivalent Widths file for GOES16 or GOES17
    eqw = pd.read_csv(
        "../data/GOES{goes}_ABI_ALLBANDS_MAR2016.eqw".format(goes=str(goes)),
        sep=r"\s+",
        skiprows=1,
        index_col="CHANNEL",
    )

    # Divide milliwats by 1000 to get watts
    scale_milliwatts_by = 1000

    # Convert units
    rad_wl = (rad_wn / scale_milliwatts_by) * (
        eqw["EQW(cm-1)"][channel] / eqw["EQW(um)"][channel]
    )

    return rad_wl


def makeABIrgb_fromReflectance(
    R_ref,
    G_ref,
    B_ref,
    gamma=2.2,
    green_coefficients=None,
):
    """Create RGB images given GOES-R ABI Channel 01, 02, and 03 datasets. Adapted from https://github.com/daniellelosos/True-Color-Image_GOES-R-ABI-L1b-Radiances

    Parameters
    ------------
    R_ref : np.ndarray
        Red band data from GOES ABI (Channel 2)
    G_ref : np.ndarray
        "Green" band data from GOES ABI (Channel 3) (ABI does not have a true green band, instead we can use the NIR "Veggie" band)
    B_ref : np.ndarray
        Blue band data from GOES ABI (Channel 1)
    gamma : float
        Gamma correction to adjust brightness of output image, defaults to gamma=2.2
    green_coefficients : dict
        Dictionary of multipliers for the red, nir, and blue bands to create a synthetic green band, defaults to green_coefficients=None, which will then set these to {'red': 0.45, 'nir': 0.1, 'blue': 0.45}

    Returns
    ------------
    RGB :  np.ndarray
        "True Color" RGB
    RGB_veggie : np.ndarray
        "False Color" RGB

    Examples
    ------------

    """

    if green_coefficients is None:
        green_coefficients = {"red": 0.45, "nir": 0.1, "blue": 0.45}

    # Apply range limits for each channel. Reflectance values must be between 0 and 1.
    R_ref = np.clip(R_ref, 0, 1)
    G_ref = np.clip(G_ref, 0, 1)
    B_ref = np.clip(B_ref, 0, 1)

    # Apply a gamma correction to the image to correct ABI detector brightness
    Red = np.power(R_ref, 1 / gamma)
    Green = np.power(G_ref, 1 / gamma)
    Blue = np.power(B_ref, 1 / gamma)

    # GOES-R Series satellites do not have a channel in the visible green range. Band 3 is a NIR channel typically used to monitor vegetation.
    # Calculate the "True" Green Band to serve as a green proxy for the RGB True Color image, using a fractional combination.
    # Source: "Generation of GOES‐16 True Color Imagery without a Green Band" - https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018EA000379
    Green_true = (
        green_coefficients["red"] * Red
        + green_coefficients["nir"] * Green
        + green_coefficients["blue"] * Blue
    )
    Green_true = np.clip(Green_true, 0, 1)  # Apply band limits again, just in case.

    # Combine three RGB channels with a stacked array, then display the resulting images.

    # The RGB array with the raw veggie band
    RGB_veggie = np.dstack([Red, Green, Blue])

    # The RGB array for the true color image
    RGB = np.dstack([Red, Green_true, Blue])

    return RGB, RGB_veggie
