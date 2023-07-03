"""
Functions for converting GOES ABI Radiance values
"""

"""Steven Pestana, March 2021 (spestana@uw.edu)

For background information, see: 
 - http://cimss.ssec.wisc.edu/goes/calibration/#eqw
 - http://cimss.ssec.wisc.edu/goes/calibration/Converting_AHI_RadianceUnits_24Feb2015.pdf
"""

import pandas as pd
import numpy as np


def goesBrightnessTemp(rad, fk1, fk2, bc1, bc2): 
    """
    Convert Radiance to Brightness Temperature for GOES-R ABI emissive bands (7-16)
    Parameters
    ------------
    rad: float, np.array, or xarray.DataArray
        radiance [mW / m^2 sr cm^-1]
    fk1: float
        Planck function coefficient 1, from GOES ABI product metadata
    fk2: float
        Planck function coefficient 2, from GOES ABI product metadata
    bc1: float
        spectral response function offset correction term, from GOES ABI product metadata
    bc2: float
        spectral response function scale correction term, from GOES ABI product metadata
    Returns
    ------------
    Tb: float, np.array, or xarray.DataArray
        brightness temperature [K]
    """
    Tb = ( fk2 / (np.log((fk1 / rad) + 1)) - bc1 ) / bc2
    return Tb


def goesReflectance(rad, kappa): 
    """
    Convert Radiance to Reflectance for GOES-R ABI reflective bands (1-6)
    
    Parameters
    ------------
    rad: float, np.array, or xarray.DataArray
        radiance [mW / m^2 sr cm^-1]
    kappa: float
        incident Lambertian equivalent radiance, from GOES ABI product metadata
    Returns
    ------------
    ref: float, np.array, or xarray.DataArray
        reflectance factor
    """
    ref = kappa * rad
    return ref

def abi_radiance_wavenumber_to_wavelength(goes, channel, rad_wn):
    """
    Convert GOES ABI Radiance units from [mW / m^2 sr cm^-1] to [W / m^2 sr um]
    
    Parameters
    ------------
    goes: int
        16 or 17 to select GOES-16 or GOES-17
    channel: int
        1-16 to select GOES ABI channel/band
    rad_wn: float, np.array, or xarray.DataArray
        GOES ABI Radiance in "wavenumber" units [mW / m^2 sr cm^-1]
    Returns
    ------------
    rad_wl:  float, np.array, or xarray.DataArray
        GOES ABI Radiance in "wavelength" units [W / m^2 sr um]
    """
    
    # Read in Band Equivalent Widths file for GOES16 or GOES17
    eqw = pd.read_csv('../data/GOES{goes}_ABI_ALLBANDS_MAR2016.eqw'.format(goes=str(goes)), sep=r'\s+', skiprows=1, index_col='CHANNEL')    
    
    # Divide milliwats by 1000 to get watts
    scale_milliwatts_by = 1000 
    
    # Convert units
    rad_wl = (rad_wn / scale_milliwatts_by) * (eqw['EQW(cm-1)'][channel]/eqw['EQW(um)'][channel])
    
    return rad_wl