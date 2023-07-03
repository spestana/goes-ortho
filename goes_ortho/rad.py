"""
Functions for converting GOES ABI Radiance values
Steven Pestana, March 2021 (spestana@uw.edu)

For background information, see: 
 - http://cimss.ssec.wisc.edu/goes/calibration/#eqw
 - http://cimss.ssec.wisc.edu/goes/calibration/Converting_AHI_RadianceUnits_24Feb2015.pdf
"""

import pandas as pd
import numpy as np


def goesBrightnessTemp(rad, fk1, fk2, bc1, bc2): 
    ''' Convert Radiance to Brightness Temperature for GOES-R ABI emissive bands (7-16)
    
    :param rad: radiance [mW / m^2 sr cm^-1]
    :type rad: float, np.array, xarray.DataArray
    :param fk1: Planck function coefficient 1, from GOES ABI product metadata
    :type fk1: float
    :param fk2: Planck function coefficient 2, from GOES ABI product metadata
    :type fk2: float
    :param bc1: spectral response function offset correction term, from GOES ABI product metadata
    :type bc1: float
    :param bc2: spectral response function scale correction term, from GOES ABI product metadata
    :type bc2: float
    
    :return: Tb, brightness temperature [K]
    :rtype: float, np.array, xarray.DataArray
    '''
    Tb = ( fk2 / (np.log((fk1 / rad) + 1)) - bc1 ) / bc2
    return Tb


def goesReflectance(rad, kappa): 
    ''' Convert Radiance to Reflectance for GOES-R ABI reflective bands (1-6)
    
    :param rad: radiance [mW / m^2 sr cm^-1]
    :type rad: float, np.array, xarray.DataArray
    :param kappa: incident Lambertian equivalent radiance, from GOES ABI product metadata
    :type kappa: float
    
    :return: ref, reflectance factor
    :rtype: float, np.array, xarray.DataArray
    '''
    ref = kappa * rad
    return ref

def abi_radiance_wavenumber_to_wavelength(goes, channel, rad_wn):
    ''' Convert GOES ABI Radiance units from [mW / m^2 sr cm^-1] to [W / m^2 sr um]
    :param goes: 16 or 17 to select GOES-16 or GOES-17
    :type goes: int
    :param channel: 1-16 to select GOES ABI channel/band
    :type channel: int
    :param rad_wn: GOES ABI Radiance in "wavenumber" units [mW / m^2 sr cm^-1]
    :type rad_wn:  float, np.array, xarray.DataArray
    
    :return: rad_wl, GOES ABI Radiance in "wavelength" units [W / m^2 sr um]
    :rtype: float, np.array, xarray.DataArray
    '''
    
    # Read in Band Equivalent Widths file for GOES16 or GOES17
    eqw = pd.read_csv('../data/GOES{goes}_ABI_ALLBANDS_MAR2016.eqw'.format(goes=str(goes)), sep=r'\s+', skiprows=1, index_col='CHANNEL')    
    
    # Divide milliwats by 1000 to get watts
    scale_milliwatts_by = 1000 
    
    # Convert units
    rad_wl = (rad_wn / scale_milliwatts_by) * (eqw['EQW(cm-1)'][channel]/eqw['EQW(um)'][channel])
    
    return rad_wl