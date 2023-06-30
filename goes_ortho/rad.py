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
    ''' Convert Radiance to Brightness Temperature for GOES-R ABI emissive bands (7-16)'''
    Tb = ( fk2 / (np.log((fk1 / rad) + 1)) - bc1 ) / bc2
    return Tb


def goesReflectance(rad, kappa): 
    ''' Convert Radiance to Reflectance for GOES-R ABI reflective bands (1-6)'''
    ref = kappa * rad
    return ref

def abi_radiance_wavenumber_to_wavelength(goes, channel, rad_wn):
    ''' Convert GOES ABI Radiance units from
                                                mW / m^2 sr cm^-1
                                        to
                                                W / m^2 sr um
    Inputs
     - goes = 16 or 17 (int); GOES-16 or GOES-17
     - channel = 1-16 (int); GOES ABI channel/band
     - rad_wn = GOES ABI Radiance in "wavenumber" units [mW / m^2 sr cm^-1]
    Outputs
     - rad_wl = GOES ABI Radiance in "wavelength" units [W / m^2 sr um]
    '''
    
    # Read in Band Equivalent Widths file for GOES16 or GOES17
    eqw = pd.read_csv('GOES{goes}_ABI_ALLBANDS_MAR2016.eqw'.format(goes=str(goes)), sep=r'\s+', skiprows=1, index_col='CHANNEL')    
    
    # Divide milliwats by 1000 to get watts
    scale_milliwatts_by = 1000 
    
    # Convert units
    rad_wl = (rad_wn / scale_milliwatts_by) * (eqw['EQW(cm-1)'][channel]/eqw['EQW(um)'][channel])
    
    return rad_wl