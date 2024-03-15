'''Functions for testing rad module.'''

import glob
import xarray as xr
import goes_ortho as go

def test_goesBrightnessTemp(setup_session):
    # get Planck function coefficients (fk1, fk2) and spectral response function correction terms (bc1, bc2) from an example GOES-R L1b Radiance product (bands 7-16)
    filepaths = glob.glob(f'./tests/resources/spestana-goes-ortho-data-*/data/F/OR_ABI-L1b-RadF-M6C16_G18_s20230370050213_e20230370059532_c20230370059568.nc')
    with xr.open_dataset(filepaths[0]) as ds:
        fk1 = ds.planck_fk1.values
        fk2 = ds.planck_fk2.values
        bc1 = ds.planck_bc1.values
        bc2 = ds.planck_bc2.values
    rad = 83.64517 # W m-2 sr-1 um-1
    # run goesBrightnessTemp()
    Tb = go.rad.goesBrightnessTemp(rad, fk1, fk2, bc1, bc2)
    assert Tb == 262.8010077455753 # K


def test_goesReflectance(setup_session):
    # incident Lambertian equivalent radiance term (kappa) from an example GOES-R L1b Radiance product (bands 1-6)
    filepaths = glob.glob(f'./tests/resources/spestana-goes-ortho-data-*/data/C/OR_ABI-L1b-RadC-M6C03_G18_s20231821906187_e20231821908560_c20231821909002.nc')
    print(filepaths)
    with xr.open_dataset(filepaths[0]) as ds:
        kappa = ds.kappa0.values
    rad = 72.73628 # W m-2 sr-1 um-1
    # run goesReflectance()
    Ref = go.rad.goesReflectance(rad, kappa)
    assert Ref == 0.24755792572951874


def test_abi_radiance_wavenumber_to_wavelength(setup_session):
    for goes in [16,17,18]:
        for channel in range(16):
            channel += 1 # add 1 so we get channel numbers 1 through 16 (rather than 0 through 15)
            rad_wl = radiance_wavenumber_to_wavelength(goes, channel, rad_wn)
            print(rad_wl)
    assert 0


