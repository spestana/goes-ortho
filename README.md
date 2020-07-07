# goes-view

---

### activate the environment:

```conda activate goes-linux```

---

### download-goes.py

Downloads GOES-16 or GOES-17 products/bands, requires command line arguments:

#### Usage:

```python ./download-goes.py --bucket <S3-BUCKET> --year <YEAR> --month <MONTH> --days <START DAY> <END DAY> --product <ABI PRODUCT CODE> --channel <ABI CHANNEL> --bounds <MIN_LAT> <MAX_LAT> <MIN_LON> <MAX_LON> --dir <DESTINATION DIRECTORY>```

#### Examples:

This will download the GOES-16 ABI Level-1b Radiance (CONUS) product for channel/band 14, for January 1-2 2020. The NetCDF files will be cropped to within latitudes 30 - 50 and longitudes -125 - -105, and saved in /storage/spestana/scratchspace.

```python ./download-goes.py --bucket noaa-goes16 --year 2020 --month 2 --days 1 2 --product ABI-L1b-RadC --channel C14 --bounds 30 50 -125 -105 --dir /storage/spestana/scratchspace```

We can do the same command with short flag names:

```python ./download-goes.py -B noaa-goes16 -Y 2020 -M 1 -D 1 2 -p ABI-L1b-RadC -c C14 -b 30 50 -125 -105 -d /storage/spestana/scratchspace```

---

### goes-timeseries.py

Creates a time series of GOES ABI radiance values for a specified point location. This takes into account the point's elevation (in meters) to correct for terrain parallax from off-nadir view angles of GOES.

#### Usage:

```python ./goes-timeseries.py -d /storage/GOES/goes16/2017/03 -l <LATITUDE> <LONGITUDE> <ELEVATION>```

#### Examples:

Gaylor Pit @ lat=37.88175, lon=-119.31212, elev=2811:

```python ./goes-timeseries.py -d /storage/GOES/goes16/2017/03 -l 37.88175 -119.31212 2811```


Grand Mesa West @ lat=39.0339, lon=-108.2140, elev=3033:

```python ./goes-timeseries.py -d /storage/GOES/goes16/2017/03 -l 39.0339 -108.2140 3033```


CUES site @  lat=37.643103, lon=-119.029146, elev=2940:

```python ./goes-timeseries.py -d /storage/GOES/goes16/2017/03 -l 37.643103 -119.029146 2940```

---

### goes_ortho.py

Functions for orthorectifying GOES-R ABI imagery using a DEM.

#### Usage:

```python
# import to use these functions
import goes_ortho

# specify filepaths for inputs
abi_filepath = '.\OR_ABI-L1b-RadC-M4C14_G16_s20171111750224_e20171111755027_c20171111755074.nc'
dem_filepath = '.\dem.tif'

# generate the pixel mapping
pixel_map = goes_ortho.make_ortho_map(abi_filepath, dem_filepath)

# orthorectify the image
goes_ortho.orthorectify_abi_rad(abi_filepath, pixel_map, out_filename='test_ortho.nc')
```

#### Example:

see the [goes-orthorectify](https://github.com/spestana/goes-view/blob/master/goes-orthorectify.ipynb) notebook for an example

---
---