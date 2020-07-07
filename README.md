# goes-view

---

### activate the environment:

```conda activate goes-linux```

---

### download-goes.py

Downloads GOES-16 or GOES-17 products/bands (specified within the script itself, to do: use command line arguments or a text file to specify what we want to download)

#### Usage:

```python ./download-goes.py```

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