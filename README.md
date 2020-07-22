# "The Mountains Are Calling, and I Must GOES" 
 
- [John Muir](https://en.wikipedia.org/wiki/John_Muir) if he had access to satellite remote sensing probably

---
# goes-ortho

The latest generation of geostationary-orbiting weather satellites make frequent observations (5-60 min) in the visible and IR wavelengths at moderate resolutions (500-~2000m) (GOES-16/17 ABI, Himawari 8/9 AHI). This lends themselves to be used to fill temporal gaps between ~12-hour repeat observations by polar-orbiting spacecraft (Aqua/Terra MODIS, SNPP VIIRS, etc). 
However, their geostationary orbits mean that outside of their sub-satellite-point on the equator, all other view angles are off-nadir, and due to the Earth's curvature in view, actual pixel sizes increase to >6 km towards the planet's limb.
Additionally, when viewing complex terrain such as the mountains of western CONUS, parallax affects the apparent position of the variable topography. Some portions of the ground suface may even become obscured from view completely by surrounding steep terrain with certain aspects.

Before using observations from these instruments for observing the land surface over mountains, orthorectification is needed to try and account for the off-nadir view angles and topographic effects in the geostationary satellite imagery. 

The terrain parallax is especially visually apparent when flipping between GOES-East and GOES-West views of a mountain range like the Sierra Nevada here:

<img src="https://github.com/spestana/goes-view/blob/master/images/GOES_east-west_vis.gif" width="600" >

The sub-pixel orthorectification method applied here uses the GOES satellite's known orbital position (from ABI product NetCDF metadata) to compute the intersection of line of sight (LOS) vectors with a DEM surface. This method is "sub-pixel" because the DEM spatial resolution is much finer (~30 m) than the GOES ABI image resolution (> 2 km). This effectively drapes ABI pixels (and their respective radiance or brightness tempreature values) over the terrain at the DEM's finer resolution.

The figure below (from the GOES ABI ATBD) illustrates the satellite's viewing geometry. The orthorectification method developed here modifies point P on the Earth's surface using information from a DEM about its elevation relative to the reference ellipsoid.

<img src="https://github.com/spestana/goes-view/blob/master/images/ABIgrid.png" width="600" >

These python scripts and jupyter notebooks help with downloading GOES ABI data (wrapper around goespy) from AWS, creating timeseries of GOES ABI brightness temperature for point locations, and orthorectifying (terrain/parallax correction) GOES ABI imagery using a DEM (here specifically for part of the Sierra Nevada in California). 

---

### Setting up the conda environment:

```bash
conda env create -f environment.yml
conda activate goes-linux; ipython kernel install --user --name goes-linux
```

---

### download-goes.py

Downloads GOES-16 or GOES-17 products/bands, requires command line arguments (wrapper around goespy.Downloader.ABI_Downloader()):

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

Functions for orthorectifying GOES-R ABI imagery using a DEM. Produces an orthorectified NetCDF at the spatial resolution of the input DEM.

This method uses the GOES satellite's known orbital position (from ABI product NetCDF metadata) to compute the intersection of line of sight (LOS) vectors with a DEM surface.

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

See the [goes-orthorectify](https://github.com/spestana/goes-view/blob/master/goes-orthorectify.ipynb) notebook for an example of orthorectifying a single GOES ABI image.

See the [goes-orthorectify-aster.py](https://github.com/spestana/goes-view/blob/master/goes-orthorectify-aster.py) script for an example of orthorectifying a batch of GOES ABI images.

#### Flowchart:

(flowchart diagram is a work in progress)

![goes-ortho-flowchart](/images/goes-ortho-flowchart.png)


---
---