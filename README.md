# "The Mountains Are Calling, and I Must GOES" 
> - [John Muir](https://en.wikipedia.org/wiki/John_Muir) if he had access to satellite remote sensing probably

##### Geospatial Data Analysis Project - Steven Pestana, 2019

## Background:
The latest generation of geostationary-orbiting weather satellites make frequent observations (5-60 min) in the visible and IR wavelengths at moderate resolutions (500-~2000m) (GOES-16/17 ABI, Himawari 8/9 AHI). This lends themselves to be used to fill temporal gaps between ~12-hour repeat observations by polar-orbiting spacecraft (Aqua/Terra MODIS, SNPP VIIRS, etc).
However, their geostationary orbits mean that outside of their sub-satellite-point on the equator, all other view angles are off-nadir, and due to the Earth's curvature in view, actual pixel sizes increase to >6 km towards the planet's limb.
Additionally, when viewing complex terrain such as the mountains of western CONUS, parallax affects the apparent position of the variable topography. Some portions of the ground suface become obscured from view completely by surrounding terrain.

Before using observations from these instruments for observing the land surface for applications where fine-scale variability is important, corrections could be applied to the geostationary satellite imagery to try and account for the off-nadir view angle and topographic effects.

## Objective(s):
* Use python to access and work with GOES-16/17 ABI data from AWS
* Overlay the coarse resolution satellite data on top of a higher resolution DEM of the Tuolumne River Watershed in the Sierra Nevada
* Mask out areas not visible to the GOES satellites as determined by their orbital position and local terrain

## Data:
* GOES-16/17 ABI (visible through longwaveIR)
* SRTM CONUS DEM at 30 m resolution
  
## Packages:
* [goes-py](https://github.com/palexandremello/goes-py) to access GOES data on AWS 
* [elevation](https://pypi.org/project/elevation/) to access SRTM DEMs
* [rasterio](https://pypi.org/project/rasterio/) for analysis and plotting
* perhaps also [gdal](https://www.gdal.org/) utilities

## Approach:
* use goes-py to access recent wintertime GOES-16 and -17 observations of CONUS
* use elevation to access SRTM3 DEMs of the Tuolumne River Watershed in the Sierra Nevada
* clip GOES imagery to the DEM areas
* determine slope angles/aspects hidden from view
* mask out the terrain hidden from view from the GOES imagery
* export a GeoTiff of this new "product"

## Expected Outcomes:
* Produce a masked GOES product that shows which areas were NOT visible to the satellite


