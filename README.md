# "The Mountains Are Calling, and I must GOES" 
> - [John Muir](https://en.wikipedia.org/wiki/John_Muir) with satellite remote sensing

##### Geospatial Data Analysis Project - Steven Pestana, 2019

## Background:
The latest generation of geostationary-orbiting weather satellites make frequent observations (5-60 min) in the visible and IR wavelengths at moderate resolutions (500-~2000m) (GOES-16/17 ABI, Himawari 8/9 AHI). This lends themselves to be used to fill temporal gaps between ~12-hour repeat observations by polar-orbiting spacecraft (Aqua/Terra MODIS, SNPP VIIRS, etc).
However, their geostationary orbits mean that outside of their sub-satellite-point on the equator, all other view angles are off-nadir, and due to the Earth's curvature in view, actual pixel sizes increase to >6 km towards the planet's limb.
Additionally, when viewing complex terrain such as the mountains of western CONUS, parallax affects the apparent position of the variable topography. Some portions of the ground suface become obscured from view completely by surrounding terrain.

Before using observations from these instruments for observing the land surface for applications where fine-scale variability is important, corrections could be applied to the geostationary satellite imagery to try and account for the off-nadir view angle and topographic effects.

## Objective(s):
* Use python to access and work with GOES-16/17 ABI data from AWS
* Overlay the coarse resolution satellite data on top of a higher resolution DEM of the Sierra Nevada
* Mask out areas not visible to the GOES satellites as determined by their orbital position and local terrain


## Data:
* GOES-16/17 ABI (visible through longwaveIR)
* SRTM CONUS DEM at 30 m resolution
  
## Packages:
* ([goes-py](https://github.com/palexandremello/goes-py)) to access GOES data on AWS 
* rasterio


## Approach:

## Expected Outcomes:
