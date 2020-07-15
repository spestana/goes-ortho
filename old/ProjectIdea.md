# Project Idea

The latest generation of geostationary-orbiting weather satellites make frequent observations (5-60 min) in the visible and IR wavelengths at moderate resolutions (500-~2000m) (GOES-16/17 ABI, Himawari 8/9 AHI). This lends themselves to be used to fill temporal gaps between ~12-hour repeat observations by polar-orbiting spacecraft (Aqua/Terra MODIS, SNPP VIIRS, etc).
However, their geostationary orbits mean that outside of their sub-satellite-point on the equator, all other view angles are off-nadir, and due to the Earth's curvature in view, actual pixel sizes increase to >6 km towards the planet's limb.
Additionally, when viewing complex terrain such as the mountains of western CONUS, parallax affects the apparent position of the variable topography. Some portions of the ground suface become obscured from view completely by surrounding terrain.

Before using observations from these instruments for observing the land surface for applications where fine-scale variability is important, corrections could be applied to the geostationary satellite imagery to try and account for the off-nadir view angle and topographic effects.

For this project I want to learn how to work with GOES-16 ABI data in python, try to better understand the above caveats to using its data over mountains, and test applying something like a sub-pixel orthorectification to the imagery over part of the Sierra Nevada.

## Packages:
* xarray and dask to work with large files
* cartopy (or another plotting tool for displaying/projecting geospatial rasters)
* something to interface with GOES data on AWS ([I've been using this so far](https://github.com/palexandremello/goes-py)), or other data sources (NOAA servers?)

## Data:
* GOES-16 ABI (visible through longwaveIR)
  * _could extend to GOES-17 if more of its data becomes available_
* high resolution DEMs
  * 30m SRTM CONUS DEM (2000)
  * ~5m Airborne Snow Observatory DEM of Tuolumne Basin, Sierra Nevada (2015)
