GOES ABI Terrain Correction
===========================

.. image:: https://zenodo.org/badge/281728618.svg
   :target: https://zenodo.org/badge/latestdoi/281728618
   :alt: DOI


.. image:: https://github.com/spestana/goes-ortho/actions/workflows/tests.yml/badge.svg
   :target: https://github.com/spestana/goes-ortho/actions/workflows/tests.yml
   :alt: tests


.. image:: https://github.com/spestana/goes-ortho/actions/workflows/build-docs.yml/badge.svg
   :target: https://github.com/spestana/goes-ortho/actions/workflows/build-docs.yml
   :alt: docs

.. image:: https://img.shields.io/conda/vn/conda-forge/goes-ortho.svg
   :target: https://anaconda.org/conda-forge/goes-ortho
   :alt: conda-forge

.. image:: https://img.shields.io/pypi/v/goes-ortho.svg
   :target: https://pypi.python.org/pypi/goes-ortho
   :alt: pypi


Orthorectifying GOES ABI imagery at sub-pixel resolution

.. image:: https://raw.githubusercontent.com/spestana/goes-ortho/main/docs/images/GOES-terrain-correction.gif
   :width: 600px

----

The latest generation of geostationary-orbiting weather satellites make frequent observations (5-60 min) in the visible and IR wavelengths at moderate resolutions (500-~2000m) (GOES-16/17 ABI, Himawari 8/9 AHI). This lends themselves to be used to fill temporal gaps between ~12-hour repeat observations by polar-orbiting spacecraft (Aqua/Terra MODIS, SNPP VIIRS, etc).
However, their geostationary orbits mean that outside of their sub-satellite-point on the equator, all other view angles are off-nadir, and due to the Earth's curvature in view, actual pixel sizes increase to >6 km towards the planet's limb.


Additionally, when viewing complex terrain such as the mountains of western CONUS, parallax affects the apparent position of the variable topography. Some portions of the ground surface may even become obscured from view completely by surrounding steep terrain with poleward-facing aspects (north-facing aspects in the Northern Hemisphere, south-facing aspects in the Southern Hemisphere).

Before using observations from these instruments for observing the land surface over mountains, orthorectification is needed to try and account for the off-nadir view angles and topographic effects in the geostationary satellite imagery.

The terrain parallax is especially visually apparent when flipping between GOES-East and GOES-West views of a mountain range like the Sierra Nevada here:

.. image:: https://raw.githubusercontent.com/spestana/goes-ortho/main/docs/images/GOES_east-west_vis.gif
   :width: 600px

The sub-pixel orthorectification method applied here uses the GOES satellite's known orbital position (from ABI product NetCDF metadata) to compute the intersection of line of sight (LOS) vectors with a DEM surface. This method is **"sub-pixel"** because the DEM spatial resolution can be much finer (here I've used ~30 m, 1 arc-second SRTM DEM) than the GOES ABI image resolution (> 2 km). This effectively drapes ABI pixels (and their respective radiance or brightness temperature values) over the terrain at the DEM's finer resolution.

The figure below (from the GOES ABI ATBD) illustrates the satellite's viewing geometry. The orthorectification method developed here modifies point P on the Earth's surface using information from a DEM about its elevation relative to the reference ellipsoid.

.. image:: https://raw.githubusercontent.com/spestana/goes-ortho/main/docs/images/ABIgrid.png
   :width: 600px

These python scripts and jupyter notebooks help with downloading GOES ABI data from AWS (wrapper around the `goespy <https://github.com/palexandremello/goes-py>`_ library), creating timeseries of GOES ABI brightness temperature for point locations, and orthorectifying (terrain correction) GOES ABI imagery using a DEM (here specifically for part of the Sierra Nevada in California).


.. image:: https://raw.githubusercontent.com/spestana/goes-ortho/main/docs/examples/make_abi_timeseries_example_plot.png
   :width: 600px


----


download-goes.py
----------------

Downloads GOES-16 or GOES-17 products/bands, requires command line arguments (wrapper around ``goespy.Downloader.ABI_Downloader()``):

Usage:
~~~~~~

.. code-block:: bash

   python ./download-goes.py --bucket <S3-BUCKET> --year <YEAR> --month <MONTH> --days <START DAY> <END DAY> --product <ABI PRODUCT CODE> --channel <ABI CHANNEL> --bounds <MIN_LON> <MIN_LAT> <MAX_LON> <MAX_LAT> --dir <DESTINATION DIRECTORY>

Examples:
~~~~~~~~~

This will download the GOES-16 ABI Level-1b Radiance (CONUS) product for channel/band 14, for January 1-2 2020. The NetCDF files will be cropped to within latitudes 30 - 50 and longitudes -125 - -105, and saved in /storage/spestana/scratchspace.

.. code-block:: bash

   python ./download-goes.py --bucket noaa-goes16 --year 2020 --month 2 --days 1 2 --product ABI-L1b-RadC --channel C14 --bounds -125 30 -105 50 --dir /storage/spestana/scratchspace

We can do the same command with short flag names:

.. code-block:: bash

   python ./download-goes.py -B noaa-goes16 -Y 2020 -M 1 -D 1 2 -p ABI-L1b-RadC -c C14 -b -125 30 -105 50 -d /storage/spestana/scratchspace

----


Flowchart:
~~~~~~~~~~

(flowchart diagram is a work in progress)

.. image:: https://raw.githubusercontent.com/spestana/goes-ortho/main/docs/images/goes-ortho-flowchart.png
   :width: 600px

