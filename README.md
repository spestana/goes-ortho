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
---