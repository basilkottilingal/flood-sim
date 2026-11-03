# Data 

  1. Topography contour can be downloaded from opentopography.org
     (File is put in .gitignore due to heavy size of the Geotiff topographiy data).
     You may use API to download raster from opentopography.
```python
import os, requests
params = {
    "demtype": "COP30",
    "south": 8.29, "north": 12.79,
    "west": 74.85,  "east": 77.40,
    "outputFormat": "GTiff",
    "API_Key": os.environ["OT_API_KEY"],
}
r = requests.get("https://portal.opentopography.org/API/globaldem", params=params, timeout=120)
r.raise_for_status()
with open("patch.tif", "wb") as f:
    f.write(r.content)
```
  2. (Missing rainfall data)
  3. [Kerala state boundary](https://github.com/opendatakerala/map.opendatakerala.org/blob/main/data/kerala_lsg_mapdata.geojson)

# References
  1. ["tiff 6.0 specification" format](https://www.itu.int/itudoc/itu-t/com16/tiff-fx/docs/tiff6.pdf)
  2. [GeoTIFF library](https://github.com/OSGeo/libgeotiff)
  3. [libtiff](https://libtiff.gitlab.io/libtiff/)
  4. Lempel–Ziv–Welch lossless compression
<a id="1">[1]</a>
Welch, T. "A technique for high-performance data compression." Computer 17.6 (1984): 8-19.
  5. [An LZW code in libtiff library](https://gitlab.com/libtiff/libtiff/-/blob/master/libtiff/tif_lzw.c)
  6. [wsServer](https://github.com/Theldus/wsServer/tree/master)
  7. [Basilisk](basilisk.fr)
  8. [GeoJSON format](https://www.rfc-editor.org/info/rfc7946/)
  9. [photogrammetry](https://github.com/awesome-photogrammetry/awesome-photogrammetry)


# Multilayer FVM solver for flood simulation

  (Planned only. Not yet implemented)
  Open source code [Basilisk](basilisk.fr) is used for simulation. 
 

# requirements (or limits) for running the code (which relies onin-built `cwrap`)
  1. gcc [fixme : update min version requirement]
  2. c99/c11 (not strict, `cwrap` translator can even pass forward gnu grammars which is a superset of std-c99)
(K&R funtion definition of C99 won't be parsed)
  3. `cwrap` : single build or static linking (with possible multiple stages of compilation). Warning : won't support creation of shared libraries
  4. Linux only. (MacOS can't run yet. 
  5. Cannot translate to a single native c file if there are multiple compilation stages.
