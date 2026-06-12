import rasterio
import numpy as np
import matplotlib.pyplot as plt
from rasterio.enums import Resampling

# ----------------------------
# CONFIG
# ----------------------------
tif_path = "topography_ESA_Copernicus_30m_resolution.tif"
downsample_factor = 10   # increase this if file is large (e.g., 20 or 50)
num_contours = 15

# ----------------------------
# READ + DOWNSAMPLE DEM
# ----------------------------

h = 0
w = 0

with rasterio.open (tif_path) as src:
  h = src.height
  w = src.width
  dem = src.read (
    1,
    out_shape = ( int (h / downsample_factor), int (w / downsample_factor) ),
    resampling = Resampling.bilinear
  )

  nodata = src.nodata

# ----------------------------
# CLEAN DATA
# ----------------------------
dem = dem.astype (float)

if nodata is not None:
  dem[dem == nodata] = np.nan

# ----------------------------
# PLOT CONTOURS
# ----------------------------
base_width = 10  # inches
base_height = base_width * h / w
plt.figure (figsize = (base_width, base_height))

contour = plt.contour (
  dem,
  levels = num_contours,
  cmap="terrain"
)

plt.clabel (contour, inline=True, fontsize=8)

plt.title ("DEM Contour Map (Downsampled)")
plt.xlabel ("X (downsampled pixels)")
plt.ylabel ("Y (downsampled pixels)")
plt.gca ().invert_yaxis ()

plt.show ()
