from __future__ import print_function

import time
import pickle

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


# Create figure.
fig = plt.figure()

# Create Basemap instance:
# - Use 'full' resolution coastlines.
# - Make sure that countries and rivers are loaded.
t0 = time.time()
bmap1 = Basemap(width=920000, height=1100000, resolution="f",
                projection="tmerc", lon_0=-4.2, lat_0=54.6)
bmap1.drawcountries()
bmap1.drawrivers()
t1 = time.time()
print("{0:.3f} secs to plot with a Basemap instance created at runtime".format(t1 - t0))

# Clear the figure.
plt.clf()

# Pickle the class instance.
with open("map.pickle", "wb") as fd:
    pickle.dump(bmap1, fd, protocol=-1)

# Read pickle back in and plot it again (should be much faster):
# - Draw coastlines and fill continents and lakes.
# - Draw political boundaries and rivers.
# - Draw parallels and meridians.
# - Draw map boundary and fill map background.
t0 = time.time()
with open("map.pickle", "rb") as fd:
    bmap2 = pickle.load(fd)
bmap2.drawcoastlines()
bmap2.fillcontinents(color="coral", lake_color="aqua")
bmap2.drawcountries(linewidth=1)
bmap2.drawrivers(color="b")
bmap2.drawparallels(np.arange(48, 65, 2), labels=[1, 1, 0, 0])
bmap2.drawmeridians(np.arange(-12, 13, 2), labels=[0, 0, 1, 1])
bmap2.drawmapboundary(fill_color="aqua")
t1 = time.time()
print("{0:.3f} secs to plot with a pickled Basemap instance".format(t1 - t0))

plt.title("High-Res British Isles", y=1.04)
plt.show()
