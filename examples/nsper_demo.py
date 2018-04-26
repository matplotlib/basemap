from __future__ import (absolute_import, division, print_function)

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import sys

def get_input(prompt):
    if sys.hexversion > 0x03000000:
        return input(prompt)
    else:
        return raw_input(prompt)

# create Basemap instance for Near-Sided Perspective (satellite view) projection.
lon_0 = float(get_input('enter reference longitude (lon_0):'))
lat_0 = float(get_input('enter reference latitude (lat_0):'))
h = float(get_input('enter altitude of camera in km (h):'))
h=h*1000.

# map with continents drawn and filled.
fig = plt.figure()
m = Basemap(projection='nsper',lon_0=lon_0,lat_0=lat_0,satellite_height=h,resolution='l')
m.drawcoastlines()
m.drawmapboundary(fill_color='aqua')
m.fillcontinents(color='coral',lake_color='aqua')
m.drawcountries()
m.drawstates()
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,10.))
m.drawmeridians(np.arange(0.,420.,20.))
m.drawmapboundary(fill_color='aqua')
plt.title('Near-Sided Perspective Map Centered on Lon=%s, Lat=%s, H=%g' %\
    (lon_0,lat_0,h/1000.),fontsize=10)

fig = plt.figure()
m1 = Basemap(projection='nsper',lon_0=lon_0,lat_0=lat_0,satellite_height=h,resolution=None)
ax = fig.add_axes([0.1,0.1,0.8,0.8], facecolor='k')
# plot just upper right quadrant (coordinates determined from global map).
m = Basemap(projection='nsper',lon_0=lon_0,lat_0=lat_0,satellite_height=h,resolution='l',llcrnrx=0.,llcrnry=0.,urcrnrx=m1.urcrnrx/2.,urcrnry=m1.urcrnry/2.)
m.drawcoastlines()
m.drawmapboundary(fill_color='aqua')
m.fillcontinents(color='coral',lake_color='aqua')
m.drawcountries()
m.drawstates()
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
plt.title('Near-Sided Perspective Map Centered on Lon=%s, Lat=%s, H=%g' %\
    (lon_0,lat_0,h/1000.),fontsize=10)

plt.show()
