from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

# create Basemap instance for Near-Sided Perspective (satellite view) projection.
lon_0 = float(raw_input('enter reference longitude (lon_0):'))
lat_0 = float(raw_input('enter reference longitude (lat_0):'))
h = float(raw_input('enter altitude of camera in km (h):'))
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
m.drawmapboundary()
plt.title('Near-Sided Perspective Map Centered on Lon=%s, Lat=%s, H=%g' %\
    (lon_0,lat_0,h),fontsize=10)
plt.show()
