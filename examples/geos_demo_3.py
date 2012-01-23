from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
lon_0=57
# global geos map
m1 = Basemap(projection='geos',lon_0=lon_0,rsphere=(6378137.00,6356752.3142),resolution=None)
ax = fig.add_axes([0.1,0.1,0.8,0.8],axisbg='k')
# plot just upper right quadrant (coordinates determined from global map).
m = Basemap(projection='geos',lon_0=lon_0,rsphere=(6378137.00,6356752.3142),resolution='l',llcrnrx=0.,llcrnry=0.,urcrnrx=m1.urcrnrx/2.,urcrnry=m1.urcrnry/2.)
m.drawcoastlines()
m.drawmapboundary(fill_color='aqua')
m.fillcontinents(color='coral',lake_color='aqua')
m.drawcountries()
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
plt.title('Geostationary Map Centered on Lon=%s' % lon_0)

fig = plt.figure()
# global ortho map
lat_0=10.
m1 = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution=None)
ax = fig.add_axes([0.1,0.1,0.8,0.8],axisbg='k')
# plot just upper right quadrant (corners determined from global map).
m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution='l',llcrnrx=0.,llcrnry=0.,urcrnrx=m1.urcrnrx/2.,urcrnry=m1.urcrnry/2.)
m.drawcoastlines()
m.drawmapboundary(fill_color='aqua')
m.fillcontinents(color='coral',lake_color='aqua')
m.drawcountries()
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
plt.title('Orthographic Map Centered on Lon=%s, Lat=%s' % (lon_0,lat_0))

plt.show()
