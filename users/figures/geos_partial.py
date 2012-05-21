from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
fig = plt.figure()
# global geostationary map centered on lon_0
lon_0=57.
# resolution = None means don't process the boundary datasets.
m1 = Basemap(projection='geos',lon_0=lon_0,resolution=None)
# add an axes with a black background
ax = fig.add_axes([0.1,0.1,0.8,0.8],axisbg='k')
# plot just upper right quadrant (corners determined from global map).
# keywords llcrnrx,llcrnry,urcrnrx,urcrnry used to define the lower
# left and upper right corners in map projection coordinates.
# llcrnrlat,llcrnrlon,ucrnrlon,urcrnrlat could be used to define
# lat/lon values of corners - but this won't work in cases such as this
# where one of the corners does not lie on the earth.
m = Basemap(projection='geos',lon_0=lon_0,resolution='l',\
    llcrnrx=0.,llcrnry=0.,urcrnrx=m1.urcrnrx/2.,urcrnry=m1.urcrnry/2.)
m.drawcoastlines()
m.drawmapboundary(fill_color='aqua')
m.fillcontinents(color='coral',lake_color='aqua')
m.drawcountries()
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,360.,60.))
m.drawmapboundary()
plt.title('Geostationary Map Showing A Quadrant of the Globe')
plt.show()
