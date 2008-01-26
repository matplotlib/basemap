from mpl_toolkits.basemap import Basemap
from pylab import title, show, arange, figure

# map with continents drawn and filled.
fig = figure()
lon_0=57
m1 = Basemap(projection='geos',lon_0=lon_0,rsphere=(6378137.00,6356752.3142),resolution=None)
ax = fig.add_axes([0.1,0.1,0.8,0.8],axisbg='k')
# plot just upper right quadrant.
m = Basemap(projection='geos',lon_0=lon_0,rsphere=(6378137.00,6356752.3142),resolution='l',llcrnrx=m1.urcrnrx/2.,llcrnry=m1.urcrnry/2.,urcrnrx=m1.urcrnrx,urcrnry=m1.urcrnry)
print m.projparams
m.drawcoastlines()
m.drawmapboundary(fill_color='aqua')
m.fillcontinents(color='coral',lake_color='aqua')
m.drawcountries()
# draw parallels and meridians.
m.drawparallels(arange(-90.,120.,30.))
m.drawmeridians(arange(0.,420.,60.))
m.drawmapboundary()
title('Geostationary Map Centered on Lon=%s' % lon_0)

fig = figure()
m1 = Basemap(projection='ortho',lon_0=lon_0,lat_0=10,resolution=None)
ax = fig.add_axes([0.1,0.1,0.8,0.8],axisbg='k')
# plot just upper right quadrant.
m = Basemap(projection='ortho',lon_0=lon_0,lat_0=10,resolution='l',llcrnrx=m1.urcrnrx/2.,llcrnry=m1.urcrnry/2.,urcrnrx=m1.urcrnrx,urcrnry=m1.urcrnry)
print m.projparams
m.drawcoastlines()
m.drawmapboundary(fill_color='aqua')
m.fillcontinents(color='coral',lake_color='aqua')
m.drawcountries()
# draw parallels and meridians.
m.drawparallels(arange(-90.,120.,30.))
m.drawmeridians(arange(0.,420.,60.))
m.drawmapboundary()
title('Orthographic Map Centered on Lon=%s' % lon_0)

show()
