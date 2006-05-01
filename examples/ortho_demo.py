from matplotlib.toolkits.basemap import Basemap
from pylab import *
# create Basemap instance for Orthographic (satellite view) projection.
lon_0 = float(raw_input('enter reference longitude (lon_0):'))
lat_0 = float(raw_input('enter reference latitude (lat_0):'))
m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0)
# plot land-sea mask.
rgba_land = (0,1,0,1)
rgba_ocean = (0,0,1,1)
# lakes=True means plot inland lakes with ocean color.
m.drawlsmask(rgba_land, rgba_ocean, lakes=True)
# draw parallels and meridians.
m.drawparallels(arange(-90.,120.,30.))
m.drawmeridians(arange(0.,420.,60.))
m.drawmapboundary()
title('Orthographic Map Centered on Lon=%s, Lat=%s' % (lon_0,lat_0))
show()
