from matplotlib.toolkits.basemap import Basemap
import cPickle
from pylab import *
# read in topo data from pickle (on a regular lat/lon grid)
topodict = cPickle.load(open('etopo20.pickle','rb'))
etopo = topodict['data']; lons = topodict['lons']; lats = topodict['lats']
# create Basemap instance for Orthographic (satellite view) projection.
lon_0 = float(raw_input('enter reference longitude (lon_0):'))
lat_0 = float(raw_input('enter reference latitude (lat_0):'))
fillcont = int(raw_input('fill continents? (1 for yes, 0 for no):'))
m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0)
# compute native map projection coordinates for lat/lon grid.
lons, lats = meshgrid(lons,lats)
x,y = m(lons,lats)
# create figure with same aspect ratio as map.
fig=m.createfigure().add_axes([0.05,0.05,0.9,0.9])
# make filled contour plot.
cs = m.contourf(x,y,etopo,30,cmap=cm.jet)
# draw coastlines.
m.drawcoastlines()
# draw a line around the map region.
m.drawmapboundary()
if fillcont:
    m.fillcontinents()
# draw parallels and meridians.
m.drawparallels(arange(-90.,120.,30.))
m.drawmeridians(arange(0.,420.,60.))
title('Orthographic Map Centered on Lon=%s, Lat=%s' % (lon_0,lat_0))
show()
