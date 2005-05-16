from matplotlib.toolkits.basemap import Basemap
import cPickle
from pylab import *
# read in topo data from pickle (on a regular lat/lon grid)
topodict = cPickle.load(open('etopo20.pickle','rb'))
etopo = topodict['data']; lons = topodict['lons']; lats = topodict['lats']
# create Basemap instance for Robinson projection.
m = Basemap(projection='robin',lon_0=0.5*(lons[0]+lons[-1]))
# compute native map projection coordinates for lat/lon grid.
lons, lats = meshgrid(lons,lats)
x,y = m(lons,lats)
# create figure with same aspect ratio as map.
fig=figure(figsize=(10,m.aspect*10))
fig.add_axes([0.1,0.1,0.8,0.8],frameon=False)
# make filled contour plot.
levels, colls = m.contourf(x,y,etopo,30,cmap=cm.jet,colors=None)
# draw coastlines.
m.drawcoastlines()
# draw a line around the map region.
m.drawmapboundary()
# draw parallels and meridians.
m.drawparallels(arange(-90.,120.,30.))
m.drawmeridians(arange(0.,420.,60.))
title('Robinson Projection')
show()

