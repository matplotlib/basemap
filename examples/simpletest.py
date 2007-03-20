from matplotlib.toolkits.basemap import Basemap
from pylab import load, show, title, meshgrid, cm, arange
# read in topo data (on a regular lat/lon grid)
etopo=load('etopo20data.gz')
lons=load('etopo20lons.gz')
lats=load('etopo20lats.gz')
# create Basemap instance for Robinson projection.
m = Basemap(projection='robin',lon_0=0.5*(lons[0]+lons[-1]))
# make filled contour plot.
x, y = m(*meshgrid(lons, lats))
cs = m.contourf(x,y,etopo,30,cmap=cm.jet)
# draw coastlines.
m.drawcoastlines()
# draw a line around the map region.
m.drawmapboundary()
# draw parallels and meridians.
m.drawparallels(arange(-60.,90.,30.),labels=[1,0,0,0])
m.drawmeridians(arange(0.,420.,60.),labels=[0,0,0,1])
# add a title.
title('Robinson Projection')
show()

