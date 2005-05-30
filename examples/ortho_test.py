from matplotlib.toolkits.basemap import Basemap
from pylab import *
import cPickle

# example of filled contour plot on orthographic projection.

# read in data on lat/lon grid.
datadict = cPickle.load(open('500hgt.pickle','rb'))
hgt = datadict['data']; lons = datadict['lons']; lats = datadict['lats']
lons, lats = meshgrid(lons, lats)

# setup of orthographic basemap
m = Basemap(lat_0=50.,lon_0=-120.,projection='ortho')
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.8,0.8],frameon=False)
# make a filled contour plot.
x, y = m(lons, lats)
levels, colls = m.contour(x,y,hgt,15,linewidths=0.5,colors='k')
levels, colls = m.contourf(x,y,hgt,15,cmap=cm.jet,colors=None)
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawmapboundary()
# draw parallels and meridians.
m.drawparallels(arange(-80.,90,20.))
m.drawmeridians(arange(0.,360.,20.))
title('Orthographic Filled Contour Demo')
show()

