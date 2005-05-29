from matplotlib.toolkits.basemap import Basemap, shiftgrid
from pylab import *
import cPickle, sys
from matplotlib.numerix import ma

# examples of filled contour plots on map projections.

# read in data on lat/lon grid.
datadict = cPickle.load(open('500hgt.pickle','rb'))
hgt = datadict['data']; lons = datadict['lons']; lats = datadict['lats']
# shift data so lons go from -180 to 180 instead of 0 to 360.
hgt,lons = shiftgrid(180.,hgt,lons,start=False)
lons, lats = meshgrid(lons, lats)

# setup of orthographic basemap
m = Basemap(resolution='c',area_thresh=10000.,projection='ortho',\
            lat_0=50.,lon_0=-120.)
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.7,0.7],frameon=False)
# make a filled contour plot.
x, y = m(lons, lats)
levels, colls = m.contour(x,y,hgt,15,linewidths=0.5,colors='k')
levels, colls = m.contourf(x,y,hgt,15,cmap=cm.jet,colors=None)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawmapboundary()
# draw parallels and meridians.
parallels = arange(-80.,90,20.)
m.drawparallels(parallels)
meridians = arange(0.,360.,20.)
m.drawmeridians(meridians)
title('Orthographic Filled Contour Demo')
show()

