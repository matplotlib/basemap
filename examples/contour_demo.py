from matplotlib.toolkits.basemap import Basemap, shiftgrid
from pylab import *
import cPickle
from matplotlib.numerix import ma

# same on contour_demo.py, but data is not interpolated to native
# projection grid.

# read in data on lat/lon grid.
datadict = cPickle.load(open('500hgt.pickle','rb'))
hgt = datadict['data']; lons = datadict['lons']; lats = datadict['lats']
# shift data so lons go from -180 to 180 instead of 0 to 360.
hgt,lons = shiftgrid(180.,hgt,lons,start=False)
lons, lats = meshgrid(lons, lats)

# set up map projection (azimuthal equidistant).
m = Basemap(llcrnrlon=-135.,llcrnrlat=-20.,urcrnrlon=45.,urcrnrlat=-20.,
             resolution='c',area_thresh=10000.,projection='aeqd',
             lat_0=90.,lon_0=-90.)

# compute x,y of lat/lon grid.
x, y = m(lons, lats)
# make for points over projection horizon.
mask = logical_or(greater(x,1.e20),greater(y,1.e20))

xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.7,0.7],frameon=False)
# make a filled contour plot.
levels, colls = m.contour(x,y,hgt,15,linewidths=0.5,colors='k',badmask=mask)
levels, colls = m.contourf(x,y,hgt,15,cmap=cm.jet,colors=None,badmask=mask)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
# draw coastlines and political boundaries.
m.drawcoastlines()
# draw parallels and meridians.
parallels = arange(0.,80,20.)
m.drawparallels(parallels,labels=[1,1,1,1])
meridians = arange(10.,360.,20.)
m.drawmeridians(meridians,labels=[0,0,1,1])

title('Azimuthal Equidistant Filled Contour Demo',y=1.075)
show()

# setup of orthographic basemap
m = Basemap(resolution='c',area_thresh=10000.,projection='ortho',\
            lat_0=50.,lon_0=-120.)

x,y = m(lons,lats)
# make for points over projection horizon.
mask = logical_or(greater(x,1.e20),greater(y,1.e20))

xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.7,0.7],frameon=False)
# make a filled contour plot.
levels, colls = m.contour(x,y,hgt,15,linewidths=0.5,colors='k',badmask=mask)
levels, colls = m.contourf(x,y,hgt,15,cmap=cm.jet,colors=None,badmask=mask)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawmapboundary()
# draw parallels and meridians.
parallels = arange(0.,80,20.)
m.drawparallels(parallels)
meridians = arange(0.,360.,20.)
m.drawmeridians(meridians)

title('Orthographic Filled Contour Demo')
show()
