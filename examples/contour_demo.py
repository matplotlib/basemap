from matplotlib.toolkits.basemap import Basemap, shiftgrid
from pylab import *

# examples of filled contour plots on map projections.

# read in data on lat/lon grid.
hgt = array(load('500hgtdata.gz'),'d')
lons = array(load('500hgtlons.gz'),'d')
lats = array(load('500hgtlats.gz'),'d')
# shift data so lons go from -180 to 180 instead of 0 to 360.
hgt,lons = shiftgrid(180.,hgt,lons,start=False)
lons, lats = meshgrid(lons, lats)

# setup of mollweide basemap
m = Basemap(resolution='c',projection='moll',lon_0=0)
fig = m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# make a filled contour plot.
x, y = m(lons, lats)
CS = m.contour(x,y,hgt,15,linewidths=0.5,colors='k')
CS = m.contourf(x,y,hgt,15,cmap=cm.jet)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawmapboundary()
m.fillcontinents()
# draw parallels and meridians.
parallels = arange(-60.,90,30.)
m.drawparallels(parallels,labels=[1,0,0,0])
meridians = arange(-360.,360.,30.)
m.drawmeridians(meridians)
title('Mollweide Filled Contour Demo')
show()

# set up Robinson map projection.
m = Basemap(resolution='c',projection='robin',lon_0=0)
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# make a filled contour plot.
x, y = m(lons, lats)
CS = m.contour(x,y,hgt,15,linewidths=0.5,colors='k')
CS = m.contourf(x,y,hgt,15,cmap=cm.jet)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawmapboundary()
m.fillcontinents()
# draw parallels and meridians.
parallels = arange(-60.,90,30.)
m.drawparallels(parallels,labels=[1,0,0,0])
meridians = arange(-360.,360.,60.)
m.drawmeridians(meridians,labels=[0,0,0,1])
title('Robinson Filled Contour Demo')
show()

# set up map projection (azimuthal equidistant).
m = Basemap(llcrnrlon=-135.,llcrnrlat=-20.,urcrnrlon=45.,urcrnrlat=-20.,
             resolution='c',projection='aeqd',
             lat_0=90.,lon_0=-90.)
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# make a filled contour plot.
x, y = m(lons, lats)
CS = m.contour(x,y,hgt,15,linewidths=0.5,colors='k')
CS = m.contourf(x,y,hgt,15,cmap=cm.jet)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawmapboundary()
m.fillcontinents()
# draw parallels and meridians.
parallels = arange(0.,80,20.)
m.drawparallels(parallels,labels=[0,0,1,1])
meridians = arange(10.,360.,20.)
m.drawmeridians(meridians,labels=[1,1,1,1])
title('Azimuthal Equidistant Filled Contour Demo',y=1.075)
show()

# setup of orthographic basemap
m = Basemap(resolution='c',projection='ortho',\
            lat_0=50.,lon_0=-120.)
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# make a filled contour plot.
x, y = m(lons, lats)
CS = m.contour(x,y,hgt,15,linewidths=0.5,colors='k')
CS = m.contourf(x,y,hgt,15,cmap=cm.jet)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
# draw coastlines and political boundaries.
m.drawcoastlines()
m.fillcontinents()
m.drawmapboundary()
# draw parallels and meridians.
parallels = arange(-80.,90,20.)
m.drawparallels(parallels)
meridians = arange(0.,360.,20.)
m.drawmeridians(meridians)
title('Orthographic Filled Contour Demo')
show()

