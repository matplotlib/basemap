from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import sys

# examples of filled contour plots on map projections.

# read in data on lat/lon grid.
hgt = np.loadtxt('500hgtdata.gz')
lons = np.loadtxt('500hgtlons.gz')
lats = np.loadtxt('500hgtlats.gz')
lons, lats = np.meshgrid(lons, lats)

# create new figure
fig=plt.figure()
# setup of sinusoidal basemap
m = Basemap(resolution='c',projection='sinu',lon_0=0)
# make a filled contour plot.
# create contour lines
CS1 = m.contour(lons,lats,hgt,15,linewidths=0.5,colors='k',latlon=True)
# fill between contour lines.
CS2 =\
m.contourf(lons,lats,hgt,CS1.levels,cmap=plt.cm.jet,extend='both',latlon=True)
m.colorbar(CS2) # draw colorbar
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawmapboundary()
m.fillcontinents()
# draw parallels and meridians.
parallels = np.arange(-60.,90,30.)
m.drawparallels(parallels,labels=[1,0,0,0])
meridians = np.arange(-360.,360.,30.)
m.drawmeridians(meridians)
plt.title('Sinusoidal Filled Contour Demo')
sys.stdout.write('plotting with sinusoidal basemap ...\n')

# create new figure
fig=plt.figure()
# setup of mollweide basemap
m = Basemap(resolution='c',projection='moll',lon_0=0)
# make a filled contour plot.
CS1 = m.contour(lons,lats,hgt,15,linewidths=0.5,colors='k',latlon=True)
CS2 =\
m.contourf(lons,lats,hgt,CS1.levels,cmap=plt.cm.jet,extend='both',latlon=True)
m.colorbar(CS2) # draw colorbar
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawmapboundary()
m.fillcontinents()
# draw parallels and meridians.
parallels = np.arange(-60.,90,30.)
m.drawparallels(parallels,labels=[1,0,0,0])
meridians = np.arange(-360.,360.,30.)
m.drawmeridians(meridians)
plt.title('Mollweide Filled Contour Demo')
sys.stdout.write('plotting with mollweide basemap ...\n')

# create new figure
fig=plt.figure()
# set up Robinson map projection.
m = Basemap(resolution='c',projection='robin',lon_0=0)
# make a filled contour plot.
CS1 = m.contour(lons,lats,hgt,15,linewidths=0.5,colors='k',latlon=True)
CS2 = m.contourf(lons,lats,hgt,CS1.levels,cmap=plt.cm.jet,extend='both',latlon=True)
m.colorbar(CS2) # draw colorbar
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawmapboundary()
m.fillcontinents()
# draw parallels and meridians.
parallels = np.arange(-60.,90,30.)
m.drawparallels(parallels,labels=[1,0,0,0])
meridians = np.arange(-360.,360.,60.)
m.drawmeridians(meridians,labels=[0,0,0,1])
plt.title('Robinson Filled Contour Demo')
sys.stdout.write('plotting with robinson basemap ...\n')

# create new figure
fig=plt.figure()
# set up map projection (azimuthal equidistant).
m = Basemap(projection='npaeqd',lon_0=-90,boundinglat=15.,resolution='c')
# make a filled contour plot.
x, y = m(lons, lats)
CS1 = m.contour(x,y,hgt,15,linewidths=0.5,colors='k')
CS2 = m.contourf(x,y,hgt,CS2.levels,cmap=plt.cm.jet,extend='both')
m.colorbar(CS2,pad='12%') # draw colorbar
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawmapboundary()
m.fillcontinents()
# draw parallels and meridians.
parallels = np.arange(0.,80,20.)
m.drawparallels(parallels,labels=[0,0,1,1])
meridians = np.arange(10.,360.,20.)
m.drawmeridians(meridians,labels=[1,1,1,1])
plt.title('Azimuthal Equidistant Filled Contour Demo',y=1.075)
sys.stdout.write('plotting with azimuthal equidistant basemap ...\n')

# create new figure
fig=plt.figure()
# setup of orthographic basemap
m = Basemap(resolution='c',projection='ortho',\
            lat_0=45.,lon_0=-120.)
# make a filled contour plot.
x, y = m(lons, lats)
CS1 = m.contour(x,y,hgt,15,linewidths=0.5,colors='k')
CS2 = m.contourf(x,y,hgt,CS1.levels,cmap=plt.cm.jet,extend='both')
m.colorbar(CS2) # draw colorbar
# draw coastlines and political boundaries.
m.drawcoastlines()
m.fillcontinents()
m.drawmapboundary()
# draw parallels and meridians.
parallels = np.arange(-80.,90,20.)
m.drawparallels(parallels)
meridians = np.arange(-360.,360.,20.)
m.drawmeridians(meridians)
plt.title('Orthographic Filled Contour Demo')
sys.stdout.write('plotting with orthographic basemap ..\n')
plt.show()
