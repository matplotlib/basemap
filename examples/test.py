# make plots of etopo bathymetry/topography data on
# various map projections, drawing coastlines, state and
# country boundaries, filling continents and drawing
# parallels/meridians

from matplotlib.toolkits.basemap import Basemap, shiftgrid
from pylab import *
import cPickle

# read in topo data from pickle (on a regular lat/lon grid)
# longitudes go from 20 to 380.
topodatin = load('etopo20data.gz')
lonsin = load('etopo20lons.gz'); latsin = load('etopo20lats.gz')

# shift data so lons go from -180 to 180 instead of 20 to 380.
topoin,lons = shiftgrid(180.,topodatin,lonsin,start=False)
lats = latsin

print 'min/max etopo20 data:'
print min(ravel(topoin)),max(ravel(topoin))

# setup cylindrical equidistant map projection (global domain).
m = Basemap(llcrnrlon=-180.,llcrnrlat=-90,urcrnrlon=180.,urcrnrlat=90.,\
            resolution='c',area_thresh=10000.,projection='cyl')
# setup figure with same aspect ratio as map.
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.75,0.75])
# plot image over map.
im = m.imshow(topoin,cm.jet)
cax = axes([0.875, 0.1, 0.05, 0.75]) # setup colorbar axes.
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines()
#m.drawcountries()
#m.drawstates()
#m.fillcontinents()
# draw parallels
delat = 30.
circles = arange(0.,90.+delat,delat).tolist()+\
          arange(-delat,-90.-delat,-delat).tolist()
m.drawparallels(circles,labels=[1,0,0,1])
# draw meridians
delon = 60.
meridians = arange(-180,180,delon)
m.drawmeridians(meridians,labels=[1,0,0,1])
title('Cylindrical Equidistant')
print 'plotting Cylindrical Equidistant example, close plot window to proceed ...'
show()

# setup miller cylindrical map projection.
m = Basemap(llcrnrlon=-180.,llcrnrlat=-90,urcrnrlon=180.,urcrnrlat=90.,\
            resolution='c',area_thresh=10000.,projection='mill')
# transform to nx x ny regularly spaced native projection grid
nx = len(lons); ny = len(lats)
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# setup figure with same aspect ratio as map.
fig=m.createfigure()
fig.add_axes([0.1,0.1,0.75,0.75])
# plot image over map.
im = m.imshow(topodat,cm.jet)
m.drawcoastlines()
# draw parallels
m.drawparallels(circles,labels=[1,1,1,1])
# draw meridians
m.drawmeridians(meridians,labels=[1,1,1,1])
title('Miller Cylindrical',y=1.1)
print 'plotting Miller Cylindrical example, close plot window to proceed ...'
show()

# setup mercator map projection (-80 to +80).
m = Basemap(llcrnrlon=-180.,llcrnrlat=-80,urcrnrlon=180.,urcrnrlat=80.,\
            resolution='c',area_thresh=10000.,projection='merc',\
            lon_0=0.5*(lons[0]+lons[-1]),lat_ts=20.)
# transform to nx x ny regularly spaced native projection grid
nx = len(lons); ny = int(80.*len(lats)/90.)
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# setup figure with same aspect ratio as map.
fig=m.createfigure()
fig.add_axes([0.1,0.1,0.75,0.75])
# plot image over map.
im = m.imshow(topodat,cm.jet)
m.drawcoastlines()
m.drawcountries()
m.drawstates()
m.fillcontinents()
# draw parallels
m.drawparallels(circles,labels=[1,1,1,1])
# draw meridians
m.drawmeridians(meridians,labels=[1,1,1,1])
title('Mercator',y=1.1)
print 'plotting Mercator example, close plot window to proceed ...'
show()

# setup cassini-soldner basemap.
m = Basemap(llcrnrlon=-6,llcrnrlat=49,urcrnrlon=4,urcrnrlat=59,\
            resolution='l',area_thresh=1000.,projection='cass',\
            lat_0=54.,lon_0=-2.)
fig=m.createfigure()
fig.add_axes([0.125,0.2,0.6,0.6])
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/20000.)+1; ny = int((m.ymax-m.ymin)/20000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# plot image over map.
im = m.imshow(topodat,cm.jet)
# get current axis instance.
ax = gca()
cax = axes([0.8, 0.2, 0.075, 0.6]) # setup colorbar axes.
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines()
# draw parallels
delat = 2.
circles = arange(40.,70.,delat)
m.drawparallels(circles,labels=[1,0,0,1],fontsize=10)
# draw meridians
delon = 2.
meridians = arange(-10,10,delon)
m.drawmeridians(meridians,labels=[1,0,0,1],fontsize=10)
title('Cassini-Soldner Projection')
print 'plotting Cassini-Soldner example, close plot window to proceed ...'
show()

# setup gnomonic basemap.
m = Basemap(llcrnrlon=-95.,llcrnrlat=-52,urcrnrlon=-35.,urcrnrlat=15.,\
            resolution='c',area_thresh=10000.,projection='gnom',\
            lat_0=-10.,lon_0=-60.)
fig=m.createfigure()
fig.add_axes([0.125,0.2,0.6,0.6])
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# plot image over map.
im = m.imshow(topodat,cm.jet)
# get current axis instance.
ax = gca()
cax = axes([0.8, 0.2, 0.075, 0.6]) # setup colorbar axes.
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines()
m.drawcountries()
# draw parallels
delat = 20.
circles = arange(-80.,100.,delat)
m.drawparallels(circles,labels=[1,0,0,1],fontsize=10)
# draw meridians
delon = 20.
meridians = arange(-180,180,delon)
m.drawmeridians(meridians,labels=[1,0,0,1],fontsize=10)
title('Gnomonic Projection')
print 'plotting Gnomonic example, close plot window to proceed ...'
show()

# setup transverse mercator basemap.
m = Basemap(llcrnrlon=170.,llcrnrlat=-45,urcrnrlon=10.,urcrnrlat=45.,\
            resolution='c',area_thresh=10000.,projection='tmerc',\
            lat_0=0.,lon_0=-90.)
fig=m.createfigure()
fig.add_axes([0.125,0.2,0.6,0.6])
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# plot image over map.
im = m.imshow(topodat,cm.jet)
# get current axis instance.
ax = gca()
cax = axes([0.8, 0.2, 0.075, 0.6]) # setup colorbar axes.
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines()
# draw parallels
delat = 20.
circles = arange(-80.,100.,delat)
m.drawparallels(circles,labels=[1,0,0,0],fontsize=10)
# draw meridians
delon = 20.
meridians = arange(-180,180,delon)
m.drawmeridians(meridians,labels=[1,0,0,0],fontsize=10)
title('Transverse Mercator Projection')
print 'plotting Transverse Mercator example, close plot window to proceed ...'
show()

# setup oblique mercator basemap.
m = Basemap(llcrnrlon=-130.,llcrnrlat=39,urcrnrlon=-124.,urcrnrlat=60.,\
            resolution='l',area_thresh=1000.,projection='omerc',\
            lon_2=-140,lat_2=55,lon_1=-120,lat_1=40)
fig=m.createfigure()
fig.add_axes([0.125,0.2,0.6,0.6])
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/20000.)+1; ny = int((m.ymax-m.ymin)/20000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# plot image over map.
im = m.imshow(topodat,cm.jet)
# get current axis instance.
ax = gca()
cax = axes([0.8, 0.2, 0.075, 0.6]) # setup colorbar axes.
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines()
m.drawcountries()
m.drawstates()
# draw parallels
delat = 3.
circles = arange(40,60,delat)
m.drawparallels(circles,labels=[1,0,0,0],fontsize=10)
# draw meridians
delon = 3.
meridians = arange(-140,-120,delon)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
title('Oblique Mercator Projection')
print 'plotting Oblique Mercator example, close plot window to proceed ...'
show()

# setup polyconic basemap.
m = Basemap(llcrnrlon=-35.,llcrnrlat=-30,urcrnrlon=80.,urcrnrlat=50.,\
            resolution='c',area_thresh=1000.,projection='poly',\
            lat_0=0.,lon_0=20.)
fig=m.createfigure()
fig.add_axes([0.125,0.2,0.6,0.6])
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# plot image over map.
im = m.imshow(topodat,cm.jet)
# get current axis instance.
ax = gca()
cax = axes([0.8, 0.2, 0.075, 0.6]) # setup colorbar axes.
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines()
m.drawcountries()
# draw parallels
delat = 20.
circles = arange(-80.,100.,delat)
m.drawparallels(circles,labels=[1,0,0,0],fontsize=10)
# draw meridians
delon = 20.
meridians = arange(-180,180,delon)
m.drawmeridians(meridians,labels=[1,0,0,1],fontsize=10)
title('Polyconic Projection')
print 'plotting Polyconic example, close plot window to proceed ...'
show()

# setup equidistant conic
m = Basemap(llcrnrlon=-90.,llcrnrlat=18,urcrnrlon=-70.,urcrnrlat=26.,\
            resolution='l',area_thresh=1000.,projection='eqdc',\
            lat_1=21.,lat_2=23.,lon_0=-80.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# setup figure with same aspect ratio as map.
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# plot image over map.
im = m.imshow(topodat,cm.jet)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes.
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines()
m.drawcountries()
m.drawstates()
m.fillcontinents(color='olive')
# draw parallels
delat = 2.
circles = arange(17,27,delat)
m.drawparallels(circles,labels=[1,0,0,0])
# draw meridians
delon = 5.
meridians = arange(-100,-60,delon)
m.drawmeridians(meridians,labels=[0,0,0,1])
title('Equidistant Conic')
print 'plotting Equidistant Conic example, close plot window to proceed ...'
show()

# setup lambert conformal map projection (North America).
m = Basemap(llcrnrlon=-145.5,llcrnrlat=1,urcrnrlon=-2.566,urcrnrlat=46.352,\
            resolution='c',area_thresh=10000.,projection='lcc',\
            lat_1=50.,lon_0=-107.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# setup figure with same aspect ratio as map.
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# plot image over map.
im = m.imshow(topodat,cm.jet)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes.
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines()
m.drawcountries()
m.drawstates()
#m.fillcontinents()
# draw parallels
delat = 20.
circles = arange(0.,90.+delat,delat).tolist()+\
          arange(-delat,-90.-delat,-delat).tolist()
m.drawparallels(circles,labels=[1,1,0,1])
# draw meridians
delon = 30.
meridians = arange(10.,360.,delon)
m.drawmeridians(meridians,labels=[1,1,0,1])
title('Lambert Conformal Conic')
print 'plotting Lambert Conformal example, close plot window to proceed ...'
show()

# setup albers equal area map projection (Europe).
m = Basemap(llcrnrlon=-10.,llcrnrlat=20,urcrnrlon=55.,urcrnrlat=75,\
            resolution='l',projection='aea',\
            lat_1=40.,lat_2=60,lon_0=35.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# setup figure with same aspect ratio as map.
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# plot image over map.
im = m.imshow(topodat,cm.jet)
im.set_clim(-4000.,3000.) # adjust range of colors.
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes.
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines()
m.drawcountries()
# draw parallels
delat = 20.
circles = arange(0.,90.+delat,delat).tolist()+\
          arange(-delat,-90.-delat,-delat).tolist()
m.drawparallels(circles,labels=[1,1,1,1])
# draw meridians
delon = 30.
meridians = arange(10.,360.,delon)
m.drawmeridians(meridians,labels=[1,1,1,1])
title('Albers Equal Area Conic',y=1.075)
print 'plotting Albers Equal Area example, close plot window to proceed ...'
show()

# setup stereographic map projection (Southern Hemisphere).
m = Basemap(llcrnrlon=-150.,llcrnrlat=20.826,urcrnrlon=30.,urcrnrlat=20.826,\
            resolution='c',area_thresh=10000.,projection='stere',\
            lat_0=-90.,lon_0=-105.,lat_ts=-90.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# setup figure with same aspect ratio as map.
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# plot image over map.
im = m.imshow(topodat,cm.jet)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes.
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines()
m.drawcountries()
#m.fillcontinents()
# draw parallels
m.drawparallels(circles)
# draw meridians
m.drawmeridians(meridians,labels=[1,1,1,1])
title('Polar Stereographic',y=1.075)
print 'plotting Stereographic example, close plot window to proceed ...'
show()

# setup lambert azimuthal map projection (Northern Hemisphere).
m = Basemap(llcrnrlon=-150.,llcrnrlat=-20.826,urcrnrlon=30.,urcrnrlat=-20.826,\
            resolution='c',area_thresh=10000.,projection='laea',\
            lat_0=90.,lon_0=-105.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# setup figure with same aspect ratio as map.
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# plot image over map.
im = m.imshow(topodat,cm.jet)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes.
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines()
m.drawcountries()
m.drawstates()
#m.fillcontinents()
# draw parallels
m.drawparallels(circles)
# draw meridians
m.drawmeridians(meridians,labels=[1,1,1,1])
title('Lambert Azimuthal Equal Area',y=1.075)
print 'plotting Lambert Azimuthal example, close plot window to proceed ...'
show()

# setup azimuthal equidistant map projection (Northern Hemisphere).
m = Basemap(llcrnrlon=-150.,llcrnrlat=40.,urcrnrlon=30.,urcrnrlat=40.,\
            resolution='c',area_thresh=10000.,projection='aeqd',\
            lat_0=90.,lon_0=-105.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# setup figure with same aspect ratio as map.
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# plot image over map.
im = m.imshow(topodat,cm.jet)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes.
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines()
m.drawcountries()
m.drawstates()
#m.fillcontinents()
# draw parallels
m.drawparallels(circles)
# draw meridians
m.drawmeridians(meridians,labels=[1,1,1,1])
title('Azimuthal Equidistant',y=1.075)
print 'plotting Azimuthal Equidistant example, close plot window to proceed ...'
show()

# projections with elliptical boundaries (orthographic, mollweide and robinson)

# can't use imshow, since it expects a rectangular image.
# instead use pcolor or contourf, and plot data 
# directly on lat/lon grid without interpolation.

# shift lons and lats by 1/2 grid increment (so values represent the vertices
# of the grid box surrounding the data value, as pcolor expects).
delon = lonsin[1]-lonsin[0]
delat = latsin[1]-latsin[0]
lons = zeros(len(lonsin)+1,'d')
lats = zeros(len(latsin)+1,'d')
lons[0:len(lonsin)] = lonsin-0.5*delon
lons[-1] = lonsin[-1]+0.5*delon
lats[0:len(latsin)] = latsin-0.5*delat
lats[-1] = latsin[-1]+0.5*delat

# setup of basemap ('ortho' = orthographic projection)
m = Basemap(projection='ortho',
            resolution='c',area_thresh=10000.,lat_0=30,lon_0=-60)
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# pcolor plot (slow)
#x,y = m(*meshgrid(lons,lats))
#p = m.pcolor(x,y,topodatin,shading='flat')
# filled contours (faster)
x,y = m(*meshgrid(lonsin,latsin))
cs = m.contourf(x,y,topodatin,20,cmap=cm.jet)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
# draw coastlines and political boundaries.
m.drawcoastlines()
# draw parallels and meridians (labelling is 
# not implemented for orthographic).
parallels = arange(-80.,90,20.)
m.drawparallels(parallels)
meridians = arange(0.,360.,20.)
m.drawmeridians(meridians)
# draw boundary around map region.
m.drawmapboundary()
title('Orthographic')
print 'plotting Orthographic example, close plot window to proceed ...'
show()

# setup of basemap ('moll' = mollweide projection)
m = Basemap(projection='moll',
            resolution='c',area_thresh=10000.,lon_0=0.5*(lonsin[0]+lonsin[-1]))
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# pcolor plot (slow)
#x,y = m(*meshgrid(lons,lats))
#p = m.pcolor(x,y,topodatin,shading='flat')
# filled contours (faster)
x,y = m(*meshgrid(lonsin,latsin))
cs = m.contourf(x,y,topodatin,20,cmap=cm.jet)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
# draw coastlines and political boundaries.
m.drawcoastlines()
# draw parallels and meridians
parallels = arange(-60.,90,30.)
m.drawparallels(parallels,labels=[1,0,0,0])
meridians = arange(0.,360.,30.)
m.drawmeridians(meridians)
# draw boundary around map region.
m.drawmapboundary()
title('Mollweide')
print 'plotting Mollweide example, close plot window to proceed ...'
show()

# setup of basemap ('robin' = robinson projection)
m = Basemap(projection='robin',
            resolution='c',area_thresh=10000.,lon_0=0.5*(lonsin[0]+lonsin[-1]))
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# pcolor plot (slow)
#x,y = m(*meshgrid(lons,lats))
#p = m.pcolor(x,y,topodatin,shading='flat')
# filled contours (faster)
x,y = m(*meshgrid(lonsin,latsin))
cs = m.contourf(x,y,topodatin,20,cmap=cm.jet)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
# draw coastlines and political boundaries.
m.drawcoastlines()
# draw parallels and meridians
parallels = arange(-60.,90,30.)
m.drawparallels(parallels,labels=[1,0,0,0])
meridians = arange(0.,360.,60.)
m.drawmeridians(meridians,labels=[0,0,0,1])
# draw boundary around map region.
m.drawmapboundary()
title('Robinson')
print 'plotting Robinson example, close plot window to proceed ...'
show()

print 'done'
