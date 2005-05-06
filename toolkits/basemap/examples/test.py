# make plots of etopo bathymetry/topography data on
# various map projections, drawing coastlines, state and
# country boundaries, filling continents and drawing
# parallels/meridians

from matplotlib.toolkits.basemap import Basemap, shiftgrid
from pylab import *
import cPickle

# read in topo data from pickle (on a regular lat/lon grid)
# longitudes go from 20 to 380.
topodict = cPickle.load(open('etopo20.pickle','rb'))
topoin = topodict['data']; lons = topodict['lons']; lats = topodict['lats']

# shift data so lons go from -180 to 180 instead of 20 to 380.
topoin,lons = shiftgrid(180.,topoin,lons,start=False)

print 'min/max etopo20 data:'
print min(ravel(topoin)),max(ravel(topoin))

# setup cylindrical equidistant map projection (global domain).
m = Basemap(-180.,-90,180.,90.,\
            resolution='c',area_thresh=10000.,projection='cyl')
# setup figure with same aspect ratio as map.
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
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
m = Basemap(-180.,-90.,180.,90.,\
            resolution='c',area_thresh=10000.,projection='mill')
# transform to nx x ny regularly spaced native projection grid
nx = len(lons); ny = len(lats)
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# setup figure with same aspect ratio as map.
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
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
m = Basemap(-180.,-80.,180.,80.,\
            resolution='c',area_thresh=10000.,projection='merc',\
            lon_0=0.5*(lons[0]+lons[-1]),lat_ts=20.)
# transform to nx x ny regularly spaced native projection grid
nx = len(lons); ny = int(80.*len(lats)/90.)
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# setup figure with same aspect ratio as map.
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
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

# setup transverse mercator basemap.
m = Basemap(170.,-45,10.,45.,\
            resolution='c',area_thresh=10000.,projection='tmerc',\
            lat_0=0.,lon_0=-90.)
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize/m.aspect,xsize))
fig.add_axes([0.125,0.2,0.6,0.6])
# transform to nx x ny regularly spaced native projection grid
ny = len(lats)
nx = ny/m.aspect
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

# setup equidistant conic
m = Basemap(-90,18.,-70,26.,\
            resolution='l',area_thresh=1000.,projection='eqdc',\
            lat_1=21.,lat_2=23.,lon_0=-80.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# setup figure with same aspect ratio as map.
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
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
m = Basemap(-145.5,1.,-2.566,46.352,\
            resolution='c',area_thresh=10000.,projection='lcc',\
            lat_1=50.,lon_0=-107.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# setup figure with same aspect ratio as map.
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
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
m = Basemap(-10.,20.,55.,75.,
            resolution='l',projection='aea',\
            lat_1=40.,lat_2=60,lon_0=35.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# setup figure with same aspect ratio as map.
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
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
m = Basemap(-150.,20.826,30.,20.826,
            resolution='c',area_thresh=10000.,projection='stere',\
            lat_0=-90.,lon_0=-105.,lat_ts=-90.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# setup figure with same aspect ratio as map.
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
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
m = Basemap(-150.,-20.826,30.,-20.826,
            resolution='c',area_thresh=10000.,projection='laea',\
            lat_0=90.,lon_0=-105.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# setup figure with same aspect ratio as map.
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
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
m = Basemap(-150.,-20.826,30.,-20.826,
            resolution='c',area_thresh=10000.,projection='aeqd',\
            lat_0=90.,lon_0=-105.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# setup figure with same aspect ratio as map.
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
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
print 'done'
