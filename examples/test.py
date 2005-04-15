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

m = Basemap(-180.,-90,180.,90.,\
            resolution='c',area_thresh=10000.,projection='cyl')
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.75,0.75])
im = imshow(topoin,cm.jet,extent=(m.llcrnrx,m.urcrnrx,m.llcrnry,m.urcrnry),origin='lower')
cax = axes([0.875, 0.1, 0.05, 0.75])
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines(ax)
#m.drawcountries(ax)
#m.drawstates(ax)
#m.fillcontinents(ax)
# draw parallels
delat = 30.
circles = arange(0.,90.+delat,delat).tolist()+\
          arange(-delat,-90.-delat,-delat).tolist()
m.drawparallels(ax,circles,labels=[1,0,0,1])
# draw meridians
delon = 60.
meridians = arange(-180,180,delon)
m.drawmeridians(ax,meridians,labels=[1,0,0,1])
ax.set_xticks([]) # no ticks
ax.set_yticks([])
title('Cylindrical Equidistant')
print 'plotting Cylindrical Equidistant example, close plot window to proceed ...'
show()

m = Basemap(-180.,-80.,180.,80.,\
            resolution='c',area_thresh=10000.,projection='merc',\
            lon_0=0.5*(lons[0]+lons[-1]),lat_ts=20.)
# transform to nx x ny regularly spaced native projection grid
nx = len(lons); ny = int(80.*len(lats)/90.)
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
fig.add_axes([0.1,0.1,0.75,0.75])
im = imshow(topodat,cm.jet,extent=(m.llcrnrx,m.urcrnrx,m.llcrnry,m.urcrnry),origin='lower')
ax = gca() # get current axis instance
m.drawcoastlines(ax)
m.drawcountries(ax)
m.drawstates(ax)
m.fillcontinents(ax)
# draw parallels
m.drawparallels(ax,circles,labels=[1,1,1,1])
# draw meridians
m.drawmeridians(ax,meridians,labels=[1,1,1,1])
ax.set_xticks([]) # no ticks
ax.set_yticks([])
title('Mercator',y=1.1)
print 'plotting Mercator example, close plot window to proceed ...'
show()

m = Basemap(-145.5,1.,-2.566,46.352,\
            resolution='c',area_thresh=10000.,projection='lcc',\
            lat_1=50.,lon_0=-107.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.7,0.7])
im = imshow(topodat,cm.jet,extent=(m.llcrnrx,m.urcrnrx,m.llcrnry,m.urcrnry),origin='lower')
cax = axes([0.875, 0.1, 0.05, 0.7])
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines(ax)
m.drawcountries(ax)
m.drawstates(ax)
#m.fillcontinents(ax)
# draw parallels
delat = 20.
circles = arange(0.,90.+delat,delat).tolist()+\
          arange(-delat,-90.-delat,-delat).tolist()
m.drawparallels(ax,circles,labels=[1,1,0,1])
# draw meridians
delon = 30.
meridians = arange(10.,360.,delon)
m.drawmeridians(ax,meridians,labels=[1,1,0,1])
ax.set_xticks([]) # no ticks
ax.set_yticks([])
title('Lambert Conformal Conic')
print 'plotting Lambert Conformal example, close plot window to proceed ...'
show()

m = Basemap(-10.,20.,55.,75.,
            resolution='l',projection='aea',\
            lat_1=40.,lat_2=60,lon_0=35.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.7,0.7])
im = imshow(topodat,cm.jet,extent=(m.llcrnrx,m.urcrnrx,m.llcrnry,m.urcrnry),origin='lower')
im.set_clim(-4000.,3000.)
cax = axes([0.875, 0.1, 0.05, 0.7])
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines(ax)
m.drawcountries(ax)
# draw parallels
delat = 20.
circles = arange(0.,90.+delat,delat).tolist()+\
          arange(-delat,-90.-delat,-delat).tolist()
m.drawparallels(ax,circles,labels=[1,1,1,1])
# draw meridians
delon = 30.
meridians = arange(10.,360.,delon)
m.drawmeridians(ax,meridians,labels=[1,1,1,1])
ax.set_xticks([]) # no ticks
ax.set_yticks([])
title('Albers Equal Area Conic',y=1.075)
print 'plotting Albers Equal Area example, close plot window to proceed ...'
show()

# north polar projection.
#m = Basemap(-150.,-20.826,30.,-20.826,
#            resolution='c',area_thresh=10000.,projection='stere',\
#            lat_0=90.,lon_0=-105.,lat_ts=90.)
# south polar projection.
m = Basemap(-150.,20.826,30.,20.826,
            resolution='c',area_thresh=10000.,projection='stere',\
            lat_0=-90.,lon_0=-105.,lat_ts=-90.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.7,0.7])
im = imshow(topodat,cm.jet,extent=(m.llcrnrx,m.urcrnrx,m.llcrnry,m.urcrnry),origin='lower')
cax = axes([0.875, 0.1, 0.05, 0.7])
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines(ax)
m.drawcountries(ax)
#m.fillcontinents(ax)
# draw parallels
m.drawparallels(ax,circles)
# draw meridians
m.drawmeridians(ax,meridians,labels=[1,1,1,1])
ax.set_xticks([]) # no ticks
ax.set_yticks([])
title('Polar Stereographic',y=1.075)
print 'plotting Stereographic example, close plot window to proceed ...'
show()

# lambert azimuthal north polar projection.
m = Basemap(-150.,-20.826,30.,-20.826,
            resolution='c',area_thresh=10000.,projection='laea',\
            lat_0=90.,lon_0=-105.,lat_ts=90.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.7,0.7])
im = imshow(topodat,cm.jet,extent=(m.llcrnrx,m.urcrnrx,m.llcrnry,m.urcrnry),origin='lower')
cax = axes([0.875, 0.1, 0.05, 0.7])
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines(ax)
m.drawcountries(ax)
m.drawstates(ax)
#m.fillcontinents(ax)
# draw parallels
m.drawparallels(ax,circles)
# draw meridians
m.drawmeridians(ax,meridians,labels=[1,1,1,1])
ax.set_xticks([]) # no ticks
ax.set_yticks([])
title('Lambert Azimuthal Equal Area',y=1.075)
print 'plotting Lambert Azimuthal example, close plot window to proceed ...'
show()
print 'done'
