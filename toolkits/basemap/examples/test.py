# make plots of etopo bathymetry/topography data on
# various map projections, drawing coastlines, state and
# country boundaries, filling continents and drawing
# parallels/meridians

from matplotlib.toolkits.basemap import Basemap, interp
from pylab import *
import cPickle

def shiftgrid(lon0,datain,lonsin,start=True):
    """ 
 shift global lat/lon grid east or west.
 assumes wraparound (or cyclic point) is included.

 lon0:  starting longitude for shifted grid 
        (ending longitude if start=False). lon0 must be on
        input grid (with the range of lonsin).
 datain:  original data.
 lonsin:  original longitudes.
 start[True]: if True, lon0 represents he starting longitude 
 of the new grid. if False, lon0 is the ending longitude.

 returns dataout,lonsout (data and longitudes on shifted grid).
    """
    if fabs(lonsin[-1]-lonsin[0]-360.) > 1.e-4:
        raise ValueError, 'cyclic point not included'
    if lon0 < lonsin[0] or lon0 > lonsin[-1]:
        raise ValueError, 'lon0 outside of range of lonsin'
    i0 = argsort(fabs(lonsin-lon0))[0]
    dataout = zeros(datain.shape,datain.typecode())
    lonsout = zeros(lonsin.shape,lonsin.typecode())
    if start:
        lonsout[0:len(lonsin)-i0] = lonsin[i0:]
    else:
        lonsout[0:len(lonsin)-i0] = lonsin[i0:]-360.
    dataout[:,0:len(lonsin)-i0] = datain[:,i0:]
    if start:
        lonsout[len(lonsin)-i0:] = lonsin[1:i0+1]+360.
    else:
        lonsout[len(lonsin)-i0:] = lonsin[1:i0+1]
    dataout[:,len(lonsin)-i0:] = datain[:,1:i0+1]
    return dataout,lonsout

def addcyclic(arrin,lonsin):
   """
 add cyclic (wraparound) point in longitude.
   """
   nlats = arrin.shape[0]
   nlons = arrin.shape[1]
   arrout  = zeros((nlats,nlons+1),arrin.typecode())
   arrout[:,0:nlons] = arrin[:,:]
   arrout[:,nlons] = arrin[:,0]
   lonsout = zeros(nlons+1,lonsin.typecode())
   lonsout[0:nlons] = lonsin[:]
   lonsout[nlons]  = lonsin[-1] + lonsin[1]-lonsin[0]
   return arrout,lonsout

# read in topo data from pickle (on a regular lat/lon grid)
# longitudes go from 20 to 380.
topodict = cPickle.load(open('etopo20.pickle','rb'))
topoin = topodict['data']; lons = topodict['lons']; lats = topodict['lats']

# shift data so lons go from -180 to 180 instead of 20 to 380.
topoin,lons = shiftgrid(180.,topoin,lons,start=False)

print 'min/max etopo20 data:'
print min(ravel(topoin)),max(ravel(topoin))

m = Basemap(lons[0],lats[0],lons[-1],lats[-1],\
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
m.drawparallels(ax,circles)
# draw meridians
delon = 60.
lon1 = int(lons[0]/delon)*delon
lon2 = (int(lons[-1]/delon)+1)*delon
meridians = arange(lon1,lon2,delon)
m.drawmeridians(ax,meridians)
ax.set_xticks([]) # no ticks
ax.set_yticks([])
title('Cylindrical Equidistant')
print 'plotting Cylindrical Equidistant example, close plot window to proceed ...'
show()

m = Basemap(lons[0],-80.,lons[-1],80.,\
            resolution='c',area_thresh=10000.,projection='merc',\
            lon_0=0.5*(lons[0]+lons[-1]),lat_ts=20.)
# define grid (nx x ny regularly spaced native projection grid)
nx = len(lons); ny = int(80.*len(lats)/90.)
lonsout, latsout = m.makegrid(nx,ny)
topodat = interp(topoin,lons,lats,lonsout,latsout)
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
fig.add_axes([0.1,0.1,0.8,0.8])
im = imshow(topodat,cm.jet,extent=(m.llcrnrx,m.urcrnrx,m.llcrnry,m.urcrnry),origin='lower')
ax = gca() # get current axis instance
m.drawcoastlines(ax)
m.drawcountries(ax)
m.drawstates(ax)
m.fillcontinents(ax)
# draw parallels
m.drawparallels(ax,circles)
# draw meridians
m.drawmeridians(ax,meridians)
ax.set_xticks([]) # no ticks
ax.set_yticks([])
title('Mercator')
print 'plotting Mercator example, close plot window to proceed ...'
show()

m = Basemap(-145.5,1.,-2.566,46.352,\
            resolution='c',area_thresh=10000.,projection='lcc',\
            lat_1=50.,lon_0=-107.)
# define grid (nx x ny regularly spaced native projection grid)
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
lonsout, latsout = m.makegrid(nx,ny)
# adjust longitudes to be consistent with etopo data
# (which goes from 20 to 380, instead of 0 to 360).
topodat = interp(topoin,lons,lats,lonsout,latsout)
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.75,0.75])
im = imshow(topodat,cm.jet,extent=(m.llcrnrx,m.urcrnrx,m.llcrnry,m.urcrnry),origin='lower')
cax = axes([0.875, 0.1, 0.05, 0.75])
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
m.drawparallels(ax,circles)
# draw meridians
delon = 30.
meridians = arange(0.,360.,delon)
m.drawmeridians(ax,meridians)
ax.set_xticks([]) # no ticks
ax.set_yticks([])
title('Lambert Conformal Conic')
print 'plotting Lambert Conformal example, close plot window to proceed ...'
show()

m = Basemap(-10.,20.,55.,75.,
            resolution='l',projection='aea',\
            lat_1=40.,lat_2=60,lon_0=35.)
# define grid (nx x ny regularly spaced native projection grid)
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
lonsout, latsout = m.makegrid(nx,ny)
# adjust longitudes to be consistent with etopo data
# (which goes from 20 to 380, instead of 0 to 360).
lonsout = where(lonsout < lons[0], lonsout+360., lonsout)
topodat = interp(topoin,lons,lats,lonsout,latsout)
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.75,0.75])
im = imshow(topodat,cm.jet,extent=(m.llcrnrx,m.urcrnrx,m.llcrnry,m.urcrnry),origin='lower')
im.set_clim(-4000.,3000.)
cax = axes([0.875, 0.1, 0.05, 0.75])
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines(ax)
m.drawcountries(ax)
# draw parallels
delat = 20.
circles = arange(0.,90.+delat,delat).tolist()+\
          arange(-delat,-90.-delat,-delat).tolist()
m.drawparallels(ax,circles)
# draw meridians
delon = 30.
meridians = arange(0.,360.,delon)
m.drawmeridians(ax,meridians)
ax.set_xticks([]) # no ticks
ax.set_yticks([])
title('Albers Equal Area Conic')
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
# define grid (nx x ny regularly spaced native projection grid)
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
lonsout, latsout = m.makegrid(nx,ny)
# adjust longitudes to be consistent with etopo data
lonsout = where(lonsout < lons[0], lonsout+360., lonsout)
topodat = interp(topoin,lons,lats,lonsout,latsout)
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.75,0.75])
im = imshow(topodat,cm.jet,extent=(m.llcrnrx,m.urcrnrx,m.llcrnry,m.urcrnry),origin='lower')
cax = axes([0.875, 0.1, 0.05, 0.75])
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines(ax)
m.drawcountries(ax)
#m.fillcontinents(ax)
# draw parallels
m.drawparallels(ax,circles)
# draw meridians
m.drawmeridians(ax,meridians)
ax.set_xticks([]) # no ticks
ax.set_yticks([])
title('Polar Stereographic')
print 'plotting Stereographic example, close plot window to proceed ...'
show()

# lambert azimuthal north polar projection.
m = Basemap(-150.,-20.826,30.,-20.826,
            resolution='c',area_thresh=10000.,projection='laea',\
            lat_0=90.,lon_0=-105.,lat_ts=90.)
# define grid (nx x ny regularly spaced native projection grid)
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
lonsout, latsout = m.makegrid(nx,ny)
# adjust longitudes to be consistent with etopo data
lonsout = where(lonsout < lons[0], lonsout+360., lonsout)
topodat = interp(topoin,lons,lats,lonsout,latsout)
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.75,0.75])
im = imshow(topodat,cm.jet,extent=(m.llcrnrx,m.urcrnrx,m.llcrnry,m.urcrnry),origin='lower')
cax = axes([0.875, 0.1, 0.05, 0.75])
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines(ax)
m.drawcountries(ax)
m.drawstates(ax)
#m.fillcontinents(ax)
# draw parallels
m.drawparallels(ax,circles)
# draw meridians
m.drawmeridians(ax,meridians)
ax.set_xticks([]) # no ticks
ax.set_yticks([])
title('Lambert Azimuthal Equal Area')
print 'plotting Lambert Azimuthal example, close plot window to proceed ...'
show()
print 'done'
