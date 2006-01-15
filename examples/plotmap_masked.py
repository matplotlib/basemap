# make plot of etopo bathymetry/topography data on
# lambert conformal conic map projection, drawing coastlines, state and
# country boundaries, and parallels/meridians.

# demonstrates use of masked arrays to mask out certain regions
# (in this case the oceans)


from matplotlib import rcParams, use
rcParams['numerix'] = 'Numeric'  # make sure Numeric is used (to read pickle)
from matplotlib.toolkits.basemap import Basemap, shiftgrid
from pylab import *
from matplotlib.numerix import ma
import matplotlib.colors as colors

# read in topo data from pickle (on a regular lat/lon grid)
# longitudes go from 20 to 380.
topoin = load('etopo20data.gz')
lonsin = load('etopo20lons.gz')
latsin = load('etopo20lats.gz')
# shift data so lons go from -180 to 180 instead of 20 to 380.
topoin,lonsin = shiftgrid(180.,topoin,lonsin,start=False)

# imshow version.

# setup of basemap ('lcc' = lambert conformal conic).
# use major and minor sphere radii from WGS84 ellipsoid.
m = Basemap(llcrnrlon=-145.5,llcrnrlat=1.,urcrnrlon=-2.566,urcrnrlat=46.352,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=50.,lon_0=-107.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat,x,y = m.transform_scalar(topoin,lonsin,latsin,nx,ny,returnxy=True)
# set up figure with same aspect ratio as map.
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7],axisbg='aqua')
# make topodat a masked array, masking values lower than sea level.
topodat = where(topodat < 0.,1.e10,topodat)
topodat = ma.masked_values(topodat, 1.e10)
# plot image over map with imshow.
im = m.imshow(topodat,cm.jet,norm=colors.normalize(vmin=-4000.0,vmax=3000.0,clip=False))
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
# plot blue dot on boulder, colorado and label it as such.
xpt,ypt = m(-104.237,40.125) 
m.plot([xpt],[ypt],'bo') 
text(xpt+100000,ypt+100000,'Boulder')
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawcountries()
m.drawstates()
# draw parallels and meridians.
# label on left, right and bottom of map.
parallels = arange(0.,80,20.)
m.drawparallels(parallels,labels=[1,1,0,1])
meridians = arange(10.,360.,30.)
m.drawmeridians(meridians,labels=[1,1,0,1])
# set title.
title('Masked ETOPO Topography - via imshow')
show()

# pcolor version.

# shift lons and lats by 1/2 grid increment (so values represent the vertices
# of the grid box surrounding the data value, as pcolor expects).
delon = lonsin[1]-lonsin[0]
delat = latsin[1]-latsin[0]
lons = zeros(len(lonsin)+1,'d')
lats = zeros(len(latsin)+1,'d')
lons[0:len(lonsin)] = lonsin-0.5*delon
lons[-1] = lonsin[-1]+0.5*delon
lats[0:len(latsin)] = latsin-0.5*delon
lats[-1] = latsin[-1]+0.5*delon
lons, lats = meshgrid(lons, lats)

# setup of basemap ('lcc' = lambert conformal conic).
# use major and minor sphere radii from WGS84 ellipsoid.
m = Basemap(llcrnrlon=-145.5,llcrnrlat=1.,urcrnrlon=-2.566,urcrnrlat=46.352,\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=50.,lon_0=-107.)
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7],axisbg='aqua')
# make topodat a masked array, masking values lower than sea level.
topodat = where(topoin < 0.,1.e10,topoin)
topodat = ma.masked_values(topodat, 1.e10)
# make a pcolor plot.
x, y = m(lons, lats)
p = m.pcolor(x,y,topodat,shading='flat',cmap=cm.jet)
clim(-4000.,3000.) # adjust colormap.
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
# plot blue dot on boulder, colorado and label it as such.
xpt,ypt = m(-104.237,40.125) 
m.plot([xpt],[ypt],'bo') 
text(xpt+100000,ypt+100000,'Boulder')
# draw coastlines and political boundaries.
m.drawcoastlines()
#m.drawcountries()
#m.drawstates()
# draw parallels and meridians.
# label on left, right and bottom of map.
parallels = arange(0.,80,20.)
m.drawparallels(parallels,labels=[1,1,0,1])
meridians = arange(10.,360.,30.)
m.drawmeridians(meridians,labels=[1,1,0,1])
title('Masked ETOPO Topography - via pcolor')
show()
