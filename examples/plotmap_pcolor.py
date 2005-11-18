# make plot of etopo bathymetry/topography data on
# orthographic map projection, drawing coastlines,
# and parallels/meridians.

# this version does not interpolate the data to the projection grid.
# instead, the x,y coordinates of the lat/lon grid are passed to pcolor.
from matplotlib import rcParams, use
rcParams['numerix'] = 'Numeric'  # make sure Numeric is used (to read pickle)

from matplotlib.toolkits.basemap import Basemap, shiftgrid
from pylab import *
import cPickle
from matplotlib.numerix import ma

# read in topo data from pickle (on a regular lat/lon grid)
# longitudes go from 20 to 380.
topodict = cPickle.load(open('etopo20.pickle','rb'))
topodatin = topodict['data']; lonsin = topodict['lons']; latsin = topodict['lats']
# shift lons and lats by 1/2 grid increment (so values represent the vertices
# of the grid box surrounding the data value, as pcolor expects).
delon = lonsin[1]-lonsin[0]
delat = latsin[1]-latsin[0]
lons = zeros(len(lonsin)+1,'f')
lats = zeros(len(latsin)+1,'f')
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
topodat = where(topodatin < 0.,1.e10,topodatin)
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
title('ETOPO Topography - Lambert Conformal Conic')
show()

# setup of basemap ('ortho' = orthographic projection)
m = Basemap(projection='ortho',
            resolution='c',area_thresh=10000.,lat_0=30,lon_0=-60)
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7],frameon=False)
# make a pcolor plot.
x, y = m(lons, lats)
p = m.pcolor(x,y,topodatin,shading='flat')
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
title('ETOPO Topography - Orthographic')
savefig('ortho')
show()
