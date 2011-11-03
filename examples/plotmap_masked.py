# make plot of etopo bathymetry/topography data on
# lambert conformal conic map projection, drawing coastlines, state and
# country boundaries, and parallels/meridians.

# demonstrates use of masked arrays to mask out certain regions
# (in this case the oceans)

from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# read in topo data (on a regular lat/lon grid)
# longitudes go from 20 to 380.
topoin = np.loadtxt('etopo20data.gz')
lonsin = np.loadtxt('etopo20lons.gz')
latsin = np.loadtxt('etopo20lats.gz')
# shift data so lons go from -180 to 180 instead of 20 to 380.
topoin,lonsin = shiftgrid(180.,topoin,lonsin,start=False)

# setup of basemap ('lcc' = lambert conformal conic).
# use major and minor sphere radii from WGS84 ellipsoid.
m = Basemap(llcrnrlon=-145.5,llcrnrlat=1.,urcrnrlon=-2.566,urcrnrlat=46.352,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=50.,lon_0=-107.)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat,x,y = m.transform_scalar(topoin,lonsin,latsin,nx,ny,returnxy=True)
# create the figure.
fig=plt.figure(figsize=(8,8))
# add an axes.
ax = fig.add_axes([0.1,0.1,0.8,0.8])
# associate this axes with the Basemap instance.
m.ax = ax
# make topodat a masked array, masking values lower than sea level.
topodat = np.where(topodat < 0.,1.e10,topodat)
topodatm = ma.masked_values(topodat, 1.e10)
palette = plt.cm.YlOrRd
palette.set_bad('aqua', 1.0)
# plot image over map with imshow.
im = m.imshow(topodatm,palette,norm=colors.normalize(vmin=0.0,vmax=3000.0,clip=False))
m.colorbar(im,pad='12%') # draw colorbar
# plot blue dot on boulder, colorado and label it as such.
xpt,ypt = m(-104.237,40.125) 
m.plot([xpt],[ypt],'bo') 
ax.text(xpt+100000,ypt+100000,'Boulder')
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawcountries()
m.drawstates()
# draw parallels and meridians.
# label on left, right and bottom of map.
parallels = np.arange(0.,80,20.)
m.drawparallels(parallels,labels=[1,1,0,1])
meridians = np.arange(10.,360.,30.)
m.drawmeridians(meridians,labels=[1,1,0,1])
# set title.
ax.set_title('Masked ETOPO Topography - via imshow')
plt.show()
