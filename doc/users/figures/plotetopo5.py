from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset

# read in etopo5 topography/bathymetry.
etopodata =\
Dataset('http://ferret.pmel.noaa.gov/thredds/dodsC/data/PMEL/etopo5.nc')
topoin = etopodata.variables['ROSE'][:]
lons = etopodata.variables['ETOPO05_X'][:]
lats = etopodata.variables['ETOPO05_Y'][:]
# shift data so lons go from -180 to 180 instead of 20 to 380.
topoin,lons = shiftgrid(180.,topoin,lons,start=False)

# plot topography/bathymetry as an image.

# create the figure and axes instances.
fig = plt.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
# setup of basemap ('lcc' = lambert conformal conic).
# use major and minor sphere radii from WGS84 ellipsoid.
m = Basemap(llcrnrlon=-145.5,llcrnrlat=1.,urcrnrlon=-2.566,urcrnrlat=46.352,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=50.,lon_0=-107.,ax=ax)
# transform to nx x ny regularly spaced 5km native projection grid
nx = int((m.xmax-m.xmin)/5000.)+1; ny = int((m.ymax-m.ymin)/5000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
# plot image over map with imshow.
im = m.imshow(topodat,cm.GMT_haxby)
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawcountries()
m.drawstates()
# draw parallels and meridians.
# label on left and bottom of map.
parallels = np.arange(0.,80,20.)
m.drawparallels(parallels,labels=[1,0,0,1])
meridians = np.arange(10.,360.,30.)
m.drawmeridians(meridians,labels=[1,0,0,1])
# use axes_grid toolkit to make colorbar axes.
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
cb = fig.colorbar(im,cax=cax)
ax.set_title('ETOPO5 Topography - Lambert Conformal Conic')
fig.savefig('etopo5.png')

# make a shaded relief plot.

# create new figure, axes instance.
fig = plt.figure()
ax = fig.add_axes([0.1,0.1,0.8,0.8])
# attach new axes image to existing Basemap instance.
m.ax = ax
# create light source object.
from matplotlib.colors import LightSource
ls = LightSource(azdeg = 90, altdeg = 20)
# convert data to rgb array including shading from light source.
# (must specify color map)
rgb = ls.shade(topodat, cm.GMT_haxby)
im = m.imshow(rgb)
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawcountries()
m.drawstates()
ax.set_title('Shaded ETOPO5 Topography - Lambert Conformal Conic')
fig.savefig('etopo5_shaded.png')
