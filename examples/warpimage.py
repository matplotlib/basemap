import pylab as P
from matplotlib.toolkits.basemap import Basemap
from matplotlib.numerix import ma
from matplotlib.image import pil_to_array
from PIL import Image

# shows how to warp an image from one map projection to another.
# image from http://visibleearth.nasa.gov/

# read in jpeg image to rgba array of normalized floats.
pilImage = Image.open('land_shallow_topo_2048.jpg')
rgba = pil_to_array(pilImage)
rgba = rgba.astype(P.Float32)/255. # convert to normalized floats.

# define lat/lon grid that image spans (projection='cyl').
nlons = rgba.shape[1]; nlats = rgba.shape[0]
delta = 360./float(nlons)
lons = P.arange(-180.+0.5*delta,180.,delta)
lats = P.arange(-90.+0.5*delta,90.,delta)

# create new figure
fig=P.figure()
# define cylindrical equidistant projection.
m = Basemap(projection='cyl',llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,resolution='l')
# plot (unwarped) rgba image.
im = m.imshow(rgba)
# draw coastlines.
m.drawcoastlines(linewidth=0.5,color='0.5')
# draw lat/lon grid lines.
m.drawmeridians(P.arange(-180,180,60),labels=[0,0,0,1],color='0.5')
m.drawparallels(P.arange(-90,90,30),labels=[1,0,0,0],color='0.5')
P.title("Blue Marble image - native 'cyl' projection",fontsize=12)
print 'plot cylindrical map (no warping needed) ...'

# create new figure
fig=P.figure()
# define orthographic projection centered on North America.
m = Basemap(projection='ortho',lat_0=40,lon_0=40,resolution='l')
# transform to nx x ny regularly spaced native projection grid
# nx and ny chosen to have roughly the same horizontal res as original image.
dx = 2.*P.pi*m.rmajor/float(nlons)
nx = int((m.xmax-m.xmin)/dx)+1; ny = int((m.ymax-m.ymin)/dx)+1
rgba_warped = ma.zeros((ny,nx,4),P.Float64)
# interpolate rgba values from proj='cyl' (geographic coords) to 'lcc'
# values outside of projection limb will be masked.
for k in range(4):
    rgba_warped[:,:,k] = m.transform_scalar(rgba[:,:,k],lons,lats,nx,ny,masked=True)
# make points outside projection limb transparent.
rgba_warped = rgba_warped.filled(0.)
# plot warped rgba image.
im = m.imshow(rgba_warped)
# draw coastlines.
m.drawcoastlines(linewidth=0.5,color='0.5')
# draw lat/lon grid lines every 30 degrees.
m.drawmeridians(P.arange(0,360,30),color='0.5')
m.drawparallels(P.arange(-90,90,30),color='0.5')
P.title("Blue Marble image warped from 'cyl' to 'ortho' projection",fontsize=12)
print 'warp to orthographic map ...'

# create new figure
fig=P.figure()
# define Lambert Conformal basemap for North America.
m = Basemap(llcrnrlon=-145.5,llcrnrlat=1.,urcrnrlon=-2.566,urcrnrlat=46.352,\
            rsphere=(6378137.00,6356752.3142),lat_1=50.,lon_0=-107.,\
            resolution='i',area_thresh=1000.,projection='lcc')
# transform to nx x ny regularly spaced native projection grid
# nx and ny chosen to have roughly the same horizontal res as original image.
dx = 2.*P.pi*m.rmajor/float(nlons)
nx = int((m.xmax-m.xmin)/dx)+1; ny = int((m.ymax-m.ymin)/dx)+1
rgba_warped = P.zeros((ny,nx,4),P.Float64)
# interpolate rgba values from proj='cyl' (geographic coords) to 'lcc'
for k in range(4):
    rgba_warped[:,:,k] = m.transform_scalar(rgba[:,:,k],lons,lats,nx,ny)
# plot warped rgba image.
im = m.imshow(rgba_warped)
# draw coastlines.
m.drawcoastlines(linewidth=0.5,color='0.5')
# draw parallels and meridians.
# label on left, right and bottom of map.
parallels = P.arange(0.,80,20.)
m.drawparallels(parallels,labels=[1,1,0,1],color='0.5')
meridians = P.arange(10.,360.,30.)
m.drawmeridians(meridians,labels=[1,1,0,1],color='0.5')
P.title("Blue Marble image warped from 'cyl' to 'lcc' projection",fontsize=12)
print 'warp to lambert conformal map ...'
P.show()
