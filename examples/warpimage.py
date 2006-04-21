import pylab as P
import Image as I
from matplotlib.toolkits.basemap import Basemap
# shows how to warp an image from one map projection to another.
# Uses PIL.
# Download image from
# http://www.space-graphics.com/earth_topo-bathy.htm, 
# convert from jpg to png.
# read in png image to rgba array of normalized floats.
try:
    rgba = P.imread('e_topo_bathy_4k.png')
except:
    msg  = """
please download image from http://www.space-graphics.com/earth_topo-bathy.htm,
convert e_topo_bathy_4k.jpg to e_topo_bathy_4k.png is present working directory."""
    raise IOError, msg
# reverse lats
rgba = rgba[::-1,:,:]
# define lat/lon grid that image spans (projection='cyl').
nlons = rgba.shape[1]; nlats = rgba.shape[0]
delta = 360./float(nlons)
lons = P.arange(-180.+0.5*delta,180.,delta)
lats = P.arange(-90.+0.5*delta,90.,delta)
# define Lambert Conformal basemap for North America.
# resolution=None means don't bother with boundary data.
m = Basemap(llcrnrlon=-145.5,llcrnrlat=1.,urcrnrlon=-2.566,urcrnrlat=46.352,\
            rsphere=(6378137.00,6356752.3142),lat_1=50.,lon_0=-107.,\
            resolution=None,projection='lcc')
# transform to nx x ny regularly spaced native projection grid
# nx and ny chosen to have roughly the same horizontal res as original image.
dx = 2.*P.pi*m.rmajor/float(nlons)
nx = int((m.xmax-m.xmin)/dx)+1; ny = int((m.ymax-m.ymin)/dx)+1
rgba_warped = P.zeros((ny,nx,4),'d')
# interpolate rgba values from proj='cyl' (geographic coords) to 'lcc'
for k in range(4):
    rgba_warped[:,:,k] = m.transform_scalar(rgba[:,:,k],lons,lats,nx,ny)
# convert normalized floats to integer RGBA values between 0 and 255
rgba_warped = (255.*rgba_warped).astype(P.UInt8)
# convert rgba values to pil image.
pilimage = I.fromstring('RGBA',(nx,ny),rgba_warped[::-1,:,:].tostring())
# plot pil image.
im = m.imshow(pilimage)
# draw parallels and meridians.
# label on left, right and bottom of map.
parallels = P.arange(0.,80,20.)
m.drawparallels(parallels,labels=[1,1,0,1])
meridians = P.arange(10.,360.,30.)
m.drawmeridians(meridians,labels=[1,1,0,1])
P.title("Global earth topo-bathy image warped from 'cyl' to 'lcc' projection",fontsize=12)
P.show()
