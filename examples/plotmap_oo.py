# make plot of etopo bathymetry/topography data on
# lambert conformal conic map projection, drawing coastlines, state and
# country boundaries, and parallels/meridians.

# the data is interpolated to the native projection grid.

##################################
# pyplot/pylab-free version of plotmap.py
##################################
# set backend to Agg.
import matplotlib, sys
matplotlib.use('Agg')

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.figure import Figure
import numpy as np
import matplotlib.cm as cm

# read in topo data (on a regular lat/lon grid)
# longitudes go from 20 to 380.
topoin = np.loadtxt('etopo20data.gz')
lons = np.loadtxt('etopo20lons.gz')
lats = np.loadtxt('etopo20lats.gz')
# shift data so lons go from -180 to 180 instead of 20 to 380.
topoin,lons = shiftgrid(180.,topoin,lons,start=False)

# create figure, axes instance.
fig = Figure()
canvas = FigureCanvas(fig)
ax = fig.add_axes([0.1,0.1,0.7,0.7])

# setup of basemap ('lcc' = lambert conformal conic).
# use major and minor sphere radii from WGS84 ellipsoid.
# pass axes instance to Basemap constructor so pylab won't
# be imported.
m = Basemap(llcrnrlon=-145.5,llcrnrlat=1.,urcrnrlon=-2.566,urcrnrlat=46.352,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=50.,lon_0=-107.,ax=ax)
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat,x,y = m.transform_scalar(topoin,lons,lats,nx,ny,returnxy=True)
# plot image over map with imshow.
im = m.imshow(topodat,cm.jet)
# setup colorbar axes instance.
pos = ax.get_position()
l, b, w, h = pos.bounds
cax = fig.add_axes([l+w+0.075, b, 0.05, h],frameon=False) # setup colorbar axes
fig.colorbar(im, cax=cax) # draw colorbar
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
ax.set_title('ETOPO Topography - Lambert Conformal Conic')
# save image (width 800 pixels with dpi=100 and fig width 8 inches).
canvas.print_figure('plotmap',dpi=100)
# done.
sys.stdout.write('image saved in plotmap.png\n')
