# make plot of etopo bathymetry/topography data on
# lambert conformal conic map projection, drawing coastlines, state and
# country boundaries, and parallels/meridians.

# the data is interpolated to the native projection grid.

##################################
# pylab-free version of plotmap.py
##################################
# set backend to Agg.
import matplotlib
matplotlib.use('Agg')

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.toolkits.basemap import Basemap, shiftgrid
from matplotlib.figure import Figure
import matplotlib.numerix as NX
import matplotlib.cm as cm

def load(fname,comments='%',delimiter=None):
    """
    Load ASCII data from fname into an array and return the array (from pylab).
    """
    if fname.endswith('.gz'):
        import gzip
        fh = gzip.open(fname)
    else:
        fh = file(fname)
    X = []
    numCols = None
    for line in fh:
        line = line[:line.find(comments)].strip()
        if not len(line): continue
        row = [float(val) for val in line.split(delimiter)]
        thisLen = len(row)
        if numCols is not None and thisLen != numCols:
            raise ValueError('All rows must have the same number of columns')
        X.append(row)

    X = NX.array(X)
    r,c = X.shape
    if r==1 or c==1:
        X.shape = max([r,c]),
    return X

# read in topo data (on a regular lat/lon grid)
# longitudes go from 20 to 380.
topoin = NX.array(load('etopo20data.gz'),'d')
lons = NX.array(load('etopo20lons.gz'),'d')
lats = NX.array(load('etopo20lats.gz'),'d')
# shift data so lons go from -180 to 180 instead of 20 to 380.
topoin,lons = shiftgrid(180.,topoin,lons,start=False)

# setup of basemap ('lcc' = lambert conformal conic).
# use major and minor sphere radii from WGS84 ellipsoid.
m = Basemap(llcrnrlon=-145.5,llcrnrlat=1.,urcrnrlon=-2.566,urcrnrlat=46.352,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',area_thresh=1000.,projection='lcc',\
            lat_1=50.,lon_0=-107.)
# create figure.
fig = Figure()
canvas = FigureCanvas(fig)
# set 'ax' instance variable so pylab won't be imported.
m.ax = fig.add_axes([0.1,0.1,0.7,0.7])
# reset figure size to have same aspect ratio as map.
# fig will be 8 inches wide.
# (don't use createfigure, since that imports pylab).
if m.aspect <= 1.:
    fig.set_figsize_inches((8,m.aspect*8.))
else:
    fig.set_figsize_inches((8/m.aspect,8.))
# transform to nx x ny regularly spaced native projection grid
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
topodat,x,y = m.transform_scalar(topoin,lons,lats,nx,ny,returnxy=True)
# plot image over map with imshow.
im = m.imshow(topodat,cm.jet)
cax = fig.add_axes([0.875, 0.1, 0.05, 0.7],frameon=False) # setup colorbar axes
fig.colorbar(im,tickfmt='%d', cax=cax) # draw colorbar
# plot blue dot on boulder, colorado and label it as such.
xpt,ypt = m(-104.237,40.125) 
m.plot([xpt],[ypt],'bo') 
m.ax.text(xpt+100000,ypt+100000,'Boulder')
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawcountries()
m.drawstates()
# draw parallels and meridians.
# label on left, right and bottom of map.
parallels = NX.arange(0.,80,20.)
m.drawparallels(parallels,labels=[1,1,0,1])
meridians = NX.arange(10.,360.,30.)
m.drawmeridians(meridians,labels=[1,1,0,1])
# set title.
m.ax.set_title('ETOPO Topography - Lambert Conformal Conic')
# save image (width 800 pixels with dpi=100 and fig width 8 inches).
canvas.print_figure('plotmap',dpi=100)
# done.
print 'image saved in plotmap.png'
