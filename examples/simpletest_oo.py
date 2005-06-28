#####################################
# pylab-free version of simpletest.py
#####################################
# set backend to Agg.
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.toolkits.basemap import Basemap
from matplotlib.figure import Figure
from matplotlib.mlab import meshgrid
import matplotlib.numerix as nx
import matplotlib.cm as cm
import cPickle
# read in topo data from pickle (on a regular lat/lon grid)
topodict = cPickle.load(open('etopo20.pickle','rb'))
etopo = topodict['data']; lons = topodict['lons']; lats = topodict['lats']
# create figure.
fig = Figure()
canvas = FigureCanvas(fig)
ax = fig.add_axes([0.125,0.175,0.75,0.75],frameon=False)
# create Basemap instance for Robinson projection.
# set 'ax' keyword so pylab won't be imported.
m = Basemap(projection='robin',lon_0=0.5*(lons[0]+lons[-1]),ax=ax)
# reset figure size to have same aspect ratio as map.
fig.set_figsize_inches((8,m.aspect*8.))
# make filled contour plot.
x, y = m(*meshgrid(lons, lats))
levels, colls = m.contourf(x,y,etopo,30,cmap=cm.jet,colors=None)
# draw coastlines.
m.drawcoastlines()
# draw a line around the map region.
m.drawmapboundary()
# draw parallels and meridians.
m.drawparallels(nx.arange(-60.,90.,30.),labels=[1,0,0,0],fontsize=10)
m.drawmeridians(nx.arange(0.,420.,60.),labels=[0,0,0,1],fontsize=10)
# add a title.
ax.set_title('Robinson Projection')
# add a colorbar.
cax = fig.add_axes([0.25, 0.05, 0.5, 0.05],frameon=False)
fig.colorbar(colls.mappable, cax=cax, tickfmt='%d', orientation='horizontal') 
# save image.
canvas.print_figure('simpletest')
# done.
print 'image saved in simpletest.png'
