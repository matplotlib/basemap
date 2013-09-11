from __future__ import print_function
######################################
# pyplot-free version of simpletest.py
######################################

import numpy as np

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.cm as cm

from mpl_toolkits.basemap import Basemap

# read in topo data (on a regular lat/lon grid)
# longitudes go from 20 to 380.
etopo = np.loadtxt('etopo20data.gz')
lons = np.loadtxt('etopo20lons.gz')
lats = np.loadtxt('etopo20lats.gz')

# create figure.
fig = Figure()
canvas = FigureCanvas(fig)

# create axes instance
ax = fig.add_axes([0.1,0.1,0.8,0.8])

# create Basemap instance for Robinson projection.
# set 'ax' keyword so pylab won't be imported.
m = Basemap(projection='robin',lon_0=0.5*(lons[0]+lons[-1]),ax=ax)

# make filled contour plot.
x, y = m(*np.meshgrid(lons, lats))
cs = m.contourf(x,y,etopo,30,cmap=cm.jet)
# draw coastlines.
m.drawcoastlines()
# draw a line around the map region.
m.drawmapboundary()
# draw parallels and meridians.
m.drawparallels(np.arange(-60.,90.,30.),labels=[1,0,0,0],fontsize=10)
m.drawmeridians(np.arange(0.,420.,60.),labels=[0,0,0,1],fontsize=10)
# add a title.
ax.set_title('Robinson Projection')

# add a colorbar (must specify mappable and fig keywords, of pyplot will be
# invoked) .
cb=m.colorbar(mappable=cs,location='right',size='5%',pad='2%',ticks=cs.levels[::3],fig=fig)

# save image (width 800 pixels with dpi=100 and fig width 8 inches).
canvas.print_figure('simpletest',dpi=100)
# done.
print('image saved in simpletest.png')

