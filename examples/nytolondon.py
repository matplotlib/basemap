# example demonstrating how to draw a great circle on a map.
import matplotlib
from matplotlib.toolkits.basemap import Basemap
from pylab import *

# setup a mercator projection.
m = Basemap(-90.,30.,30.,60.,\
            resolution='c',area_thresh=10000.,projection='merc',\
            lat_ts=20.)
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
fig.add_axes([0.1,0.1,0.8,0.8])
ax = gca() # get current axis instance
# compute 100 points along a great circle from NY to London
# (in map projection coordinates).
nylat = 40.78
nylon = -73.98
lonlat = 51.53
lonlon = 0.08
x,y = m.gcpoints(nylon,nylat,lonlon,lonlat,100)
# plot the great circle.
ax.plot(x,y,linewidth=2)
ax.update_datalim(((m.llcrnrx, m.llcrnry),(m.urcrnrx,m.urcrnry)))
ax.set_xlim((m.llcrnrx, m.urcrnrx))
ax.set_ylim((m.llcrnry, m.urcrnry))
m.drawcoastlines(ax)
m.fillcontinents(ax)
# draw parallels
circles = [35,45,55]
m.drawparallels(ax,circles,labels=[1,1,0,1])
# draw meridians
delon = 30.
meridians = arange(-180,180,delon)
m.drawmeridians(ax,meridians,labels=[1,1,0,1])
title('Great Circle from New York to London')
show()
