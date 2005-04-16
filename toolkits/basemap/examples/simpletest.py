from matplotlib.toolkits.basemap import Basemap
import cPickle
from pylab import *
# read in topo data from pickle (on a regular lat/lon grid)
topodict = cPickle.load(open('etopo20.pickle','rb'))
etopo = topodict['data']; lons = topodict['lons']; lats = topodict['lats']
m = Basemap(lons[0],lats[0],lons[-1],lats[-1])
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
fig.add_axes([0.1,0.1,0.8,0.8])
im = m.imshow(etopo)
# draw coastlines and fill continents.
m.drawcoastlines()
m.fillcontinents()
# draw parallels, label on bottom.
circles = arange(-90.,120.,30.)
m.drawparallels(circles,labels=[1,0,0,0])
# draw meridians, label on left.
meridians = arange(0.,390.,60.)
m.drawmeridians(meridians,labels=[0,0,0,1])
title('Cylindrical Equidistant')
show()

