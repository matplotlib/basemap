from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
# read in topo data (on a regular lat/lon grid)
etopo=np.loadtxt('etopo20data.gz')
lons=np.loadtxt('etopo20lons.gz')
lats=np.loadtxt('etopo20lats.gz')
# create Basemap instance for Robinson projection.
m = Basemap(projection='robin',lon_0=0.5*(lons[0]+lons[-1]))
# make filled contour plot.
x, y = m(*np.meshgrid(lons, lats))
cs = m.contourf(x,y,etopo,30,cmap=plt.cm.jet)
# draw coastlines.
m.drawcoastlines()
# draw a line around the map region.
m.drawmapboundary()
# draw parallels and meridians.
m.drawparallels(np.arange(-60.,90.,30.),labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,420.,60.),labels=[0,0,0,1])
# add a title.
plt.title('Robinson Projection')
plt.show()
