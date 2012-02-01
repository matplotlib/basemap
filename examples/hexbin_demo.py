# example showing how to plot scattered data with hexbin.
from numpy.random import uniform
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
# create north polar stereographic basemap
m = Basemap(lon_0=270, boundinglat=20, projection='npstere',round=True)
# number of points to plot.
npts = 20000
# generate random points on a sphere,
# so that every small area on the sphere is expected
# to have the same number of points.
# http://mathworld.wolfram.com/SpherePointPicking.html
u = uniform(0.,1.,size=npts)
v = uniform(0.,1.,size=npts)
lons1 = 360.*u
lats1 = (180./np.pi)*np.arccos(2*v-1) - 90.
# toss points outside of map region.
lats = np.compress(lats1 > 20, lats1)
lons = np.compress(lats1 > 20, lons1)
# convert to map projection coordinates.
x, y = m(lons, lats)
# function to plot at those points.
xscaled = 4.*(x-0.5*(m.xmax-m.xmin))/m.xmax
yscaled = 4.*(y-0.5*(m.ymax-m.ymin))/m.ymax
z = xscaled*np.exp(-xscaled**2-yscaled**2)
# make plot
CS = plt.hexbin(x,y,C=z,gridsize=50,cmap=plt.cm.jet)
# draw coastlines, lat/lon lines.
m.drawcoastlines()
m.drawparallels(np.arange(0,81,20))
m.drawmeridians(np.arange(-180,181,60))
m.colorbar() # draw colorbar
plt.title('hexbin demo')
plt.show()
