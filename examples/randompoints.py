from pylab import *
import math, random
from matplotlib.toolkits.basemap import Basemap
from matplotlib.numerix.random_array import uniform

# Plot a bunch of randomly distributed points on the earth.

# set up stereographic map centered on N. Pole.
m = Basemap(lon_0=-105,boundinglat=30.,
            resolution='l',area_thresh=10000.,projection='npstere')
# number of points to plot.
npts = 750
# generate random points on a sphere,
# so that every small area on the sphere is expected
# to have the same number of points.
# http://mathworld.wolfram.com/SpherePointPicking.html
try: # this works for numpy
    u = uniform(0.,1.,size=npts)
    v = uniform(0.,1.,size=npts)
except: # this works for Numeric/numarray
    u = uniform(0.,1.,shape=npts)
    v = uniform(0.,1.,shape=npts)
lons = 360.*u
lats = (180./math.pi)*arccos(2*v-1) - 90.
# transform lons and lats to map coordinates.
x,y = m(lons,lats)
# plot them as filled circles on the map.
# first, create figure with same aspect ratio as map.
fig=m.createfigure()
# background color will be used for 'wet' areas.
fig.add_axes([0.1,0.1,0.8,0.8],axisbg='aqua')
# use zorder=10 to make sure markers are drawn last.
# (otherwise they are covered up when continents are filled)
m.scatter(x,y,marker='o',c='k',s=25,zorder=10)
# draw coasts and fill continents.
m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color='coral')
# draw parallels and meridians.
delat = 20.
circles = arange(0.,90.,delat).tolist()+\
          arange(-delat,-90,-delat).tolist()
m.drawparallels(circles)
delon = 45.
meridians = arange(0,360,delon)
m.drawmeridians(meridians,labels=[1,1,1,1])
title('Random Points',y=1.075)
show()
