from matplotlib.toolkits.basemap import Basemap
import pyproj
from pylab import title, show, arange, pi

# the shortest route from the center of the map
# to any other point is a straight line in the azimuthal
# equidistant projection. Such lines show the true scale
# on the earth's surface.
# So, for the specified point, this script draws a map that shows
# in which direction to depart for other points on earth and how far
# it will be to reach that destination.
# The specified point shows up as a red dot in the center of the map.

# This example shows how to use the pyproj module included in basemap
# to choose the lon/lat values for the corners of the map.

# user enters the lon/lat of the point, and it's name
lon_0 = float(raw_input('input reference lon (degrees):'))
lat_0 = float(raw_input('input reference lat (degrees):'))
location = raw_input('name of location:')

# use the pyproj module to find 4 corners of the map.
rsphere = 6370997
p = pyproj.Proj(proj='aeqd',lon_0=lon_0,lat_0=lat_0,R=rsphere)
# upper right corner is 14000 km from center point.
xurcrnr,yurcrnr = 14000000,14000000
# find the lon/lat for the lower left and upper right corners.
llcrnrlon,llcrnrlat = p(-xurcrnr,-yurcrnr,inverse=True)
urcrnrlon,urcrnrlat = p(xurcrnr,yurcrnr,inverse=True)
# use these values to setup Basemap instance.
m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,\
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            resolution='c',projection='aeqd',\
            lat_0=lat_0,lon_0=lon_0,rsphere=rsphere)
# draw coasts and fill continents.
m.drawcoastlines(linewidth=0.5)
m.fillcontinents()
# 20 degree graticule.
m.drawparallels(arange(-80,81,20))
m.drawmeridians(arange(-180,180,20))
# draw a red dot at the center.
xpt, ypt = m(lon_0, lat_0)
m.plot([xpt],[ypt],'ro') 
# draw the title.
title('The World According to Garp in '+location)
show()
