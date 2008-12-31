from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

# the shortest route from the center of the map
# to any other point is a straight line in the azimuthal
# equidistant projection. Such lines show the true scale
# on the earth's surface.
# So, for the specified point, this script draws a map that shows
# in which direction to depart for other points on earth and how far
# it will be to reach that destination.
# The specified point shows up as a red dot in the center of the map.

# user enters the lon/lat of the point, and it's name
lon_0 = float(raw_input('input reference lon (degrees):'))
lat_0 = float(raw_input('input reference lat (degrees):'))
location = raw_input('name of location:')

m = Basemap(resolution='c',projection='aeqd',lat_0=lat_0,lon_0=lon_0)
# fill background.
m.drawmapboundary(fill_color='aqua')
# draw coasts and fill continents.
m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color='coral',lake_color='aqua')
# 20 degree graticule.
m.drawparallels(np.arange(-80,81,20))
m.drawmeridians(np.arange(-180,180,20))
# draw a black dot at the center.
xpt, ypt = m(lon_0, lat_0)
m.plot([xpt],[ypt],'ko') 
# draw the title.
plt.title('The World According to Garp in '+location)
plt.show()
