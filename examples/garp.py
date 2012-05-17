from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import sys

def get_input(prompt):
    if sys.hexversion > 0x03000000:
        return input(prompt)
    else:
        return raw_input(prompt)

# the shortest route from the center of the map
# to any other point is a straight line in the azimuthal
# equidistant projection. Such lines show the true scale
# on the earth's surface.
# So, for the specified point, this script draws a map that shows
# in which direction to depart for other points on earth and how far
# it will be to reach that destination.
# The specified point shows up as a red dot in the center of the map.

# user enters the lon/lat of the point, and it's name
lon_0 = float(get_input('input reference lon (degrees):'))
lat_0 = float(get_input('input reference lat (degrees):'))
location = get_input('name of location:')

# no width/height or lat/lon corners specified, so whole world
# is plotted in a circle.
m = Basemap(resolution='c',projection='aeqd',lat_0=lat_0,lon_0=lon_0)

# draw coastlines and fill continents.
# **it's easy to make this fail with global aeqd plots.
#  For example, if the center point is at the North Pole,
#  the continent filling routines get confused and fills
#  the outside of Antartica instead of the inside**

#m.drawmapboundary(fill_color='white')
#m.drawcoastlines(linewidth=0.5)
#m.fillcontinents(color='black',lake_color='white')
#m.drawparallels(np.arange(-80,81,20),color='0.7')
#m.drawmeridians(np.arange(-180,180,20),color='0.7')

# draw lsmask instead of drawing continents (slower, but more robust).

m.drawlsmask(land_color='black',ocean_color='white',lakes=True,resolution='l',grid=5)
m.drawparallels(np.arange(-80,81,20),color='0.7')
m.drawmeridians(np.arange(-180,180,20),color='0.7')
m.drawmapboundary()

# blue marble background (pretty, but slow).

#m.bluemarble(scale=0.5)
#m.drawparallels(np.arange(-80,81,20),color='0.5')
#m.drawmeridians(np.arange(-180,180,20),color='0.5')
#m.drawmapboundary(color='0.5')

# draw a red dot at the center.
xpt, ypt = m(lon_0, lat_0)
m.plot([xpt],[ypt],'ro') 

# draw the title.
plt.title('The World According to Garp in '+location)

plt.show()
