from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import sys

def get_input(prompt):
    if sys.hexversion > 0x03000000:
        return input(prompt)
    else:
        return raw_input(prompt)

# create Basemap instance for Orthographic (satellite view) projection.
lon_0 = float(get_input('enter reference longitude (lon_0):'))
lat_0 = float(get_input('enter reference latitude (lat_0):'))

# map with land/sea mask plotted
fig = plt.figure()
resolution = 'l'; grid = 5
m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution=resolution)
# land coral, oceans aqua.
# lakes=True means plot inland lakes with ocean color.
# resolution = 5 (default) means use 5 min dataset (can use 2.5)
m.drawcoastlines()
m.drawlsmask(land_color='coral',ocean_color='aqua', lakes=True,\
        resolution=resolution,grid=grid)
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
m.drawmapboundary()
plt.title('Orthographic Map Centered on Lon=%s, Lat=%s' % (lon_0,lat_0))

# map with continents drawn and filled (continent filling fails for
# lon=-120,lat=60).
fig = plt.figure()
m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,resolution=resolution)
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
m.drawcountries()
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
m.drawmapboundary(fill_color='aqua')
# add a map scale.
length = 5000 
x1,y1 = 0.3*m.xmax, 0.25*m.ymax
lon1,lat1 = m(x1,y1,inverse=True)
m.drawmapscale(lon1,lat1,lon1,lat1,length,fontsize=8,barstyle='fancy',\
               labelstyle='fancy',units='km')
plt.title('Orthographic Map Centered on Lon=%s, Lat=%s' % (lon_0,lat_0))
plt.show()
