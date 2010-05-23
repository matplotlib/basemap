from mpl_toolkits.basemap import Basemap, NetCDFFile
import matplotlib.pyplot as plt
import numpy as np
# read in orography of icosahedral global grid.
f = NetCDFFile('C02562.orog.nc')
lons = (180./np.pi)*f.variables['grid_center_lon'][:]
lats = (180./np.pi)*f.variables['grid_center_lat'][:]
z = f.variables['zs'][:]
map = Basemap(projection='ortho',lon_0=-105,lat_0=40)
x,y = map(lons, lats)
map.drawcoastlines()
map.drawmapboundary()
# tri=True forwards to axes.tripcolor
map.pcolor(x,y,z,tri=True,shading='faceted')
plt.show()
