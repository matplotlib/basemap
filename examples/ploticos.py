from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from netCDF4 import Dataset as NetCDFFile
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
#z = ma.masked_where(z < 1.e-5, z) # for testing masked arrays.
map.pcolor(x,y,z,tri=True,shading='flat',edgecolors='k',cmap=plt.cm.hot_r,vmin=0,vmax=3000)
#map.contourf(x,y,z,np.arange(0,3000,150),tri=True)
#map.contour(x,y,z,np.arange(0,3000,150),tri=True)
plt.title('pcolor plot on a global icosahedral mesh')
plt.show()
