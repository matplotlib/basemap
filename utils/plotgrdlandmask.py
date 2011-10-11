from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import sys
from netCDF4 import Dataset

# 5 min lsmask, resolution = 'i': 
# grdlandmask -Ggrdlandmask5min_i.nc -I5m -R-180/180/-90/90 -Di -N0/1/2/1/2
# -A100+l
# 2.5 min lsmask, resolution = 'h': 
# grdlandmask -Ggrdlandmask2p5min_h.nc -I2.5m -R-180/180/-90/90 -Dh -N0/1/2/1/2
# -A10+l
filename = sys.argv[1]

nc = Dataset(filename)
lons = nc.variables['x'][:]
nlons = len(lons)
lats = nc.variables['y'][:]
nlats = len(lats)
lsmask = nc.variables['z'][:].astype(np.uint8)
resolution = 'i'

m =\
Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,resolution=resolution,projection='cyl')
m.drawcoastlines()
m.drawlsmask(land_color='coral',ocean_color='aqua',lsmask=lsmask,lsmask_lons=lons,lsmask_lats=lats,lakes=True)
plt.title('%s by %s land-sea mask from grdlandmask' % (nlons,nlats))
#lsmask.tofile('5minlsmask_gshhs_i.dat')
import gzip
f = gzip.open('2p5minlsmask_gshhs_i.gz','wb')
f.write(lsmask.tostring())
f.close()

plt.show()
