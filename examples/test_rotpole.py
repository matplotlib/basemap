from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

nc = Dataset('wm201_Arctic_JJA_1990-2008_moyenneDesMoyennes.nc')
lats = nc.variables['lat'][:]
lons = nc.variables['lon'][:]
rlats = nc.variables['rlat'][:]
rlons = nc.variables['rlon'][:]
rlons, rlats = np.meshgrid(rlons, rlats)
data = nc.variables['air'][0,0,:,:].squeeze()
data = np.ma.masked_values(data,-999.)
rotpole = nc.variables['rotated_pole']

m = Basemap(projection='npstere',lon_0=10,boundinglat=30,resolution='c')
x,y = m(lons,lats)
m.drawcoastlines()
m.contourf(x,y,data,20)
m.drawmeridians(np.arange(-180,180,20))
m.drawparallels(np.arange(20,80,20))
m.colorbar()
plt.title('rotated pole data in polar stere map')

plt.figure()
m = Basemap(projection='rotpole',lon_0=rotpole.grid_north_pole_longitude-180.,\
            o_lon_p=rotpole.north_pole_grid_longitude,\
            o_lat_p=rotpole.grid_north_pole_latitude,\
            llcrnrlat = lats[0,0], urcrnrlat = lats[-1,-1],\
            llcrnrlon = lons[0,0], urcrnrlon = lons[-1,-1],resolution='c')
x,y = m(lons,lats)
m.drawcoastlines()
m.contourf(x,y,data,20)
m.drawmeridians(np.arange(-180,180,20))
m.drawparallels(np.arange(20,80,20))
m.colorbar()
plt.title('rotated pole data in native map using real sphere corner lat/lons' )

plt.figure()
m = Basemap(projection='rotpole',lon_0=rotpole.grid_north_pole_longitude-180.,\
            o_lon_p=rotpole.north_pole_grid_longitude,\
            o_lat_p=rotpole.grid_north_pole_latitude,\
            llcrnry = rlats[0,0], urcrnry = rlats[-1,-1],\
            llcrnrx = rlons[0,0], urcrnrx = rlons[-1,-1],resolution='c')
x,y = m(lons,lats)
m.drawcoastlines()
m.contourf(x,y,data,20)
m.drawmeridians(np.arange(-180,180,20))
m.drawparallels(np.arange(20,80,20))
m.colorbar()
plt.title('rotated pole data in native map using rotated sphere corner lat/lons' )

plt.show()
