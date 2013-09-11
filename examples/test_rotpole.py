from __future__ import print_function
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

# o_lon_p, o_lat_p: true lat/lon of pole in rotated coordinate system

# mapping to CF metadata convention:
# grid_north_pole_longitude = normalize180(180 + lon_0), where normalize180
#  is a function that maps to interval [-180,180).
# grid_north_pole_latitude = o_lat_p
# north_pole_grid_longitude = o_lon_p (optional, assumed zero if not present)

def normalize180(lon):
    """Normalize lon to range [180, 180)"""
    lower = -180.; upper = 180.
    if lon > upper or lon == lower:
        lon = lower + abs(lon + upper) % (abs(lower) + abs(upper))
    if lon < lower or lon == upper:
        lon = upper - abs(lon - lower) % (abs(lower) + abs(upper))
    return lower if lon == upper else lon

lon_0 = normalize180(rotpole.grid_north_pole_longitude-180.)
o_lon_p = rotpole.north_pole_grid_longitude
o_lat_p = rotpole.grid_north_pole_latitude
print( rotpole )
print( 'lon_0,o_lon_p,o_lat_p=',lon_0,o_lon_p,o_lat_p)

m= Basemap(projection='rotpole',lon_0=lon_0,o_lon_p=o_lon_p,o_lat_p=o_lat_p,\
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
m= Basemap(projection='rotpole',lon_0=lon_0,o_lon_p=o_lon_p,o_lat_p=o_lat_p,\
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
