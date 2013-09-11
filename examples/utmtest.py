from __future__ import print_function
from mpl_toolkits.basemap import pyproj
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
# proj4 definition of UTM zone 17.
# coordinates agree to the 3rd decimal place (less than 1mm)
# when compared to result from
# http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html
p1 = pyproj.Proj(proj='utm',zone='17')
miamilon=-80.28; miamilat=25.92
x,y = p1(miamilon,miamilat)
print('utm coordinates from pyproj: ',x,y)
lon_0=-81; lat_0=0.
# mimic this with Basemap (WGS84 ellipsoid, with scaling factor of 0.9996 at
# central meridian).  latitude bounds are given for the MGRS band N, which covers Florida.
m = Basemap(projection='tmerc',lon_0=lon_0,lat_0=lat_0,
            k_0=0.9996,rsphere=(6378137.00,6356752.314245179),
            llcrnrlon=lon_0-3,llcrnrlat=24,urcrnrlon=lon_0+3.,urcrnrlat=32,resolution='i')
x0,y0 = m(lon_0,lat_0) 
x1,y1 = m(miamilon,miamilat)
# add coordinate offsets (Basemap coordinate system has x=0,y=0 at lower
# left corner of domain)
x = x1-x0+5.e5
y = y1-y0
print('utm coordinates from Basemap:',x,y)
m.drawcoastlines()
m.drawstates()
m.fillcontinents(color='coral',lake_color='aqua')
m.drawmapboundary(fill_color='aqua')
m.drawparallels(range(24,33),labels=[1,0,0,0])
m.drawmeridians(range(-84,-77),labels=[0,0,0,1])
plt.title('UTM zone 17N')
plt.show()
