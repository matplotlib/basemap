from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

# create Basemap instance for Geostationary (satellite view) projection.
lon_0 = float(raw_input('enter reference longitude (lon_0):'))

# map with land/sea mask plotted
fig=plt.figure()
m = Basemap(projection='geos',lon_0=lon_0,rsphere=(6378137.00,6356752.3142),resolution=None)
# plot land-sea mask.
# land red, oceans blue.
# lakes=True means plot inland lakes with ocean color.
m.drawlsmask(land_color='red',ocean_color='blue',lakes=True)
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
m.drawmapboundary()
plt.title('Geostationary Map Centered on Lon=%s' % (lon_0))

# map with continents drawn and filled.
fig = plt.figure()
m = Basemap(projection='geos',lon_0=lon_0,rsphere=(6378137.00,6356752.3142),resolution='l')
m.drawcoastlines()
m.drawmapboundary(fill_color='aqua')
m.fillcontinents(color='coral',lake_color='aqua')
m.drawcountries()
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
m.drawmapboundary()
plt.title('Geostationary Map Centered on Lon=%s' % (lon_0))
plt.show()
