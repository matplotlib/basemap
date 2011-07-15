from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
# make land-sea-lake mask from built-in coastline dataset using is_land method.
# this is very slow!
lons = np.arange(-180,181,2.5)
lats = np.arange(-90,91,2.5)
lsmask = np.zeros((len(lats),len(lons)),dtype=np.uint8)
m = Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,resolution='c',projection='cyl')
for j,lat in enumerate(lats):
    for i,lon in enumerate(lons):
        lsmask[j,i] = m.is_land(lon,lat,lsmask=True)
m.drawlsmask(land_color='coral',ocean_color='aqua',lsmask=lsmask,lsmask_lons=lons,lsmask_lats=lats,lakes=True)
plt.title('2.5 degree land-sea mask')
plt.show()
