from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
# make land-sea-lake mask from built-in coastline dataset using is_land method.
# this is very slow!
minutes = 150
resolution = 'c'

delon = minutes/60.
nlons = int(360./delon); nlons = nlons + 1
nlats = int(180./delon); nlats = nlats + 1
print minutes,nlons,nlats,resolution
lons = np.linspace(-180,180,nlons)
lats = np.linspace(-90,90,nlats)
lsmask = np.zeros((len(lats),len(lons)),dtype=np.uint8)
print lsmask.shape, lsmask.dtype

m =\
Basemap(llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,resolution=resolution,projection='cyl')
for j,lat in enumerate(lats):
    #print j
    for i,lon in enumerate(lons):
        lsmask[j,i] = m.is_land(lon,lat)
m.drawlsmask(land_color='coral',ocean_color='aqua',lsmask=lsmask,lsmask_lons=lons,lsmask_lats=lats,lakes=True)
plt.title('%s minute degree land-sea mask' % minutes)
#lsmask.tofile('%sminlsmask_gshhs_%s.dat' % (minutes, resolution))

plt.show()
