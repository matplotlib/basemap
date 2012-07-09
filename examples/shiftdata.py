from mpl_toolkits.basemap import Basemap,shiftgrid
import numpy as np
import matplotlib.pyplot as plt
# read in topo data (on a regular lat/lon grid)
# (skip last column of data)
etopo=np.loadtxt('etopo20data.gz')[:,:-1]
lons=np.loadtxt('etopo20lons.gz')[:-1]
lats=np.loadtxt('etopo20lats.gz')
lons, lats = np.meshgrid(lons, lats)
# create Basemap instance.
m = Basemap(projection='kav7',lon_0=0)
# shift data and longitudes to fit map region
etopo, lons = m.shiftdata(etopo, lons)
# transform lats/lons to map projection coords
x, y = m(lons, lats)
# make filled contour plot
cs = m.contourf(x,y,etopo,30,cmap=plt.cm.jet)
# draw coastlines.
m.drawcoastlines()
# draw parallels and meridians.
m.drawparallels(np.arange(-60.,90.,30.),labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.),labels=[0,0,0,1],fontsize=12)
plt.title('test shiftdata method')
plt.show()
