import mpl_toolkits.basemap as bm
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
# change default value of latlon kwarg to True.
bm.latlon_default=True
# read in topo data (on a regular lat/lon grid)
etopo=np.loadtxt('etopo20data.gz')
lons=np.loadtxt('etopo20lons.gz')
lats=np.loadtxt('etopo20lats.gz')
# mask land regions.
etopo = ma.masked_where(etopo > 0, etopo)
lons, lats = np.meshgrid(lons, lats)
# create Basemap instance.
m = bm.Basemap(projection='kav7',lon_0=0)
# latlon_default=True, so contourf expecting lons,lats (not x,y)
# data automatically shifted in longitude to fit map projection region.
cs = m.contourf(lons,lats,etopo,30,cmap=plt.cm.jet)
# draw coastlines.
m.drawcoastlines()
# draw parallels and meridians.
m.drawparallels(np.arange(-60.,90.,30.),labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.),labels=[0,0,0,1],fontsize=12)
plt.title('test latlon=True')
plt.show()
