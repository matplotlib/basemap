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
etopo = ma.masked_where(etopo > 0, etopo)
lons, lats = np.meshgrid(lons, lats)
# create Basemap instance.
m = bm.Basemap(projection='kav7',lon_0=0)
# can either shift data manually...
# shift data and longitudes to fit map region
#lons, etopo = m.shiftdata(lons, etopo)
# transform lats/lons to map projection coords
#x, y = m(lons, lats)
# make filled contour plot
#cs = m.contourf(x,y,etopo,30,cmap=plt.cm.jet)
# or let contourf to it for you by passing lons/lats and
# setting latlon=True.
cs = m.contourf(lons,lats,etopo,30,cmap=plt.cm.jet)
#cs = m.pcolormesh(lons,lats,etopo,shading='flat')
#cs = m.pcolor(lons,lats,etopo,shading='flat')
# draw coastlines.
m.drawcoastlines()
# draw parallels and meridians.
m.drawparallels(np.arange(-60.,90.,30.),labels=[1,0,0,0])
m.drawmeridians(np.arange(0.,360.,60.),labels=[0,0,0,1],fontsize=12)
plt.title('test shiftdata method')
plt.show()
