from mpl_toolkits.basemap import Basemap, shiftgrid, maskoceans, interp
import numpy as np 
import matplotlib.pyplot as plt

# example showing how to mask out 'wet' areas on a contour or pcolor plot.

topodatin = np.loadtxt('etopo20data.gz')
lonsin = np.loadtxt('etopo20lons.gz')
latsin = np.loadtxt('etopo20lats.gz')

# shift data so lons go from -180 to 180 instead of 20 to 380.
topoin,lons1 = shiftgrid(180.,topodatin,lonsin,start=False)
lats1 = latsin

fig=plt.figure()
# setup basemap
m=Basemap(resolution='l',projection='lcc',lon_0=-100,lat_0=40,width=8.e6,height=6.e6)
lons, lats = np.meshgrid(lons1,lats1)
x, y = m(lons, lats)
# interpolate land/sea mask to topo grid, mask ocean values.
# output may look 'blocky' near coastlines, since data is at much
# lower resolution than land/sea mask.
topo = maskoceans(lons, lats, topoin)
# make contour plot (ocean values will be masked)
CS=m.contourf(x,y,topo,np.arange(-300,3001,50),cmap=plt.cm.jet,extend='both')
#im=m.pcolormesh(x,y,topo,cmap=plt.cm.jet,vmin=-300,vmax=3000)
# draw coastlines.
m.drawcoastlines()
plt.title('ETOPO data with marine areas masked (original grid)')

fig=plt.figure()
# interpolate topo data to higher resolution grid (to better match
# the land/sea mask). Output looks less 'blocky' near coastlines.
nlats = 3*topoin.shape[0]
nlons = 3*topoin.shape[1]
lons = np.linspace(-180,180,nlons)
lats = np.linspace(-90,90,nlats)
lons, lats = np.meshgrid(lons, lats)
x, y = m(lons, lats)
topo = interp(topoin,lons1,lats1,lons,lats,order=1)
# interpolate land/sea mask to topo grid, mask ocean values.
topo = maskoceans(lons, lats, topo)
# make contour plot (ocean values will be masked)
CS=m.contourf(x,y,topo,np.arange(-300,3001,50),cmap=plt.cm.jet,extend='both')
#im=m.pcolormesh(x,y,topo,cmap=plt.cm.jet,vmin=-300,vmax=3000)
# draw coastlines.
m.drawcoastlines()
plt.title('ETOPO data with marine areas masked (data on finer grid)')
plt.show()
