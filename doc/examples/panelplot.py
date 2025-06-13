from __future__ import (absolute_import, division, print_function)

from mpl_toolkits.basemap import Basemap
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator
import numpy as np
import matplotlib.pyplot as plt


# read in data on lat/lon grid.
hgt  = np.loadtxt('500hgtdata.gz')
lons = np.loadtxt('500hgtlons.gz')
lats = np.loadtxt('500hgtlats.gz')
lons, lats = np.meshgrid(lons, lats)

# Example to show how to make multi-panel plots.

# 2-panel plot, oriented vertically, colorbar on bottom.

rcParams['figure.subplot.hspace'] = 0.4 # more height between subplots
rcParams['figure.subplot.wspace'] = 0.5 # more width between subplots

# create new figure
fig=plt.figure()
# panel 1
mnh = Basemap(lon_0=-105,boundinglat=20.,
             resolution='c',area_thresh=10000.,projection='nplaea')
xnh,ynh = mnh(lons,lats)
ax = fig.add_subplot(211)
CS = mnh.contour(xnh,ynh,hgt,15,linewidths=0.5,colors='k')
CS = mnh.contourf(xnh,ynh,hgt,15,cmap=plt.cm.Spectral)
# colorbar on bottom.
mnh.colorbar(location='bottom',pad='12%',ticks=CS.levels[0::4])
mnh.drawcoastlines(linewidth=0.5)
delat = 30.
circles = np.arange(0.,90.,delat).tolist()+\
          np.arange(-delat,-90,-delat).tolist()
mnh.drawparallels(circles,labels=[1,0,0,0])
delon = 45.
meridians = np.arange(0,360,delon)
mnh.drawmeridians(meridians,labels=[1,0,0,1])
plt.title('NH 500 hPa Height (cm.Spectral)')

# panel 2
msh = Basemap(lon_0=-105,boundinglat=-20.,
             resolution='c',area_thresh=10000.,projection='splaea')
xsh,ysh = msh(lons,lats)
ax = fig.add_subplot(212)
CS = msh.contour(xsh,ysh,hgt,15,linewidths=0.5,colors='k')
CS = msh.contourf(xsh,ysh,hgt,15,cmap=plt.cm.Spectral)
# colorbar on bottom.
msh.colorbar(location='bottom',pad='12%',ticks=MultipleLocator(320))
msh.drawcoastlines(linewidth=0.5)
msh.drawparallels(circles,labels=[1,0,0,0])
msh.drawmeridians(meridians,labels=[1,0,0,1])
plt.title('SH 500 hPa Height (cm.Spectral)')

# 2-panel plot, oriented horizontally, colorbar on right.

# adjust default subplot parameters a bit
rcParams['figure.subplot.left'] = 0.1   # move left edge of subplot over a bit
rcParams['figure.subplot.right'] = 0.85
rcParams['figure.subplot.top'] = 0.85

# panel 1
fig = plt.figure()
ax = fig.add_subplot(121)
CS = mnh.contour(xnh,ynh,hgt,15,linewidths=0.5,colors='k')
CS = mnh.contourf(xnh,ynh,hgt,15,cmap=plt.cm.RdBu)
# colorbar on right
mnh.colorbar(location='right',ticks=MultipleLocator(160))
mnh.drawcoastlines(linewidth=0.5)
mnh.drawparallels(circles,labels=[1,0,0,0])
mnh.drawmeridians(meridians,labels=[1,0,0,1])
plt.title('NH 500 hPa Height (cm.RdBu)')

# panel 2
ax = fig.add_subplot(122)
CS = msh.contour(xsh,ysh,hgt,15,linewidths=0.5,colors='k')
CS = msh.contourf(xsh,ysh,hgt,15,cmap=plt.cm.RdBu)
# colorbar on right.
msh.colorbar(location='right',ticks=MultipleLocator(160))
msh.drawcoastlines(linewidth=0.5)
msh.drawparallels(circles,labels=[1,0,0,0])
msh.drawmeridians(meridians,labels=[1,0,0,1])
plt.title('SH 500 hPa Height (cm.RdBu)')
plt.show()
