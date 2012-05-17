# example demonstrating how to draw a great circle on a map.
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import sys

# setup lambert azimuthal map projection.
# create new figure
fig=plt.figure()
m = Basemap(llcrnrlon=-100.,llcrnrlat=20.,urcrnrlon=20.,urcrnrlat=60.,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='c',area_thresh=10000.,projection='merc',\
            lat_0=40.,lon_0=-20.,lat_ts=20.)
# nylat, nylon are lat/lon of New York
nylat = 40.78
nylon = -73.98
# lonlat, lonlon are lat/lon of London.
lonlat = 51.53
lonlon = 0.08
# find 1000 points along the great circle.
#x,y = m.gcpoints(nylon,nylat,lonlon,lonlat,1000)
# draw the great circle.
#m.plot(x,y,linewidth=2)
# drawgreatcircle performs the previous 2 steps in one call.
m.drawgreatcircle(nylon,nylat,lonlon,lonlat,linewidth=2,color='b')
m.drawcoastlines()
m.fillcontinents()
# draw parallels
circles = np.arange(10,90,20)
m.drawparallels(circles,labels=[1,1,0,1])
# draw meridians
meridians = np.arange(-180,180,30)
m.drawmeridians(meridians,labels=[1,1,0,1])
plt.title('Great Circle from New York to London (Mercator)')
sys.stdout.write('plotting Great Circle from New York to London (Mercator)\n')

# create new figure
fig=plt.figure()
# setup a gnomonic projection.
m = Basemap(llcrnrlon=-100.,llcrnrlat=20.,urcrnrlon=20.,urcrnrlat=60.,\
            resolution='c',area_thresh=10000.,projection='gnom',\
            lat_0=40.,lon_0=-45.)
# nylat, nylon are lat/lon of New York
nylat = 40.78
nylon = -73.98
# lonlat, lonlon are lat/lon of London.
lonlat = 51.53
lonlon = 0.08
# find 1000 points along the great circle.
#x,y = m.gcpoints(nylon,nylat,lonlon,lonlat,1000)
# draw the great circle.
#m.plot(x,y,linewidth=2)
# drawgreatcircle performs the previous 2 steps in one call.
m.drawgreatcircle(nylon,nylat,lonlon,lonlat,linewidth=2,color='b')
m.drawcoastlines()
m.fillcontinents()
# draw parallels
circles = np.arange(10,90,20)
m.drawparallels(circles,labels=[0,1,0,0])
# draw meridians
meridians = np.arange(-180,180,30)
m.drawmeridians(meridians,labels=[1,1,0,1])
plt.title('Great Circle from New York to London (Gnomonic)')
sys.stdout.write('plotting Great Circle from New York to London (Gnomonic)\n')
plt.show()
