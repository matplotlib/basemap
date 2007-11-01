from matplotlib.toolkits.basemap import Basemap
from pylab import *
import time

# create new figure
fig=figure()
# create Basemap instance. Use 'high' resolution coastlines.
t1 = time.clock()
m = Basemap(llcrnrlon=-11.,llcrnrlat=49.,urcrnrlon=5.,urcrnrlat=59.,
            resolution='h',projection='tmerc',lon_0=-8.,lat_0=0.)
print 'time to create instance with high-res boundaries = ',time.clock()-t1
# draw coastlines and fill continents.
m.drawcoastlines()
m.fillcontinents()
# draw political boundaries.
m.drawcountries(linewidth=1)
# draw major rivers.
m.drawrivers(color='b')
# draw parallels
circles = arange(48,65,2).tolist()
m.drawparallels(circles,labels=[1,1,0,0])
# draw meridians
meridians = arange(-12,13,2)
m.drawmeridians(meridians,labels=[0,0,1,1])
title("High-Res British Isles",y=1.075)
show()
