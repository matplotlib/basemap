from matplotlib.toolkits.basemap import Basemap
from pylab import *
import cPickle

# create Basemap instance. Use 'crude' resolution coastlines.
m = Basemap(llcrnrlon=-11.,llcrnrlat=50.5,urcrnrlon=-5.,urcrnrlat=56.,
            resolution='c',area_thresh=1000.,projection='tmerc',lon_0=-8.,lat_0=0.)
# create figure with same aspect ratio as map.
fig=m.createfigure()
# draw coastlines and fill continents.
m.drawcoastlines()
m.fillcontinents()
# draw political boundaries.
m.drawcountries()
# draw parallels
circles = arange(50,60,1).tolist()
m.drawparallels(circles,labels=[1,1,0,0])
# draw meridians
meridians = arange(-12,0,1)
m.drawmeridians(meridians,labels=[0,0,1,1])
title("Crude Res Boundaries ('c')",y=1.075)
show()

# create Basemap instance. Use 'low' resolution coastlines.
m = Basemap(llcrnrlon=-11.,llcrnrlat=50.5,urcrnrlon=-5.,urcrnrlat=56.,
            resolution='l',area_thresh=1000.,projection='tmerc',lon_0=-8.,lat_0=0.)
# create figure with same aspect ratio as map.
fig=m.createfigure()
# draw coastlines and fill continents.
m.drawcoastlines()
m.fillcontinents()
# draw political boundaries.
m.drawcountries()
# draw parallels
m.drawparallels(circles,labels=[1,1,0,0])
# draw meridians
m.drawmeridians(meridians,labels=[0,0,1,1])
title("Low Res Boundaries ('l')",y=1.075)
show()

import time
t1 = time.clock()
# create Basemap instance. Use 'intermediate' resolution coastlines.
m = Basemap(llcrnrlon=-11.,llcrnrlat=50.5,urcrnrlon=-5.,urcrnrlat=56.,
            resolution='i',area_thresh=1000.,projection='tmerc',lon_0=-8.,lat_0=0.)
print time.clock()-t1,'seconds to create class instance with intermediate res coastlines'
# cPickle the class instance.
cPickle.dump(m,open('map.pickle','wb'),-1)
# create figure with same aspect ratio as map.
fig=m.createfigure()
# draw coastlines and fill continents.
m.drawcoastlines()
m.fillcontinents()
# draw political boundaries.
m.drawcountries()
# draw parallels
m.drawparallels(circles,labels=[1,1,0,0])
# draw meridians
m.drawmeridians(meridians,labels=[0,0,1,1])
title("Intermediate Res Boundaries ('i')",y=1.075)
show()

# read cPickle back in and plot it again (should be much faster).
t1 = time.clock()
m2 = cPickle.load(open('map.pickle','rb'))
print time.clock()-t1,' to read the intermediate res coastline class instance back in from a cPickle'
# create figure with same aspect ratio as map.
fig=m2.createfigure()
# draw coastlines and fill continents.
m2.drawcoastlines()
m2.fillcontinents()
# draw political boundaries.
m2.drawcountries(linewidth=1.0)
# draw major rivers.
m2.drawrivers(color='b')
# draw parallels
m2.drawparallels(circles,labels=[1,1,0,0])
# draw meridians
m2.drawmeridians(meridians,labels=[0,0,1,1])
title("Intermediate Res Boundaries ('i') with major rivers",y=1.075)
show()
