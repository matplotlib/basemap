from matplotlib.toolkits.basemap import Basemap, basemap_datadir
from pylab import show, title, arange, figure, draw, ion, ioff, clf
import cPickle, time, sys, os

# turn interactive mode on.
ion()
# create new figure
fig=figure()
# create Basemap instance. Use 'crude' resolution coastlines.
m = Basemap(llcrnrlon=-11.,llcrnrlat=50.5,urcrnrlon=-5.,urcrnrlat=56.,
            resolution='c',projection='tmerc',lon_0=-8.,lat_0=0.)
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
print 'plotting with crude res boundaries ...'
title("Crude Res Boundaries ('c')",y=1.05)
draw()
time.sleep(5)

# create new figure
#fig=figure()
clf()
# create Basemap instance. Use 'low' resolution coastlines.
m = Basemap(llcrnrlon=-11.,llcrnrlat=50.5,urcrnrlon=-5.,urcrnrlat=56.,
            resolution='l',projection='tmerc',lon_0=-8.,lat_0=0.)
# draw coastlines and fill continents.
m.drawcoastlines()
m.fillcontinents()
# draw political boundaries.
m.drawcountries()
# draw parallels
m.drawparallels(circles,labels=[1,1,0,0])
# draw meridians
m.drawmeridians(meridians,labels=[0,0,1,1])
print 'plotting with low res boundaries ...'
title("Low Res Boundaries ('l')",y=1.05)
draw()
time.sleep(5)

t1 = time.clock()
# create new figure
#fig=figure()
clf()
print 'plotting with intermediate res boundaries ...'
# create Basemap instance. Use 'intermediate' resolution coastlines.
m = Basemap(llcrnrlon=-11.,llcrnrlat=50.5,urcrnrlon=-5.,urcrnrlat=56.,
            resolution='i',projection='tmerc',lon_0=-8.,lat_0=0.)
print time.clock()-t1,'seconds to create class instance with intermediate res coastlines'
# cPickle the class instance.
cPickle.dump(m,open('map.pickle','wb'),-1)
# draw coastlines and fill continents.
m.drawcoastlines()
m.fillcontinents()
# draw political boundaries.
m.drawcountries()
# draw parallels
m.drawparallels(circles,labels=[1,1,0,0])
# draw meridians
m.drawmeridians(meridians,labels=[0,0,1,1])
title("Intermediate Res Boundaries ('i')",y=1.05)
draw()
time.sleep(5)

# create new figure
#fig=figure()
clf()
# read cPickle back in and plot it again (should be much faster).
t1 = time.clock()
m2 = cPickle.load(open('map.pickle','rb'))
print time.clock()-t1,' to read the intermediate res coastline class instance back in from a cPickle'
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
print 'plotting with intermediate res boundaries from a saved pickle ...'
title("Intermediate Res Boundaries ('i') with major rivers",y=1.05)
draw()
time.sleep(5)
# only do high-res coastlines if they are installed, otherwise stop here.
if not os.path.isfile(os.path.join(basemap_datadir,'gshhs_h.txt')):
    sys.exit(0)

import time
t1 = time.clock()
# create new figure
#fig=figure()
clf()
print 'plotting with high res boundaries ...'
# create Basemap instance. Use 'high' resolution coastlines.
m = Basemap(llcrnrlon=-11.,llcrnrlat=50.5,urcrnrlon=-5.,urcrnrlat=56.,
            resolution='h',projection='tmerc',lon_0=-8.,lat_0=0.)
print time.clock()-t1,'seconds to create class instance with high res coastlines'
# cPickle the class instance.
cPickle.dump(m,open('map.pickle','wb'),-1)
# draw coastlines and fill continents.
m.drawcoastlines()
m.fillcontinents()
# draw political boundaries.
m.drawcountries()
# draw parallels
m.drawparallels(circles,labels=[1,1,0,0])
# draw meridians
m.drawmeridians(meridians,labels=[0,0,1,1])
title("High Res Boundaries ('i')",y=1.05)
draw()
time.sleep(5)

# create new figure
#fig=figure()
clf()
# read cPickle back in and plot it again (should be much faster).
t1 = time.clock()
m2 = cPickle.load(open('map.pickle','rb'))
print time.clock()-t1,' to read the high res coastline class instance back in from a cPickle'
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
print 'plotting with high res boundaries from a saved pickle ...'
title("High Res Boundaries ('i') with major rivers",y=1.05)
ioff()
show()
