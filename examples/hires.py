from matplotlib.toolkits.basemap import Basemap
from pylab import show, title, arange, figure, clf
import cPickle, time

# create figure with aqua background (will be oceans)
fig = figure()

# create Basemap instance. Use 'high' resolution coastlines.
t1 = time.clock()
#m = Basemap(llcrnrlon=-10.5,llcrnrlat=49.5,urcrnrlon=3.5,urcrnrlat=59.5,
#            resolution='h',projection='tmerc',lon_0=-4,lat_0=0)
m = Basemap(width=920000,height=1100000,
            resolution='h',projection='tmerc',lon_0=-4.2,lat_0=54.6)
# make sure countries and rivers are loaded
m.drawcountries()
m.drawrivers()
print time.clock()-t1,' secs to create original Basemap instance'

# cPickle the class instance.
cPickle.dump(m,open('map.pickle','wb'),-1)

# clear the figure
clf()
# read cPickle back in and plot it again (should be much faster).
t1 = time.clock()
m2 = cPickle.load(open('map.pickle','rb'))
# draw coastlines and fill continents.
m.drawcoastlines()
# fill continents and lakes
m.fillcontinents(color='coral',lake_color='aqua')
# draw political boundaries.
m.drawcountries(linewidth=1)
# fill map projection region light blue (this will
# paint ocean areas same color as lakes).
m.drawmapboundary(fill_color='aqua')
# draw major rivers.
m.drawrivers(color='b')
print time.clock()-t1,' secs to plot using using a pickled Basemap instance'
# draw parallels
circles = arange(48,65,2).tolist()
m.drawparallels(circles,labels=[1,1,0,0])
# draw meridians
meridians = arange(-12,13,2)
m.drawmeridians(meridians,labels=[0,0,1,1])
title("High-Res British Isles",y=1.075)
show()
