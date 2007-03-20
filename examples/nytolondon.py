# example demonstrating how to draw a great circle on a map.
from matplotlib.toolkits.basemap import Basemap
from pylab import title, arange, show, figure

# setup lambert azimuthal map projection.
# create new figure
fig=figure()
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
circles = arange(10,90,20)
m.drawparallels(circles,labels=[1,1,0,1])
# draw meridians
meridians = arange(-180,180,30)
m.drawmeridians(meridians,labels=[1,1,0,1])
title('Great Circle from New York to London (Mercator)')
print 'plotting Great Circle from New York to London (Mercator)'

# create new figure
fig=figure()
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
circles = arange(10,90,20)
m.drawparallels(circles,labels=[0,1,0,0])
# draw meridians
meridians = arange(-180,180,30)
m.drawmeridians(meridians,labels=[1,1,0,1])
title('Great Circle from New York to London (Gnomonic)')
print 'plotting Great Circle from New York to London (Gnomonic)'
show()
