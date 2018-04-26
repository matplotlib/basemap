from __future__ import (absolute_import, division, print_function)

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

fig=plt.figure()


# Example of Great Circles which exit and reenter the map
m = Basemap(projection='robin', lat_0=0, lon_0=0,resolution='c')

m.drawmapboundary(fill_color='#ffffff',color='#ffffff')
m.fillcontinents(color='#f2f2f2',lake_color='#ffffff')
m.drawcoastlines(color='#e0e0e0')
m.drawcountries(color='#e0e0e0')

parallels = np.arange(-90,90,10.)
meridians = np.arange(10.,351.,20.)
m.drawparallels(parallels,labels=[False,False,False,False], color='#d3d3d3',dashes=[1,3])
m.drawmeridians(meridians,labels=[False,False,False,False], color='#d3d3d3',dashes=[1,3])

# Choose the lat and longtitude of two points

p1 = (45.27,-75.42) # roughly Ottawa
p2 = (44.05,-4.28)  # roughly Paris
p3 = (-38.58,145.05) # roughly Victoria

la1, lo1 = p1
la2, lo2 = p2
la3, lo3 = p3

# Drawing points; you need to convert from lon-lat to xy-cartesian
x1,y1 = m(lo1,la1)
x2,y2 = m(lo2,la2)
x3,y3 = m(lo3,la3)

# Convert back by setting inverse=True
lon,lat = m(x1,y1,inverse=True)

# Plot pionts using markers
m.plot(x1,y1, marker='.', markersize=8, color='#000000')  
m.plot(x2,y2, marker='.', markersize=8, color='#000000')  
m.plot(x3,y3, marker='.', markersize=8, color='#000000')  

# Of note, the map uses metres as it's smallest distance, so 1000*x is x km
# The offset is in projection cooordinates
plt.text(x1+100000,y1+100000,'Ottawa')
plt.text(x2+100000,y2+100000,'Paris')
plt.text(x3+100000,y3+100000,'Victoria')

# Draw a great circle line joining them
m.drawgreatcircle(lo1,la1,lo2,la2,linewidth=1,color='#000000',alpha=1,del_s=100)
m.drawgreatcircle(lo2,la2,lo3,la3,linewidth=1,color='#000000',alpha=1,del_s=100)

# Drawing a great circle which exits and reenters the map works on certain projections
m.drawgreatcircle(lo1,la1,lo3,la3,linewidth=1,color='#FF0000',alpha=1,del_s=100)

plt.show()
