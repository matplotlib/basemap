from pylab import *
import math, sys, random
from matplotlib.toolkits.basemap import Basemap

def get_dist(lon1,lons,lat1,lats):
    """compute great circle distance between lat1,lon1 and lats,lons"""
    arg = sin(lat1)*sin(lats)+cos(lat1)*cos(lats)*cos(lon1-lons)
    arg = clip(arg,-1.,1.)
    return arccos(arg)

d2r = math.pi/180.
# define a regular 2.5x2.5 degree grid
delat= 2.5
nlons = 144
nlats = 37
lons = delat*arange(nlons)
lats = 90.-delat*arange(nlats)
lons, lats = meshgrid(lons, lats)
lons = (d2r*lons.flat).tolist()
lats = (d2r*lats.flat).tolist()
# randomly shuffle locations.
random.shuffle(lons)
random.shuffle(lats)
lons = array(lons,'f')
lats = array(lats,'f')

# minimum separation distance in km.
rcrit = 500.

# set up lambert azimuthal map centered on N. Pole.
m = Basemap(-150.,0.,30.,0.,
            resolution='c',area_thresh=10000.,projection='laea',\
            lat_0=90.,lon_0=-105.,lat_ts=90.)

print len(lons), ' obs before thinning'

# calculate distance between each ob and all preceding obs in list.
# throw out those that are closer than rcrit.
import time
time1 = time.clock()
nob = 0
lats_out = []
lons_out = []
for lon,lat in zip(lons,lats):
   if nob:
      r = (m.rsphere/1000.)*get_dist(lon,lons[0:nob],lat,lats[0:nob])
      if min(r) > rcrit: 
          lats_out.append(lat)
          lons_out.append(lon)
   nob = nob + 1

print time.clock()-time1
print len(lons_out), ' obs after thinning'

# transform lons and lats to map coordinates.
x,y = m(array(lons_out)/d2r, array(lats_out)/d2r)
# find just those points in map region.
xx=[]
yy=[]
for xi,yi in zip(x,y):
    if (xi>m.llcrnrx and xi<m.urcrnrx) and (yi>m.llcrnry and yi<m.urcrnry):
        xx.append(xi)
        yy.append(yi)
# plot them as filled circles on the map.
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
fig.add_axes([0.1,0.1,0.8,0.8])
# use zorder=10 to make sure markers are drawn last.
# (otherwise they are covered up when continents are filled)
m.scatter(xx,yy,marker='o',c='k',s=25,zorder=10)
# draw coasts and fill continents.
m.drawcoastlines()
m.fillcontinents()
# draw parallels and meridians.
delat = 30.
circles = arange(0.,90.,delat).tolist()+\
          arange(-delat,-90,-delat).tolist()
m.drawparallels(circles)
delon = 45.
meridians = arange(0,360,delon)
m.drawmeridians(meridians,labels=[1,1,1,1])
title('Randomly Spaced Locations (Min Dist = %g km, %g points)'% (rcrit,len(lons_out)),y=1.075)
show()
