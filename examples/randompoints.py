from pylab import show, title, arange, figure, title, arccos, pi, cm, text
from matplotlib.colors import rgb2hex
from matplotlib.toolkits.basemap import Basemap
from numpy.random import uniform

# Plot a bunch of randomly distributed points on the earth.

# set up stereographic map centered on N. Pole.
m = Basemap(lon_0=-105,boundinglat=30.,
            resolution='l',area_thresh=10000.,projection='npstere')
# number of points to plot.
npts = 750
# generate random points on a sphere,
# so that every small area on the sphere is expected
# to have the same number of points.
# http://mathworld.wolfram.com/SpherePointPicking.html
u = uniform(0.,1.,size=npts)
v = uniform(0.,1.,size=npts)
z = uniform(0,100,size=npts)
lons = 360.*u
lats = (180./pi)*arccos(2*v-1) - 90.
# transform lons and lats to map coordinates.
x,y = m(lons,lats)
# plot them as filled circles on the map.
# first, create a figure.
fig=figure()
# background color will be used for 'wet' areas.
fig.add_axes([0.1,0.1,0.8,0.8],axisbg='aqua')
# draw colored markers.
# use zorder=10 to make sure markers are drawn last.
# (otherwise they are covered up when continents are filled)
#m.scatter(x,y,25,z,cmap=cm.jet,marker='o',faceted=False,zorder=10) 
# create a list of strings containing z values
# or, plot actual numbers as color-coded text strings.
zn = [ '%2i' % zz for zz in z ]
# plot numbers on map, colored by value.
for numstr,zval,xpt,ypt in zip(zn,z,x,y):
    # only plot values inside map region.
    if xpt > m.xmin and xpt < m.xmax and ypt > m.ymin and ypt < m.ymax:
        hexcolor = rgb2hex(cm.jet(zval/100.)[:3])
        text(xpt,ypt,numstr,fontsize=9,weight='bold',color=hexcolor)
# draw coasts and fill continents.
m.drawcoastlines(linewidth=0.5)
m.fillcontinents(color='coral')
# draw parallels and meridians.
delat = 20.
circles = arange(0.,90.,delat).tolist()+\
          arange(-delat,-90,-delat).tolist()
m.drawparallels(circles)
delon = 45.
meridians = arange(0,360,delon)
m.drawmeridians(meridians,labels=[1,1,1,1])
title('Random Points',y=1.075)
show()
