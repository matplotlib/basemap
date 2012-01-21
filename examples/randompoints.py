import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import rgb2hex
from mpl_toolkits.basemap import Basemap
from numpy.random import uniform

# Plot a bunch of randomly distributed points on the earth.

# set up stereographic map centered on N. Pole.
m = Basemap(lon_0=-105,boundinglat=20.,round=True,
            resolution='l',area_thresh=10000.,projection='npstere')
# number of points to plot.
npts = 1000
# generate random points on a sphere,
# so that every small area on the sphere is expected
# to have the same number of points.
# http://mathworld.wolfram.com/SpherePointPicking.html
u = uniform(0.,1.,size=npts)
v = uniform(0.,1.,size=npts)
lons = 360.*u
lats = (180./np.pi)*np.arccos(2*v-1) - 90.
z = uniform(0,100,size=npts) # this field controls color of dots.
# transform lons and lats to map coordinates.
x,y = m(lons,lats)
# plot them as filled circles on the map.
# first, create a figure.
fig=plt.figure()
# draw colored markers.
# use zorder=10 to make sure markers are drawn last.
# (otherwise they are covered up when continents are filled)
m.scatter(x,y,25,z,cmap=plt.cm.jet,marker='o',edgecolors='none',zorder=10) 
# plot colorbar for markers.
m.colorbar(pad='8%')
# create a list of strings containing z values
# or, plot actual numbers as color-coded text strings.
#zn = [ '%2i' % zz for zz in z ]
## plot numbers on map, colored by value.
#for numstr,zval,xpt,ypt in zip(zn,z,x,y):
#    # only plot values inside map region.
#    if xpt > m.xmin and xpt < m.xmax and ypt > m.ymin and ypt < m.ymax:
#        hexcolor = rgb2hex(plt.cm.jet(zval/100.)[:3])
#        plt.text(xpt,ypt,numstr,fontsize=9,weight='bold',color=hexcolor)
# draw coasts and fill continents/lakes.
m.drawcoastlines(linewidth=0.5,color='y')
m.drawcountries(color='y')
m.drawstates(color='y')
m.fillcontinents(color='grey',lake_color='black')
# color ocean areas 
m.drawmapboundary(fill_color='black')
# draw parallels and meridians.
m.drawparallels(np.arange(-80,81,20),color='y')
m.drawmeridians(np.arange(-180,181,45),labels=[1,1,0,0],color='y')
plt.title('Colored Random Points',y=1.05)
plt.show()
