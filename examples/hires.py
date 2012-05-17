from __future__ import print_function
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import pickle, time

# create figure with aqua background (will be oceans)
fig = plt.figure()

# create Basemap instance. Use 'high' resolution coastlines.
t1 = time.clock()
#m = Basemap(llcrnrlon=-10.5,llcrnrlat=49.5,urcrnrlon=3.5,urcrnrlat=59.5,
#            resolution='h',projection='tmerc',lon_0=-4,lat_0=0)
m = Basemap(width=920000,height=1100000,
            resolution='f',projection='tmerc',lon_0=-4.2,lat_0=54.6)
# make sure countries and rivers are loaded
m.drawcountries()
m.drawrivers()
print(time.clock()-t1,' secs to create original Basemap instance')

# pickle the class instance.
pickle.dump(m,open('map.pickle','wb'),-1)

# clear the figure
plt.clf()
# read pickle back in and plot it again (should be much faster).
t1 = time.clock()
m2 = pickle.load(open('map.pickle','rb'))
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
print(time.clock()-t1,' secs to plot using using a pickled Basemap instance')
# draw parallels
circles = np.arange(48,65,2).tolist()
m.drawparallels(circles,labels=[1,1,0,0])
# draw meridians
meridians = np.arange(-12,13,2)
m.drawmeridians(meridians,labels=[0,0,1,1])
plt.title("High-Res British Isles",y=1.04)
plt.show()
