from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
# create new figure, axes instances.
fig=plt.figure()
ax=fig.add_axes([0.1,0.1,0.8,0.8])
# setup Mercator map projection.
m = Basemap(llcrnrlat=-36.,llcrnrlon=15.,urcrnrlat=-19.,urcrnrlon=36.,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='i',projection='merc')
# jhblat, jhblon are lat/lon of Johannesburg
jhblat = -26.1367; jhblon = 28.2411
# draw equidistant circle 500 and 1000 km around Johannesburg
m.drawequidistantcircle(jhblon, jhblat, 500, color='b')
m.drawequidistantcircle(jhblon, jhblat, 1000, color='r')
# Draw in map details
m.drawcoastlines()
m.drawcountries()
m.fillcontinents()
# draw parallels
m.drawparallels(np.arange(-90,0,10),labels=[1,1,0,1])
# draw meridians
m.drawmeridians(np.arange(-180,180,10),labels=[1,1,0,1])
ax.set_title('Equidistant circles around Johannesburg')
plt.show()
