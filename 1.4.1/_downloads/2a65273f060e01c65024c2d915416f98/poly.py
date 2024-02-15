from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
# setup polyconic basemap 
# by specifying lat/lon corners and central point.
# area_thresh=1000 means don't plot coastline features less
# than 1000 km^2 in area.
m = Basemap(llcrnrlon=-35.,llcrnrlat=-30,urcrnrlon=80.,urcrnrlat=50.,\
            resolution='l',area_thresh=1000.,projection='poly',\
            lat_0=0.,lon_0=20.)
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
m.drawmapboundary(fill_color='aqua') 
plt.title("Polyconic Projection")
plt.show()
