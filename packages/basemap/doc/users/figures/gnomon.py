from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
m = Basemap(width=15.e6,height=15.e6,\
            projection='gnom',lat_0=60.,lon_0=-30.)
m.drawmapboundary(fill_color='aqua')
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
m.drawparallels(np.arange(10,90,20))
m.drawmeridians(np.arange(-180,180,30))
plt.title('Gnomonic Projection')
plt.show()
