from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
# setup lambert azimuthal equal area basemap.
# lat_ts is latitude of true scale.
# lon_0,lat_0 is central point.
m = Basemap(width=12000000,height=8000000,
            resolution='l',projection='laea',\
            lat_ts=50,lat_0=50,lon_0=-107.)
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
m.drawmapboundary(fill_color='aqua') 
plt.title("Lambert Azimuthal Equal Area Projection")
plt.savefig('laea.png')
