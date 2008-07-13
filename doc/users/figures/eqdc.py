from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
# setup equidistant conic basemap.
# lat_1 is first standard parallel.
# lat_2 is second standard parallel.
# lon_0,lat_0 is central point.
# resolution = 'l' for low-resolution coastlines.
m = Basemap(width=12000000,height=9000000,
            resolution='l',projection='eqdc',\
            lat_1=45.,lat_2=55,lat_0=50,lon_0=-107.)
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
m.drawmapboundary(fill_color='aqua') 
plt.title("Equidistant Conic Projection")
plt.savefig('eqdc.png')
