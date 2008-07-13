from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
# setup albers equal area conic basemap
# lat_1 is first standard parallel.
# lat_2 is second standard parallel.
# lon_0,lat_0 is central point.
m = Basemap(width=8000000,height=7000000,
            resolution='l',projection='aea',\
            lat_1=40.,lat_2=60,lon_0=35,lat_0=50)
m.drawcoastlines()
m.drawcountries()
m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
m.drawmapboundary(fill_color='aqua') 
# draw tissot's indicatrix to show distortion.
ax = plt.gca()
for y in np.linspace(m.ymax/20,19*m.ymax/20,10):
    for x in np.linspace(m.xmax/20,19*m.xmax/20,12):
        lon, lat = m(x,y,inverse=True)
        seg = m.tissot(lon,lat,1.,100)
        poly = Polygon(seg,facecolor='green',zorder=10,alpha=0.5)
        ax.add_patch(poly)
plt.title("Albers Equal Area Projection")
plt.savefig('aea.png')
