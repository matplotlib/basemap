from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
# lon_0, lat_0 are the center point of the projection.
# satellite_height is the altitude of the camera.
# resolution = 'l' means use low resolution coastlines.
h = 3000.
m = Basemap(projection='nsper',lon_0=-105,lat_0=40,
        satellite_height=h*1000.,resolution='l')
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
m.drawmapboundary(fill_color='aqua') 
plt.title("Full Disk Near-Sided Perspective Projection %d km above earth" %
        h,fontsize=10)
plt.show()
