from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
# setup oblique mercator basemap.
# width is width of map projection region in km (xmax-xmin_
# height is height of map projection region in km (ymax-ymin)
# lon_0, lat_0 are the central longitude and latitude of the projection.
# lat_1,lon_1 and lat_2,lon_2  are two pairs of points that define
# the projection centerline.
# Map projection coordinates are automatically rotated to true north.
# To avoid this, set no_rot=True.
# area_thresh=1000 means don't plot coastline features less
# than 1000 km^2 in area.
m = Basemap(height=16700000,width=12000000,
            resolution='l',area_thresh=1000.,projection='omerc',\
            lon_0=-100,lat_0=15,lon_2=-120,lat_2=65,lon_1=-50,lat_1=-55)
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
m.drawmapboundary(fill_color='aqua') 
plt.title("Oblique Mercator Projection")
plt.show()
