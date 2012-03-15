from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# Set up the primary map
fig = plt.figure()
ax = fig.add_subplot(111)
bmap =\
Basemap(projection='lcc',width=6000.e3,height=4000.e3,lon_0=-90,lat_0=40,resolution='l',ax=ax)
bmap.fillcontinents(color='coral', lake_color='aqua')
bmap.drawcountries()
bmap.drawstates()
bmap.drawmapboundary(fill_color='aqua')
bmap.drawcoastlines()
plt.title('map with an inset showing where the map is')

# axes for inset map.
axin = inset_axes(bmap.ax,width="30%",height="30%",loc=4)

# inset map is global, with primary map projection region drawn on it.
omap = Basemap(projection='ortho',lon_0=-105,lat_0=40,ax=axin,anchor='NE')
omap.drawcountries(color='white')
omap.fillcontinents(color='gray') #color = 'coral'               
bx, by = omap(bmap.boundarylons, bmap.boundarylats)
xy = list(zip(bx,by))
mapboundary = Polygon(xy,edgecolor='red',linewidth=2,fill=False)
omap.ax.add_patch(mapboundary)

plt.show()
