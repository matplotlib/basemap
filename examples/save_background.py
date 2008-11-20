import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

# this example shows how to save a map background and
# reuse it in another figure.

# make sure we have all the same properties on all figs
figprops = dict(figsize=(8,6), dpi=100, facecolor='white')

# generate the first figure.
fig1 = plt.figure(1,**figprops)
ax1 = fig1.add_subplot(111)
# create basemap instance, plot coastlines.
map = Basemap(projection='moll',lon_0=0)
map.drawcoastlines()
map.drawmapboundary(fill_color='aqua')
map.fillcontinents(color='coral',lake_color='aqua')
fig1.canvas.draw()
background = fig1.canvas.copy_from_bbox(fig1.bbox)
fig1.savefig('figure1.png', dpi=100)


# generate the second figure, re-using the background
# from figure 1.
fig2 = plt.figure(2,frameon=False,**figprops)
ax2 = fig2.add_subplot(111, frameon=False, xticks=[], yticks=[])
# restore previous background.
fig2.canvas.restore_region(background)
# draw parallels and meridians on existing background.
map.drawparallels(range(-90,90,30))
map.drawmeridians(range(-180,180,60))
fig2.savefig('figure2.png', dpi=100)

print 'images saved in figure1.png and figure2.png'
