# hack to draw a round polar stereographic plot.
from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# read in topo data (on a regular lat/lon grid)
# longitudes go from 20 to 380.
etopo = np.loadtxt('etopo20data.gz')
lons = np.loadtxt('etopo20lons.gz')
lats = np.loadtxt('etopo20lats.gz')

print 'min/max etopo20 data:'
print etopo.min(),etopo.max()

# create projection.
#m = Basemap(boundinglat=-20,lon_0=90,projection='spstere')
m = Basemap(boundinglat=20,lon_0=270,projection='npstere)
print m.xmin, m.xmax,m.ymin,m.ymax
# compute native map projection coordinates for lat/lon grid.
x,y = m(*np.meshgrid(lons,lats))
# make filled contour plot.
cs =\
m.contourf(x,y,etopo,np.linspace(-7500,4500,41),cmap=cm.GMT_haxby,extend='both')
# colorbar
cb = m.colorbar()
# draw coastlines.
coasts = m.drawcoastlines()
# draw parallels and meridians.
parallels = m.drawparallels(np.arange(20.,90,20.))
merids = m.drawmeridians(np.arange(0.,360.,60.))
plt.box(on=False) # don't draw axes frame.
# create clip path.
clipit =\
patches.Circle((0.5*(m.xmax+m.xmin),0.5*(m.ymax+m.ymin)),radius=0.5*(m.xmax-m.xmin),fc='none')
ax = plt.gca()
p = ax.add_patch(clipit)
p.set_clip_on(False)
# clip coastlines.
coasts.set_clip_path(clipit)
# clip contours.
for cntr in cs.collections:
    cntr.set_clip_path(clipit)
# clip meridian lines.
for merid in merids:
    lines,labels = merids[merid]
    for l in lines:
        l.set_clip_path(clipit)
plt.title('a round polar plot')
plt.show()
