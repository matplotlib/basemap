from __future__ import unicode_literals
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# example showing how to create custom tick labels for a cylindrical
# projection.

def lat2str(deg):
    min = 60 * (deg - np.floor(deg))
    deg = np.floor(deg)
    dir = 'N'
    if deg < 0:
        if min != 0.0:
            deg += 1.0
            min -= 60.0
        dir = 'S'
    return ("%d\N{DEGREE SIGN} %g' %s") % (np.abs(deg),np.abs(min),dir)

def lon2str(deg):
    min = 60 * (deg - np.floor(deg))
    deg = np.floor(deg)
    dir = 'E'
    if deg < 0:
        if min != 0.0:
            deg += 1.0
            min -= 60.0
        dir = 'W'
    return ("%d\N{DEGREE SIGN} %g' %s") % (np.abs(deg),np.abs(min),dir)

# (1) use matplotlib custom tick formatter
#     instead of Basemap labelling methods.

# create figure.
fig=plt.figure()
# create Basemap instance (regular lat/lon projection).
# suppress_ticks=False allows custom axes ticks to be used
# Ticks are suppressed by default, so Basemap methods
# drawparallels and drawmeridians used to draw labelled lat/lon grid.
m = Basemap(llcrnrlon=-156.5,llcrnrlat=18.75,urcrnrlon=-154.5,urcrnrlat=20.5,
            resolution='h',projection='cyl',suppress_ticks=False)
# draw coastlines, fill land and lake areas.
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
# background color will be used for oceans.
m.drawmapboundary(fill_color='aqua')
# get axes instance.
ax = plt.gca()
# add custom ticks.
# This only works for projection='cyl'.
def xformat(x, pos=None): return lon2str(x)
xformatter = FuncFormatter(xformat)
ax.xaxis.set_major_formatter(xformatter)
def yformat(y, pos=None): return lat2str(y)
yformatter = FuncFormatter(yformat)
ax.yaxis.set_major_formatter(yformatter)
ax.fmt_xdata = lambda x: lon2str(x)
ax.fmt_ydata = lambda y: lat2str(y)
ax.grid()
ax.set_title('Hawaii')

# (2) use Basemap labelling methods, but pass a
#     custom formatting function with the 'fmt' keyword.

# create figure.
fig = plt.figure()
# create Basemap instance.
m = Basemap(llcrnrlon=-156.5,llcrnrlat=18.75,urcrnrlon=-154.5,urcrnrlat=20.5,
            resolution='h',projection='cyl')
# draw coastlines, fill land and lake areas.
m.drawcoastlines()
m.fillcontinents(color='coral',lake_color='aqua')
# background color will be used for oceans.
m.drawmapboundary(fill_color='aqua')
# label meridians and parallels, passing string formatting function
# with 'fmt' keyword.
m.drawparallels(np.linspace(18,21,7),labels=[1,0,0,0],fmt=lat2str,dashes=[2,2])
m.drawmeridians(np.linspace(-157,-154,7),labels=[0,0,0,1],fmt=lon2str,dashes=[2,2])
plt.title('Hawaii')

plt.show()
