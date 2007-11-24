from matplotlib.toolkits.basemap import Basemap
import pylab, numpy
from matplotlib.ticker import FuncFormatter

# example showing how to create custom tick labels for a cylindrical
# projection.

def deg2str(deg, dir='E', fmt="%3.1f"):
    min = 60 * (deg - numpy.floor(deg))
    deg = numpy.floor(deg)
    if deg < 0:
        if min != 0.0:
            deg += 1.0
            min -= 60.0
        if dir=='E':
            dir='W'
        if dir=='N':
            dir='S'
    return (u"%d\N{DEGREE SIGN}" + fmt + "' %s") % (numpy.abs(deg), numpy.abs(min), dir)

# create figure.
fig=pylab.figure()
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
ax = pylab.gca()
# add custom ticks.
# This only works for projection='cyl'.
def xformat(x, pos=None): return deg2str(x, 'E', fmt="%2.0f")
xformatter = FuncFormatter(xformat)
ax.xaxis.set_major_formatter(xformatter)
def yformat(y, pos=None): return deg2str(y, 'N', fmt="%2.0f")
yformatter = FuncFormatter(yformat)
ax.yaxis.set_major_formatter(yformatter)
ax.fmt_xdata = lambda x: deg2str(x, 'E', fmt="%5.3f")
ax.fmt_ydata = lambda y: deg2str(y, 'N', fmt="%5.3f")
ax.grid()
ax.set_title('Hawaii')
pylab.show()
