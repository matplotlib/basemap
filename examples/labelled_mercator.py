# This example shows how to set up an axis tick Formatter object so that
# Mercator projection Y-axis ticks are labelled with latitude.
# The Mercator projection Y-axis is not latitude, so if you don't do this,
# your Y-axis will be labelled in native map projection coordinates in meters.
#
# The MercYAxisFormatter class translates between native coordinates and
# latitude.  We tell an Axes to use it by calling ax.yaxis.set_major_formatter.

# Courtesy of Michael Brady of NASA/JPL.

import matplotlib
from matplotlib.toolkits.basemap import Basemap
from pylab import *

class MercYAxisFormatter(matplotlib.ticker.Formatter):
   """The format function for Mercur projection Y-axis.  Translates plot y
   in meters to latitude.
   """

   def __init__(self, baseMap):
      """baseMap is the Basemap object that will be used to translate between
      lon/lat and native map coordinates.
      """
      self.baseMap = baseMap

   def __call__(self, y, pos=1):
      """Return the label for tick value y at position pos.
      """
      lon, lat = self.baseMap(0.0, y, inverse=True)
      return "%g" % lat

m = Basemap(-180.,-80.,180.,80.,\
            resolution='c',area_thresh=10000.,projection='merc',\
            lat_ts=20.)
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
fig.add_axes([0.1,0.1,0.8,0.8])
ax = gca() # get current axis instance
ax.update_datalim(((m.llcrnrx, m.llcrnry),(m.urcrnrx,m.urcrnry)))
ax.set_xlim((m.llcrnrx, m.urcrnrx))
ax.set_ylim((m.llcrnry, m.urcrnry))
m.drawcoastlines(ax)
m.drawcountries(ax)
m.drawstates(ax)
m.fillcontinents(ax)
# draw parallels
delat = 30.
circles = arange(0.,90.,delat).tolist()+\
          arange(-delat,-90,-delat).tolist()
m.drawparallels(ax,circles)
# convert parallels to native map projection coordinates
nativeCoordCircles = m([0]*len(circles), circles)[1]
ax.set_yticks(nativeCoordCircles)
ax.yaxis.set_major_formatter(MercYAxisFormatter(m))
# draw meridians
delon = 60.
lon1 = int(m.llcrnrlon/delon)*delon
lon2 = (int(m.urcrnrlon/delon)+1)*delon
meridians = arange(lon1,lon2,delon)
m.drawmeridians(ax,meridians)
ax.set_xticks(meridians)
title('Mercator with Longitude and Latitude Labels')
print 'plotting Mercator Y-axis latitude ticks example, close plot window to proceed ...'
show()
