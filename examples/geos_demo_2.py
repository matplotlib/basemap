"""geos_demo_2.py

This script shows how to plot data onto the Geostationary Satellite projection
when the data is from a portion of the full Earth image. The script assumes that
the data is already contained in a regular grid in the geos projection and that
the corner points of the data grid are known in lat-long.

Dependencies: Matplotlib, Basemap toolkit, Python Imaging Library

"""
from PIL import Image
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import pil_to_array

plot_name = 'geos_demo.png'
overlay_color = 'black'

# read in jpeg image to rgb array
pilImage = Image.open('200706041200-msg-ch01-SAfrica.jpg')
#data = asarray(pilImage)
data = pil_to_array(pilImage)
data = data[:, :, 0] # get data from first channel in the image

# define data region and projection parameters
ll_lon = 9.74
ll_lat = -35.55
ur_lon = 48.45
ur_lat = 0.2
lon_0 = 0.0
satellite_height = 35785831.0

fig = plt.figure(figsize=(7,7))
ax = fig.add_axes((0.1,0.1,0.8,0.8))
# create Basemap instance for a Geostationary projection.
m = Basemap(projection='geos', lon_0=lon_0, satellite_height=satellite_height,
            resolution='l', llcrnrlon=ll_lon, llcrnrlat=ll_lat, urcrnrlon=ur_lon, urcrnrlat=ur_lat)
# add data
m.imshow(data, cmap=plt.cm.gray, interpolation='nearest')
plt.clim(0, 255)
# draw coastlines.
m.drawcoastlines(linewidth=0.5, color=overlay_color)
m.drawcountries(linewidth=0.5, color=overlay_color)
# can't label meridians on bottom, because labels would
# be outside map projection region.
m.drawmeridians(np.arange(10,76,5), labels=[0,0,1,0], color=overlay_color)
m.drawparallels(np.arange(-90,90,5), labels=[1,0,0,0], color=overlay_color)
# add timestamp and save
fig = plt.gcf()
fig.text(x=0.275, y=0.025, s=u'Meteosat-9 VIS 0.6 channel - 12:00 UTC 04/06/2007\n    \N{COPYRIGHT SIGN} EUMETSAT 2007',
            horizontalalignment='left',
            verticalalignment='bottom',
            fontsize=10,
            fontweight='bold',
            bbox=dict(facecolor='gray', alpha=0.25, pad=15))
fig.set_size_inches((8, 6))
plt.title('Meteosat Geostationary Satellite Image - Portion of Full Earth',y=1.05,fontsize=12)

plt.show()
#fig.savefig(plot_name)
#print 'Plot saved to %s' % (plot_name)
