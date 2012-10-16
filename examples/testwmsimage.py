"""
example showing how to use OWSlib to retrieve an image
from a WMS server and display it on a map
"""
from mpl_toolkits.basemap import Basemap
from owslib.wms import WebMapService
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import imread
import urllib2
lon_min = -118.8; lon_max = -108.6
lat_min = 22.15;  lat_max = 32.34
m = Basemap(llcrnrlon=lon_min, urcrnrlat=lat_max, urcrnrlon=lon_max,
            llcrnrlat=lat_min,resolution='i',epsg=4326)
urlbase='http://motherlode.ucar.edu:8080/thredds/wms/fmrc/NCEP/NAM/CONUS_12km/NCEP-NAM-CONUS_12km-noaaport_best.ncd?'
wms = WebMapService(urlbase)
img = wms.getmap(layers=['Temperature_height_above_ground'],
      styles=['boxfill/rainbow'],srs='EPSG:4326',
      time=datetime.utcnow().strftime('%Y-%m-%dT00:00:00.000Z'),
      bbox=(lon_min,lat_min,lon_max,lat_max),
      size=(500,500),format='image/png',elevation='2',
      colorscalerange='271.2,308',numcolorbands='20',logscale=False)
cs=m.imshow(imread(urllib2.urlopen(img.url)),origin='upper')
m.drawcoastlines(linewidth=0.25)
parallels = np.arange(20,36,2.)
a=m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
meridians = np.arange(-120,-100,2.)
b=m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
plt.show()
