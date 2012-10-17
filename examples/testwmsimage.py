"""
example showing how to use OWSlib to retrieve an image
from a WMS server and display it on a map (using the
wmsimage convienence method)
"""
from mpl_toolkits.basemap import Basemap, pyproj
from owslib.wms import WebMapService
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import imread
import urllib2

serverurl='http://motherlode.ucar.edu:8080/thredds/wms/fmrc/NCEP/NAM/CONUS_12km/NCEP-NAM-CONUS_12km-noaaport_best.ncd?'

lon_min = -118.8; lon_max = -108.6
lat_min = 22.15;  lat_max = 32.34
m = Basemap(llcrnrlon=lon_min, urcrnrlat=lat_max,
            urcrnrlon=lon_max, llcrnrlat=lat_min,resolution='i',epsg=4326)
plt.figure()
m.wmsimage(serverurl,xpixels=500,verbose=True,
      layers=['Temperature_height_above_ground'],
      styles=['boxfill/rainbow'],
      time=datetime.utcnow().strftime('%Y-%m-%dT00:00:00.000Z'),
      elevation='2',
      colorscalerange='271.2,308',numcolorbands='20',logscale=False)
m.drawcoastlines(linewidth=0.25)
parallels = np.arange(20,36,2.)
a=m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
meridians = np.arange(-120,-100,2.)
b=m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

epsg = 32661
m = Basemap(epsg=epsg, resolution='l',width=20000.e3,height=20000.e3)
plt.figure()
m.wmsimage(serverurl,xpixels=500,
      layers=['Temperature_height_above_ground'],
      styles=['boxfill/rainbow'],
      time=datetime.utcnow().strftime('%Y-%m-%dT00:00:00.000Z'),
      elevation='2',
      colorscalerange='271.2,308',numcolorbands='20',logscale=False)
m.drawcoastlines(linewidth=0.5)

plt.show()
