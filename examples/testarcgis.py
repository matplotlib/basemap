from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

# test arcgisimage method for retrieving images from web map servers.

plt.figure(figsize=(8,8))
epsg=4326
m=Basemap(projection='cyl',llcrnrlon=-90,llcrnrlat=30,urcrnrlon=-60,urcrnrlat=60,resolution='i')
# default 'blue marble' image.
m.arcgisimage(verbose=True)
m.drawmeridians(np.arange(-180,180,10),labels=[0,0,0,1],color='y')
m.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0],color='y')
m.drawcoastlines(linewidth=0.5,color='0.5')
plt.title('test WMS map background EPSG=%s'%epsg)

plt.figure(figsize=(8,8))
epsg = 3413; width = 18000.e3; height = 18000.e3
m=Basemap(epsg=epsg,resolution='l',width=width,height=height)
# default 'blue marble' image.
m.arcgisimage()
m.drawmeridians(np.arange(-180,180,30),labels=[0,0,0,1],color='y')
m.drawparallels(np.arange(0,80,10),labels=[1,0,0,0],color='y')
m.drawcoastlines(linewidth=0.5,color='0.5')
plt.title('test WMS map background EPSG=%s' % epsg)

plt.figure(figsize=(10,6.6666))
epsg = 2263; width=600.e3; height = 400.e3
m=Basemap(epsg=epsg,resolution='h',width=width,height=height)
# specify a different server.
m.arcgisimage(server='http://maps.ngdc.noaa.gov',service='etopo1',verbose=True)
m.drawmeridians(np.arange(-180,180,2),labels=[0,0,0,1])
m.drawparallels(np.arange(0,80,1),labels=[1,0,0,0])
m.drawcoastlines(linewidth=0.25)
plt.title('test WMS map background EPSG=%s' % epsg)

plt.figure(figsize=(6,8))
epsg = 3086; lon1 = -85; lat1 = 24; lon2 = -79; lat2 =32
m=Basemap(epsg=epsg,resolution='i',llcrnrlon=lon1,llcrnrlat=lat1,\
           urcrnrlon=lon2,urcrnrlat=lat2)
# specify a different service
m.arcgisimage(service='NatGeo_World_Map',verbose=True,xpixels=600,interpolation='bicubic')
m.drawmeridians(np.arange(-180,180,2),labels=[0,0,0,1])
m.drawparallels(np.arange(0,80,1),labels=[1,0,0,0])
m.drawcoastlines(linewidth=0.25)
plt.title('test WMS map background EPSG=%s' % epsg)

plt.show()
