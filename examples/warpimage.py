from __future__ import print_function
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

# illustrate use of warpimage method to display an image background
# on the map projection region.  Default background is the 'blue
# marble' image from NASA (http://visibleearth.nasa.gov).

# create new figure
fig=plt.figure()
# define orthographic projection centered on North America.
m = Basemap(projection='ortho',lat_0=40,lon_0=-100,resolution='l')
# display a non-default image.
m.warpimage(image='earth_lights_lrg.jpg')
# draw coastlines.
m.drawcoastlines(linewidth=0.5,color='0.5')
# draw lat/lon grid lines every 30 degrees.
m.drawmeridians(np.arange(0,360,30),color='0.5')
m.drawparallels(np.arange(-90,90,30),color='0.5')
plt.title("Lights at Night image warped from 'cyl' to 'ortho' projection",fontsize=12)
print('warp to orthographic map ...')

# create new figure
fig=plt.figure()
# define projection centered on North America.
m = Basemap(projection='mbtfpq',lon_0=-100,resolution='l')
m.bluemarble(scale=0.5)
# draw coastlines.
m.drawcoastlines(linewidth=0.5,color='0.5')
# draw lat/lon grid lines every 30 degrees.
m.drawmeridians(np.arange(0,360,60),color='0.5')
m.drawparallels(np.arange(-90,90,30),color='0.5')
plt.title("Blue Marble image warped from 'cyl' to 'mbtfpq' projection",fontsize=12)
print('warp to McBryde-Thomas Flat-Polar Quartic map ...')

# create new figure
fig=plt.figure()
# define projection centered on North America.
m = Basemap(projection='hammer',lon_0=-100,resolution='l')
m.bluemarble(scale=0.5)
# draw coastlines.
m.drawcoastlines(linewidth=0.5,color='0.5')
# draw lat/lon grid lines every 30 degrees.
m.drawmeridians(np.arange(0,360,60),color='0.5')
m.drawparallels(np.arange(-90,90,30),color='0.5')
plt.title("Blue Marble image warped from 'cyl' to 'hammer' projection",fontsize=12)
print('warp to Hammer map ...')

# create new figure
fig=plt.figure()
# define cylindrical equidistant projection.
m = Basemap(projection='cyl',llcrnrlon=-180,llcrnrlat=-90,urcrnrlon=180,urcrnrlat=90,resolution='l')
# plot (unwarped) rgba image.
im = m.bluemarble(scale=0.5)
# draw coastlines.
m.drawcoastlines(linewidth=0.5,color='0.5')
# draw lat/lon grid lines.
m.drawmeridians(np.arange(-180,180,60),labels=[0,0,0,1],color='0.5')
m.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0],color='0.5')
plt.title("Blue Marble image - native 'cyl' projection",fontsize=12)
print('plot native cylindrical map (no warping needed) ...')

# create new figure
fig=plt.figure()
# define cylindrical equidistant projection.
m = Basemap(projection='cyl',llcrnrlon=0,llcrnrlat=-60,urcrnrlon=360,urcrnrlat=60,resolution='l')
# plot (unwarped) rgba image.
im = m.bluemarble(scale=0.5)
# draw coastlines.
m.drawcoastlines(linewidth=0.5,color='0.5')
# draw lat/lon grid lines.
m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1],color='0.5')
m.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0],color='0.5')
plt.title("Blue Marble image - non-native 'cyl' projection",fontsize=12)
print('plot non-native cylindrical map (warping needed) ...')

# create new figure
fig=plt.figure()
# define orthographic projection centered on Europe.
m = Basemap(projection='ortho',lat_0=40,lon_0=40,resolution='l')
# plot a gray-scale image specified from a URL.
im = m.warpimage("http://earthobservatory.nasa.gov/Features/BlueMarble/Images/gebco_bathy.5400x2700.jpg")
# draw coastlines.
m.drawcoastlines(linewidth=0.5,color='0.5')
# draw lat/lon grid lines every 30 degrees.
m.drawmeridians(np.arange(0,360,30),color='0.5')
m.drawparallels(np.arange(-90,90,30),color='0.5')
plt.title("Blue Marble image warped from 'cyl' to 'ortho' projection",fontsize=12)
print('warp to orthographic map ...')

# create new figure
fig=plt.figure()
# define Lambert Conformal basemap for North America.
m = Basemap(llcrnrlon=-145.5,llcrnrlat=1.,urcrnrlon=-2.566,urcrnrlat=46.352,\
            rsphere=(6378137.00,6356752.3142),lat_1=50.,lon_0=-107.,\
            resolution='i',area_thresh=1000.,projection='lcc')
im = m.bluemarble(scale=0.5)
# draw coastlines.
m.drawcoastlines(linewidth=0.5,color='0.5')
# draw parallels and meridians.
# label on left, right and bottom of map.
parallels = np.arange(0.,80,20.)
m.drawparallels(parallels,labels=[1,1,0,1],color='0.5')
meridians = np.arange(10.,360.,30.)
m.drawmeridians(meridians,labels=[1,1,0,1],color='0.5')
plt.title("Blue Marble image warped from 'cyl' to 'lcc' projection",fontsize=12)
print('warp to lambert conformal map ...')

# create new figure
fig=plt.figure()
# define oblique mercator map.
m = Basemap(height=24000000,width=12000000,
            resolution=None,projection='omerc',\
            lon_0=-100,lat_0=15,lon_2=-120,lat_2=65,lon_1=-50,lat_1=-55)
# plot warped rgba image.
im = m.bluemarble(scale=0.5)
# draw lat/lon grid lines every 20 degrees.
m.drawmeridians(np.arange(0,360,20),color='0.5')
m.drawparallels(np.arange(-80,81,20),color='0.5')
plt.title("Blue Marble image warped from 'cyl' to 'omerc' projection",fontsize=12)
print('warp to oblique mercator map ...')

plt.show()
