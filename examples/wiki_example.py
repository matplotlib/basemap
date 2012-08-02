from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
# set up orthographic map projection with
# perspective of satellite looking down at 50N, 100W.
# use low resolution coastlines.
bmap = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
# draw coastlines, country boundaries, fill continents.
bmap.drawcoastlines(linewidth=0.25)
bmap.drawcountries(linewidth=0.25)
bmap.fillcontinents(color='coral',lake_color='aqua')
# draw the edge of the map projection region (the projection limb)
bmap.drawmapboundary(fill_color='aqua')
# draw lat/lon grid lines every 30 degrees.
bmap.drawmeridians(np.arange(0,360,30))
bmap.drawparallels(np.arange(-90,90,30))
# lat/lon coordinates of five cities.
lats=[40.02,32.73,38.55,48.25,17.29]
lons=[-105.16,-117.16,-77.00,-114.21,-88.10]
cities=['Boulder, CO','San Diego, CA',
        'Washington, DC','Whitefish, MT','Belize City, Belize']
# compute the native map projection coordinates for cities.
xc,yc = bmap(lons,lats)
# plot filled circles at the locations of the cities.
bmap.plot(xc,yc,'bo')
# plot the names of those five cities.
for name,xpt,ypt in zip(cities,xc,yc):
    plt.text(xpt+50000,ypt+50000,name,fontsize=9)
# make up some data on a regular lat/lon grid.
nlats = 73; nlons = 145; delta = 2.*np.pi/(nlons-1)
lats = (0.5*np.pi-delta*np.indices((nlats,nlons))[0,:,:])
lons = (delta*np.indices((nlats,nlons))[1,:,:])
wave = 0.75*(np.sin(2.*lats)**8*np.cos(4.*lons))
mean = 0.5*np.cos(2.*lats)*((np.sin(2.*lats))**2 + 2.)
# compute native map projection coordinates of lat/lon grid.
x, y = bmap(lons*180./np.pi, lats*180./np.pi)
# contour data over the map.
cs = bmap.contour(x,y,wave+mean,15,linewidths=1.5)
plt.title('filled continent background')

# as above, but use land-sea mask image as map background.
fig = plt.figure()
bmap.drawmapboundary()
bmap.drawmeridians(np.arange(0,360,30))
bmap.drawparallels(np.arange(-90,90,30))
# plot filled circles at the locations of the cities.
bmap.plot(xc,yc,'wo')
# plot the names of five cities.
for name,xpt,ypt in zip(cities,xc,yc):
    plt.text(xpt+50000,ypt+50000,name,fontsize=9,color='w')
# contour data over the map.
cs = bmap.contour(x,y,wave+mean,15,linewidths=1.5)
plt.title('land-sea mask background')
bmap.drawlsmask(ocean_color='aqua',land_color='coral')

# as above, but use blue marble image as map background.
fig = plt.figure()
bmap.drawmapboundary()
bmap.drawmeridians(np.arange(0,360,30))
bmap.drawparallels(np.arange(-90,90,30))
# plot filled circles at the locations of the cities.
bmap.plot(xc,yc,'wo')
# plot the names of five cities.
for name,xpt,ypt in zip(cities,xc,yc):
    plt.text(xpt+50000,ypt+50000,name,fontsize=9,color='w')
# contour data over the map.
cs = bmap.contour(x,y,wave+mean,15,linewidths=1.5)
plt.title('blue marble background')
bmap.bluemarble()

# as above, but use shaded relief image as map background.
fig = plt.figure()
bmap.drawmapboundary()
bmap.drawmeridians(np.arange(0,360,30))
bmap.drawparallels(np.arange(-90,90,30))
# plot filled circles at the locations of the cities.
bmap.plot(xc,yc,'wo')
# plot the names of five cities.
for name,xpt,ypt in zip(cities,xc,yc):
    plt.text(xpt+50000,ypt+50000,name,fontsize=9,color='w')
# contour data over the map.
cs = bmap.contour(x,y,wave+mean,15,linewidths=1.5)
plt.title('shaded relief background')
bmap.shadedrelief()

# as above, but use etopo image as map background.
fig = plt.figure()
bmap.drawmapboundary()
bmap.drawmeridians(np.arange(0,360,30))
bmap.drawparallels(np.arange(-90,90,30))
# plot filled circles at the locations of the cities.
bmap.plot(xc,yc,'wo')
# plot the names of five cities.
for name,xpt,ypt in zip(cities,xc,yc):
    plt.text(xpt+50000,ypt+50000,name,fontsize=9,color='w')
# contour data over the map.
cs = bmap.contour(x,y,wave+mean,15,linewidths=1.5)
plt.title('etopo background')
bmap.etopo()

# as above, but use etopo image as map background overlaid with
# land-sea mask image where land areas are transparent (so etopo
# image shows through over land).
fig = plt.figure()
bmap.drawmapboundary()
bmap.drawmeridians(np.arange(0,360,30))
bmap.drawparallels(np.arange(-90,90,30))
# plot filled circles at the locations of the cities.
bmap.plot(xc,yc,'wo')
# plot the names of five cities.
for name,xpt,ypt in zip(cities,xc,yc):
    plt.text(xpt+50000,ypt+50000,name,fontsize=9,color='w')
# contour data over the map.
cs = bmap.contour(x,y,wave+mean,15,linewidths=1.5)
plt.title('etopo background with oceans masked')
bmap.etopo()
bmap.drawlsmask(ocean_color='DarkBlue',land_color=(255,255,255,1))

plt.show()
