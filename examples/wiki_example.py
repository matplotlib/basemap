from matplotlib.toolkits.basemap import Basemap
import pylab as p
# set up orthographic map projection with
# perspective of satellite looking down at 50N, 100W.
# use low resolution coastlines.
map = Basemap(projection='ortho',lat_0=50,lon_0=-100,resolution='l')
# draw coastlines, country boundaries, fill continents.
map.drawcoastlines()
map.drawcountries()
map.fillcontinents(color='coral')
# draw the edge of the map projection region (the projection limb)
map.drawmapboundary()
# draw lat/lon grid lines every 30 degrees.
map.drawmeridians(p.arange(0,360,30))
map.drawparallels(p.arange(-90,90,30))
# lat/lon coordinates of five cities.
lats=[40.02,34.00,38.55,48.25,17.29]
lons=[-105.16,-119.40,-77.00,-114.21,-88.10]
cities=['Boulder, CO','Santa Cruz, CA',
        'Washington, DC','Whitefish, MT','Belize City, Belize']
# compute the native map projection coordinates for cities.
x,y = map(lons,lats)
# plot filled circles at the locations of the cities.
map.plot(x,y,'bo')
# plot the names of those five cities.
for name,xpt,ypt in zip(cities,x,y):
    p.text(xpt+50000,ypt+50000,name)
# make up some data on a regular lat/lon grid.
nlats = 73; nlons = 145; delta = 2.*p.pi/(nlons-1)
lats = (0.5*p.pi-delta*p.indices((nlats,nlons))[0,:,:])
lons = (delta*p.indices((nlats,nlons))[1,:,:])
wave = 0.75*(p.sin(2.*lats)**8*p.cos(4.*lons))
mean = 0.5*p.cos(2.*lats)*((p.sin(2.*lats))**2 + 2.)
# compute native map projection coordinates of lat/lon grid.
x, y = map(lons*180./p.pi, lats*180./p.pi)
# contour data over the map.
cs = map.contour(x,y,wave+mean,15,linewidths=1.5)
p.show()

