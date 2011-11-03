"""
example showing how to plot data from a DEM file and an ESRI shape file using 
gdal (http://pypi.python.org/pypi/GDAL).
"""
from osgeo import gdal, ogr
from mpl_toolkits.basemap import Basemap, cm
import numpy as np
import matplotlib.pyplot as plt
from numpy import ma

# read 2.5 minute U.S. DEM file using gdal.
# (http://www.prism.oregonstate.edu/docs/meta/dem_25m.htm)
gd = gdal.Open('us_25m.dem')
array = gd.ReadAsArray()
# get lat/lon coordinates from DEM file.
coords = gd.GetGeoTransform()
nlons = array.shape[1]; nlats = array.shape[0]
delon = coords[1]
delat = coords[5]
lons = coords[0] + delon*np.arange(nlons)
lats = coords[3] + delat*np.arange(nlats)[::-1] # reverse lats
# setup figure.
fig = plt.figure(figsize=(11,6))
# setup basemap instance.
m = Basemap(llcrnrlon=-119,llcrnrlat=22,urcrnrlon=-64,urcrnrlat=49,
            projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
# create masked array, reversing data in latitude direction
# (so that data is oriented in increasing latitude, as transform_scalar requires).
topoin = ma.masked_values(array[::-1,:],-999.)
# transform DEM data to a 4 km native projection grid
nx = int((m.xmax-m.xmin)/4000.)+1; ny = int((m.ymax-m.ymin)/4000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny,masked=True)
# plot DEM image on map.
im = m.imshow(topodat,cmap=cm.GMT_haxby_r)
# draw meridians and parallels.
m.drawparallels(np.arange(20,71,10),labels=[1,0,0,0])
m.drawmeridians(np.arange(-120,-40,10),labels=[0,0,0,1])
# plot state boundaries from shapefile using ogr.
g = ogr.Open ("st99_d00.shp")
L = g.GetLayer(0) # data is in 1st layer.
for feat in L: # iterate over features in layer
    geo = feat.GetGeometryRef()
    # iterate over geometries. 
    for count in range(geo.GetGeometryCount()):
        geom = geo.GetGeometryRef(count)
        if not geom.GetGeometryCount(): # just one geometry.
            # get lon,lat points
            lons = [geom.GetX(i) for i in range(geom.GetPointCount())]
            lats = [geom.GetY(i) for i in range(geom.GetPointCount())]
            # convert to map projection coords.
            x, y = m(lons,lats)
            # plot on map.
            m.plot(x,y,'k')
        else: # iterate over nested geometries.
            for cnt in range( geom.GetGeometryCount()):
                g = geom.GetGeometryRef( cnt )
                lons = [g.GetX(i) for i in range(g.GetPointCount())]
                lats = [g.GetY(i) for i in range(g.GetPointCount())]
                x, y = m(lons,lats)
                m.plot(x,y,'k')
# draw colorbar.
m.colorbar(im)
plt.title(gd.GetDescription()+' with state boundaries from '+g.GetName(),y=1.05)
plt.show()
