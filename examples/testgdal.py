"""
example showing how to plot a USGS DEM file using 
gdal (http://gdal.maptools.org).

Data files must be downloaded manually from USGS:
http://edcftp.cr.usgs.gov/pub/data/DEM/250/D/denver-w.gz
http://edcftp.cr.usgs.gov/pub/data/nationalatlas/countyp020.tar.gz
"""
import gdal
from matplotlib.toolkits.basemap import Basemap
from gdalconst import *
import pylab as p

# download from 
# http://edcftp.cr.usgs.gov/pub/data/DEM/250/D/denver-w.gz
gd = gdal.Open('denver-w')
# get data from DEM file
array = gd.ReadAsArray()
# get lat/lon coordinates from DEM file.
coords = gd.GetGeoTransform()
llcrnrlon = coords[0]
urcrnrlon = llcrnrlon+(array.shape[1]-1)*coords[1]
urcrnrlat = coords[3]
llcrnrlat = urcrnrlat+(array.shape[0]-1)*coords[5]
# create Basemap instance.
m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,projection='cyl')
# create a figure, add an axes
# (leaving room for a colorbar).
fig = p.figure()
ax = fig.add_axes([0.1,0.1,0.75,0.75])
# plot image from DEM over map.
im = m.imshow(array,origin='upper')
# make a colorbar.
cax = p.axes([0.875, 0.1, 0.05, 0.75]) # setup colorbar axes.
p.colorbar(cax=cax) # draw colorbar
p.axes(ax)  # make the original axes current again
# draw meridians and parallels.
m.drawmeridians(p.linspace(llcrnrlon+0.1,urcrnrlon-0.1,5),labels=[0,0,0,1],fmt='%4.2f')
m.drawparallels(p.linspace(llcrnrlat+0.1,urcrnrlat-0.1,5),labels=[1,0,0,0],fmt='%4.2f')
# plot county boundaries from
# http://edcftp.cr.usgs.gov/pub/data/nationalatlas/countyp020.tar.gz
shp_info = m.readshapefile('countyp020','counties',drawbounds=True,linewidth=1.0)
# plot some cities.
lons = [-105.22,-105.513,-105.316,-105.47]; lats = [39.76,39.801,39.633,39.41]
names =  ['Golden','Central City','Evergreen','Bailey']
x,y = m(lons,lats)
m.plot(x,y,'ko')
for name,xx,yy in zip(names,x,y):
    p.text(xx+0.01,yy+0.01,name)
p.title(gd.GetDescription()+' USGS DEM with county boundaries')
p.show()
