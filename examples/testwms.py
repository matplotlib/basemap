from mpl_toolkits.basemap import Basemap, pyproj
import matplotlib.pyplot as plt
import urllib2
import numpy as np

# blue marble image on cylindrical projection.
plt.figure()
lon1 = -10; lon2 = 30; lat1 = 30; lat2 = 60
basemap_url =\
"http://server.arcgisonline.com/ArcGIS/rest/services/ESRI_Imagery_World_2D/MapServer/export?\
bbox=%d,%d,%d,%d&\
bboxSR=4326&\
imageSR=4326&\
size=800,600&\
dpi=128&\
format=png32&\
f=image" % (lon1, lat1, lon2, lat2)
m = Basemap(projection='cyl',llcrnrlon=lon1,urcrnrlon=lon2,\
            llcrnrlat=lat1,urcrnrlat=lat2,resolution='i')
m.drawcoastlines()
m.drawmeridians(np.arange(-5,26,5),labels=[0,0,0,1],color='y')
m.drawparallels(np.arange(25,56,5),labels=[1,0,0,0],color='y')
m.imshow(plt.imread(urllib2.urlopen(basemap_url)),origin='upper')

# blue marble image on south polar stereographic projection.
width = 12000.e3
plt.figure()
basemap_url =\
"http://server.arcgisonline.com/ArcGIS/rest/services/ESRI_Imagery_World_2D/MapServer/export?\
bbox=%d,%d,%d,%d&\
bboxSR=3412&\
imageSR=3412&\
size=800,800&\
dpi=128&\
format=png32&\
f=image" % (-width/2,-width/2,width/2,width/2)
m =\
Basemap(projection='stere',resolution='i',lon_0=0,lat_0=-90,lat_ts=-70,\
        width=width,height=width,rsphere=(6378273,6356889.449))
m.imshow(plt.imread(urllib2.urlopen(basemap_url)),origin='upper')
m.drawmeridians(np.arange(-180,180,30),labels=[0,0,0,1],color='y')
m.drawparallels(np.arange(-80,-0,10),labels=[1,0,0,0],color='y')
m.drawcoastlines()

# blue marble image on north american lambert conformal projection.
width = 8000.e3; height=8000.e3
plt.figure()
basemap_url =\
"http://server.arcgisonline.com/ArcGIS/rest/services/ESRI_Imagery_World_2D/MapServer/export?\
bbox=%d,%d,%d,%d&\
bboxSR=102009&\
imageSR=102009&\
size=800,800&\
dpi=128&\
format=png32&\
f=image" % (-width/2,-height/2,width/2,height/2)
m =\
Basemap(projection='lcc',resolution='i',lat_1=20,lat_2=60,\
        lat_0=40.,lon_0=-96,\
        width=width,height=height,rsphere=(6378137,6356752.3141))
m.imshow(plt.imread(urllib2.urlopen(basemap_url)),origin='upper')
m.drawmeridians(np.arange(-180,180,10),labels=[0,0,0,1],color='y')
m.drawparallels(np.arange(0,80,10),labels=[1,0,0,0],color='y')
m.drawcoastlines()

# world physical map on europe lambert equal area projection.
lon1 = -8.9067; lat1 = 33.2307; lon2 = 72.9617; lat2 = 58.9174
lat_0=52; lon_0=10
plt.figure()
m =\
Basemap(projection='laea',resolution='i',\
        llcrnrlat=lat1,llcrnrlon=lon1,urcrnrlon=lon2,urcrnrlat=lat2,\
        lat_0=lat_0,lon_0=lon_0,\
        rsphere=(6378137,6356752.3141))
p = pyproj.Proj(
    proj='laea',lat_0=lat_0,lon_0=lon_0,x_0=4321000,y_0=3210000,ellps='GRS80')
x1,y1 = p(lon1,lat1)
x2,y2 = p(lon2,lat2)
basemap_url =\
"http://server.arcgisonline.com/ArcGIS/rest/services/World_Physical_Map/MapServer/export?\
bbox=%d,%d,%d,%d&\
bboxSR=3035&\
imageSR=3035&\
size=898,768&\
dpi=128&\
format=png32&\
f=image" % (x1, y1, x2, y2)
m.imshow(plt.imread(urllib2.urlopen(basemap_url)),origin='upper')
m.drawmeridians(np.arange(-180,180,10),labels=[0,0,0,1])
m.drawparallels(np.arange(0,80,10),labels=[1,0,0,0])
m.drawcoastlines(linewidth=0.25)

plt.show()
