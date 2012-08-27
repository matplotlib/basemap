from mpl_toolkits.basemap import Basemap, pyproj
import matplotlib.pyplot as plt
from matplotlib.image import imread
from matplotlib.cbook import dedent
import urllib2
import numpy as np

class Basemap2(Basemap):
    def wmsmap(self,server='http://server.arcgisonline.com/ArcGIS',layer='World_Physical_Map',\
               xpixels=800,dpi=96,format='png',verbose=False,**kwargs):
        # constructs a URL using the ArcGIS Server REST API, retrieves
        # an image and displays it on the map. In order to use this method,
        # the Basemap instance must be created using the epsg keyword to
        # specify the map projection.
        if not hasattr(self,'epsg'):
            msg = dedent("""
            Basemap instance must be creating using an EPSG code
            (http://spatialreference.org) in order to use the wmsmap method""")
            raise ValueError('no epsg code')
        p = pyproj.Proj(init="epsg:%s" % self.epsg, preserve_units=True)
        x1,y1 = p(self.llcrnrlon,self.llcrnrlat)
        x2,y2 = p(self.urcrnrlon,self.urcrnrlat)
        ypixels = int(self.aspect*xpixels)
        basemap_url = \
"%s/rest/services/%s/MapServer/export?\
bbox=%d,%d,%d,%d&\
bboxSR=%s&\
imageSR=%s&\
size=%s,%s&\
dpi=%s&\
format=%s&\
f=image" %\
(server,layer,x1,y1,x2,y2,self.epsg,self.epsg,xpixels,ypixels,dpi,format)
        if verbose: print basemap_url
        return m.imshow(imread(urllib2.urlopen(basemap_url)),origin='upper')

plt.figure(figsize=(10,5))
epsg = 2263; width=600.e3; height = 400.e3
m=Basemap2(epsg=epsg,resolution='h',width=width,height=height)
# default
#m.wmsmap()
# specify a different layer, pixel resolution.
#m.wmsmap(layer='NatGeo_World_Map',xpixels=1200)
# specify a different server.
m.wmsmap(server='http://maps.ngdc.noaa.gov',layer='etopo1',verbose=True)
m.drawmeridians(np.arange(-180,180,2),labels=[0,0,0,1])
m.drawparallels(np.arange(0,80,1),labels=[1,0,0,0])
m.drawcoastlines(linewidth=0.25)
plt.title('test WMS map background')

plt.show()
