from mpl_toolkits.basemap import Basemap, pyproj
import matplotlib.pyplot as plt
import urllib2
import numpy as np

class Basemap2(Basemap):
    def wmsmap(self,layer='World_Physical_Map',xpixels=800,**kwargs):
        from matplotlib.image import imread
        if not hasattr(self,'epsg'):
            raise ValueError('no epsg code')
        p = pyproj.Proj(init="epsg:%s" % self.epsg, preserve_units=True)
        x1,y1 = p(self.llcrnrlon,self.llcrnrlat)
        x2,y2 = p(self.urcrnrlon,self.urcrnrlat)
        ypixels = int(self.aspect*xpixels)
        basemap_url =\
"http://server.arcgisonline.com/ArcGIS/rest/services/%s/MapServer/export?\
bbox=%d,%d,%d,%d&\
bboxSR=%s&\
imageSR=%s&\
size=%s,%s&\
dpi=128&\
format=png32&\
f=image" % (layer,x1,y1,x2,y2,self.epsg,self.epsg,xpixels,ypixels)
        return m.imshow(imread(urllib2.urlopen(basemap_url)),origin='upper')

plt.figure(figsize=(10,5))
epsg = 2263; width=1000.e3; height = 500.e3
m=Basemap2(epsg=epsg,resolution='h',width=width,height=height)
# default
#m.wmsmap()
# specify a different layer, pixel resolution.
m.wmsmap(layer='NatGeo_World_Map',xpixels=1200)
m.drawmeridians(np.arange(-180,180,2),labels=[0,0,0,1])
m.drawparallels(np.arange(0,80,2),labels=[1,0,0,0])
m.drawcoastlines(linewidth=0.25)
plt.title('test WMS map background')

plt.show()
