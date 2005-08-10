import pylab as p
import matplotlib.numerix as nx
from matplotlib.toolkits.basemap import Basemap as Basemap_base
from matplotlib.collections import LineCollection
from matplotlib.colors import rgb2hex
import random
from shapelib import ShapeFile

# add a "drawclimdivs" method for Basemap, which draws climate
# division boundaries.

# requires pyshapelib from Thuban (http://thuban.intevation.org/).
# cd to libraries/pyshapelib in Thuban source distribution, run
# 'python setup.py install'.

class Basemap(Basemap_base):

    def drawclimdivs(self,linewidth=0.5,color='k',antialiased=1,ax=None):
        """read in NCDC climate division boundaries and plot on map"""
        # open shapefile, read vertices for each object, convert
        # to map projection coordinates.
        shp = ShapeFile('divisions')
        self.climdivsegs = []
        for npoly in range(shp.info()[0]):
            shp_object = shp.read_object(npoly)
            verts = shp_object.vertices()[0]
            lons, lats = zip(*verts)
            x, y = self(lons,lats)
            self.climdivsegs.append(zip(x,y))
        # get current axes instance (if none specified).
        if ax is None and self.ax is None:
            try: 
                ax = p.gca()
            except:
                import pylab as p
                ax = p.gca()
        elif ax is None and self.ax is not None:
            ax = self.ax
        # make LineCollections for each polygon.
        lines = LineCollection(self.climdivsegs,antialiaseds=(antialiased,))
        lines.set_color(color)
        lines.set_linewidth(linewidth)
        ax.add_collection(lines)
        # make sure axis ticks are turned off
        if self.noticks == True:
            ax.set_xticks([])
            ax.set_yticks([])
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)

# Lambert Conformal map of lower 48 states.
m = Basemap(llcrnrlon=-119,llcrnrlat=22,urcrnrlon=-64,urcrnrlat=49,
            projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
fig=p.figure(figsize=(8,m.aspect*8))
fig.add_axes([0.1,0.1,0.8,0.8])
# draw climate division boundaries.
m.drawclimdivs(linewidth=1.0)
# cycle through climate divisions, fill each with a random color.
for seg in m.climdivsegs:
    xx,yy = zip(*seg)
    color = (random.uniform(0,1),random.uniform(0,1),random.uniform(0,1))
    p.fill(xx,yy,rgb2hex(color))
# draw meridians and parallels.
m.drawparallels(nx.arange(25,65,20),labels=[1,0,0,0])
m.drawmeridians(nx.arange(-120,-40,20),labels=[0,0,0,1])
p.title('NCDC Climate Divisions')
p.show()
