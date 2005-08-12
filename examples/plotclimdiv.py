import pylab as p
import matplotlib.numerix as nx
from matplotlib.toolkits.basemap import Basemap as Basemap
from matplotlib.collections import LineCollection
from matplotlib.colors import rgb2hex
import random

# requires pyshapelib from Thuban (http://thuban.intevation.org/).
# cd to libraries/pyshapelib in Thuban source distribution, run
# 'python setup.py install'.

# Lambert Conformal map of lower 48 states.
m = Basemap(llcrnrlon=-119,llcrnrlat=22,urcrnrlon=-64,urcrnrlat=49,
            projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
fig=p.figure(figsize=(8,m.aspect*8))
fig.add_axes([0.1,0.1,0.8,0.8])
# draw climate division boundaries.
shp_info = m.readshapefile('divisions','climdivs',drawbounds=True)
print shp_info
# make sure the shapefile has polygons (and not just lines).
if shp_info[1] != 5:
    print 'warning: shapefile does not contain polygons'
# choose a color for each climate division (randomly).
colors={}
divnames=[]
print m.climdivs_info[0].keys()
for shapedict in m.climdivs_info:
    divname = shapedict['ST']+repr(shapedict['DIV'])
    colors[divname] = (random.uniform(0,1),random.uniform(0,1),random.uniform(0,1))
    divnames.append(divname)
# cycle through climate divnames, color each one.
for nshape,seg in enumerate(m.climdivs):
    xx,yy = zip(*seg)
    color = rgb2hex(colors[divnames[nshape]]) 
    p.fill(xx,yy,color,edgecolor=color)
# draw meridians and parallels.
m.drawparallels(nx.arange(25,65,20),labels=[1,0,0,0])
m.drawmeridians(nx.arange(-120,-40,20),labels=[0,0,0,1])
p.title('NCDC Climate Divisions')
p.show()
