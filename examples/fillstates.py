import pylab as p
import matplotlib.numerix as nx
from matplotlib.toolkits.basemap import Basemap as Basemap
from matplotlib.collections import LineCollection
from matplotlib.colors import rgb2hex
import random

# Lambert Conformal map of lower 48 states.
m = Basemap(llcrnrlon=-119,llcrnrlat=22,urcrnrlon=-64,urcrnrlat=49,
            projection='lcc',lat_1=33,lat_2=45,lon_0=-95)
fig=m.createfigure()
# draw state boundaries.
# data from U.S Census Bureau
# http://www.census.gov/geo/www/cob/st1990.html
shp_info = m.readshapefile('st99_d90','states',drawbounds=True)
print shp_info
# choose a color for each state (randomly).
colors={}
statenames=[]
print m.states_info[0].keys()
for shapedict in m.states_info:
    statename = shapedict['NAME']
    colors[statename] = (random.uniform(0,1),random.uniform(0,1),random.uniform(0,1))
    statenames.append(statename)
# cycle through state names, color each one.
for nshape,seg in enumerate(m.states):
    xx,yy = zip(*seg)
    color = rgb2hex(colors[statenames[nshape]]) 
    p.fill(xx,yy,color,edgecolor=color)
# draw meridians and parallels.
m.drawparallels(nx.arange(25,65,20),labels=[1,0,0,0])
m.drawmeridians(nx.arange(-120,-40,20),labels=[0,0,0,1])
p.title('Filling State Polygons')
p.show()
