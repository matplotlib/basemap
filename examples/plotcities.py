import pylab as p
import numpy
from matplotlib.toolkits.basemap import Basemap as Basemap

# cities colored by population rank.
m = Basemap()
shp_info = m.readshapefile('cities','cities')
x, y = zip(*m.cities)
pop = []
for item in m.cities_info:
    pop.append(int(item['POPULATION']))
pop = numpy.array(pop)
poprank = numpy.argsort(pop)
m.drawcoastlines()
m.fillcontinents()
m.scatter(x,y,25,poprank,cmap=p.cm.jet_r,marker='o',faceted=False,zorder=10)
p.title('City Locations colored by Population Rank')
p.show()
