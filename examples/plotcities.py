import pylab as p
from matplotlib.mlab import prctile_rank
import numpy
from mpl_toolkits.basemap import Basemap as Basemap

# cities colored by population rank.

m = Basemap()
shp_info = m.readshapefile('cities','cities')
x, y = zip(*m.cities)
pop = []
for item in m.cities_info:
    population = item['POPULATION']
    if population < 0: continue # population missing
    pop.append(population)
popranks = prctile_rank(pop,100)
colors = []
for rank in popranks:
    colors.append(p.cm.jet(float(rank)/100.))
m.drawcoastlines()
m.fillcontinents()
m.scatter(x,y,25,colors,marker='o',faceted=False,zorder=10)
p.title('City Locations colored by Population Rank')
p.show()
