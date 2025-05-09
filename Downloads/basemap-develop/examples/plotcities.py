from __future__ import (absolute_import, division, print_function)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap as Basemap


def prctile_rank(x, p): return np.searchsorted(
    np.percentile(x, np.atleast_1d(p)), x)

# cities colored by population rank.


m = Basemap()
shp_info = m.readshapefile('cities', 'cities')
x, y, pop = [], [], []
for item, (x_i, y_i) in zip(m.cities_info, m.cities):
    population = item['POPULATION']
    if population >= 0:
        pop.append(population)
        x.append(x_i)
        y.append(y_i)
popranks = prctile_rank(pop, np.linspace(0, 100, 101))
colors = []
for rank in popranks:
    colors.append(plt.cm.jet(float(rank) / 100.))
m.drawcoastlines()
m.fillcontinents()
m.scatter(x, y, 25, colors, marker='o', edgecolors='none', zorder=10)
plt.title('City Locations colored by Population Rank')
plt.show()
