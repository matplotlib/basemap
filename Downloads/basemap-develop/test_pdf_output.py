import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

fig = plt.figure(figsize=(10, 6), dpi=100)
m = Basemap(projection='cyl',
            llcrnrlat=40, urcrnrlat=47,
            llcrnrlon=-72, urcrnrlon=-67,
            resolution='l',
            fix_aspect=False)

m.drawcoastlines()

lat_ticks = np.arange(41, 47)
lon_ticks = np.arange(-72, -67)

m.drawparallels(lat_ticks, labels=[1, 0, 0, 0], linewidth=0)
m.drawmeridians(lon_ticks, labels=[0, 0, 0, 1], linewidth=0)

plt.savefig("test_output.pdf")
print("PDF saved as test_output.pdf")
