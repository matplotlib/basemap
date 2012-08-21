import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def draw_map_background(m, ax):
    ax.set_axis_bgcolor('#729FCF')
    m.fillcontinents(color='#FAFAFA', ax=ax, zorder=0)
    m.drawcounties(ax=ax)
    m.drawstates(ax=ax)
    m.drawcountries(ax=ax)
    m.drawcoastlines(ax=ax)
KM = 1000.
clat = 39.3
clon = -94.7333
wid = 5500 * KM
hgt = 3500 * KM
m = Basemap(width=wid, height=hgt, rsphere=(6378137.00,6356752.3142),
            resolution='i', area_thresh=2500., projection='lcc',
            lat_1=38.5, lat_2=38.5, lat_0=clat, lon_0=clon)
fig = plt.figure()
ax = fig.add_subplot(111)
draw_map_background(m, ax)
plt.show()