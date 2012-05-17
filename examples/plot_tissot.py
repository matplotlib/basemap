import numpy as np
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap 
from mpl_toolkits.basemap import __version__ as basemap_version

# Tissot's Indicatrix (http://en.wikipedia.org/wiki/Tissot's_Indicatrix). 
# These diagrams illustrate the distortion inherent in all map projections.
# In conformal projections, where angles are conserved around every location, 
# the Tissot's indicatrix are all circles, with varying sizes. In equal-area 
# projections, where area proportions between objects are conserved, the 
# Tissot's indicatrix have all unit area, although their shapes and 
# orientations vary with location.

# requires Basemap version 0.99.1
if basemap_version < '0.99.1':
    raise SystemExit("this example requires Basemap version 0.99.1 or higher")

# create Basemap instances with several different projections
m1 = Basemap(llcrnrlon=-180,llcrnrlat=-80,urcrnrlon=180,urcrnrlat=80,
              projection='cyl')
m2 = Basemap(llcrnrlon=-180,llcrnrlat=-80,urcrnrlon=180,urcrnrlat=80,
              projection='mill')
m3 = Basemap(lon_0=-60,lat_0=45,projection='ortho')
m4 = Basemap(llcrnrlon=-180,llcrnrlat=-80,urcrnrlon=180,urcrnrlat=80,
             projection='merc',lat_ts=20)
m5 = Basemap(lon_0=270,lat_0=90,boundinglat=5,projection='npstere')
m6 = Basemap(lon_0=270,lat_0=90,boundinglat=5,projection='nplaea')
m7 = Basemap(lon_0=0,projection='moll')
m8 = Basemap(lon_0=0,projection='robin')
m9 = Basemap(lon_0=0,projection='hammer')
m10 = Basemap(lon_0=0,projection='mbtfpq')
m11 = Basemap(lon_0=0,projection='eck4')
m12 = Basemap(lon_0=0,projection='kav7')
m13 = Basemap(lon_0=270,lat_0=90,boundinglat=5,projection='npaeqd')

for m in [m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13]:
    # make a new figure.
    fig = plt.figure()
    # draw "circles" at specified longitudes and latitudes.
    for parallel in range(-70,71,20):
        for meridian in range(-150,151,60):
            poly = m.tissot(meridian,parallel,6,100,facecolor='green',zorder=10,alpha=0.5)
    # draw meridians and parallels.
    m.drawparallels(np.arange(-60,61,30),labels=[1,0,0,0])
    m.drawmeridians(np.arange(-180,180,60),labels=[0,0,0,1])
    # draw coastlines, fill continents, plot title.
    m.drawcoastlines()
    m.drawmapboundary(fill_color='aqua') 
    m.fillcontinents(color='coral',lake_color='aqua')
    title = 'Tissot Diagram: projection = %s' % m.projection
    sys.stdout.write(title+'\n')
    plt.title(title)

plt.show()
