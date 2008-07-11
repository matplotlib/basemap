import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap 
from matplotlib.patches import Polygon

# Tissot's Indicatrix (http://en.wikipedia.org/wiki/Tissot's_Indicatrix). 
# These diagrams illustrate the distortion inherent in all map projections.
# In conformal projections, where angles are conserved around every location, 
# the Tissot's indicatrix are all circles, with varying sizes. In equal-area 
# projections, where area proportions between objects are conserved, the 
# Tissot's indicatrix have all unit area, although their shapes and 
# orientations vary with location.

# create Basemap instances with several different projections
m1 = Basemap(llcrnrlon=-180,llcrnrlat=-80,urcrnrlon=180,urcrnrlat=80,
              projection='cyl')
m2 = Basemap(lon_0=-60,lat_0=45,projection='ortho')
# use WGS84 ellipsoid in this one.
m3 = Basemap(llcrnrlon=-180,llcrnrlat=-70,urcrnrlon=180,urcrnrlat=70,
              projection='merc',lat_ts=20,rsphere=(6378137.0,6356752.3142))
m4 = Basemap(lon_0=270,lat_0=90,boundinglat=10,projection='npstere')
m5 = Basemap(lon_0=270,lat_0=90,boundinglat=10,projection='nplaea')

for m in [m1,m2,m3,m4,m5]:
    # make a new figure.
    fig = plt.figure()
    ax = plt.gca()
    # draw "circles" at specified longitudes and latitudes.
    for parallel in range(-60,61,30):
        for meridian in range(-165,166,30):
            seg = m.tissot(meridian,parallel,6,100)
            poly = Polygon(seg,facecolor='green',zorder=10)
            ax.add_patch(poly)
    # draw meridians and parallels.
    m.drawparallels(np.arange(-60,61,30))
    m.drawmeridians(np.arange(-180,180,60))
    # draw coastlines, fill continents, plot title.
    m.drawcoastlines()
    m.fillcontinents()
    title = 'Tissot Diagram: projection = %s' % m.projection
    print title
    plt.title(title)

plt.show()
