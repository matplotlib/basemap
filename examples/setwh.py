# examples of using the 'width' and 'height' keywords
# to the Basemap constructor.

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

# setup projection parameters
lat_0 = 40.
lon_0 = -100.
width = 6000000.
height = 2.*width/3.
delat = 25.
circles = np.arange(0.,90.+delat,delat).tolist()+\
          np.arange(-delat,-90.-delat,-delat).tolist()
delon = 30.
meridians = np.arange(10.,360.,delon)
npanel = 0
# plots of the US.
projs = ['lcc','aeqd','aea','laea','eqdc','stere']
fig = plt.figure(figsize=(7,7))
for proj in projs:
    m = Basemap(width=width,height=height,
                resolution='c',projection=proj,\
                lat_0=lat_0,lon_0=lon_0)
    npanel = npanel + 1
    fig.add_subplot(3,2,npanel)
    m.drawcoastlines()
    m.drawcountries()
    m.fillcontinents()
    m.drawstates()
    m.drawparallels(circles)
    m.drawmeridians(meridians)
    plt.title('proj = '+proj+' centered on %sW, %sN' % (lon_0,lat_0),fontsize=10)

plt.show()
