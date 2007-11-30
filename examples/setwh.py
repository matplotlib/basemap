# examples of using the 'width' and 'height' keywords
# to the Basemap constructor.

from matplotlib.toolkits.basemap import Basemap
from pylab import arange, show, title, figure

# setup projection parameters
lat_0 = 40.
lon_0 = -100.
width = 6000000.
height = 2.*width/3.
delat = 25.
circles = arange(0.,90.+delat,delat).tolist()+\
          arange(-delat,-90.-delat,-delat).tolist()
delon = 30.
meridians = arange(10.,360.,delon)
npanel = 0
# plots of the US.
projs = ['lcc','aeqd','aea','laea','eqdc','stere']
fig = figure(figsize=(7,7))
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
    title('proj = '+proj+' centered on %sW, %sN' % (lon_0,lat_0),fontsize=10)

show()
