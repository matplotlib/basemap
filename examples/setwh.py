# examples of using the 'width' and 'height' keywords
# to the Basemap constructor.

from matplotlib.toolkits.basemap import Basemap
from pylab import arange, show, title, figure

# setup projection centered on lon_0,lat_0
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
fig = figure(figsize=(8,12))
for proj in projs:
    m = Basemap(width=width,height=height,
                resolution='c',projection=proj,\
                lat_0=lat_0,lon_0=lon_0)
    npanel = npanel + 1
    fig.add_subplot(3,2,npanel)
    # setup figure with same aspect ratio as map.
    m.drawcoastlines()
    m.drawcountries()
    m.fillcontinents()
    m.drawstates()
    m.drawparallels(circles)
    m.drawmeridians(meridians)
    title('proj = '+proj+' centered on %sW, %sN' % (lon_0,lat_0),fontsize=10)

proj = 'omerc'
delat = 10.
circles = arange(0.,90.+delat,delat).tolist()+\
          arange(-delat,-90.-delat,-delat).tolist()
delon = 10.
meridians = arange(10.,360.,delon)
lat_1 = 40; lat_2 = 55
lon_1 = -120; lon_2 = -140
lat_0 = 47.5 ; lon_0 = -130
fig = figure()
width = 1500000.
height = 2.5*width
m = Basemap(width=width,height=height,
            resolution='l',projection=proj,\
            lon_1=lon_1,lon_2=lon_2,\
            lat_1=lat_1,lat_2=lat_2,\
            lat_0=lat_0,lon_0=lon_0)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.fillcontinents()
m.drawstates(linewidth=0.5)
m.drawparallels(circles)
m.drawmeridians(meridians)
title('proj = '+proj+' centered on %sW, %sN' % (lon_0,lat_0),fontsize=10)


#lon_0 = -8
#lat_0 = 53.3
#width=350000.
#height= 1.33*width
#fig=figure()
#proj='tmerc'
#circles = arange(50,60,1)
#meridians = arange(-12,2,1)
#m = Basemap(width=width,height=height,
#            resolution='i',projection=proj,\
#            lat_0=lat_0,lon_0=lon_0)
#m.drawcoastlines(linewidth=0.5)
#m.drawcountries(linewidth=0.5)
#m.fillcontinents()
#m.drawstates(linewidth=0.5)
#m.drawparallels(circles)
#m.drawmeridians(meridians)
#title('proj = '+proj+' centered on %sW, %sN' % (lon_0,lat_0),fontsize=10)

show()
