# make plot of etopo bathymetry/topography data on
# lambert conformal conic map projection, drawing coastlines, state and
# country boundaries, and parallels/meridians.

from matplotlib.toolkits.basemap import Basemap, interp
from pylab import *
import cPickle

# read in topo data from pickle (on a regular lat/lon grid)
# longitudes go from 20 to 380.
topodict = cPickle.load(open('etopo20.pickle','rb'))
topoin = topodict['data']; lons = topodict['lons']; lats = topodict['lats']

# setup of basemap ('lcc' = lambert conformal conic).
m = Basemap(-145.5,1.,-2.566,46.352,\
            resolution='c',area_thresh=10000.,projection='lcc',\
            lat_1=50.,lon_0=-107.)
# define grid (nx x ny regularly spaced native projection grid)
nx = int((m.xmax-m.xmin)/40000.)+1; ny = int((m.ymax-m.ymin)/40000.)+1
lonsout, latsout = m.makegrid(nx,ny)
# make sure there are no neg lons.
lonsout = where(lonsout < 0., lonsout + 360., lonsout)
topodat = interp(topoin,lons,lats,lonsout,latsout)
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.7,0.7])
im = imshow(topodat,cm.jet,extent=(m.xmin, m.xmax, m.ymin, m.ymax),origin='lower')

cax = axes([0.875, 0.1, 0.05, 0.7])
colorbar(tickfmt='%d', cax=cax) # draw colorbar

axes(ax)  # make the original axes current again
m.drawcoastlines(ax)
m.drawcountries(ax)
m.drawstates(ax)
#m.fillcontinents(ax)
# draw parallels
delat = 20.
parallels = arange(0.,90.+delat,delat).tolist()+\
          arange(-delat,-90.-delat,-delat).tolist()
m.drawparallels(ax,parallels,labels=[1,1,0,1])
# draw meridians
delon = 30.
meridians = arange(10.,360.,delon)
m.drawmeridians(ax,meridians,labels=[1,1,0,1])

title('ETOPO Topography - Lambert Conformal Conic')
show()
