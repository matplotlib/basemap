from matplotlib.toolkits.basemap import Basemap, interp
from pylab import *
import cPickle
from matplotlib import __version__

if __version__ <= '0.70':
    print 'warning - contour example will not work with this version of matplotlib if axes.patch is not applied'

# read in data on lat/lon grid.
datadict = cPickle.load(open('500hgt.pickle','rb'))
hgt = datadict['data']; lons = datadict['lons']; lats = datadict['lats']

# set up map projection (lambert azimuthal equal area).
m = Basemap(-150.,-20.,30.,-20.,
             resolution='c',area_thresh=10000.,projection='laea',
             lat_0=90.,lon_0=-105.)
# interpolate to map projection grid.
nx = 101
ny = 101
lonsout, latsout = m.makegrid(nx,ny)
hgt = interp(hgt,lons,lats,lonsout,latsout)
dx = (m.xmax-m.xmin)/(nx-1)
dy = (m.ymax-m.ymin)/(ny-1)  
x = m.xmin+dx*indices((ny,nx))[1,:,:]
y = m.ymin+dy*indices((ny,nx))[0,:,:]

#m = Basemap(lons[0],lats[0],lons[-1],lats[-1],\
#              resolution='c',area_thresh=10000.,projection='cyl')
#x, y = meshgrid(lons, lats)

fig = figure(figsize=(8,8))

plots = ['contour','pcolor']

for np,plot in enumerate(plots):

    fig.add_subplot(1,2,np+1)
    #fig.add_subplot(2,1,np+1)

    # plot data.
    print plot+' plot ...'
    if plot == 'pcolor':
        pcolor(x,y,hgt,shading='flat')
    elif plot == 'imshow':
        im = imshow(hgt,extent=(m.xmin, m.xmax, m.ymin, m.ymax),origin='lower')
    elif plot == 'contour':
        levels, colls = contour(hgt,x=x,y=y,levels=20,linewidths=1.,colors='k')

    # set size of plot to match aspect ratio of map.
    ax = gca()
    l,b,w,h = ax.get_position()
    b = 0.5 - 0.5*w*m.aspect; h = w*m.aspect
    #l = 0.5 - 0.5*h/m.aspect; w = h/m.aspect
    ax.set_position([l,b,w,h])
    corners = ((m.xmin,m.ymin), (m.xmax,m.ymax))
    ax.update_datalim( corners )                                          
    axis([m.xmin, m.xmax, m.ymin, m.ymax])  
    
    # draw map.
    if plot in ['pcolor','imshow']:                                         
        m.drawcoastlines(ax)
    else:
        m.fillcontinents(ax,color='gray')
	
    # draw parallels
    delat = 30.
    delon = 30.
    circles = arange(0.,90.+delat,delat).tolist()+\
	   arange(-delat,-90.-delat,-delat).tolist()
    m.drawparallels(ax,circles)

    # draw meridians
    meridians = arange(0.,360.,delon)
    m.drawmeridians(ax,meridians)

    ax.set_xticks([]) # no ticks
    ax.set_yticks([])
    title('500 hPa Height - '+plot)

show()
