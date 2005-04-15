from matplotlib.toolkits.basemap import Basemap, interp, shiftgrid
from pylab import *
import cPickle

# read in data on lat/lon grid.
datadict = cPickle.load(open('500hgt.pickle','rb'))
hgt = datadict['data']; lons = datadict['lons']; lats = datadict['lats']
# shift data so lons go from -180 to 180 instead of 0 to 360.
hgt,lons = shiftgrid(180.,hgt,lons,start=False)

# set up map projection (lambert azimuthal equal area).
m = Basemap(-135.,-20.,45.,-20.,
             resolution='c',area_thresh=10000.,projection='laea',
             lat_0=90.,lon_0=-90.)

# transform to map projection coordinates.
nx = 101; ny = 101
hgt,x,y = m.transform_scalar(hgt,lons,lats,nx,ny,returnxy=True)

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
        im = imshow(hgt,extent=(m.llcrnrx,m.urcrnrx,m.llcrnry,m.urcrnry),origin='lower')
    elif plot == 'contour':
        levels, colls = contour(x,y,hgt,15,linewidths=0.5,colors='k')
        levels, colls = contourf(x,y,hgt,15,cmap=cm.jet,colors=None)

    # set size of plot to match aspect ratio of map.
    ax = gca()
    l,b,w,h = ax.get_position()
    b = 0.5 - 0.5*w*m.aspect; h = w*m.aspect
    #l = 0.5 - 0.5*h/m.aspect; w = h/m.aspect
    ax.set_position([l,b,w,h])
    corners = ((m.llcrnrx,m.llcrnry), (m.urcrnrx,m.urcrnry))
    ax.update_datalim( corners )                                          
    axis([m.llcrnrx, m.urcrnrx, m.llcrnry, m.urcrnry])  
    
    # draw map.
    m.drawcoastlines(ax)
	
    # draw parallels
    delat = 30.
    delon = 90.
    circles = arange(10.,90.+delat,delat).tolist()
    m.drawparallels(ax,circles,labels=[0,0,1,1])

    # draw meridians
    meridians = arange(0.,360.,delon)
    m.drawmeridians(ax,meridians,labels=[1,1,1,1],fontsize=10)

    title('500 hPa Height - '+plot,y=1.075)

show()
