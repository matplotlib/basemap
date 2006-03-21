from matplotlib.toolkits.basemap import Basemap, interp
from pylab import *

# read in data.
file = open('fcover.dat','r')
ul=[];vl=[];pl=[]
nlons=73; nlats=73
dellat = 2.5; dellon = 5.
for line in file.readlines():
   l = line.replace('\n','').split()
   ul.append(float(l[0]))
   vl.append(float(l[1]))
   pl.append(float(l[2]))
u = reshape(array(ul,'f'),(nlats,nlons))
v = reshape(array(vl,'f'),(nlats,nlons))
p = reshape(array(pl,'f'),(nlats,nlons))
lats1 = -90.+dellat*arange(nlats)
lons1 = -180.+dellon*arange(nlons)
lons, lats = meshgrid(lons1, lats1)

# plot vectors in geographical (lat/lon) coordinates.

# north polar projection.
m = Basemap(lon_0=-135,boundinglat=25,
            resolution='c',area_thresh=10000.,projection='npstere')
# setup figure with same aspect ratio as map.
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# rotate wind vectors to map projection coordinates.
# (also compute native map projections coordinates of lat/lon grid)
# only do Northern Hemisphere.
urot,vrot,x,y = m.rotate_vector(u[36:,:],v[36:,:],lons[36:,:],lats[36:,:],returnxy=True)
# plot filled contours over map.
cs = m.contourf(x,y,p[36:,:],15,cmap=cm.jet)
# plot wind vectors over map.
m.quiver(x,y,urot,vrot)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes.
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines()
m.drawcountries()
# draw parallels
delat = 20.
circles = arange(0.,90.+delat,delat).tolist()+\
          arange(-delat,-90.-delat,-delat).tolist()
m.drawparallels(circles,labels=[1,1,1,1])
# draw meridians
delon = 45.
meridians = arange(-180,180,delon)
m.drawmeridians(meridians,labels=[1,1,1,1])
title('Surface Winds Winds and Pressure',y=1.075)
show()

# plot vectors in map projection coordinates.

# north polar projection.
m = Basemap(lon_0=-135,boundinglat=25,
            resolution='c',area_thresh=10000.,projection='npstere')
# transform from spherical to map projection coordinates (rotation
# and interpolation).
nxv = 41; nyv = 41
nxp = 101; nyp = 101
spd = sqrt(u**2+v**2)
udat, vdat, xv, yv = m.transform_vector(u,v,lons1,lats1,nxv,nyv,returnxy=True)
pdat, xp, yp = m.transform_scalar(p,lons1,lats1,nxp,nyp,returnxy=True)
# setup figure with same aspect ratio as map.
fig=m.createfigure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# plot image over map
im = m.imshow(pdat,cm.jet)
# plot wind vectors over map.
m.quiver(xv,yv,udat,vdat)
cax = axes([0.875, 0.1, 0.05, 0.7]) # setup colorbar axes.
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines()
m.drawcountries()
# draw parallels
delat = 20.
circles = arange(0.,90.+delat,delat).tolist()+\
          arange(-delat,-90.-delat,-delat).tolist()
m.drawparallels(circles,labels=[1,1,1,1])
# draw meridians
delon = 45.
meridians = arange(-180,180,delon)
m.drawmeridians(meridians,labels=[1,1,1,1])
title('Surface Winds Winds and Pressure',y=1.075)
show()
