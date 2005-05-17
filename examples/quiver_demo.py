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
lats = -90.+dellat*arange(nlats)
lons = -180.+dellon*arange(nlons)

# plot vectors in geographical (lat/lon) coordinates.

# north polar projection.
m = Basemap(-180.,10.,0.,10.,
            resolution='c',area_thresh=10000.,projection='stere',\
            lat_0=90.,lon_0=-135.,lat_ts=90.)
# setup figure with same aspect ratio as map.
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# rotate wind vectors to map projection coordinates.
# (also compute native map projections coordinates of lat/lon grid)
# only do Northern Hemisphere.
urot,vrot,x,y = m.rotate_vector(u[36:,:],v[36:,:],lons,lats[36:],returnxy=True)
# plot filled contours over map.
levels, colls = m.contourf(x,y,p[36:,:],15,cmap=cm.jet,colors=None)
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
m = Basemap(-180.,10.,0.,10.,
            resolution='c',area_thresh=10000.,projection='stere',\
            lat_0=90.,lon_0=-135.,lat_ts=90.)
# transform from spherical to map projection coordinates (rotation
# and interpolation).
nxv = 41; nyv = 41
nxp = 101; nyp = 101
spd = sqrt(u**2+v**2)
print max(ravel(spd))
udat, vdat, xv, yv = m.transform_vector(u,v,lons,lats,nxv,nyv,returnxy=True)
pdat, xp, yp = m.transform_scalar(p,lons,lats,nxp,nyp,returnxy=True)
# setup figure with same aspect ratio as map.
xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# plot filled contours over map
levels, colls = m.contourf(xp,yp,pdat,15,cmap=cm.jet,colors=None)
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
