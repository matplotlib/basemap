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

# north polar projection.
m = Basemap(-180.,10.,0.,10.,
            resolution='c',area_thresh=10000.,projection='stere',\
            lat_0=90.,lon_0=-135.,lat_ts=90.)

# transform from spherical to map projection coordinates.
nxv = 41; nyv = 41
nxp = 101; nyp = 101
spd = sqrt(u**2+v**2)
print max(ravel(spd))
udat, vdat, xv, yv = m.transform_vector(u,v,lons,lats,nxv,nyv,returnxy=True)
pdat, xp, yp = m.transform_scalar(p,lons,lats,nxp,nyp,returnxy=True)

print min(ravel(udat)),max(ravel(udat))
print min(ravel(vdat)),max(ravel(vdat))
print min(ravel(pdat)),max(ravel(pdat))
spd = sqrt(udat**2+vdat**2)
print max(ravel(spd))

xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.7,0.7])
levels, colls = contourf(xp,yp,pdat,15,cmap=cm.jet,colors=None)

scale = max([(m.xmax-m.xmin)/(nxv-1),(m.ymax-m.ymin)/(nyv-1)])
quiver(xv,yv,udat,vdat,scale)
ax=gca()
corners = ((m.llcrnrx,m.llcrnry), (m.urcrnrx,m.urcrnry))
ax.update_datalim( corners )                                          
axis([m.llcrnrx, m.urcrnrx, m.llcrnry, m.urcrnry])  
cax = axes([0.875, 0.1, 0.05, 0.7])
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines(ax)
m.drawcountries(ax)
# draw parallels
delat = 20.
circles = arange(0.,90.+delat,delat).tolist()+\
          arange(-delat,-90.-delat,-delat).tolist()
m.drawparallels(ax,circles,labels=[1,1,1,1])
# draw meridians
delon = 45.
meridians = arange(-180,180,delon)
m.drawmeridians(ax,meridians,labels=[1,1,1,1])
ax.set_xticks([]) # no ticks
ax.set_yticks([])
title('Surface Winds Winds and Pressure',y=1.075)
show()

# mercator projection
m = Basemap(-180.,-80.,180.,80.,\
            resolution='c',area_thresh=10000.,projection='merc',\
            lon_0=0.,lat_ts=20.)

# transform from spherical to map projection coordinates.
nxv = 41; nyv = 41
nxp = 101; nyp = 101
udat, vdat, xv, yv = m.transform_vector(u,v,lons,lats,nxv,nyv,returnxy=True)
pdat, xp, yp = m.transform_scalar(p,lons,lats,nxp,nyp,returnxy=True)

print min(ravel(udat)),max(ravel(udat))
print min(ravel(vdat)),max(ravel(vdat))
print min(ravel(pdat)),max(ravel(pdat))
spd = sqrt(udat**2+vdat**2)
print max(ravel(spd))

xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.7,0.7])
levels, colls = contourf(xp,yp,pdat,15,cmap=cm.jet,colors=None)

scale = max([(m.xmax-m.xmin)/(nxv-1),(m.ymax-m.ymin)/(nyv-1)])
quiver(xv,yv,udat,vdat,scale)
ax=gca()
corners = ((m.llcrnrx,m.llcrnry), (m.urcrnrx,m.urcrnry))
ax.update_datalim( corners )                                          
axis([m.llcrnrx, m.urcrnrx, m.llcrnry, m.urcrnry])  
cax = axes([0.875, 0.1, 0.05, 0.7])
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines(ax)
m.drawcountries(ax)
# draw parallels
delat = 30.
circles = arange(0.,90.+delat,delat).tolist()+\
          arange(-delat,-90.-delat,-delat).tolist()
m.drawparallels(ax,circles,labels=[1,0,0,1])
# draw meridians
delon = 60.
meridians = arange(-180,180,delon)
m.drawmeridians(ax,meridians,labels=[1,0,0,1])
ax.set_xticks([]) # no ticks
ax.set_yticks([])
title('Surface Winds Winds and Pressure',y=1.075)
show()

#cylindrical equidistant
m = Basemap(-180,-90.,180.,90.,\
            resolution='c',area_thresh=10000.,projection='cyl')

# transform from spherical to map projection coordinates.
# (for 'cyl' projection this does nothing).
nxv = 41; nyv = 41
nxp = 101; nyp = 101
udat, vdat, xv, yv = m.transform_vector(u,v,lons,lats,nxv,nyv,returnxy=True)
pdat, xp, yp = m.transform_scalar(p,lons,lats,nxp,nyp,returnxy=True)

print min(ravel(udat)),max(ravel(udat))
print min(ravel(vdat)),max(ravel(vdat))
print min(ravel(pdat)),max(ravel(pdat))
spd = sqrt(udat**2+vdat**2)
print max(ravel(spd))

xsize = rcParams['figure.figsize'][0]
fig=figure(figsize=(xsize,m.aspect*xsize))
ax = fig.add_axes([0.1,0.1,0.7,0.7])
levels, colls = contourf(xp,yp,pdat,15,cmap=cm.jet,colors=None)

scale = max([(m.xmax-m.xmin)/(nxv-1),(m.ymax-m.ymin)/(nyv-1)])
quiver(xv,yv,udat,vdat,scale)
ax=gca()
corners = ((m.llcrnrx,m.llcrnry), (m.urcrnrx,m.urcrnry))
ax.update_datalim( corners )                                          
axis([m.llcrnrx, m.urcrnrx, m.llcrnry, m.urcrnry])  
cax = axes([0.875, 0.1, 0.05, 0.7])
colorbar(tickfmt='%d', cax=cax) # draw colorbar
axes(ax)  # make the original axes current again
m.drawcoastlines(ax)
m.drawcountries(ax)
# draw parallels
delat = 30.
circles = arange(0.,90.+delat,delat).tolist()+\
          arange(-delat,-90.-delat,-delat).tolist()
m.drawparallels(ax,circles,labels=[1,0,0,1])
# draw meridians
delon = 60.
meridians = arange(-180,180,delon)
m.drawmeridians(ax,meridians,labels=[1,0,0,1])
ax.set_xticks([]) # no ticks
ax.set_yticks([])
title('Surface Winds Winds and Pressure',y=1.075)
show()
