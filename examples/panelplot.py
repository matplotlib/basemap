from matplotlib.toolkits.basemap import Basemap
from matplotlib import rcParams
import pylab as P


# read in data on lat/lon grid.
hgt  = P.array(P.load('500hgtdata.gz'),'d')
lons = P.array(P.load('500hgtlons.gz'),'d')
lats = P.array(P.load('500hgtlats.gz'),'d')
lons, lats = P.meshgrid(lons, lats)

# Example to show how to make panel plots, taking care 
# to preserve aspect ratio of map (so distances on map
# are not distorted).

# If you're making a single panel plot, the createfigure method
# will create a figure with the same aspect ratio as the map
# projection region (so that the map looks right when you plot it).
# However, when using subplot to make multiple panels this won't work.

# This example shows how to use axis('scaled') to fix this.
# It makes two 2-panel plots, one oriented vertically and one
# horizontally.  The map should be perfected square in both cases.
# If you comment out the axis('scaled'), you'll see that the map
# become stretched horizontally (in plot 1) or vertically 
# (in plot 2).

# 2-panel plot, oriented vertically, colorbar on bottom.

rcParams['figure.subplot.hspace'] = 0.4 # more height between subplots
rcParams['figure.subplot.wspace'] = 0.5 # more width between subplots

# panel 1
mnh = Basemap(llcrnrlon=-150,llcrnrlat=-20.826,urcrnrlon=30,urcrnrlat=-20.826,\
             resolution='c',area_thresh=10000.,projection='laea',\
             lat_0=90.,lon_0=-105.)
xnh,ynh = mnh(lons,lats)
fig = P.figure(figsize=(8,8))
ax = fig.add_subplot(211)
P.axis('scaled')
# recenter plot
l,b,w,h = ax.get_position()
l = 0.5*(1.-w)
ax.set_position([l,b,w,h])
CS = mnh.contour(xnh,ynh,hgt,15,linewidths=0.5,colors='k')
CS = mnh.contourf(xnh,ynh,hgt,15,cmap=P.cm.Spectral)
# colorbar on bottom.
cax = P.axes([l, b-0.05, w, 0.025]) # setup colorbar axes
P.colorbar(tickfmt='%d', cax=cax,orientation='horizontal',clabels=CS.levels[0::3]) # draw colorbar
P.axes(ax)  # make the original axes current again
mnh.drawcoastlines(linewidth=0.5)
delat = 30.
circles = P.arange(0.,90.,delat).tolist()+\
          P.arange(-delat,-90,-delat).tolist()
mnh.drawparallels(circles,labels=[1,0,0,0])
delon = 60.
meridians = P.arange(0,360,delon)
mnh.drawmeridians(meridians,labels=[1,0,0,1])
P.title('NH 500 hPa Height (cm.Spectral)')

# panel 2
msh = Basemap(llcrnrlon=-150,llcrnrlat=20.826,urcrnrlon=30,urcrnrlat=20.826,\
              resolution='c',area_thresh=10000.,projection='laea',\
              lat_0=-90.,lon_0=-105.)
xsh,ysh = msh(lons,lats)
ax = fig.add_subplot(212)
P.axis('scaled') # watch what happens if you comment this out!
# recenter plot
l,b,w,h = ax.get_position()
l = 0.5*(1.-w)
ax.set_position([l,b,w,h])
CS = msh.contour(xsh,ysh,hgt,15,linewidths=0.5,colors='k')
CS = msh.contourf(xsh,ysh,hgt,15,cmap=P.cm.Spectral)
# colorbar on bottom.
l,b,w,h = ax.get_position()
cax = P.axes([l, b-0.05, w, 0.025]) # setup colorbar axes
P.colorbar(tickfmt='%d', cax=cax,orientation='horizontal',clabels=CS.levels[0::3]) # draw colorbar
P.axes(ax)  # make the original axes current again
msh.drawcoastlines(linewidth=0.5)
msh.drawparallels(circles,labels=[1,0,0,0])
msh.drawmeridians(meridians,labels=[1,0,0,1])
P.title('SH 500 hPa Height (cm.Spectral)')
P.show()

# 2-panel plot, oriented horizontally, colorbar on right.

# panel 1
fig = P.figure(figsize=(8,8))
ax = fig.add_subplot(121)
P.axis('scaled')
# center the plot
l,b,w,h = ax.get_position()
l = l-0.025; b = 0.5*(1.-h)
ax.set_position([l,b,w,h])
CS = mnh.contour(xnh,ynh,hgt,15,linewidths=0.5,colors='k')
CS = mnh.contourf(xnh,ynh,hgt,15,cmap=P.cm.RdBu)
# colorbar on right.
cax = P.axes([l+w+0.025, b, 0.025, h]) # setup colorbar axes
P.colorbar(tickfmt='%d', cax=cax, clabels=CS.levels[0::2]) # draw colorbar
P.axes(ax)  # make the original axes current again
mnh.drawcoastlines(linewidth=0.5)
mnh.drawparallels(circles,labels=[1,0,0,0])
mnh.drawmeridians(meridians,labels=[1,0,0,1])
P.title('NH 500 hPa Height (cm.RdBu)')

# panel 2
ax = fig.add_subplot(122)
P.axis('scaled')
# center the plot
l,b,w,h = ax.get_position()
l = l-0.025; b = 0.5*(1.-h)
ax.set_position([l,b,w,h])
CS = msh.contour(xsh,ysh,hgt,15,linewidths=0.5,colors='k')
CS = msh.contourf(xsh,ysh,hgt,15,cmap=P.cm.RdBu)
# colorbar on right.
cax = P.axes([l+w+0.025, b, 0.025, h]) # setup colorbar axes
P.colorbar(tickfmt='%d', cax=cax, clabels=CS.levels[0::2]) # draw colorbar
P.axes(ax)  # make the original axes current again
msh.drawcoastlines(linewidth=0.5)
msh.drawparallels(circles,labels=[1,0,0,0])
msh.drawmeridians(meridians,labels=[1,0,0,1])
P.title('SH 500 hPa Height (cm.RdBu)')
P.show()
