from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

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
u = np.reshape(np.array(ul,np.float32),(nlats,nlons))
v = np.reshape(np.array(vl,np.float32),(nlats,nlons))
p = np.reshape(np.array(pl,np.float32),(nlats,nlons))
lats1 = -90.+dellat*np.arange(nlats)
lons1 = -180.+dellon*np.arange(nlons)
lons, lats = np.meshgrid(lons1, lats1)
# convert from mps to knots.
u = 1.944*u; v = 1.944*v

# plot barbs in map projection coordinates.

# stereogrpaphic projection.
m = Basemap(width=10000000,height=10000000,lon_0=-90,lat_0=45.,lat_ts=45,
            resolution='l',projection='stere')
x,y = m(lons,lats)
# transform from spherical to map projection coordinates (rotation
# and interpolation).
nxv = 25; nyv = 25
udat, vdat, xv, yv = m.transform_vector(u,v,lons1,lats1,nxv,nyv,returnxy=True)
# create a figure, add an axes.
fig=plt.figure(figsize=(8,6))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
# plot color-filled contours over map
levs = np.arange(960,1051,4)
cs1 = m.contour(x,y,p,levs,colors='k',linewidths=0.5)
cs2 = m.contourf(x,y,p,levs)
# plot barbs.
m.barbs(xv,yv,udat,vdat,length=6,barbcolor='k',flagcolor='r',linewidth=0.5)
# plot colorbar for pressure
m.colorbar(pad='12%') # draw colorbar
# draw coastlines
m.drawcoastlines()
# draw parallels
m.drawparallels(np.arange(0,81,20),labels=[1,1,0,0])
# draw meridians
m.drawmeridians(np.arange(-180,0,20),labels=[0,0,0,1])
plt.title('Surface Wind Barbs and Pressure (NH)')

# stereogrpaphic projection (SH).
# 'flip_barb' flag is automatically set for SH data, so that
# barbs point toward lower pressure (in both Hemisphere).
m = Basemap(width=10000000,height=10000000,lon_0=-90,lat_0=-45.,lat_ts=-45,
            resolution='l',projection='stere')
x,y = m(lons,lats)
# transform from spherical to map projection coordinates (rotation
# and interpolation).
nxv = 25; nyv = 25
udat, vdat, xv, yv = m.transform_vector(u,v,lons1,lats1,nxv,nyv,returnxy=True)
# create a figure, add an axes.
fig=plt.figure(figsize=(8,6))
ax = fig.add_axes([0.1,0.1,0.8,0.8])
# plot color-filled contours over map
levs = np.arange(960,1051,4)
cs1 = m.contour(x,y,p,levs,colors='k',linewidths=0.5)
cs2 = m.contourf(x,y,p,levs)
# plot barbs.
m.barbs(xv,yv,udat,vdat,length=6,barbcolor='k',flagcolor='r',linewidth=0.5)
# plot colorbar for pressure
m.colorbar(pad='12%') # draw colorbar
# draw coastlines
m.drawcoastlines()
# draw parallels
m.drawparallels(np.arange(-80,-19,20),labels=[1,1,0,0])
# draw meridians
m.drawmeridians(np.arange(-180,0,20),labels=[0,0,1,0])
plt.title('Surface Wind Barbs and Pressure (SH)',y=1.04)
plt.show()
