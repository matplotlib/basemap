# example using matplotlib.animation to create a movie
# reads data over http - needs an active internet connection.

import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
import datetime, time
from mpl_toolkits.basemap import Basemap, shiftgrid
from netCDF4 import Dataset as NetCDFFile, date2index, num2date
import matplotlib.animation as animation

# times for March 1993 'storm of the century'
date1 = datetime.datetime(1993,3,10,0)
date2 = datetime.datetime(1993,3,17,0)

# set OpenDAP server URL.
URL="http://nomad2.ncep.noaa.gov:9090/dods/reanalyses/reanalysis-2/6hr/pgb/pgb"
try:
    data = NetCDFFile(URL)
except:
    raise IOError('opendap server not providing the requested data')

# read lats,lons,times.
latitudes = data.variables['lat'][:]
longitudes = data.variables['lon'][:].tolist()
times = data.variables['time']
ntime1 = date2index(date1,times,calendar='standard')
ntime2 = date2index(date2,times,calendar='standard')
# get sea level pressure and 10-m wind data.
slpdata = data.variables['presmsl']
udata = data.variables['ugrdprs']
vdata = data.variables['vgrdprs']
# mult slp by 0.01 to put in units of millibars.
slpin = 0.01*slpdata[ntime1:ntime2+1,:,:]
uin = udata[ntime1:ntime2+1,0,:,:] 
vin = vdata[ntime1:ntime2+1,0,:,:] 
dates = num2date(times[ntime1:ntime2+1], times.units, calendar='standard')
# add cyclic points
slp = np.zeros((slpin.shape[0],slpin.shape[1],slpin.shape[2]+1),np.float64)
slp[:,:,0:-1] = slpin; slp[:,:,-1] = slpin[:,:,0]
u = np.zeros((uin.shape[0],uin.shape[1],uin.shape[2]+1),np.float64)
u[:,:,0:-1] = uin; u[:,:,-1] = uin[:,:,0]
v = np.zeros((vin.shape[0],vin.shape[1],vin.shape[2]+1),np.float64)
v[:,:,0:-1] = vin; v[:,:,-1] = vin[:,:,0]
longitudes.append(360.); longitudes = np.array(longitudes)
# make 2-d grid of lons, lats
lons, lats = np.meshgrid(longitudes,latitudes)
# make orthographic basemap.
m = Basemap(resolution='c',projection='ortho',lat_0=60.,lon_0=-60.)
uin = udata[ntime1:ntime2+1,0,:,:] 
vin = vdata[ntime1:ntime2+1,0,:,:] 
# create figure, add axes (leaving room for colorbar on right)
fig = plt.figure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# set desired contour levels.
clevs = np.arange(960,1061,5)
# compute native x,y coordinates of grid.
x, y = m(lons, lats)
# define parallels and meridians to draw.
parallels = np.arange(-80.,90,20.)
meridians = np.arange(0.,360.,20.)
# number of repeated frames at beginning and end is n1.
nframe = 0; n1 = 10
pos = ax.get_position()
l, b, w, h = pos.bounds
# loop over times, make contour plots, draw coastlines, 
# parallels, meridians and title.
nt = 0; date = dates[nt]
CS1 = m.contour(x,y,slp[nt,:,:],clevs,linewidths=0.5,colors='k')
CS2 = m.contourf(x,y,slp[nt,:,:],clevs,cmap=plt.cm.RdBu_r)
# plot wind vectors on lat/lon grid.
# rotate wind vectors to map projection coordinates.
#urot,vrot = m.rotate_vector(u[nt,:,:],v[nt,:,:],lons,lats)
# plot wind vectors over map.
#Q = m.quiver(x,y,urot,vrot,scale=500) 
# plot wind vectors on projection grid (looks better).
# first, shift grid so it goes from -180 to 180 (instead of 0 to 360
# in longitude).  Otherwise, interpolation is messed up.
ugrid,newlons = shiftgrid(180.,u[nt,:,:],longitudes,start=False)
vgrid,newlons = shiftgrid(180.,v[nt,:,:],longitudes,start=False)
# transform vectors to projection grid.
urot,vrot,xx,yy = m.transform_vector(ugrid,vgrid,newlons,latitudes,51,51,returnxy=True,masked=True)
# plot wind vectors over map.
Q = m.quiver(xx,yy,urot,vrot,scale=500,zorder=10)
# make quiver key.
qk = plt.quiverkey(Q, 0.1, 0.1, 20, '20 m/s', labelpos='W')
# draw coastlines, parallels, meridians, title.
m.drawcoastlines(linewidth=1.5)
m.drawparallels(parallels)
m.drawmeridians(meridians)
txt = plt.title('SLP and Wind Vectors '+str(date))
# plot colorbar on a separate axes (only for first frame)
cax = plt.axes([l+w-0.05, b, 0.03, h]) # setup colorbar axes
fig.colorbar(CS2,drawedges=True, cax=cax) # draw colorbar
cax.text(0.0,-0.05,'mb')
plt.axes(ax) # reset current axes

def updatefig(nt):
    global CS1,CS2,Q
    date = dates[nt]
    for c in CS1.collections: c.remove()
    CS1 = m.contour(x,y,slp[nt,:,:],clevs,linewidths=0.5,colors='k')
    for c in CS2.collections: c.remove()
    CS2 = m.contourf(x,y,slp[nt,:,:],clevs,cmap=plt.cm.RdBu_r)
    ugrid,newlons = shiftgrid(180.,u[nt,:,:],longitudes,start=False)
    vgrid,newlons = shiftgrid(180.,v[nt,:,:],longitudes,start=False)
    urot,vrot,xx,yy = m.transform_vector(ugrid,vgrid,newlons,latitudes,51,51,returnxy=True,masked=True)
    txt.set_text('SLP and Wind Vectors '+str(date))
    Q.set_UVC(urot,vrot)

ani = animation.FuncAnimation(fig, updatefig, frames=len(dates))

#ani.save('movie.mp4')

plt.show()
