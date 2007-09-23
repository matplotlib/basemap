# make a sequence of png files which may then be displayed
# as an animation using a tool like imagemagick animate, or
# converted to an animate gif (using imagemagick convert).

# requires pydap module.
try:
    from dap import client
except:
    raise ImportError,"requires pyDAP module (version 2.1 or higher) from http://pydap.org"
import pylab as p
from matplotlib.numerix import ma
import datetime, sys, time, subprocess
from matplotlib.toolkits.basemap import Basemap, shiftgrid

hrsgregstart = 13865688 # hrs from 00010101 to 15821015 in Julian calendar.
# times in many datasets use mixed Gregorian/Julian calendar, datetime 
# module uses a proleptic Gregorian calendar. So, I use datetime to compute
# hours since start of Greg. calendar (15821015) and add this constant to
# get hours since 1-Jan-0001 in the mixed Gregorian/Julian calendar.
gregstart = datetime.datetime(1582,10,15) # datetime.datetime instance

def dateto_hrs_since_day1CE(curdate):
    """given datetime.datetime instance, compute hours since 1-Jan-0001"""
    if curdate < gregstart:
        msg = 'date must be after start of gregorian calendar (15821015)!'
        raise ValueError, msg
    difftime = curdate-gregstart
    hrsdiff = 24*difftime.days + difftime.seconds/3600
    return hrsdiff+hrsgregstart

def hrs_since_day1CE_todate(hrs):
    """return datetime.datetime instance given hours since 1-Jan-0001"""
    if hrs < 0.0:
        msg = "hrs must be positive!"
        raise ValueError, msg
    delta = datetime.timedelta(hours=1)
    hrs_sincegreg = hrs - hrsgregstart
    curdate = gregstart + hrs_sincegreg*delta
    return curdate

# times for March 1993 'storm of the century'
YYYYMMDDHH1 = '1993031000'
YYYYMMDDHH2 = '1993031700'

YYYY = YYYYMMDDHH1[0:4]
if YYYY != YYYYMMDDHH2[0:4]:
    raise ValueError,'dates must be in same year'

# set OpenDAP server URL.
URLbase="http://nomad3.ncep.noaa.gov:9090/dods/reanalyses/reanalysis-2/6hr/pgb/"
URL=URLbase+'pres'
URLu=URLbase+'wind'
URLv=URLbase+'wind'
print URL
print URLu
print URLv
try:
    data = client.open(URL)
    datau = client.open(URLu)
    datav = client.open(URLv)
except:
    raise IOError, 'opendap server not providing the requested data'

# read lats,lons,times.
print data.keys()
print datau.keys()
print datav.keys()
latitudes = data['lat'][:]
longitudes = data['lon'][:].tolist()
times = data['time'][:]
# put times in YYYYMMDDHH format.
dates=[]
for t in times:
    t = t*24
    fdate = hrs_since_day1CE_todate(int(t))
    dates.append(fdate.strftime('%Y%m%d%H'))
if YYYYMMDDHH1 not in dates or YYYYMMDDHH2 not in dates:
    raise ValueError, 'date1 or date2 not a valid date (must be in form YYYYMMDDHH, where HH is 00,06,12 or 18)'
# find indices bounding desired times.
ntime1 = dates.index(YYYYMMDDHH1)
ntime2 = dates.index(YYYYMMDDHH2)
print 'ntime1,ntime2:',ntime1,ntime2
if ntime1 >= ntime2:
    raise ValueError,'date2 must be greater than date1'
# get sea level pressure and 10-m wind data.
slpdata = data['presmsl']
udata = datau['ugrdprs']
vdata = datau['vgrdprs']
# mult slp by 0.01 to put in units of millibars.
slpin = 0.01*p.squeeze(slpdata[ntime1:ntime2+1,:,:])
uin = p.squeeze(udata[ntime1:ntime2+1,0,:,:]) 
vin = p.squeeze(vdata[ntime1:ntime2+1,0,:,:]) 
datelabels = dates[ntime1:ntime2+1]
# add cyclic points
slp = p.zeros((slpin.shape[0],slpin.shape[1],slpin.shape[2]+1),p.Float64)
slp[:,:,0:-1] = slpin; slp[:,:,-1] = slpin[:,:,0]
u = p.zeros((uin.shape[0],uin.shape[1],uin.shape[2]+1),p.Float64)
u[:,:,0:-1] = uin; u[:,:,-1] = uin[:,:,0]
v = p.zeros((vin.shape[0],vin.shape[1],vin.shape[2]+1),p.Float64)
v[:,:,0:-1] = vin; v[:,:,-1] = vin[:,:,0]
longitudes.append(360.); longitudes = p.array(longitudes)
# make 2-d grid of lons, lats
lons, lats = p.meshgrid(longitudes,latitudes)
print 'min/max slp,u,v'
print min(p.ravel(slp)),max(p.ravel(slp))
print min(p.ravel(uin)),max(p.ravel(uin))
print min(p.ravel(vin)),max(p.ravel(vin))
print 'dates'
print datelabels
# make orthographic basemap.
m = Basemap(resolution='c',projection='ortho',lat_0=60.,lon_0=-60.)
p.ion() # interactive mode on.
uin = p.squeeze(udata[ntime1:ntime2+1,0,:,:]) 
vin = p.squeeze(vdata[ntime1:ntime2+1,0,:,:]) 
datelabels = dates[ntime1:ntime2+1]
# make orthographic basemap.
m = Basemap(resolution='c',projection='ortho',lat_0=60.,lon_0=-60.)
p.ion() # interactive mode on.
# create figure, add axes (leaving room for colorbar on right)
fig = p.figure()
ax = fig.add_axes([0.1,0.1,0.7,0.7])
# set desired contour levels.
clevs = p.arange(960,1061,5)
# compute native x,y coordinates of grid.
x, y = m(lons, lats)
# define parallels and meridians to draw.
parallels = p.arange(-80.,90,20.)
meridians = p.arange(0.,360.,20.)
# number of repeated frames at beginning and end is n1.
nframe = 0; n1 = 10
l,b,w,h=ax.get_position()
# loop over times, make contour plots, draw coastlines, 
# parallels, meridians and title.
for nt,date in enumerate(datelabels[1:]):
    CS = m.contour(x,y,slp[nt,:,:],clevs,linewidths=0.5,colors='k',animated=True)
    CS = m.contourf(x,y,slp[nt,:,:],clevs,cmap=p.cm.RdBu_r,animated=True)
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
    Q = m.quiver(xx,yy,urot,vrot,scale=500)
    # make quiver key.
    qk = p.quiverkey(Q, 0.1, 0.1, 20, '20 m/s', labelpos='W')
    # draw coastlines, parallels, meridians, title.
    m.drawcoastlines(linewidth=1.5)
    m.drawparallels(parallels)
    m.drawmeridians(meridians)
    p.title('SLP and Wind Vectors '+date)
    if nt == 0: # plot colorbar on a separate axes (only for first frame)
        cax = p.axes([l+w-0.05, b, 0.03, h]) # setup colorbar axes
        fig.colorbar(CS,drawedges=True, cax=cax) # draw colorbar
        cax.text(0.0,-0.05,'mb')
        p.axes(ax) # reset current axes
    p.draw() # draw the plot
    # save first and last frame n1 times 
    # (so gif animation pauses at beginning and end)
    if nframe == 0 or nt == slp.shape[0]-1:
       for n in range(n1):
           p.savefig('anim%03i'%nframe+'.png')
           nframe = nframe + 1
    else:
       p.savefig('anim%03i'%nframe+'.png')
       nframe = nframe + 1
    ax.clear() # clear the axes for the next plot.

print """
Now display animation using imagemagick 'animate -delay 10 anim*png'
or, create an animated gif using 'convert -delay 10 anim*png anim.gif'
and display in a web browser (or Powerpoint).

Will try to run imagemagick animate in 5 seconds ..."""
time.sleep(5)
subprocess.call('animate -delay 10 anim*png',shell=True)
