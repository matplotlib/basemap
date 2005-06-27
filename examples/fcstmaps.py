# this example reads weather forecasts initialized on a day
# specified on the command line, and makes a multi-panel plot.
# Requires the OpenDAP module (a pure-python module with no dependencies
# other than Numeric) from http://opendap.oceanografia.org, and an active
# internet connection.

try:
    from opendap import client
except:
    raise ImportError,"requires opendap module (version 1.3.0 or higher) from http://opendap.oceanografia.org"
from pylab import *
from matplotlib.numerix import ma
import datetime, sys
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

# read in date (YYYYMMDDHH) from command line.

if len(sys.argv) == 1:
    YYYYMMDDHH = raw_input('Enter date forecast was started on (YYYYMMDDHH):')
else:
    YYYYMMDDHH = sys.argv[1]
YYYYMM = YYYYMMDDHH[0:6]
YYYYMMDD = YYYYMMDDHH[0:8]
HH = YYYYMMDDHH[8:10]

# set OpenDAP server URL.

try:
    URLbase="http://nomads.ncdc.noaa.gov:9090/dods/NCDC_NOAAPort_ETA/"
    URL=URLbase+YYYYMM+'/'+YYYYMMDD+'/meso-eta_211_'+YYYYMMDD+'_'+HH+'00_fff'
    print URL
    data = client.Dataset(URL)
except:
    raise IOError, 'nomad server not providing the requested data.'


# read levels, lats,lons,times.

levels = data.variables['lev']
latitudes = data.variables['lat']
longitudes = data.variables['lon']
fcsttimes = data.variables['time']
times = fcsttimes[:]
# put forecast times in YYYYMMDDHH format.
verifdates = []
fcsthrs=[]
print times
for time in times:
    fcsthrs.append(int((time-times[0])*24))
    fdate = hrs_since_day1CE_todate(int(time*24.0)) 
    verifdates.append(fdate.strftime('%Y%m%d%H'))
print fcsthrs
print verifdates
levs = levels[:]
lats = latitudes[:]
lons = longitudes[:]
lons, lats = meshgrid(lons,lats)

# unpack 2-meter temp forecast data.

print data.variables.keys()
t2mvar = data.variables['t2m']
missval = t2mvar.missing_value
t2m = t2mvar[:,:,:]
if missval < 0:
    t2m = ma.masked_values(where(t2m>-1.e20,t2m,1.e20), 1.e20)
else:
    t2m = ma.masked_values(where(t2m<1.e20.e-12,t2m,1.e20), 1.e20)
t2min = amin(t2m.compressed()); t2max= amax(t2m.compressed())

llcrnrlat = 22.0
urcrnrlat = 48.0
latminout = 22.0
llcrnrlon = -125.0
urcrnrlon = -60.0
standardpar = 50.0
centerlon=-105.
# create Basemap instance for Lambert Conformal Conic projection.
m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
            rsphere=6371200.,
            resolution='l',area_thresh=5000.,projection='lcc',
            lat_1=standardpar,lon_0=centerlon)
x, y = m(lons, lats)
# setup figure size so each panel will have same aspect ratio as map.
xsize = rcParams['figure.figsize'][0]
ysize = (3./2.)*m.aspect*xsize
fig=figure(figsize=(xsize,ysize))
yoffset = (m.urcrnry-m.llcrnry)/30.
for npanel,fcsthr in enumerate(arange(0,72,12)):
    nt = fcsthrs.index(fcsthr)
    ax = fig.add_subplot(320+npanel+1)
    # make a pcolor plot.
    #p = m.pcolor(x,y,t2m[nt,:,:],shading='flat',cmap=cm.jet)
    # make a filled contour plot.
    levels, colls = m.contourf(x,y,t2m[nt,:,:],20,cmap=cm.jet,colors=None)
    clim(t2min,t2max) # fix color range to be the same for all panels
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    m.drawparallels(arange(25,75,20),labels=[1,0,0,0],fontsize=8,fontstyle='oblique')
    m.drawmeridians(arange(-140,0,20),labels=[0,0,0,1],fontsize=8,yoffset=yoffset,fontstyle='oblique')
    # panel title
    title(repr(fcsthr)+'-h forecast valid '+verifdates[nt],fontsize=12)
# figure title
figtext(0.5,0.95,u"2-m temp (\N{DEGREE SIGN}K) forecasts from %s"%verifdates[0],horizontalalignment='center',fontsize=14)
# a single colorbar.
cax = axes([0.25, 0.03, 0.5, 0.025])
colorbar(tickfmt='%d', cax=cax, orientation='horizontal') 
show()
