# this example reads today's numerical weather forecasts 
# from the NOAA OpenDAP servers and makes a multi-panel plot.
# Requires the pyDAP module (a pure-python module)
# from http://pydap.org, and an active intenet connection.

try:
    from dap import client
except:
    raise ImportError,"requires pyDAP module (version 2.1 or higher) from http://pydap.org"
from pylab import title, show, figure, cm, arange, frange, figtext, \
                  meshgrid, axes, colorbar, where, amin, amax, around
import sys
from matplotlib.numerix import ma
import datetime
from matplotlib.toolkits.basemap import Basemap

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

# today's date is default.
if len(sys.argv) > 1:
    YYYYMMDD = sys.argv[1]
else:
   YYYYMMDD = datetime.datetime.today().strftime('%Y%m%d')
YYYYMM = YYYYMMDD[0:6]

# set OpenDAP server URL.
HH='09'
URLbase="http://nomad3.ncep.noaa.gov:9090/dods/sref/sref"
URL=URLbase+YYYYMMDD+"/sref_eta_ctl1_"+HH+"z"
print URL+'\n'
try:
    data = client.open(URL)
except:
    msg = """
opendap server not providing the requested data.
Try another date by providing YYYYMM on command line."""
    raise IOError, msg


# read levels, lats,lons,times.

print data.keys()
levels = data['lev']
latitudes = data['lat']
longitudes = data['lon']
fcsttimes = data['time']
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

t2mvar = data['tmp2m']
missval = t2mvar.missing_value
t2m = t2mvar[:,:,:]
if missval < 0:
    t2m = ma.masked_values(where(t2m>-1.e20,t2m,1.e20), 1.e20)
else:
    t2m = ma.masked_values(where(t2m<1.e20,t2m,1.e20), 1.e20)
t2min = amin(t2m.compressed()); t2max= amax(t2m.compressed())
print t2min,t2max
clevs = frange(around(t2min/10.)*10.-5.,around(t2max/10.)*10.+5.,4)
print clevs[0],clevs[-1]
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
# create figure.
fig=figure(figsize=(8,8))
yoffset = (m.urcrnry-m.llcrnry)/30.
for npanel,fcsthr in enumerate(arange(0,72,12)):
    nt = fcsthrs.index(fcsthr)
    ax = fig.add_subplot(320+npanel+1)
    #cs = m.contour(x,y,t2m[nt,:,:],clevs,colors='k')
    cs = m.contourf(x,y,t2m[nt,:,:],clevs,cmap=cm.jet)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    m.drawparallels(arange(25,75,20),labels=[1,0,0,0],fontsize=8,fontstyle='oblique')
    m.drawmeridians(arange(-140,0,20),labels=[0,0,0,1],fontsize=8,yoffset=yoffset,fontstyle='oblique')
    # panel title
    title(repr(fcsthr)+'-h forecast valid '+verifdates[nt],fontsize=12)
# figure title
figtext(0.5,0.95,u"2-m temp (\N{DEGREE SIGN}K) forecasts from %s"%verifdates[0],
        horizontalalignment='center',fontsize=14)
# a single colorbar.
cax = axes([0.1, 0.03, 0.8, 0.025])
colorbar(cax=cax, orientation='horizontal')
show()
