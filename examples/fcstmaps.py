# this example reads today's numerical weather forecasts 
# from the NOAA OpenDAP servers and makes a multi-panel plot.
from pylab import title, show, figure, cm,  figtext, \
                  meshgrid, axes, colorbar
import numpy
import sys
from numpy import ma
import datetime
from matplotlib.toolkits.basemap import Basemap, NetCDFFile, addcyclic, num2date


# today's date is default.
if len(sys.argv) > 1:
    YYYYMMDD = sys.argv[1]
else:
   YYYYMMDD = datetime.datetime.today().strftime('%Y%m%d')
YYYYMM = YYYYMMDD[0:6]

# set OpenDAP server URL.
URLbase="http://nomad3.ncep.noaa.gov:9090/dods/mrf/mrf"
URL=URLbase+YYYYMMDD+'/mrf'+YYYYMMDD
print URL+'\n'
try:
    data = NetCDFFile(URL)
except:
    msg = """
opendap server not providing the requested data.
Try another date by providing YYYYMM on command line."""
    raise IOError, msg


# read lats,lons,times.

print data.variables.keys()
latitudes = data.variables['lat']
longitudes = data.variables['lon']
fcsttimes = data.variables['time']
times = fcsttimes[0:6] # first 6 forecast times.
ntimes = len(times)
# put forecast times in YYYYMMDDHH format.
verifdates = []
fcsthrs=[]
for time in times:
    print time, times[0]
    fcsthrs.append(int((time-times[0])*24))
    fdate = num2date(time,fcsttimes.units)
    verifdates.append(fdate.strftime('%Y%m%d%H'))
print fcsthrs
print verifdates
lats = latitudes[:]
nlats = len(lats)
lons1 = longitudes[:]
nlons = len(lons1)

# unpack 2-meter temp forecast data.

t2mvar = data.variables['tmp2m']
t2min = t2mvar[0:ntimes,:,:]
t2m = numpy.zeros((ntimes,nlats,nlons+1),t2min.dtype)
# create Basemap instance for Orthographic projection.
m = Basemap(lon_0=-90,lat_0=60,projection='ortho')
# add wrap-around point in longitude.
for nt in range(ntimes):
    t2m[nt,:,:], lons = addcyclic(t2min[nt,:,:], lons1)
# convert to celsius.
t2m = t2m-273.15
# contour levels
clevs = numpy.arange(-30,30.1,2.)
lons, lats = meshgrid(lons, lats)
x, y = m(lons, lats)
# create figure.
fig=figure(figsize=(6,8))
# make subplots.
for nt,fcsthr in enumerate(fcsthrs):
    ax = fig.add_subplot(321+nt)
    cs = m.contourf(x,y,t2m[nt,:,:],clevs,cmap=cm.jet,extend='both')
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries()
    m.drawparallels(numpy.arange(-80,81,20))
    m.drawmeridians(numpy.arange(0,360,20))
    # panel title
    title(repr(fcsthr)+'-h forecast valid '+verifdates[nt],fontsize=9)
# figure title
figtext(0.5,0.95,u"2-m temp (\N{DEGREE SIGN}C) forecasts from %s"%verifdates[0],
        horizontalalignment='center',fontsize=14)
# a single colorbar.
cax = axes([0.1, 0.03, 0.8, 0.025])
colorbar(cax=cax, orientation='horizontal')
show()
