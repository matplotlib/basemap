from __future__ import print_function
from __future__ import unicode_literals
# this example reads today's numerical weather forecasts
# from the NOAA OpenDAP servers and makes a multi-panel plot.
# This version demonstrates the use of the AxesGrid toolkit.
import numpy as np
import matplotlib.pyplot as plt
import sys
import numpy.ma as ma
import datetime
from mpl_toolkits.basemap import Basemap, addcyclic
from mpl_toolkits.axes_grid1 import AxesGrid
from netCDF4 import Dataset as NetCDFFile, num2date


# today's date is default.
if len(sys.argv) > 1:
    YYYYMMDD = sys.argv[1]
else:
    YYYYMMDD = datetime.datetime.today().strftime('%Y%m%d')

# set OpenDAP server URL.
try:
    URLbase="http://nomads.ncep.noaa.gov:9090/dods/gfs/gfs"
    URL=URLbase+YYYYMMDD+'/gfs_00z'
    print(URL)
    data = NetCDFFile(URL)
except:
    msg = """
opendap server not providing the requested data.
Try another date by providing YYYYMMDD on command line."""
    raise IOError(msg)


# read lats,lons,times.

print(data.variables.keys())
latitudes = data.variables['lat']
longitudes = data.variables['lon']
fcsttimes = data.variables['time']
times = fcsttimes[0:6] # first 6 forecast times.
ntimes = len(times)
# convert times for datetime instances.
fdates = num2date(times,units=fcsttimes.units,calendar='standard')
# make a list of YYYYMMDDHH strings.
verifdates = [fdate.strftime('%Y%m%d%H') for fdate in fdates]
# convert times to forecast hours.
fcsthrs = []
for fdate in fdates:
    fdiff = fdate-fdates[0]
    fcsthrs.append(fdiff.days*24. + fdiff.seconds/3600.)
print(fcsthrs)
print(verifdates)
lats = latitudes[:]
nlats = len(lats)
lons1 = longitudes[:]
nlons = len(lons1)

# unpack 2-meter temp forecast data.

t2mvar = data.variables['tmp2m']

# create figure, set up AxesGrid.
fig=plt.figure(figsize=(6,8))
grid = AxesGrid(fig, [0.05,0.01,0.9,0.9],
                nrows_ncols=(3, 2),
		axes_pad=0.25,
		cbar_mode='single',
		cbar_pad=0.3,
		cbar_size=0.1,
                cbar_location='top',
                share_all=True,
		)

# create Basemap instance for Orthographic projection.
m = Basemap(lon_0=-90,lat_0=60,projection='ortho')
# add wrap-around point in longitude.
t2m = np.zeros((ntimes,nlats,nlons+1),np.float32)
for nt in range(ntimes):
    t2m[nt,:,:], lons = addcyclic(t2mvar[nt,:,:], lons1)
# convert to celsius.
t2m = t2m-273.15
# contour levels
clevs = np.arange(-30,30.1,2.)
lons, lats = np.meshgrid(lons, lats)
x, y = m(lons, lats)
# make subplots.
for nt,fcsthr in enumerate(fcsthrs):
    ax = grid[nt]
    m.ax = ax
    cs = m.contourf(x,y,t2m[nt,:,:],clevs,cmap=plt.cm.jet,extend='both')
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries()
    m.drawparallels(np.arange(-80,81,20))
    m.drawmeridians(np.arange(0,360,20))
    # panel title
    ax.set_title('%d-h forecast valid '%fcsthr+verifdates[nt],fontsize=9)
# figure title
plt.figtext(0.5,0.95,
            "2-m temp (\N{DEGREE SIGN}C) forecasts from %s"%verifdates[0],
            horizontalalignment='center',fontsize=14)
# a single colorbar.
cbar = fig.colorbar(cs, cax=grid.cbar_axes[0], orientation='horizontal')
plt.show()
