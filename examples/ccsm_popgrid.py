"""
This example shows how to plot data on rectangular 2D grids
(grids that are not rectlinear in geographic or native map projection
coordinates).

An example of such a grid is the 'POP' grid which is used in
the ocean component NCAR Community Climate System Model (CCSM).
"POP" stands for "Parallel Ocean Program", which was developed
at Los Alamos.

These grids may be thought of as rectangular arrays wrapped around the
globe in the usual way, with one subscript, call it I, associated with
longitude and the other subscript, call it J, associated with latitude,
and then deformed in such a way as to move the top edge of the array to
a circle centered somewhere other than over the North Pole (typically,
over Greenland or Canada) and the bottom edge of the array to a circle
that is centered on the South Pole, but lies entirely within Antarctica.
The lines defined by the rows and columns of the rectangular arrays are
locally orthogonal to each other.

POP grids are used extensively locally in oceanographic and ice models.

This example uses the pupynere pure-python netCDF 3 reader
from Roberto De Almeida for reading netCDF files (included).

"""
from pupynere import NetCDFFile
import pylab as pl
from matplotlib import rcParams
try:
    import numpy
    rcParams['numerix'] = 'numpy'
except:
    raise ImportError('this example requires numpy')
import matplotlib.numerix.ma as MA
import matplotlib.numerix as N
from matplotlib.toolkits.basemap import Basemap

# read in data from netCDF file.
infile    = 'ccsm_popgrid.nc'
fpin      = NetCDFFile(infile)
tlat      = fpin.variables['TLAT'][:]
tlon      = fpin.variables['TLONG'][:]
temp      = fpin.variables['TEMP'][:]
fillvalue = fpin.variables['TEMP'].attributes['_FillValue']
fpin.close()

# make longitudes monotonically increasing.
tlon = N.where(N.greater_equal(tlon,min(tlon[:,0])),tlon-360,tlon)

# create a masked array with temperature data (continents masked).
temp = MA.masked_values(temp,fillvalue)

# stack grids side-by-side (in longitiudinal direction), so
# any range of longitudes may be plotted on a world map.
tlon = N.concatenate((tlon,tlon+360),1)
tlat = N.concatenate((tlat,tlat),1)
temp = MA.concatenate((temp,temp),1)
tlon = tlon-360.

pl.figure(figsize=(8.5,11))
pl.subplot(2,1,1)
# subplot 1 just shows POP grid cells.
map = Basemap(projection='merc', lat_ts=20, llcrnrlon=-180, \
      urcrnrlon=180, llcrnrlat=-84, urcrnrlat=84, resolution='c')

map.drawcoastlines()
map.fillcontinents(color='white')

x, y = map(tlon,tlat)
im = map.pcolor(x,y,MA.masked_array(N.zeros(temp.shape,'f'), temp.mask),\
                shading='faceted',cmap=pl.cm.cool,vmin=0,vmax=0)
# disclaimer:  these are not really the grid cells because of the
# way pcolor interprets the x and y args.
pl.title('(A) CCSM POP Grid Cells')

# subplot 2 is a contour plot of surface temperature from the
# CCSM ocean model.
pl.subplot(2,1,2)
map.drawcoastlines()
map.fillcontinents(color='white')

CS1 = map.contourf(x,y,temp,15)
CS2 = map.contour(x,y,temp,15,colors='black',linewidths=0.5)
pl.title('(B) Surface Temp contours on POP Grid')

pl.show()
#pl.savefig('ccsm_popgrid.ps')

