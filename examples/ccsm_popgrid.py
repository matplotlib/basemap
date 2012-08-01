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
"""
from matplotlib import rcParams
import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset as NetCDFFile

# read in data from netCDF file.
infile    = 'ccsm_popgrid.nc'
fpin      = NetCDFFile(infile)
tlat      = fpin.variables['TLAT'][:]
tlon      = fpin.variables['TLONG'][:]
# masked array returned, masked where data == _FillValue
temp      = fpin.variables['TEMP'][:]
fpin.close()

# make longitudes monotonically increasing.
tlon = np.where(np.greater_equal(tlon,min(tlon[:,0])),tlon-360,tlon)

# stack grids side-by-side (in longitiudinal direction), so
# any range of longitudes may be plotted on a world map.
tlon = np.concatenate((tlon,tlon+360),1)
tlat = np.concatenate((tlat,tlat),1)
temp = ma.concatenate((temp,temp),1)
tlon = tlon-360.

plt.figure(figsize=(6,8))
plt.subplot(2,1,1)
# subplot 1 just shows POP grid cells.
m = Basemap(projection='merc', lat_ts=20, llcrnrlon=-180, \
      urcrnrlon=180, llcrnrlat=-84, urcrnrlat=84, resolution='c')

m.drawcoastlines()
m.fillcontinents(color='white')

x, y = m(tlon,tlat)
im = m.pcolor(x,y,ma.masked_array(np.zeros(temp.shape,'f'), temp.mask),
                shading='faceted', antialiased=True, cmap=plt.cm.cool,
                vmin=0, vmax=0)
# disclaimer:  these are not really the grid cells because of the
# way pcolor interprets the x and y args.
plt.title('(A) CCSM POP Grid Cells')

# subplot 2 is a contour plot of surface temperature from the
# CCSM ocean model.
plt.subplot(2,1,2)
m.drawcoastlines()
m.fillcontinents(color='white')

CS1 = m.contourf(x,y,temp,15)
CS2 = m.contour(x,y,temp,15,colors='black',linewidths=0.5)
plt.title('(B) Surface Temp contours on POP Grid')

plt.show()
#plt.savefig('ccsm_popgrid.ps')

