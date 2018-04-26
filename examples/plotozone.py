from __future__ import (absolute_import, division, print_function)

# make plot of ozone concentration data on
# lambert conformal conic map projection, drawing coastlines, state and
# country boundaries, and parallels/meridians.

# the data is interpolated to the native projection grid.
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
plt.rcParams['text.usetex'] = False

# read in netCDF4 file. Results from CAMx v6
# test case, converted to netcdf by PseudoNetCDF
# pseudonetcdf.googlecode.com
camx = netCDF4.Dataset('camx.sample.nc')

#alternatively read directly from CAMx uamiv file
#if available
#
# from PseudoNetCDF.camxfiles.Memmaps import uamiv
# camx = uamiv('camx.bin')

# Get Ozone Variable
o3 = camx.variables['O3']

# Get projection space
llcrnrx = camx.XORIG
llcrnry = camx.YORIG
urcrnrx = llcrnrx + (o3[:].shape[-1] + 1) * camx.XCELL
urcrnry = llcrnry + (o3[:].shape[-2] + 1) * camx.XCELL

# Get edge values for pcolor
xedge = np.linspace(0, urcrnrx - llcrnrx, camx.NCOLS + 1)
yedge = np.linspace(0, urcrnry - llcrnry, camx.NCOLS + 1)
X, Y = np.meshgrid(xedge, yedge)


# setup of basemap ('lcc' = lambert conformal conic).
# projection parameters from CAMx file
m = Basemap(projection = 'lcc',
            lon_0=camx.P_GAM, lat_0 = 40.,
            lat_1 = camx.P_ALP, lat_2 = camx.P_BET,
            llcrnrx = llcrnrx, llcrnry = llcrnry,
            urcrnry = urcrnry, urcrnrx = urcrnrx)

# create the figure.
fig=plt.figure(figsize=(8,8))

# add an axes.
ax = fig.add_axes([0.1,0.1,0.8,0.8])
ax.set_facecolor('lightgrey')
# associate this axes with the Basemap instance.
m.ax = ax

# plot tile plot with pcolor
# Use first time and first layer (i.e., o3[0, 0] (time, layer, row, col))
# Edge cells have precisely 0 value, and are masked
# to avoid an unnecessary color range.
# Each color bin contains precisely 10% of values
# which makes for a pretty plot.
from matplotlib.colors import ListedColormap
WhGrYlBu = ListedColormap(['#ffffff', '#b7f6ff', '#70edff', '#29e4ff', '#00e1fb', '#0fffc6', '#3bffa4', '#68ff82', '#94ff60', '#c0ff3e', '#edff1c', '#fff400', '#ffc700', '#ff9b00', '#ff6e00', '#ff4200', '#ff1500', '#e80000', '#bb0000', '#8f0000'])
#.from_list('WhGrYlBu', ['white', 'white', 'cyan', 'lightblue', 'lightgreen', 'green', 'yellow', 'orange', 'red', 'red'])

toplot = np.ma.masked_values(o3[0, 0], 0.) * 1000.
bounds = np.percentile(toplot.compressed().ravel(), np.linspace(5, 95, 9).tolist())
ptch = m.pcolor(X, Y, toplot, cmap = WhGrYlBu, norm = plt.matplotlib.colors.BoundaryNorm(bounds, 20), vmin = bounds[0], vmax = bounds[-1])

# Add a colorbar using proportional spacing, but
# colors based on 10 distinct bins
cb = m.colorbar(ptch, location='right',pad='10%', boundaries = bounds, spacing = 'proportional', format = '%.3f', extend = 'both') # draw colorbar

# Add units to the colorbar
cb.ax.set_xlabel('%s*1000.' % o3.units.strip())


# plot blue dot on Houston, Baton Rouge, and Atlanta
def add_dot(lon, lat, label):
    xpt,ypt = m(lon,lat) 
    m.plot([xpt],[ypt],'bo') 
    ax.annotate(label, xy=(xpt, ypt), xytext=(xpt+1e5, ypt+1e5),
            bbox=dict(boxstyle="round4", fc="w"),
            arrowprops=dict(facecolor='black'),
            )
add_dot(-95.361328,29.754505, 'Houston') 
add_dot(-91.140320, 30.458283, 'Baton Rouge')
add_dot(-84.387982, 33.748995, 'Atlanta')
# draw coastlines and political boundaries.
m.drawcoastlines()
m.drawcountries()
m.drawstates()
# draw parallels and meridians.
# label on left, right and bottom of map.
parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,labels=[1,1,0,1])
meridians = np.arange(-120., 70.,10.)
m.drawmeridians(meridians,labels=[1,1,0,1])

# set title.
ax.set_title('O$_3$ as predicted by the CAMx v6 Test-Case\neach color division has 10% of cells 5-95% and 5% in triagles')
import textwrap
histstr = 'Processing: %s' % '\n'.join(textwrap.wrap(camx.history.strip(), 140))

fig.text(0.01, 0.01, histstr, horizontalalignment = 'left', verticalalignment = 'bottom', size = 8)
plt.draw()
plt.show()
