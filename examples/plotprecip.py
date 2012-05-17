from __future__ import print_function
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset as NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
import copy
from matplotlib import rcParams

# make tick labels smaller
rcParams['xtick.labelsize']=9
rcParams['ytick.labelsize']=9

# plot rainfall from NWS using special precipitation
# colormap used by the NWS, and included in basemap.

nc = NetCDFFile('nws_precip_conus_20061222.nc')
# data from http://water.weather.gov/precip/
prcpvar = nc.variables['amountofprecip']
data = 0.01*prcpvar[:]
latcorners = nc.variables['lat'][:]
loncorners = -nc.variables['lon'][:]
plottitle = prcpvar.long_name+' for period ending '+prcpvar.dateofdata
print(data.min(), data.max())
print(latcorners)
print(loncorners)
print(plottitle)
print(data.shape)
lon_0 = -nc.variables['true_lon'].getValue()
lat_0 = nc.variables['true_lat'].getValue()
# create polar stereographic Basemap instance.
m = Basemap(projection='stere',lon_0=lon_0,lat_0=90.,lat_ts=lat_0,\
            llcrnrlat=latcorners[0],urcrnrlat=latcorners[2],\
            llcrnrlon=loncorners[0],urcrnrlon=loncorners[2],\
            rsphere=6371200.,resolution='l',area_thresh=10000)
# create figure
fig = plt.figure(figsize=(8.5,11))
plt.subplot(211)
ax = plt.gca()
# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.drawstates()
m.drawcountries()
# draw parallels.
delat = 10.0
parallels = np.arange(0.,90,delat)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
delon = 10.
meridians = np.arange(180.,360.,delon)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
ny = data.shape[0]; nx = data.shape[1]
lons, lats = m.makegrid(nx, ny) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.
# draw filled contours.
clevs = [0,1,2.5,5,7.5,10,15,20,30,40,50,70,100,150,200,250,300,400,500,600,750]
cs = m.contourf(x,y,data,clevs,cmap=cm.s3pcpn)
# draw colorbar.
cbar = m.colorbar(cs,location='bottom',pad="10%")
cbar.set_label('mm')
# plot title
plt.title(plottitle+'- contourf',fontsize=10)

plt.subplot(212)
ax = plt.gca()
# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.drawstates()
m.drawcountries()
# draw parallels.
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
# draw meridians
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
# draw image
im = m.imshow(data,cmap=cm.s3pcpn,interpolation='nearest',vmin=0,vmax=750)
# make a copy of the image object, change
# colormap to linear version of the precip colormap.
im2 = copy.copy(im)
im2.set_cmap(cm.s3pcpn_l)
# draw colorbar using im2, not im (hack to prevent colors from being
# too compressed at the low end on the colorbar - results
# from highly nonuniform colormap)
cb = m.colorbar(im2,location='bottom',pad="10%")
cb.set_label('mm')
# reset colorbar tick labels.
cb.set_ticks(np.linspace(clevs[0],clevs[-1],len(clevs)))
cb.set_ticklabels(['%g' % clev for clev in clevs])
# plot title
plt.title(plottitle+' - imshow',fontsize=10)
plt.show() # display onscreen.
