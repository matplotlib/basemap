# example showing how to use Line Integral Convolution to visualize a vector
# flow field (from Hurricane Earl).  Produces something akin to streamlines.
# Requires vectorplot scikit (http://scikits.appspot.com/vectorplot).
from netCDF4 import Dataset as NetCDFFile
from mpl_toolkits.basemap import Basemap, interp
import numpy as np
import matplotlib.pyplot as plt

# H*wind data from http://www.aoml.noaa.gov/hrd/data_sub/wind.html
ncfile = NetCDFFile('rita.nc')
udat = ncfile.variables['sfc_u'][0,:,:]
vdat = ncfile.variables['sfc_v'][0,:,:]
lons1 = ncfile.variables['longitude'][:]
lats1 = ncfile.variables['latitude'][:]
lat0 = lats1[len(lats1)/2]; lon0 = lons1[len(lons1)/2]
lons, lats = np.meshgrid(lons1,lats1)
ncfile.close()

# downsample to finer grid for nicer looking plot.
nlats = 2*udat.shape[0]; nlons = 2*udat.shape[1]
lons = np.linspace(lons1[0],lons1[-1],nlons)
lats = np.linspace(lats1[0],lats1[-1],nlats)
lons, lats = np.meshgrid(lons, lats)
udat = interp(udat,lons1,lats1,lons,lats,order=3)
vdat = interp(vdat,lons1,lats1,lons,lats,order=3)
speed = np.sqrt(udat**2+vdat**2)


fig = plt.figure(figsize=(8,8))
m = Basemap(projection='cyl',llcrnrlat=lats1[0],llcrnrlon=lons1[0],urcrnrlat=lats1[-1],urcrnrlon=lons1[-1],resolution='i')
x, y = m(lons,lats)
plt.streamplot(x,y,udat,vdat,color=speed,linewidth=2,density=2,cmap=plt.cm.spectral)
m.colorbar()
m.drawcoastlines()
m.drawmeridians(np.arange(-120,-60,2),labels=[0,0,0,1])
m.drawparallels(np.arange(0,30,2),labels=[1,0,0,0])
plt.title('Hurricane Rita flow field visualized with streamlines',\
        fontsize=13)
plt.show()
