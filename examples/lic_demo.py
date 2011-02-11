# example showing how to use Line Integral Convolution to visualize a vector
# flow field (from Hurricane Earl).  Produces something akin to streamlines.
# Requires vectorplot scikit (http://scikits.appspot.com/vectorplot).
try:
    from netCDF4 import Dataset as NetCDFFile
except ImportError:
    from mpl_toolkits.basemap import NetCDFFile
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
try:
    from vectorplot import lic_internal
except ImportError:
    raise  ImportError('need vectorplot scikit for this example')

# H*wind data from http://www.aoml.noaa.gov/hrd/data_sub/wind.html
ncfile = NetCDFFile('rita.nc')
u = ncfile.variables['sfc_u'][0,:,:]
v = ncfile.variables['sfc_v'][0,:,:]
lons1 = ncfile.variables['longitude'][:]
lats1 = ncfile.variables['latitude'][:]
lat0 = lats1[len(lats1)/2]; lon0 = lons1[len(lons1)/2]
print lat0,lon0
lons, lats = np.meshgrid(lons1,lats1)
ncfile.close()

fig = plt.figure(figsize=(8,8))
m = Basemap(projection='stere',lat_0=lat0,lon_0=lon0,width=1.e6,height=1.e6,resolution='i')
nxv = 301; nyv = 301; kernellen=41
udat, vdat, xv, yv = m.transform_vector(u,v,lons1,lats1,nxv,nyv,returnxy=True)
texture = np.random.rand(udat.shape[0],udat.shape[1]).astype(np.float32)
kernel = np.sin(np.arange(kernellen)*np.pi/kernellen).astype(np.float32)
image = lic_internal.line_integral_convolution(udat.astype(np.float32),\
        vdat.astype(np.float32), texture, kernel)
im = m.imshow(image,plt.cm.gist_stern)
m.drawcoastlines()
m.drawmeridians(np.arange(0,360,2),labels=[0,0,0,1])
m.drawparallels(np.arange(-30,30,2),labels=[1,0,0,0])
plt.title('Hurricane Rita flow field visualized with Line Integral Convolution',\
        fontsize=13)
plt.show()
