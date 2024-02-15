from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset, date2index
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve
date = datetime(2007,12,15,0) # date to plot.
# open dataset.
sstpath, sstheader = urlretrieve("https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/sst.day.mean.{0}.nc".format(date.year))
dataset = Dataset(sstpath)
timevar = dataset.variables['time']
timeindex = date2index(date,timevar) # find time index for desired date.
# read sst.  Will automatically create a masked array using
# missing_value variable attribute. 'squeeze out' singleton dimensions.
sst = dataset.variables['sst'][timeindex,:].squeeze()
# read ice.
dataset.close()
icepath, iceheader = urlretrieve("https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/icec.day.mean.{0}.nc".format(date.year))
dataset = Dataset(icepath)
ice = dataset.variables['icec'][timeindex,:].squeeze()
# read lats and lons (representing centers of grid boxes).
lats = dataset.variables['lat'][:]
lons = dataset.variables['lon'][:]
dataset.close()
latstep, lonstep = np.diff(lats[:2]), np.diff(lons[:2])
lats = np.append(lats - 0.5 * latstep, lats[-1] + 0.5 * latstep)
lons = np.append(lons - 0.5 * lonstep, lons[-1] + 0.5 * lonstep)
lons, lats = np.meshgrid(lons,lats)
# create figure, axes instances.
fig = plt.figure()
ax = fig.add_axes([0.05,0.05,0.9,0.9])
# create Basemap instance.
# coastlines not used, so resolution set to None to skip
# continent processing (this speeds things up a bit)
m = Basemap(projection='kav7',lon_0=0,resolution=None)
# draw line around map projection limb.
# color background of map projection region.
# missing values over land will show up this color.
m.drawmapboundary(fill_color='0.3')
# plot sst, then ice with pcolor
im1 = m.pcolormesh(lons,lats,sst,shading='flat',cmap=plt.cm.jet,latlon=True)
im2 = m.pcolormesh(lons,lats,ice,shading='flat',cmap=plt.cm.gist_gray,latlon=True)
# draw parallels and meridians, but don't bother labelling them.
m.drawparallels(np.arange(-90.,99.,30.))
m.drawmeridians(np.arange(-180.,180.,60.))
# add colorbar
cb = m.colorbar(im1,"bottom", size="5%", pad="2%")
# add a title.
ax.set_title('SST and ICE analysis for %s'%date)
plt.show()
