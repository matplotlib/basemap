from mpl_toolkits.basemap import Basemap, NetCDFFile
import numpy as np
import matplotlib.pyplot as plt
import sys
# read in sea-surface temperature and ice data
# can be a local file, a URL for a remote opendap dataset,
# or (if PyNIO is installed) a GRIB or HDF file.
if len(sys.argv) == 1:
    date = '20071215'
else:
    date = sys.argv[1]
if date[0:4] > '2005':
    ncfile = NetCDFFile('http://nomads.ncdc.noaa.gov:8085/thredds/dodsC/oisst/'+date[0:4]+'/AVHRR/sst4-navy-eot.'+date+'.nc')
else:
    ncfile = NetCDFFile('http://nomads.ncdc.noaa.gov:8085/thredds/dodsC/oisst/'+date[0:4]+'/AVHRR/sst4-path-eot.'+date+'.nc')
# read sst.  Will automatically create a masked array using
# missing_value variable attribute.
sst = ncfile.variables['sst'][:]
# read ice.
ice = ncfile.variables['ice'][:]
# read lats and lons (representing centers of grid boxes).
lats = ncfile.variables['lat'][:]
lons = ncfile.variables['lon'][:]
# shift lats, lons so values represent edges of grid boxes
# (as pcolor expects).
delon = lons[1]-lons[0]
delat = lats[1]-lats[0]
lons = (lons - 0.5*delon).tolist()
lons.append(lons[-1]+delon)
lons = np.array(lons,np.float64)
lats = (lats - 0.5*delat).tolist()
lats.append(lats[-1]+delat)
lats = np.array(lats,np.float64)
# create Basemap instance for mollweide projection.
# coastlines not used, so resolution set to None to skip
# continent processing (this speeds things up a bit)
#m = Basemap(projection='ortho',lon_0=-110,lat_0=20,resolution=None)
m = Basemap(projection='moll',lon_0=lons.mean(),lat_0=0,resolution=None)
# compute map projection coordinates of grid.
x, y = m(*np.meshgrid(lons, lats))
# draw line around map projection limb.
# color background of map projection region.
# missing values over land will show up this color.
m.drawmapboundary(fill_color='0.3')
# plot ice, then with pcolor
im1 = m.pcolor(x,y,sst,shading='flat',cmap=plt.cm.jet)
im2 = m.pcolor(x,y,ice,shading='flat',cmap=plt.cm.gist_gray)
# draw parallels and meridians, but don't bother labelling them.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
# draw horizontal colorbar.
plt.colorbar(im1,orientation='horizontal')
# display the plot with a title.
plt.title('SST and ICE analysis for %s'%date)
plt.show()
