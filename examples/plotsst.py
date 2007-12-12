from matplotlib.toolkits.basemap import Basemap, NetCDFFile
import pylab, numpy
# read in sea-surface temperature data
# can be a local file, a URL for a remote opendap dataset,
# or (if PyNIO is installed) a GRIB or HDF file.
ncfile = NetCDFFile('http://nomads.ncdc.noaa.gov:8085/thredds/dodsC/oisst/2007/AVHRR/sst4-navy-eot.20071201.nc')
sst = ncfile.variables['sst'][:]
lats = ncfile.variables['lat'][:]
lons = ncfile.variables['lon'][:]
# create Basemap instance for mollweide projection.
# coastlines not used, so resolution set to None to skip
# continent processing (this speeds things up a bit)
m = Basemap(projection='moll',lon_0=lons.mean(),lat_0=0,resolution=None)
# compute map projection coordinates of grid.
x, y = m(*numpy.meshgrid(lons, lats))
m.drawmapboundary(fill_color='k')
# plot with pcolor
im = m.pcolormesh(x,y,sst,shading='flat',cmap=pylab.cm.gist_ncar)
# draw parallels and meridians, but don't bother labelling them.
m.drawparallels(numpy.arange(-90.,120.,30.))
m.drawmeridians(numpy.arange(0.,420.,60.))
# draw line around map projection limb.
# color map region background black (missing values will be this color)
# draw horizontal colorbar.
pylab.colorbar(orientation='horizontal')
pylab.show()
