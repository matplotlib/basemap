from netCDF4 import Dataset
import time, calendar, datetime, numpy
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
def datetomsecs(d):
    """convert from datetime to msecs since the unix epoch began"""
    return int(calendar.timegm(time.struct_time(d.timetuple()))*1000)
# set date range
date1 = datetime.datetime(2010,1,1,0)
date2 = datetime.datetime(2010,1,8,0)
t1 = datetomsecs(date1); t2 = datetomsecs(date2)
# build constraint expression to get locations of floats in specified time
# range.
urlbase='http://dapper.pmel.noaa.gov/dapper/argo/argo_all.cdp'
sel="?location.JULD,location.LATITUDE,location.LONGITUDE&location.JULD>%s&location.JULD<%s"%(t1,t2)
# retrieve data
dset = Dataset(urlbase+sel)
lats = dset.variables['location.LATITUDE'][:]
lons = dset.variables['location.LONGITUDE'][:]
# draw map with markers for float locations
m = Basemap(projection='hammer',lon_0=180)
x, y = m(lons,lats)
m.drawmapboundary(fill_color='#99ffff')
m.fillcontinents(color='#cc9966',lake_color='#99ffff')
m.scatter(x,y,3,marker='o',color='k')
plt.title('Locations of %s ARGO floats active between %s and %s' %\
        (len(lats),date1.strftime('%Y%m%d'),date2.strftime('%Y%m%d')))
plt.show()
