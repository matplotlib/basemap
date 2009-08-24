import numpy as np
from mpl_toolkits.basemap import Basemap, netcdftime
import matplotlib.pyplot as plt
from datetime import datetime

# example showing how to compute the day/night terminator and shade nightime
# areas on a map.

def epem(date):
    """
    input: date - datetime object (assumed UTC)
    ouput: gha - Greenwich hour angle, the angle between the Greenwich
           meridian and the meridian containing the subsolar point.
           dec - solar declination.
    """
    dg2rad = np.pi/180.
    rad2dg = 1./dg2rad
    # compute julian day from UTC datetime object.
    jday = netcdftime.JulianDayFromDate(date)
    jd = np.floor(jday) # truncate to integer.
    # utc hour.
    ut = d.hour + d.minute/60. + d.second/3600.
    # calculate number of centuries from J2000
    t = (jd + (ut/24.) - 2451545.0) / 36525.
    # mean longitude corrected for aberration
    l = (280.460 + 36000.770 * t) % 360
    # mean anomaly
    g = 357.528 + 35999.050 * t
    # ecliptic longitude
    lm = l + 1.915 * np.sin(g*dg2rad) + 0.020 * np.sin(2*g*dg2rad)
    # obliquity of the ecliptic
    ep = 23.4393 - 0.01300 * t
    # equation of time
    eqtime = -1.915*np.sin(g*dg2rad) - 0.020*np.sin(2*g*dg2rad) \
            + 2.466*np.sin(2*lm*dg2rad) - 0.053*np.sin(4*lm*dg2rad)
    # Greenwich hour angle
    gha = 15*ut - 180 + eqtime
    # declination of sun
    dec = np.arcsin(np.sin(ep*dg2rad) * np.sin(lm*dg2rad)) * rad2dg
    return gha, dec

def daynightgrid(date, nlons):
    """
    date is datetime object (assumed UTC).
    nlons is # of longitudes used to compute terminator."""
    nlats = ((nlons-1)/2)+1
    dg2rad = np.pi/180.
    lons = np.linspace(-180,180,nlons)
    # compute greenwich hour angle and solar declination
    # from datetime object (assumed UTC).
    tau, dec = epem(date) 
    longitude = lons + tau
    lats = np.arctan(-np.cos(longitude*dg2rad)/np.tan(dec*dg2rad))/dg2rad
    lons2 = np.linspace(-180,180,nlons)
    lats2 = np.linspace(-90,90,nlats)
    lons2, lats2 = np.meshgrid(lons2,lats2)
    daynight = np.ones(lons2.shape, np.float)
    for nlon in range(nlons):
        daynight[:,nlon] = np.where(lats2[:,nlon]>lats[nlon],0,daynight[:,nlon])
    return lons2,lats2,daynight

# current time in UTC.
d = datetime.utcnow()

# miller projection 
map = Basemap(projection='mill',lon_0=0)
# plot coastlines, draw label meridians and parallels.
map.drawcoastlines()
map.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0])
map.drawmeridians(np.arange(-180,180,60),labels=[0,0,0,1])
# fill continents 'coral' (with zorder=0), color wet areas 'aqua'
map.drawmapboundary(fill_color='aqua')
map.fillcontinents(color='coral',lake_color='aqua',zorder=0)
# create grid of day=0, night=1
# 1441 means use a 0.25 degree grid.
lons,lats,daynight = daynightgrid(d,1441)
x,y = map(lons, lats)
# contour this grid with 1 contour level, specifying colors.
# (gray for night, white for day). Use alpha transparency so
# map shows through.
CS=map.contourf(x,y,daynight,1,colors=['white','0.7'],alpha=0.5)
plt.title('Day/Night Map for %s (UTC)' %d )
plt.show()
