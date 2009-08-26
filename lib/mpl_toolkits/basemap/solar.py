# some simple functions to calculate solar position, day-night terminator
import numpy as np
from numpy import ma
import netcdftime

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
    ut = date.hour + date.minute/60. + date.second/3600.
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

def daynight_terminator(date, delta, lonmin, lonmax):
    """
    date is datetime object (assumed UTC).
    nlons is # of longitudes used to compute terminator."""
    dg2rad = np.pi/180.
    lons = np.arange(lonmin,lonmax+0.5*delta,delta,dtype=np.float32)
    # compute greenwich hour angle and solar declination
    # from datetime object (assumed UTC).
    tau, dec = epem(date) 
    # compute day/night terminator from hour angle, declination.
    longitude = lons + tau
    lats = np.arctan(-np.cos(longitude*dg2rad)/np.tan(dec*dg2rad))/dg2rad
    return lons, lats, tau, dec

def daynight_grid(date, delta, lonmin, lonmax):
    """
    date is datetime object (assumed UTC).
    delta is the grid interval (in degrees) used to compute terminator."""
    lons, lats, tau, dec = daynight_terminator(date, delta, lonmin, lonmax)
    # create day/night grid (1 for night, 0 for day)
    lats2 = np.arange(-90,90+0.5*delta,delta,dtype=np.float32)
    nlons = len(lons); nlats = len(lats2)
    lons2, lats2 = np.meshgrid(lons,lats2)
    lats = lats[np.newaxis,:]*np.ones((nlats,nlons),dtype=np.float32)
    daynight = np.ones(lons2.shape, np.int8)
    if dec > 0: # NH summer
        daynight = np.where(lats2>lats,0,daynight)
    else: # NH winter
        daynight = np.where(lats2<lats,0,daynight)
    daynight = ma.array(daynight,mask=1-daynight) # mask day areas.
    return lons2,lats2,daynight
