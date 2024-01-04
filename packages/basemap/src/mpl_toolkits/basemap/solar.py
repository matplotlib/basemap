"""Simple functions to calculate solar position and day-night terminator."""
# pylint: disable=invalid-name
from __future__ import division

import numpy as np


def JulianDayFromDate(date, calendar="standard"):
    """Return Julian day from a :class:`datetime.datetime` object.

    Algorithm:

        Meeus, Jean (1998) Astronomical Algorithms (2nd Edition).
        Willmann-Bell, Virginia. p. 63

    Paramaters
    ---------

    date : datetime.datetime
        a :class:`datetime.datetime` object

    calendar : {'standard', 'gregorian', 'proleptic_gregorian',
                'julian'}, optional
        if 'standard' or 'gregorian', Julian day follows the Julian
        calendar on and before 1582-10-05, and the Gregorian calendar
        after 1582-10-15; if 'proleptic_gregorian', Julian day follows
        the Gregorian calendar; if 'julian', Julian day follows the
        Julian calendar

    Returns
    -------

    jd : float
        the Julian day, with resolution of 1 second
    """

    # Based on `redate.py` by David Finlayson.
    year, month, day = date.year, date.month, date.day
    hour, minute, second = date.hour, date.minute, date.second

    # Convert time to fractions of a day.
    day = day + hour / 24.0 + minute / 1440.0 + second / 86400.0

    # Start Meeus algorithm (variables are in his notation).
    if month < 3:
        month = month + 12
        year = year - 1
    A = int(year / 100)
    jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + day - 1524.5

    # Optionally adjust `jd` for the switch from Julian to Gregorian calendar.
    # Here assumed to have occurred the day after 1582 October 4
    if calendar in ["standard", "gregorian"]:
        if jd >= 2299170.5:
            # 1582 October 15 (Gregorian Calendar).
            B = 2 - A + int(A / 4)
        elif jd < 2299160.5:
            # 1582 October 5 (Julian Calendar).
            B = 0
        else:
            raise ValueError("impossible date (falls in gap between end of "
                             "Julian calendar and start of Gregorian calendar")
    elif calendar == "proleptic_gregorian":
        B = 2 - A + int(A / 4)
    elif calendar == "julian":
        B = 0
    else:
        raise ValueError("unknown calendar '{0}' (must be one of 'standard', "
                         "'gregorian', 'proleptic_gregorian' or 'julian'".format(calendar))

    # Adjust for Julian calendar if necessary.
    jd = jd + B

    return jd


def epem(date):
    """Return the Greenwich hour angle.

    Parameters
    ----------

    date : datetime.datetime
        a :class:`datetime.datetime` object (assumed UTC)

    Returns
    -------

    gha : float
        Greenwich hour angle, i.e. the angle between the Greenwich
        meridian and the meridian containing the subsolar point.

    dec : float
        solar declination
    """

    dg2rad = np.pi / 180.
    rad2dg = 1. / dg2rad

    # Compute Julian day from UTC `datetime.datetime` object (note that
    # `datetime.datetime` objects use proleptic Gregorian calendar).
    jday = JulianDayFromDate(date, calendar="proleptic_gregorian")
    jd = np.floor(jday)  # Truncate to integer.

    # Compute UTC hour.
    ut = date.hour + date.minute / 60. + date.second / 3600.

    # Calculate number of centuries from J2000.
    t = (jd + (ut / 24.) - 2451545.0) / 36525.

    # Mean longitude corrected for aberration.
    l = (280.460 + 36000.770 * t) % 360  # noqa: E741

    # Mean anomaly.
    g = 357.528 + 35999.050 * t

    # Ecliptic longitude.
    lm = l + 1.915 * np.sin(g * dg2rad) + 0.020 * np.sin(2 * g * dg2rad)

    # Obliquity of the ecliptic.
    ep = 23.4393 - 0.01300 * t

    # Equation of time.
    eqtime = (-1.915 * np.sin(g * dg2rad) - 0.020 * np.sin(2 * g * dg2rad) +
              +2.466 * np.sin(2 * lm * dg2rad) - 0.053 * np.sin(4 * lm * dg2rad))

    # Greenwich hour angle.
    gha = 15 * ut - 180 + eqtime

    # Declination of Sun.
    dec = np.arcsin(np.sin(ep * dg2rad) * np.sin(lm * dg2rad)) * rad2dg

    return gha, dec


def daynight_terminator(date, delta, lonmin, lonmax):
    """Return the day-night terminator.

    Parameters
    ----------

    date : datetime.datetime
        a :class:`datetime.datetime` object (assumed UTC)

    delta : float
        input longitude grid step in degrees used to compute the
        day-night terminator

    lonmin : float
        minimum input longitude in degrees

    lonmax : float
        maximum input longitude in degrees

    Returns
    -------

    lons : array-like
        array of longitudes of the day-night terminator

    lats : array-like
        array of latitudes of the day-night terminator

    tau : float
        Greenwich hour angle for the input datetime

    dec : float
        solar declination for the input datetime
    """

    dg2rad = np.pi / 180.
    lons = np.arange(lonmin, lonmax + 0.5 * delta, delta, dtype=np.float32)

    # Compute Greenwich hour angle and solar declination.
    tau, dec = epem(date)

    # compute day-night terminator from Greenwich hour angle and declination.
    longitude = lons + tau
    lats = np.arctan(-np.cos(longitude * dg2rad) / np.tan(dec * dg2rad)) / dg2rad

    return lons, lats, tau, dec


def daynight_grid(date, delta, lonmin, lonmax):
    """Return day-night mask array.

    Parameters
    ----------

    date : datetime.datetime
        a :class:`datetime.datetime` object (assumed UTC)

    delta : float
        input longitude grid step in degrees used to compute the
        day-night terminator

    lonmin : float
        minimum input longitude in degrees

    lonmax : float
        maximum input longitude in degrees

    Returns
    -------

    lons2 : array-like
        meshgrid of longitudes of the day-night mask array

    lats2 : array-like
        meshgrid of latitudes of the day-night mask array

    daynight : ~numpy.ma.MaskedArray
        day-night mask array (masked for day, 1 for night)
    """

    lons, lats, _, dec = daynight_terminator(date, delta, lonmin, lonmax)

    # Create meshgrid of longitudes and latitudes.
    lats2 = np.arange(-90, 90 + 0.5 * delta, delta, dtype=np.float32)
    lons2, lats2 = np.meshgrid(lons, lats2)

    # Create day-night grid (0 for day, 1 for night).
    nlats, nlons = len(lats2), len(lons)
    lats = lats[np.newaxis, :] * np.ones((nlats, nlons), dtype=np.float32)
    daynight = np.ones(lons2.shape, np.int8)
    if dec > 0:  # NH summer
        daynight = np.where(lats2 > lats, 0, daynight)
    else:        # NH winter
        daynight = np.where(lats2 < lats, 0, daynight)

    # Create day-night masked array (with day areas masked).
    daynight_mask = 1 - daynight
    daynight = np.ma.array(daynight, mask=daynight_mask)

    return lons2, lats2, daynight
