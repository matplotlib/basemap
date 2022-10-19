"""Make a multi-panel plot from numerical weather forecast in NOAA OPeNDAP.

This version demonstrates the use of the AxesGrid toolkit.
"""
from __future__ import print_function

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import addcyclic
from mpl_toolkits.axes_grid1 import AxesGrid


def main(date, verbose=True):
    """Main function."""

    # Open dataset from OPeNDAP URL.
    url = "http://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs%Y%m%d/gfs_0p25_00z"
    try:
        data = netCDF4.Dataset(date.strftime(url), "r")
        if verbose:
            print("Data variables:")
            print(sorted(data.variables))
    except OSError as err:
        err.args = (err.args[0], "date not found in OPeNDAP server")
        raise

    # Read lats, lons, and times.
    latitudes = data.variables["lat"]
    longitudes = data.variables["lon"]
    fcsttimes = data.variables["time"]
    times = fcsttimes[0:6]  # First 6 forecast times.
    ntimes = len(times)

    # Convert times for datetime instances.
    fdates = netCDF4.num2date(
        times, units=fcsttimes.units, calendar="standard")

    # Make a list of YYYYMMDDHH strings.
    verifdates = [fdate.strftime("%Y%m%d%H") for fdate in fdates]
    if verbose:
        print("Forecast datetime strings:")
        print(verifdates)

    # Convert times to forecast hours.
    fcsthrs = []
    for fdate in fdates:
        fdiff = fdate - fdates[0]
        fcsthrs.append(fdiff.days * 24. + fdiff.seconds / 3600.)
    if verbose:
        print("Forecast datetime hours:")
        print(fcsthrs)

    # Unpack 2-meter temp forecast data.
    lats = latitudes[:]
    nlats = len(lats)
    lons1 = longitudes[:]
    nlons = len(lons1)
    t2mvar = data.variables["tmp2m"]

    # Create Basemap instance for orthographic projection.
    bmap = Basemap(lon_0=-90, lat_0=60, projection="ortho")

    # Add wrap-around point in longitude.
    t2m = np.zeros((ntimes, nlats, nlons + 1), np.float32)
    for nt in range(ntimes):
        t2m[nt, :, :], lons = addcyclic(t2mvar[nt, :, :], lons1)

    # Convert to Celsius.
    t2m = t2m - 273.15

    # Define contour levels.
    clevs = np.arange(-30, 30.1, 2.0)
    lons, lats = np.meshgrid(lons, lats)
    x, y = bmap(lons, lats)

    # Create figure and AxesGrid instance.
    fig = plt.figure(figsize=(6, 8))
    grid = AxesGrid(
        fig,
        [0.05, 0.01, 0.9, 0.9],
        nrows_ncols=(3, 2),
        axes_pad=0.5,
        cbar_mode="single",
        cbar_pad=0.75,
        cbar_size=0.1,
        cbar_location="top",
        share_all=True)

    # Make subplots.
    for nt, fcsthr in enumerate(fcsthrs):
        bmap.ax = grid[nt]
        cs = bmap.contourf(x, y, t2m[nt, :, :], clevs,
                           cmap=plt.cm.jet, extend="both")
        bmap.drawcoastlines(linewidth=0.5)
        bmap.drawcountries()
        bmap.drawparallels(np.arange(-80, 81, 20))
        bmap.drawmeridians(np.arange(0, 360, 20))
        # Set panel title.
        bmap.ax.set_title(
            "%d-h forecast valid " % fcsthr + verifdates[nt], fontsize=9)

    # Set figure title.
    plt.figtext(
        0.5, 0.95,
        "2-m temp (\N{DEGREE SIGN}C) forecasts from %s" % verifdates[0],
        horizontalalignment="center", fontsize=14)

    # Draw a single colorbar.
    cax = grid.cbar_axes[0]
    fig.colorbar(cs, cax=cax, orientation="horizontal")
    plt.show()


if __name__ == "__main__":

    import sys
    import datetime as dt

    # Parse input date (default: today).
    if len(sys.argv) > 1:
        dateobj = dt.datetime.strptime(sys.argv[1], "%Y%m%d")
    else:
        dateobj = dt.datetime.today()
    main(dateobj, verbose=True)
