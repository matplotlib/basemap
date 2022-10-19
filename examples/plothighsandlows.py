"""Plot H's and L's on a sea-level pressure map."""
from __future__ import print_function

import datetime as dt
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import addcyclic
from scipy.ndimage import minimum_filter
from scipy.ndimage import maximum_filter


def extrema(mat, mode="wrap", window=10):
    """Find the indices of local extrema (min and max) in the input array."""

    minimum = minimum_filter(mat, size=window, mode=mode)
    maximum = maximum_filter(mat, size=window, mode=mode)

    # Return the indices of the maxima, minima.
    # (mat == maximum) true if pixel is equal to the local max.
    # (mat == minimum) true if pixel is equal to the local in.
    return np.nonzero(mat == minimum), np.nonzero(mat == maximum)


def main():
    """Main function."""

    # Plot 00 UTC today.
    url = "http://nomads.ncep.noaa.gov/dods/gfs_0p25/gfs%Y%m%d/gfs_0p25_00z"
    date = dt.datetime.now()

    # Open OPeNDAP dataset.
    data = netCDF4.Dataset(date.strftime(url))

    # Read lats and lons.
    lats = data.variables["lat"][:]
    lons1 = data.variables["lon"][:]

    # Read prmsl and convert to hPa (mbar).
    prmsl = 0.01 * data.variables["prmslmsl"][0]

    # The window parameter controls the number of highs and lows detected
    # (higher value, fewer highs and lows).
    local_min, local_max = extrema(prmsl, mode="wrap", window=50)

    # Create Basemap instance.
    bmap = Basemap(projection="mill",
                   llcrnrlon=0, llcrnrlat=-80,
                   urcrnrlon=360, urcrnrlat=80)

    # Add wrap-around point in longitude.
    prmsl, lons = addcyclic(prmsl, lons1)

    # Define contour levels.
    clevs = np.arange(900, 1100., 5.)

    # Find x, y of map projection grid.
    lons, lats = np.meshgrid(lons, lats)
    x, y = bmap(lons, lats)

    # Create figure.
    fig = plt.figure(figsize=(8, 4.5))
    fig.add_axes([0.05, 0.05, 0.9, 0.85])
    bmap.contour(x, y, prmsl, clevs, colors="k", linewidths=1.0)
    bmap.drawcoastlines(linewidth=1.25)
    bmap.fillcontinents(color="0.8")
    bmap.drawparallels(np.arange(-80, 81, 20), labels=[1, 1, 0, 0])
    bmap.drawmeridians(np.arange(0, 360, 60), labels=[0, 0, 0, 1])
    xlows, xhighs = x[local_min], x[local_max]
    ylows, yhighs = y[local_min], y[local_max]
    lowvals, highvals = prmsl[local_min], prmsl[local_max]

    # Plot lows as blue L's, with min pressure value underneath.
    # Do not plot if there is already a L or H within dmin meters.
    xyplotted = []
    yoffset = 0.022 * (bmap.ymax - bmap.ymin)
    dmin = yoffset
    for x, y, p in zip(xlows, ylows, lowvals):
        if bmap.xmin < x < bmap.xmax and bmap.ymin < y < bmap.ymax:
            dist = [np.sqrt((x - x0)**2 + (y - y0)**2) for x0, y0 in xyplotted]
            if not dist or min(dist) > dmin:
                bbox = dict(boxstyle="square", ec="None", fc=(1, 1, 1, 0.5))
                plt.text(x, y, "L", fontsize=14, fontweight="bold",
                         ha="center", va="center", color="b")
                plt.text(x, y - yoffset, repr(int(p)), fontsize=9,
                         ha="center", va="top", color="b", bbox=bbox)
                xyplotted.append((x, y))
    # Plot highs as red H's, with max pressure value underneath.
    xyplotted = []
    for x, y, p in zip(xhighs, yhighs, highvals):
        if bmap.xmin < x < bmap.xmax and bmap.ymin < y < bmap.ymax:
            dist = [np.sqrt((x - x0)**2 + (y - y0)**2) for x0, y0 in xyplotted]
            if not dist or min(dist) > dmin:
                bbox = dict(boxstyle="square", ec="None", fc=(1, 1, 1, 0.5))
                plt.text(x, y, "H", fontsize=14, fontweight="bold",
                         ha="center", va="center", color="r")
                plt.text(x, y - yoffset, repr(int(p)), fontsize=9,
                         ha="center", va="top", color="r", bbox=bbox)
                xyplotted.append((x, y))

    # Set plot title and show.
    plt.title("Mean Sea-Level Pressure (with Highs and Lows) %s" % date)
    plt.show()


if __name__ == "__main__":
    main()
