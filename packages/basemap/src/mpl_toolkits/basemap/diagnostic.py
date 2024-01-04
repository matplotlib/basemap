"""Diagnostic and debugging functions for :mod:`mpl_toolkits.basemap`."""


def proj4_version():
    """Return PROJ library version as string (requires :mod:`pyproj`)."""

    import pyproj

    try:
        # Case for `pyproj` versions from 1.9.6.
        return pyproj.proj_version_str
    except AttributeError:
        # Case for `pyproj` versions before 1.9.6, where PROJ version
        # string is generated from its version in floating point format.
        proj4_ver_num = pyproj.Proj(proj="latlong").proj_version

        # Reformat floating point number into string and return
        # (e.g. 4.90 becomes "4.9.0").
        return ".".join(str(int(proj4_ver_num * 100)))


def package_versions():
    """Return version information for package dependencies."""

    from collections import namedtuple
    from sys import version as sys_version

    # pylint: disable=no-name-in-module
    from numpy import __version__ as numpy_version
    from matplotlib import __version__ as matplotlib_version
    # pylint: enable=no-name-in-module

    from pyproj import __version__ as pyproj_version
    from shapefile import __version__ as pyshp_version

    from _geoslib import __geos_version__ as geos_version
    from . import __version__ as basemap_version

    try:
        # PROJ geodesic version can be read for `pyproj` 1.9.6 or newer.
        from pyproj import geodesic_version_str as geodesic_version
    except ImportError:  # pragma: no cover
        geodesic_version = "unknown"

    # Import optional dependencies.
    try:
        from OWSLib import __version__ as owslib_version
    except ImportError:  # pragma: no cover
        owslib_version = "not installed"
    try:
        from PIL import __version__ as pillow_version
    except ImportError:  # pragma: no cover
        pillow_version = "not installed"

    # pylint: disable=invalid-name
    BasemapPackageVersions = namedtuple(
        "BasemapPackageVersions",
        "Python, basemap, numpy, matplotlib, GEOS, pyproj, PROJ4, geodesic, "
        "pyshp, Pillow, OWSLib")
    # pylint: enable=invalid-name

    return BasemapPackageVersions(**{
        # Mandatory dependencies.
        "Python":
            sys_version,
        "basemap":
            basemap_version,
        "numpy":
            numpy_version,
        "matplotlib":
            matplotlib_version,
        "GEOS":
            str(geos_version.decode("ascii")),
        "pyproj":
            pyproj_version,
        "PROJ4":
            proj4_version(),
        "geodesic":
            geodesic_version,
        "pyshp":
            pyshp_version,
        # Optional dependencies.
        "Pillow":
            pillow_version,
        "OWSLib":
            owslib_version,
    })


def check_proj_inv_hammer(segfault_protection=True):
    """Return if installed PROJ supports inverse of Hammer projection.

    Parameters
    ----------

    segfault_protection : bool, optional
       if True (default), perform check while protecting from segfault;
       if False, perform check allowing Python to segfault (if segfault
       occurs, inverse of Hammer projection is not supported)

    Returns
    -------

    result : {True, False, "Unknown"}
        check result as bool or "Unknown" if check is inconclusive
    """

    import pyproj
    from packaging.version import Version

    if Version(proj4_version()) > Version("4.9.2"):
        return True

    if Version(pyproj.__version__) > Version("1.9.5.1") or not segfault_protection:
        hammer = pyproj.Proj(proj="hammer")
        xy_coordinates = hammer(-30.0, 40.0)
        try:
            hammer(xy_coordinates, inverse=True)
            return True
        except RuntimeError:
            return False

    return "Unknown"
