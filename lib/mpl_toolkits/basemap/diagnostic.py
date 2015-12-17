"""
These are diagnostic and debugging functions for basemap.
"""

def package_versions():
    """
    Gives version information for dependent packages.

    returns namedtuple BasemapPackageVersions
    """
    from collections import namedtuple
    from sys import version as sys_version

    from matplotlib import __version__ as matplotlib_version
    from numpy import __version__ as numpy_version
    import pyproj
    from shapefile import __version__ as pyshp_version

    import _geoslib
    from mpl_toolkits.basemap import __version__ as basemap_version
    
    # import optional dependencies
    try:
        from OWSLib import __version__ as OWSLib_version
    except ImportError:
        OWSLib_version = 'not installed'

    try:
        from PIL import VERSION as pil_version
        try:
            from PIL import PILLOW_VERSION as pillow_version
        except ImportError:
            pillow_version = 'not installed'
    except ImportError:
        pil_version = 'not installed'
        pillow_version = 'not installed'
    
    # Get PROJ.4 version info in a floating point number
    proj_ver_num = pyproj.Proj(init='epsg:4326').proj_version
    # reformats floating point number into string (4.90 becomes '4.9.0')
    proj4_version = '.'.join(list(str(int(proj_ver_num*100))))
    
    BasemapPackageVersions = namedtuple(
                               'BasemapPackageVersions',
                               """Python, basemap, matplotlib,
                                  numpy, pyproj, pyshp, PROJ4, GEOS,
                                  OWSLib, PIL, Pillow""")

    return BasemapPackageVersions(
                   Python = sys_version,
                   basemap = basemap_version,
                   matplotlib = matplotlib_version,
                   numpy = numpy_version,
                   pyproj = pyproj.__version__,
                   pyshp = pyshp_version,
                   PROJ4 = proj4_version,
                   GEOS = _geoslib.__geos_version__,
                   # optional dependencies below
                   OWSLib = OWSLib_version,
                   PIL = pil_version,
                   Pillow = pillow_version)
