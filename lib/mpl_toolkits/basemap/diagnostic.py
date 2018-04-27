from __future__ import (absolute_import, division, print_function)

"""
These are diagnostic and debugging functions for basemap.
"""

def proj4_version():
    """
    Gives the proj.4 library's version number. (requires pyproj to be installed)

    returns string, so proj.4 version 4.9.3 will return "4.9.3"
    """
    import pyproj
    try:
        return pyproj.proj_version_str
    except AttributeError:
        # for pyproj versions 1.9.5.1 and before, this will run
        # Get PROJ4 version in a floating point number
        proj4_ver_num = pyproj.Proj(proj='latlong').proj_version
        
        # reformats floating point number into string (4.90 becomes '4.9.0')
        # Exploits single number version numbers for proj4,
        return '.'.join( str(int(proj4_ver_num*100)) )
    
    
def package_versions():
    """
    Gives version information for dependent packages.

    returns namedtuple BasemapPackageVersions
    """
    from collections import namedtuple
    from sys import version as sys_version

    from matplotlib import __version__ as matplotlib_version
    from numpy import __version__ as numpy_version
    from pyproj import __version__ as pyproj_version
    from shapefile import __version__ as pyshp_version

    import _geoslib
    from mpl_toolkits.basemap import __version__ as basemap_version
    
    try:
        # geodesic is a part of proj.4 library
        # new variable in pyproj versions greater than 1.9.5.1
        from pyproj import geodesic_version_str as geodesic_version
    except ImportError:
        geodesic_version = 'Unknown'
    
    # import optional dependencies
    try:
        from OWSLib import __version__ as OWSLib_version
    except ImportError:
        OWSLib_version = 'not installed'

    try:
        from PIL import __version__ as pillow_version
    except ImportError:
        pillow_version = 'not installed'
    
    
    BasemapPackageVersions = namedtuple(
                               'BasemapPackageVersions',
                               """Python, basemap, matplotlib,
                                  numpy, pyproj, pyshp, PROJ4, geodesic, 
                                  GEOS, OWSLib, Pillow""")

    return BasemapPackageVersions(
                   Python = sys_version,
                   basemap = basemap_version,
                   matplotlib = matplotlib_version,
                   numpy = numpy_version,
                   pyproj = pyproj_version,
                   pyshp = pyshp_version,
                   PROJ4 = proj4_version(),
                   geodesic = geodesic_version,
                   GEOS = _geoslib.__geos_version__,
                   # optional dependencies below
                   OWSLib = OWSLib_version,
                   Pillow = pillow_version)

def check_proj_inv_hammer(segfault_protection=True):
    """
    Check if the inverse of the hammer projection is supported by installed
    version of PROJ4.
    
    segfault_protection   True (default) - test while protecting from segfault
                          False -  testing that might cause Python to segfault.
                                   BE CAREFUL setting this flag to False!
                                   If it segfaults, this the inverse hammer is not supported.

    returns True      - inverse hammer is supported
            False     - inverse hammer is not supported
            "Unknown" - support is Unknown
    """
    from distutils.version import LooseVersion
    from pyproj import __version__ as pyproj_version
    
    if LooseVersion(proj4_version()) > LooseVersion('4.9.2'):
        return True
    
    if LooseVersion(pyproj_version) > LooseVersion('1.9.5.1') \
            or segfault_protection is False:
        from pyproj import Proj
        hammer = Proj(proj='hammer')
        
        x, y = hammer(-30.0, 40.0)
        try:
            lon, lat = hammer(x, y, inverse=True)
            return True
        except RuntimeError:            
            return False
    
    return 'Unknown'
