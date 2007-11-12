from ctypes import CDLL
from ctypes.util import find_library
import os, sys

def find_geoslib():
    lgeos = False
    # if GEOS_DIR env var set, look there for shared lib.
    if os.environ.has_key('GEOS_DIR'):
        geos_dir = os.environ['GEOS_DIR']
        from ctypes import CDLL
        if sys.platform == 'win32':
            try:
                lgeos = CDLL(os.path.join(geos_dir,'libgeos_c.dll'))
            except (ImportError, WindowsError):
                # Try GEOS DLL from the Windows PostGIS installer for
                # PostgreSQL 8.2 before failing
                try:
                    lgeos = CDLL(os.path.join(geos_dir,'libgeos_c-1.dll'))
                except:
                    pass
        elif sys.platform == 'darwin':
            try:
                lgeos = CDLL(os.path.join(geos_dir,'libgeos_c.dylib'))
            except:
                pass
        else:
            # Try the major versioned name first, falling back on the unversioned name.
            try:
                lgeos = CDLL(os.path.join(geos_dir,'libgeos_c.so.1'))
            except ImportError:
                try:
                    lgeos = CDLL(os.path.join(geos_dir,'libgeos_c.so'))
                except:
                    pass
    # if GEOS_DIR env var not set, use find_library to look for shared lib.
    else:
        if sys.platform == 'win32':
            try:
                lgeos = find_library('libgeos_c.dll')
            except (ImportError, WindowsError):
                # Try GEOS DLL from the Windows PostGIS installer for
                # PostgreSQL 8.2 before failing
                lgeos = find_library('libgeos_c-1.dll')
        elif sys.platform == 'darwin':
            try:
                lgeos = find_library('libgeos_c.dylib')
            except: 
                # fink installs in /sw, but find_library doesn't look there.
                lgeos = find_library('/sw/lib/libgeos_c.dylib')
        else:
            # Try the major versioned name first, falling back on the unversioned name.
            try:
                lgeos = find_library('libgeos_c.so.1')
            except ImportError:
                lgeos = find_library('libgeos_c.so')
    if not lgeos:
        raise ImportError("""
Cannot find libgeos_c. Please install version 2.2.3.
If it is installed already, please set the GEOS_DIR environment
variable to the location in which it is installed.""")
    else:  
        # lgeos is either a string (containing the absolute
        # path to the shared lib) or a CDLL instance.
        # convert it to a CDLL instance if it is not one already.
        try:
            lgeos = CDLL(lgeos)
        except:
            pass
    return lgeos
