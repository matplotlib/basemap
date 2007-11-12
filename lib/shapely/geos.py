"""
Exports the libgeos_c shared lib, GEOS-specific exceptions, and utilities.
"""

import atexit
from ctypes import CDLL, CFUNCTYPE, c_char_p
from ctypes.util import find_library
import sys

from find_geoslib import find_geoslib
lgeos = find_geoslib()

# Exceptions

class ReadingError(Exception):
    pass

class DimensionError(Exception):
    pass

class TopologicalError(Exception):
    pass

class PredicateError(Exception):
    pass


# GEOS error handlers, which currently do nothing.

def error_handler(fmt, list):
    pass
error_h = CFUNCTYPE(None, c_char_p, c_char_p)(error_handler)

def notice_handler(fmt, list):
    pass
notice_h = CFUNCTYPE(None, c_char_p, c_char_p)(notice_handler)

# Register a cleanup function

def cleanup():
    lgeos.finishGEOS()

atexit.register(cleanup)

lgeos.initGEOS(notice_h, error_h)


