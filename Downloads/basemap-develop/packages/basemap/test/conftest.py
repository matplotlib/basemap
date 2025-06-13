"""Configuration file for :mod:`mpl_toolkits.basemap` test suite.

The :mod:`mpl_toolkits.basemap` test suite is configured so that all
tests involving :mod:`mpl_toolkits.basemap` are run using the 'Agg'
backend, because there is no guarantee that another backend will be
available at runtime.
"""
from __future__ import absolute_import

import os
import sys

if "MPLBACKEND" not in os.environ:
    os.environ["MPLBACKEND"] = "Agg"

try:
    from mpl_toolkits import basemap
except ImportError:
    import basemap
    sys.modules["mpl_toolkits.basemap"] = basemap
