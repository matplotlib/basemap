.. _installing:

**********
Installing
**********

Dependencies
============

**Requirements**

These are external packages which you will need to install before
installing Basemap.


Matplotlib 1.0.0 (or later, `download <https://matplotlib.org/users/installing.html>`__)

Python 2.6 (or later, including Python 3) (`download <http://www.python.org/download/>`__)
    Matplotlib 2.2 LTS requires Python 2.7 or later
    Matplotlib 3.0 requires Python 3.5 or later

NumPy 1.2.1 (or later)
    Array support for Python (`download <http://www.numpy.org/>`__)

`PROJ4 <https://trac.osgeo.org/proj/>`__ Cartographic Projections Library.

**Required library that ships with Basemap**

`GEOS <http://trac.osgeo.org/geos/>`__ (Geometry Engine - Open Source) library 3.1.1 or later.
    Source code is included in the geos-3.3.3 directory.
    When building from source, must be built and installed separately
    from basemap (see build instructions below).
    Included in Windows binary installers.

**Optional libraries**

Pillow
    Python Imaging Library (`download <https://python-pillow.org/>`__),
    only needed for :func:`~mpl_toolkits.basemap.Basemap.bluemarble`, :func:`~mpl_toolkits.basemap.Basemap.etopo`, :func:`~mpl_toolkits.basemap.Basemap.shadedrelief` and :func:`~mpl_toolkits.basemap.Basemap.warpimage` instance methods.

Installation
============

Download either Windows binary installers or source tarballs
`here <https://github.com/matplotlib/basemap/releases/>`__.

To install from the source, follow these steps:


* Install pre-requisite requirements.

* Untar the basemap version X.Y.Z source tar.gz file, and
  and cd to the basemap-X.Y.Z directory.

* Install the GEOS library.  If you already have it on your
  system, just set the environment variable GEOS_DIR to point to the location
  of libgeos_c and geos_c.h (if libgeos_c is in /usr/local/lib and
  geos_c.h is in /usr/local/include, set GEOS_DIR to /usr/local).
  Then go to next step.  If you don't have it, you can build it from
  the source code included with basemap by following these steps::

      cd geos-3.3.3
      export GEOS_DIR=<where you want the libs and headers to go>
      # A reasonable choice on a Unix-like system is /usr/local, or
      # if you don't have permission to write there, your home directory.
      ./configure --prefix=$GEOS_DIR
      make; make install

* cd back to the top level basemap directory (basemap-X.Y.Z) and
  run the usual ``python setup.py install``.  Check your installation
  by running ``from mpl_toolkits.basemap import Basemap`` at the Python
  prompt.

* To test, cd to the examples directory and run ``python simpletest.py``.
  To run all the examples (except those that have extra dependencies
  or require an internet connection), execute ``python run_all.py``.
