.. _installing:

**********
Installing
**********

Dependencies
============

**Requirements**

These are external packages which you will need to install before
installing basemap. 


matplotlib 0.98 (or later, `download <http://sf.net/projects/matplotlib/>`__)

python 2.4 (or later but not python3)
    matplotlib requires python 2.4 or later (`download <http://www.python.org/download/>`__)

numpy 1.1 (or later)
    array support for python (`download <http://sourceforge.net/project/showfiles.php?group_id=1369&package_id=175103>`__)

**Required libraries that ship with basemap**

`GEOS <http://trac.osgeo.org/geos/>`__ (Geometry Engine - Open Source) library 2.2.3 or later (2.2.3 recommended).
    Source code is included in the geos-2.2.3 directory. 
    When building from source, must be built and installed separately
    from basemap (see build instructions below).
    Included in Windows binary installers.

`PROJ4 <http://trac.osgeo.org/proj/>`__ Cartographic Projections Library.
    Patched version automatically built into basemap.

`pyshapelib <http://intevation.de/pipermail/thuban-devel/2004-May/000184.html>`__
    C library with python interface for reading ESRI shapefiles.
    If not present, will be installed with basemap.

`pupnyere <http://pypi.python.org/pypi/pupynere/>`__ 
    Pure python `netCDF <http://www.unidata.ucar.edu/software/netcdf/>`__
    patched version automatically installed with basemap.

`pydap <http://code.google.com/p/pydap>`__ 
    Pure python `OPeNDAP <http://opendap.org>`__ implementation.
    If not present, client (not server) will be installed with basemap.

`httplib2 <http://code.google.com/p/httplib2>`__ (needed for pydap client).
    If not present, will be installed with basemap.
  

**Optional libraries**

PIL
    Python Imaging Library (`download <http://www.pythonware.com/products/pil/>`__),
    only needed for :func:`~mpl_toolkits.basemap.Basemap.bluemarble` and :func:`~mpl_toolkits.basemap.Basemap.warpimage` instance methods.

Installation
============

Windows binary installers are available for
`download <http://sourceforge.net/project/showfiles.php?group_id=80706&package_id=142792/&abmode=1>`__. 

For other platforms, download the source release and follow these steps:


* Install pre-requisite requirements.

* Untar the basemap version X.Y.Z source tar.gz file, and
  and cd to the basemap-X.Y.Z directory.

* Install the GEOS library.  If you already have it on your
  system, just set the environment variable GEOS_DIR to point to the location 
  of libgeos_c and geos_c.h (if libgeos_c is in /usr/local/lib and
  geos_c.h is in /usr/local/include, set GEOS_DIR to /usr/local).
  Then go to next step.  If you don't have it, you can build it from
  the source code included with basemap by following these steps::

      cd geos-2.2.3
      export GEOS_DIR=<where you want the libs and headers to go>
      # A reasonable choice on a Unix-like system is /usr/local, or
      # if you don't have permission to write there, your home directory.
      ./configure --prefix=$GEOS_DIR 
      make; make install

* cd back to the top level basemap directory (basemap-X.Y.Z) and
  run the usual ``python setup.py install``.  Check your installation
  by running ``from mpl_toolkits.basemap import Basemap`` at the python
  prompt.

  Basemap includes three auxilliary packages, pydap, pyshapelib httplib2.
  By default, setup.py checks to 
  see if these are already installed, and if so does not try to overwrite 
  them. If you get import errors related to any of these packages, 
  edit setup.cfg and set the appropriate entry to True to force 
  installation of the included versions.

* To test, cd to the examples directory and run ``python simpletest.py``.
  To run all the examples (except those that have extra dependencies
  or require an internet connection), execute ``python run_all.py``.
