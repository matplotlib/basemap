#Basemap

Plot on map projections (with coastlines and political boundaries)
using matplotlib.

##Requirements

* Python 2.6 (or higher)

* matplotlib

* numpy 

* [pyproj](https://github.com/jswhit/pyproj)

* [pyshp](https://github.com/GeospatialPython/pyshp)

* The GEOS (Geometry Engine - Open Source) library (version 3.1.1 or higher).
Source code is included in the geos-3.3.3 directory.

* On linux, if your python was installed via a package management system, make
sure the corresponding "python-dev" package is also installed.  Otherwise, you
may not have the python header (Python.h), which is required to build python
C extensions.

###Optional

* [OWSLib](https://github.com/geopython/OWSLib) (optional) It is needed
for the BaseMap.wmsimage function.

* [Pillow](https://python-pillow.github.io/) (optional)  It is
needed for Basemap warpimage, bluemarble, shadedrelief, and etop methods.
PIL should work on Python 2.x.  Pillow is a maintained fork of PIL.

##Copyright

source code for the GEOS library is
included in the 'geos-3.3.3' directory under the terms given in
LICENSE_geos.

the land-sea mask, coastline, lake, river and political boundary data are extracted
from datasets provided with the [Generic Mapping Tools (GMT)](http://gmt.soest.hawaii.edu)
and are included under the terms given in LICENSE_data.

Everything else (including src/_geos.c, and src/_geos.pyx):

copyright (c) 2011 by Jeffrey Whitaker.

Permission to use, copy, modify, and distribute this software and its
documentation for any purpose and without fee is hereby granted,
provided that the above copyright notices appear in all copies and that
both the copyright notices and this permission notice appear in
supporting documentation.
THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, INDIRECT OR
CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.

##Documentation 

see http://matplotlib.github.com/basemap/

see scripts in 'examples' directory for example usage.

read the FAQ and/or email the matplotlib-users mailing list if 
you have problems or questions.

##Install

0. Install pre-requisite python modules numpy and matplotlib.

1. Then download basemap-X.Y.Z.tar.gz (approx 100 mb) from
the [GitHub Releases](https://github.com/matplotlib/basemap/releases) page, unpack and cd to basemap-X.Y.Z.

2. Install the GEOS library.  If you already have it on your
system, just set the environment variable GEOS_DIR to point to the location 
of libgeos_c and geos_c.h (if libgeos_c is in /usr/local/lib and
geos_c.h is in /usr/local/include, set GEOS_DIR to /usr/local).
Then go to step (3).  If you don't have it, you can build it from
the source code included with basemap by following these steps:

	```
	 > cd geos-3.3.3
	 > export GEOS_DIR=<where you want the libs and headers to go>
	   A reasonable choice on a Unix-like system is /usr/local, or
	   if you don't have permission to write there, your home directory.
	 > ./configure --prefix=$GEOS_DIR 
	 > make; make install
	```

3. cd back to the top level basemap directory (basemap-X.Y.Z) and
run the usual 'python setup.py install'.  Check your installation
by running "from mpl_toolkits.basemap import Basemap" at the python
prompt.

4. To test, cd to the examples directory and run 'python simpletest.py'.
To run all the examples (except those that have extra dependencies
or require an internet connection), execute 'python run_all.py'.

##Contact

Ben Root <ben.v.root@gmail.com>

##Thanks

Special thanks to John Hunter, Andrew Straw, Eric Firing, Rob Hetland, Scott Sinclair, Ivan Lima, Erik Andersen, Michael Hearne, Jesper Larsen, Ryan May, David Huard, Mauro Cavalcanti, Jonas Bluethgen, Chris Murphy, Pierre Gerard-Marchant, Christoph Gohlke, Eric Bruning, Stephane Raynaud, Tom Loredo, Patrick Marsh, Phil Elson, and Henry Hammond for valuable contributions.
