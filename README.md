# Basemap

Plot on map projections (with coastlines and political boundaries)
using matplotlib.

:warning: **Warning**: this package is being deprecated in favour of
[cartopy](https://scitools.org.uk/cartopy/docs/latest/).

## Requirements

Basic requirements are the following:

* Python 2.6 (or higher)
* matplotlib
* numpy
* [pyproj](https://github.com/pyproj4/pyproj)
* [pyshp](https://github.com/GeospatialPython/pyshp)
* The GEOS (Geometry Engine - Open Source) library (version 3.1.1 or
  higher). Source code is included in the `geos-3.3.3` directory.
* On Linux, if your Python was installed via a package management
  system, make sure the corresponding `python-dev` package is also
  installed. Otherwise, you may not have the Python header (`Python.h`),
  which is required to build Python C extensions.

Optional requirements include:

* [OWSLib](https://github.com/geopython/OWSLib). It is needed for the
  `BaseMap.wmsimage` function.

* [Pillow](https://python-pillow.github.io/). It is needed for Basemap
  warpimage, bluemarble, shadedrelief, and etop methods. PIL should work
  on Python 2.x. Pillow is a maintained fork of PIL.

## License

The source code and data assets are under the following licenses:

* `basemap`:
  * [HPND]: Python source code and Python wrapper for GEOS.
  * [LGPL-2.1-only]: GEOS source code and dynamic library.
* `basemap-data`:
  * [GPL-2.0-or-later]: land-sea mask, coastline, lake, river and
    political boundary data derived from GMT.
  * [MIT]: EPSG file.
  * [HPND]: remaining data files.
* `basemap-data-hires`:
  * [GPL-2.0-or-later]: land-sea mask, coastline, lake, river and
    political boundary data derived from GMT.

For a full description, please visit the README and LICENSE files of
each package in the corresponding package folders.

[MIT]:
https://spdx.org/licenses/MIT.html
[HPND]:
https://spdx.org/licenses/HPND.html
[LGPL-2.1-only]:
https://spdx.org/licenses/LGPL-2.1-only.html
[GPL-2.0-or-later]:
https://spdx.org/licenses/GPL-2.0-or-later.html


## Documentation

See http://matplotlib.github.io/basemap/

See scripts in `examples` directory for example usage.

Read the FAQ and/or email the matplotlib-users mailing list if you have
problems or questions.

## Install

0. Install pre-requisite Python modules numpy and matplotlib.

1. Then download `basemap-X.Y.Z.tar.gz` (approx 100 MB) from the
[GitHub Releases](https://github.com/matplotlib/basemap/releases) page,
unpack and `cd` to `basemap-X.Y.Z`.

2. Install the GEOS library. If you already have it on your system, just
   set the environment variable `GEOS_DIR` to point to the location of
   `libgeos_c` and `geos_c.h` (if `libgeos_c` is in `/usr/local/lib` and
   `geos_c.h` is in `/usr/local/include`, set `GEOS_DIR` to
   `/usr/local`). Then go to step (3). If you don't have it, you can
   build it from the source code included with basemap by following
   these steps:
   ```sh
   cd geos-3.3.3
   export GEOS_DIR=<where you want the libs and headers to go>
   # A reasonable choice on a Unix-like system is /usr/local, or
   # if you don't have permission to write there, your home directory.
   ./configure --prefix=$GEOS_DIR
   make; make install
   ```

3. `cd` back to the top level basemap directory (`basemap-X.Y.Z`) and
   run the usual `python setup.py install`. Check your installation by
   running ``"from mpl_toolkits.basemap import Basemap"`` at the Python
   prompt.

4. To test, `cd` to the examples folder and run `python simpletest.py`.
   To run all the examples (except those that have extra dependencies or
   require an internet connection), execute `python run_all.py`.

An alternative method is using `pip`:
```
python -m pip install --user git+https://github.com/matplotlib/basemap.git
```

## Contact

Ben Root <ben.v.root@gmail.com>

## Thanks

Special thanks to John Hunter, Andrew Straw, Eric Firing, Rob Hetland,
Scott Sinclair, Ivan Lima, Erik Andersen, Michael Hearne, Jesper Larsen,
Ryan May, David Huard, Mauro Cavalcanti, Jonas Bluethgen, Chris Murphy,
Pierre Gerard-Marchant, Christoph Gohlke, Eric Bruning, Stephane
Raynaud, Tom Loredo, Patrick Marsh, Phil Elson, and Henry Hammond for
valuable contributions.
