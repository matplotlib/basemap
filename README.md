# Basemap

Plot on map projections (with coastlines and political boundaries)
using matplotlib.

## Requirements

Basic requirements are the following:

* Python 2.6 (or higher)
* [matplotlib](https://github.com/matplotlib/matplotlib)
* [numpy](https://github.com/numpy/numpy)
* [pyproj](https://github.com/pyproj4/pyproj)
* [pyshp](https://github.com/GeospatialPython/pyshp)

Optional requirements include:

* [OWSLib](https://github.com/geopython/OWSLib). It is needed for the
  `Basemap.wmsimage` function.

* [Pillow](https://github.com/python-pillow/Pillow). It is needed for
  the methods `Basemap.bluemarble`, `Basemap.etopo`,
  `Basemap.shadedrelief` and `Basemap.warpimage`.

## Installation

The `basemap-data` and `basemap-data-hires` packages are available in
PyPI and can be installed with [`pip`](https:/pip.pypa.io/):
```sh
python -m pip install basemap-data
python -m pip install basemap-data-hires
```

Precompiled `basemap` binary wheels for Windows and GNU/Linux are also
available in PyPI (architectures x86 and x64, Python 2.7 and 3.5+):
```sh
python -m pip install basemap
```

Otherwise, you will need to install `basemap` from source as follows:

1. Install pre-requisite Python modules:
   - [cython](https://github.com/cython/cython)
   - [numpy](https://github.com/numpy/numpy)

2. Download the `basemap` source code and move to the `packages/basemap`
   folder:
   ```sh
   git clone --depth 1 https://github.com/matplotlib/basemap.git
   cd basemap/packages/basemap
   ```

3. Build the [GEOS](https://github.com/libgeos/geos) library. You may
   use the helper provided in `utils`, i.e.
   ```sh
   export GEOS_DIR=<your desired location>
   python -c "import utils; utils.GeosLibrary('3.5.1').build(installdir='${GEOS_DIR}')"
   ```
   or you can link directly to the system library if it is already
   installed. `GEOS_DIR` must point to the GEOS installation prefix;
   e.g. if `libgeos_c.so` is located in `/usr/lib` and `geos_c.h` is
   located in `/usr/include`, then you must set `GEOS_DIR` to `/usr`.

4. Build and install the `basemap` binary wheel:
   ```sh
   python -m pip install .
   ```
   On Linux, if your Python was installed through a package management
   system, make sure that you have the Python header `Python.h` required
   to build Cython extensions (e.g. on Debian-like systems, you should
   have the package `python-dev` installed).

5. Check that the package installed correctly by executing:
   ```sh
   python -c "from mpl_toolkits.basemap import Basemap"
   ```
   You can also test the examples available in the `examples` folder.

## License

The source code and data assets are under the following licenses:

* `basemap`: [MIT].
  * GEOS bundled dynamic library is under the [LGPL-2.1-only] license.
* `basemap-data`: [LGPL-3.0-or-later].
  * The EPSG file and the JPG images are also under the [MIT] license.
* `basemap-data-hires`: [LGPL-3.0-or-later].

For a full description, please visit the `README` and `LICENSE` files of
each package.

[MIT]:
https://spdx.org/licenses/MIT.html
[LGPL-2.1-only]:
https://spdx.org/licenses/LGPL-2.1-only.html
[LGPL-3.0-or-later]:
https://spdx.org/licenses/LGPL-3.0-or-later.html

## Documentation

See https://matplotlib.github.io/basemap/

See scripts in `examples` directory for example usage.

Read the FAQ and/or email the matplotlib-users mailing list if you have
problems or questions.

## Contact

Ben Root <ben.v.root@gmail.com>

Víctor Molina García ([@molinav](https://github.com/molinav))

## Thanks

Special thanks to John Hunter, Andrew Straw, Eric Firing, Rob Hetland,
Scott Sinclair, Ivan Lima, Erik Andersen, Michael Hearne, Jesper Larsen,
Ryan May, David Huard, Mauro Cavalcanti, Jonas Bluethgen, Chris Murphy,
Pierre Gerard-Marchant, Christoph Gohlke, Eric Bruning, Stephane
Raynaud, Tom Loredo, Patrick Marsh, Phil Elson, and Henry Hammond for
valuable contributions.

## Known bugs

The `Basemap.fillcontinents` method doesn't always do the right thing.
Matplotlib always tries to fill the inside of a polygon. Under certain
situations, what is the inside of a coastline polygon can be ambiguous,
and the outside may be filled instead of the inside. A workaround is to
change the map projection region slightly or mask the land areas with
the `Basemap.drawlsmask` method instead of filling the coastline
polygons (this is illustrated in the `ortho_demo.py` example). 
