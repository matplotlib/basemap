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
* GEOS library (Geometry Engine, Open Source) version 3.1.1 or higher.
* On Linux, if your Python was installed via a package management
  system, make sure the corresponding `python-dev` package is also
  installed. Otherwise, you may not have the Python header (`Python.h`),
  which is required to build Python C extensions.

Optional requirements include:

* [OWSLib](https://github.com/geopython/OWSLib). It is needed for the
  `BaseMap.wmsimage` function.

* [Pillow](https://python-pillow.github.io/). It is needed for Basemap
  warpimage, bluemarble, shadedrelief, and etopo methods. PIL should
  work on Python 2.x. Pillow is a maintained fork of PIL.

## Installation

The `basemap-data` and `basemap-data-hires` packages are available in
PyPI and can be installed with [`pip`](https:/pip.pypa.io/):
```sh
python -m pip install basemap-data
python -m pip install basemap-data-hires
```

Precompiled binaries for GNU/Linux are also available in PyPI:
```sh
python -m pip install basemap
```

Otherwise, you will need to install `basemap` from source as follows:

1. Install pre-requisite Python modules:
   - `cython`.
   - `numpy`.

2. Download the `basemap` source code and move to the `packages/basemap`
   folder:
   ```sh
   git clone https://github.com/matplotlib/basemap.git
   cd basemap/packages/basemap
   ```

3. Build the GEOS library. You may use the helper provided in `utils`, i.e.
   ```sh
   export GEOS_DIR=<your desired location>
   python -c "import utils; utils.GeosLibrary('3.6.5').build(installdir='${GEOS_DIR}')"
   ```
   or you can link directly to the system library if it is already installed.
   `GEOS_DIR` must point to the GEOS installation prefix; e.g. if `libgeos_c.so`
   is located in `/usr/lib` and `geos_c.h` is located in `/usr/include`, then
   you must set `GEOS_DIR` to `/usr`.

4. Build the basemap wheel from the `packages/basemap` folder and install it:
   ```sh
   python setup.py bdist_wheel
   python -m pip install dist/*.whl
   ```

5. Check that the package installed correctly by executing in the terminal:
   ```sh
   python -c "from mpl_toolkits.basemap import Basemap"
   ```
   You can also test the `basemap` examples available in the `examples` folder.

## License

The source code and data assets are under the following licenses:

* `basemap`: [MIT].
  * GEOS bundled dynamic library is under the [LGPL-2.1-only] license.
* `basemap-data`: [LGPL-3.0-or-later].
  * The EPSG file and the JPG images are also under the [MIT] license.
* `basemap-data-hires`: [LGPL-3.0-or-later].

For a full description, please visit the README and LICENSE files of
each package in the corresponding package folders.

[MIT]:
https://spdx.org/licenses/MIT.html
[LGPL-2.1-only]:
https://spdx.org/licenses/LGPL-2.1-only.html
[LGPL-3.0-or-later]:
https://spdx.org/licenses/LGPL-3.0-or-later.html

## Documentation

See http://matplotlib.github.io/basemap/

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
