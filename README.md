# basemap

Plot on map projections (with coastlines and political boundaries) using
[`matplotlib`].

## Installation

Precompiled binary wheels for Windows, GNU/Linux and MacOS are available
on PyPI and can be installed with [`pip`]:
```sh
python -m pip install basemap
```

For specific details on how to install [`basemap`] through [`conda`] or
from source, please refer to the [`basemap` installation instructions]
in the documentation.

## Requirements

This package depends on [`basemap-data`] with the basic [`basemap`]
data assets supporting the essential functionality.

This package depends optionally on [`basemap-data-hires`] with
the high-resolution data assets, which can be installed manually
with [`pip`]:
```sh
python -m pip install basemap-data-hires
```

This package depends optionally on [`OWSLib`] for the `Basemap` method
`Basemap.wmsimage`.

## License

The library is licensed under the terms of the [MIT] license (see
[`LICENSE`]). The GEOS dynamic library bundled with the package wheels
is provided under the terms of the [LGPL-2.1-only] license as given in
[`LICENSE.geos`].

## Documentation

See https://matplotlib.org/basemap/.

See scripts in the `doc/examples` directory for example usage.

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
polygons (this is illustrated in the [`ortho_demo.py`] example).


[`pip`]:
https://pip.pypa.io/
[`conda`]:
https://github.com/conda/conda

[`matplotlib`]:
https://matplotlib.org/
[`basemap`]:
https://matplotlib.org/basemap/
[`basemap-data`]:
https://pypi.org/project/basemap-data
[`basemap-data-hires`]:
https://pypi.org/project/basemap-data-hires
[`OWSLib`]:
https://pypi.org/project/OWSLib

[`basemap` installation instructions]:
https://matplotlib.org/basemap/stable/users/installation.html
[`ortho_demo.py`]:
https://github.com/matplotlib/basemap/blob/v2.0.0/doc/examples/ortho_demo.py

[MIT]:
https://spdx.org/licenses/MIT.html
[LGPL-2.1-only]:
https://spdx.org/licenses/LGPL-2.1-only.html

[`LICENSE`]:
https://github.com/matplotlib/basemap/blob/v2.0.0/LICENSE
[`LICENSE.geos`]:
https://github.com/matplotlib/basemap/blob/v2.0.0/LICENSE.geos
