# basemap-data

Plot on map projections (with coastlines and political boundaries) using
[`matplotlib`].

This is a support package for [`basemap`] with the basic data assets
required by [`basemap`] to work.

## Installation

The package is available in PyPI and can be installed with [`pip`]:
```python
python -m pip install basemap-data
```

## License

The land-sea mask, coastline, lake, river and political boundary data
are extracted from the [GSHHG] datasets (version 2.3.6) using [GMT]
(5.x series) and are included under the terms of the [LGPLv3+] license
(see [`COPYING`] and [`COPYING.LESSER`]).

The other files are included under the terms of the [MIT] license. See
[`LICENSE.epsg`] for the EPSG file (taken from the PROJ.4 package) and
[`LICENSE.mit`] for the rest.


[`matplotlib`]:
https://matplotlib.org/
[`basemap`]:
https://matplotlib.org/basemap/
[`pip`]:
https://pip.pypa.io/

[GSHHG]:
https://www.soest.hawaii.edu/pwessel/gshhg
[GMT]:
https://www.generic-mapping-tools.org/

[LGPLv3+]:
https://spdx.org/licenses/LGPL-3.0-or-later.html
[MIT]:
https://spdx.org/licenses/MIT.html

[`COPYING`]:
https://github.com/matplotlib/basemap/blob/v1.3.2/packages/basemap_data/COPYING
[`COPYING.LESSER`]:
https://github.com/matplotlib/basemap/blob/v1.3.2/packages/basemap_data/COPYING.LESSER
[`LICENSE.epsg`]:
https://github.com/matplotlib/basemap/blob/v1.3.2/packages/basemap_data/LICENSE.epsg
[`LICENSE.mit`]:
https://github.com/matplotlib/basemap/blob/v1.3.2/packages/basemap_data/LICENSE.mit
