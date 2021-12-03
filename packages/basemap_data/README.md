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
are extracted from datasets provided with the [Generic Mapping Tools
(GMT)] and are included under the terms given in [`COPYING`].

The EPSG file belongs to the PROJ.4 package and is licensed here under
the [MIT license], as stated in [`LICENSE.epsg`].

Everything else is licensed under the [Historical Permission Notice and
Disclaimer (HPND)], as stated in [`LICENSE.hpnd`].


[`matplotlib`]:
https://matplotlib.org/
[`basemap`]:
https://matplotlib.org/basemap/
[`pip`]:
https://pip.pypa.io/
[`COPYING`]:
https://github.com/molinav/basemap/blob/develop/packages/basemap_data/COPYING
[`LICENSE.epsg`]:
https://github.com/molinav/basemap/blob/develop/packages/basemap_data/LICENSE.epsg
[`LICENSE.hpnd`]:
https://github.com/molinav/basemap/blob/develop/packages/basemap_data/LICENSE.hpnd
[Generic Mapping Tools (GMT)]:
http://gmt.soest.hawaii.edu
[Historical Permission Notice and Disclaimer (HPND)]:
https://opensource.org/licenses/HPND
[MIT license]:
https://opensource.org/licenses/MIT
