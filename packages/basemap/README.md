# basemap

Plot on map projections (with coastlines and political boundaries) using
[`matplotlib`].

This package depends on the support package [`basemap-data`] with the
basic [`basemap`] data assets, and optionally on the support package
[`basemap-data-hires`] with high-resolution data assets.

## Installation

Precompiled binary wheels for Windows and GNU/Linux are available in
PyPI (architectures x86 and x64, Python 2.7 and 3.5+) and can be
installed with [`pip`]:
```sh
python -m pip install basemap
```

If you need to install from source, please visit the
[GitHub repository](https://github.com/matplotlib/basemap) for a
step-by-step description.

## License

The library is licensed under the terms of the [MIT] license (see
[`LICENSE`]). The GEOS dynamic library bundled with the package wheels
is provided under the terms of the [LGPLv2.1] license as given in
[`LICENSE.geos`].


[`matplotlib`]:
https://matplotlib.org/
[`basemap`]:
https://matplotlib.org/basemap/
[`basemap-data`]:
https://pypi.org/project/basemap-data
[`basemap-data-hires`]:
https://pypi.org/project/basemap-data-hires
[`pip`]:
https://pip.pypa.io/

[LGPLv2.1]:
https://spdx.org/licenses/LGPL-2.1-only.html
[MIT]:
https://spdx.org/licenses/MIT.html

[`LICENSE`]:
https://github.com/matplotlib/basemap/blob/v1.3.7/packages/basemap/LICENSE
[`LICENSE.geos`]:
https://github.com/matplotlib/basemap/blob/v1.3.7/packages/basemap/LICENSE.geos
