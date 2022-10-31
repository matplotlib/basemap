# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog], and the project adheres to
[Semantic Versioning] since version 1.3.0.

[Keep a Changelog]:
https://keepachangelog.com/en/1.0.0/
[Semantic Versioning]:
https://semver.org/spec/v2.0.0.html


## [1.3.6] - 2022-10-31

### Added
- Support for Python 3.11 (PR [#563], solves issue [#561]).
- Optional argument `toolset` in `GeosLibrary.build` method.

### Changed
- Upgrade `matplotlib` upper pin to 3.7.
- Upgrade `pyproj` upper pin to 3.5.

### Fixed
- Set MSVC 14.0 (VS2015) to build the precompiled Windows wheels in
  GitHub workflows (PR [#564]).

## [1.3.5] - 2022-10-25

### Fixed
- Fix broken array slicing inside `addcyclic` (PR [#559], solves issue
  [#555], thanks to @fragkoul).
- Fix `GeosLibrary` wrapper to also work with GEOS >= 3.7.0 on Windows
  and GNU/Linux.
- Fix wrong Antarctica coastline boundary with GEOS >= 3.9.0 (PR [#560],
  solves issue [#522]).

## [1.3.4] - 2022-08-10

### Changed
- Upgrade `numpy` upper pin to 1.24.
- Upgrade `pyshp` upper pin to 2.4.
- Upgrade `sphinx` upper pin to 5.0 and require at least Python 3.6 to
  build the docs.

### Fixed
- Update `numpy` build dependency to ensure that builds also work on
  MacOS (fixes issue [#547], thanks to @SongJaeIn for testing).
- Fix broken implementation of `Basemap.arcgisimage` (PR [#548], solves
  issue [#546]).
- Enforce up-to-date `numpy` dependency when possible:
  - Set `numpy >= 1.19` for Python == 3.6 due to `numpy` vulnerabilities
    [CVE-2021-41495] and [CVE-2021-41496].
  - Set `numpy >= 1.22` for Python >= 3.8 due to `numpy` vulnerability
    [CVE-2021-34141].
- Enforce up-to-date `pillow` dependency when possible:
  - Set `pillow >= 9.0.1` for Python >= 3.7 due to `pillow`
    vulnerability [CVE-2022-24303].

## [1.3.3] - 2022-05-11

### Changed
- Reformat `basemap.cm` using `flake8` and `black`.

### Fixed
- Fix issue in `drawcoastlines` with shape of vertices array
  (PR [#538] by @guziy, fixes issue [#512]).
- Fix setup to identify GEOS dylib on MacOS correctly (PR [#541],
  fixes issue [#539], thanks to @ronaldbradford and @CaffreyR for
  testing).

### Removed
- Remove dependency on `six` (PR [#537], fixes issue [#536]).

## [1.3.2] - 2022-02-10

### Added
- Add `"/usr/lib/x86_64-linux-gnu"` to list of candidate locations to
  search the GEOS library during installation.

### Changed
- Split lint and test requirements into two separate files.

### Fixed
- Fix setup encoding comment to deal with corner case under PowerShell.
- Enforce dependency `numpy >= 1.21` for Python >= 3.7 due to `numpy`
  vulnerability [CVE-2021-33430].
- Fix wrong marker for `unittest2` in development requirements.
- Fix `sdist` so that packages can be built from source distributions
  (PR [#532] by @DWesl, fixes [#533]).
- Specify Cython language level for `_geoslib` extension explicitly.
- Enforce up-to-date `pillow` dependency when possible:
  - `pillow >= 9.0.0` for Python >= 3.7 due to `pillow` vulnerabilities
    [CVE-2022-22815], [CVE-2022-22816] and [CVE-2022-22817].
  - `pillow >= 8.3.2` for Python >= 3.6 due to `pillow` vulnerabilities
    [CVE-2020-35653], [CVE-2020-35654], [CVE-2020-35655],
    [CVE-2021-23437], [CVE-2021-25287], [CVE-2021-25288],
    [CVE-2021-25290], [CVE-2021-25291], [CVE-2021-25292],
    [CVE-2021-25293], [CVE-2021-27921], [CVE-2021-27922],
    [CVE-2021-27923], [CVE-2021-28675], [CVE-2021-28676],
    [CVE-2021-28677], [CVE-2021-28678] and [CVE-2021-34552].
  - `pillow >= 7.1.0` for Python >= 3.5 due to `pillow` vulnerabilities
    [CVE-2020-10177], [CVE-2020-10378], [CVE-2020-10379],
    [CVE-2020-10994] and [CVE-2020-11538].
  - `pillow >= 6.2.2` For Python == 2.7 due to `pillow` vulnerabilities
    [CVE-2019-16865], [CVE-2019-19911], [CVE-2020-5310], [CVE-2020-5312]
    and [CVE-2020-5313].

### Removed
- Remove deprecation notices (issue [#527]).

## [1.3.1] - 2022-01-22

### Added
- Support for Python 3.10 (issues [#530] and [#531]).

### Changed
- Upgrade `numpy` upper pin to 1.23.
- Upgrade `matplotlib` upper pin to 3.6.
- Upgrade development requirements for Python 3.10.
- Move `doc` folder into `packages/basemap`.

### Fixed
- Fix error message when trying to load high- and full-resolution datasets
  without installing the `basemap-data-hires` package.

## [1.3.0] - 2021-12-28

### Added
- Precompiled binary wheels available in PyPI.
- Complete workflow to build the project wheels for Windows and GNU/Linux
  using GitHub Actions.

### Changed
- Reorganise the package structure. In summary, the former `basemap` package
  is split in three:
  - `basemap` itself contains the Python modules.
  - `basemap-data` contains the mandatory data assets required by `basemap`
    to provide minimal functionality.
  - `basemap-data-hires` contains the high-resolution data assets.

  This change together with the precompiled binary wheels in PyPI should solve
  most of the former installation problems (see issues [#403], [#405], [#422],
  [#436], [#445], [#456], [#461], [#488], [#489], [#491], [#510], [#513],
  [#525], [#526] and [#535]).
- Upgrade default GEOS library dependency to 3.5.1.
- Update and clarify licenses. In summary:
  - `basemap`: MIT license.
    - GEOS bundled dynamic library is under the LGPL-2.1-only license.
  - `basemap-data`: LGPL-3.0-or-later.
    - The EPSG file and the JPG images are also under the MIT license.
  - `basemap-data-hires`: LGPL-3.0-or-later.

### Fixed
- Fix `Basemap.pcolormesh` for `"ortho"` projection (PR [#476]).
- Fix `Basemap.arcgisimage` for cylindrical coordinates (PR [#505]).
- Force `setup.py` to cythonize `_geoslib.pyx` at compile time (issues [#487],
  [#518] and [#521]).
- Update `README` files and apply corrections and changes to outdated content
  (issue [#179]).

### Removed
- Bundled GEOS source code. The same source code can be downloaded using the
  `GeosLibrary` class in `utils` (issue [#228]).
- Precompiled `_geoslib.c` file.

## [1.2.2] - 2020-08-04

### Fixed
- Some incompatibilities with `matplotlib` v3.3+.
- Some incompatibilities with newer `libgeos` versions (tested against v1.6.1).
- Some incompatibilities with the `make.py`.

## [1.2.1] - 2019-08-08

### Added
- Some documentation updates.

### Fixed
- More compatibility bugfixes.
- Fix a bug introduced in v[1.1.0] in `addcyclic`.

## [1.2.0] - 2018-09-26

### Fixed
- Mostly compatibility bugfixes.
- Fix build using `zsh` as shell instead of `bash` (issues [#362] and [#383]).

## [1.1.0] - 2017-05-04

### Added
- `drawmapscale` supports `"feet"` as an input unit (PR #289).
- `nightshade` supports timezone-aware datetime objects as long as the
  timezone is UTC. Still assumes timezone-naive objects are UTC as well
  (issue #272).
- Add `"textcolor"` kwarg to `drawmeridians` and `drawparallels` (issue #145).
- Add `"facecolor"` keyword argument to `drawcounties` method; gives user
  ability to fill counties with specified `matplotlib` color argument.

### Changed
- Update packaged data to use GSHHG v2.3.6 (PR #311 & #320).
- Convert `pyshp` and `pyproj` into external dependencies (PR #234).
- Don't assume grid is regular when adding cyclic point in `addcyclic`.
  New kwargs `"axis"` and `"cyclic"` added. More than one array can be
  handled at a time, as long as the last one is longitude (issue #139).
- On non-Windows platforms, only link against libgeos C API (and not against
  C++ API anymore; issue #140). This is recommended by the authors/maintainers
  of GEOS. Also, this means that on e.g. Debian Linux, `basemap` can now be
  installed from source simply using `pip install basemap` when the `libgeos`
  packages are installed globally using the package management.

### Fixed
- Properly clip plots regardless of projection (issue #175).
- Properly perform input validation/casting (issue #260).
- Support newer versions of OWSLib (PR #259).
- Use integer indexes (PR #246).
- Stop `ResourceWarnings` for unclosed files in Python 3 (PR #244).
- Fix for coastline drawing glitch (issue #123).
- Fix `drawgreatcircle` bug so that lines exiting and reentering a projection
  region don't draw horizontally across the map.

## [1.0.7] - 2013-08-17

### Added
- Make `basemap` a namespace package (issue #114).
- Include `mpl_toolkits/__init__.py`, since the one installed by `matplotlib`
  is now inaccessible inside an egg (in 1.4.x).
- Support for rotated pole transformation (`projection = "rotpole"`).

### Changed
- Update `pyproj` (with fixes to geodesic calculations).

### Fixed
- Fix `drawmeridians` so meridians reach edge of plot when map projection
  region is *very* small (issue #113).
- Fix `warpimage` with `"hammer"` projection (issue #100).
- Fix tolerances for detecting jumps in meridians and parallels for
  `"cyl"` and `"rotpole"` projections (issue #108).
- Update `pyproj` (with fixes to geodesic calculations).
- Fix clipping to map projection region done in `contourf` and `contour`
  methods so it doesn't assume `x` and `y` are increasing (issue #110).

## [1.0.6] - 2013-01-13

### Added
- Add `"epsg"` keyword for defining projections.
- Add `"ellps"` keyword (`"rsphere"` ignored if `"ellps"` given).
- Add `"linestyle"` keyword to all draw* methods.
- Add support for cylindrical equal area (`"cea"`) projection.
- Add `drawcounties` method (PR #65). Thanks to Patrick Marsh.
- Add `arcgisimage` method for displaying background image retrieved from an
  ArcGIS server using the REST API (requires using `"epsg"` keyword to define
  projection).
- Add `wmsimage` method for displaying background image retrieved from
  an OGC-compliant WMS server using OWSLib (http://pypi.python.org/OWSLib)
  (requires using `"epsg"` keyword to define projection).
- Add `"latlon"` keyword to `plot` and `scatter` methods (PR #64).
- Add module variable `"latlon_default"` that can be used to switch default
  value of `"latlon"` kwarg to True, so that plotting methods can be passed
  `lats` and `lons` (geographic coordinates in degrees) instead of `x` and `y`
  (projection coordinates).

### Changed
- Update `pyproj` to version 1.9.3 (remove geographiclib Python code with
  C code from `PROJ.4`).
- Allow for latitude values slightly greater than 90 to exist in shapefiles,
  (by truncating to 90). Still raise exception if latitude exceeds 90.01.
- Make `drawcoastlines` use line segments instead of coastline polygons, to
  avoid *thickening* of lines around edges of map.

### Fixed
- Fix `drawcounties` for Python 3.3.
- Fix drawing of meridians and parallels in very small map regions (issue #79).
- Fix `shiftdata` method so it shifts mask along with data (PR #68).
- Fix bug that caused plotting to fail when `"latlon"` keyword is explicitly
  set to False (PR #66).
- Fix masking of grid cells outside the map projection in `contour` and
  `contourf`. In these methods, all points outside the map projection region
  were masked. This meant that if a grid cell was partly inside and partly
  outside the map projection region, nothing was drawn, leaving a gap along the
  edge of the map. This was particularly noticeable for very coarse resolution
  grids. This commit only masks those points more than one grid length beyond
  the edge of the map projection region (issue #88).

## [1.0.5] - 2012-08-06

### Added
- Add `"latlon"` keyword to plotting methods. If `latlon = True`, `x` and `y`
  values are assumed to longitudes and latitudes in degrees. The data and
  longitudes are shifted to the map projection region (for cylindrical and
  pseudo-cylindrical projections) using the `shiftdata` method, and lons/lats
  are converted to map projection coords. Default value is False. New example
  `shiftdata.py` added to illustrate usage (issue #54).

### Fixed
- Fix bug triggered when `drawlsmask` method was called more than once.
- Fix further corner cases with splitting of parallels that cross
  the dateline (issue #40).
- Fix error in `contour` method that caused a bogus mask to be applied
  to the data (issue #58).
- Fix `bluemarble` and `warpimage` methods to account for change in orientation
  of arrays returned to `matplotlib` function `pil_to_array` (issue #51).
- Fix glitch with drawing meridians and filling coastline polygons with
  `"omerc"` projection that includes pole.

## [1.0.4] - 2012-06-13

### Fixed
- Fix bug that caused Europe coastlines to disappear from some maps
  (evident from `nytolondon.py` example).
- Fix splitting of parallels in conic projections that cross the dateline
  (issue 40).

## [1.0.3] - 2012-05-18

### Added
- Add `"alpha"` keyword to `fillcontinents` (to set transparency).
- Add `"k_0"` keyword for `"tmerc"` projection so UTM zones can be created.
  Also can be used with `"lcc"`, `"omerc"` and `"stere"`. Add `utmtest.py`
  example.
- Add `streamplot` method, along with `streamplot_demo.py` example.
- Add `Basemap` attributes `boundarylats` abd `boundarylons` (arrays
  describing map boundaries; useful for illustrating map projection
  region on another map). Example illustrating usage `make_inset.py` added.
- Add `"round"` keyword to `Basemap.__init__` for pole-centered projections
  to make them round (clipped at `boundinglat`) instead of square.
- Add `hexbin` method, along with `hexbin_demo.py` example.
- Add `"kav7"` (Kavrayskiy VII) and `"eck4"` (Eckert IV) projections (PR #9).

### Changed
- Update GEOS from 3.3.1 to 3.3.3.
- Upgrade PROJ.4 source to version 4.8.0 and `pyproj` to version 1.9.2. New
  `pyproj` source from pyproj.googlecode.com includes more robust and accurate
  pure Python code for geodesic computations from `geographiclib`.
- Update coastlines, rivers, political boundaries to GSHHS 2.2.0/GMT 4.5.7.
  The `fillcontinents` bug (filling the outside instead of the inside of a
  coastline polygon) is now much harder to trigger.
- Make `drawmapboundary` use axes `"bgcolor"` as default `"fill_color"`. If
  no color fill is wanted, set `"fill_color"` to `'none'` (a string).

### Fixed
- Fix some more Python 3 compatibility issues (all examples now work with
  Python 3.2).
- Fix broken daynight terminator function.
- Bug in `drawparallels` that results in `KeyError` when drawing parallels
  very close together (0.1 degrees).
- Bugs in celestial projections and non-standard sphere radii (PR #6).
- Fix constant in Robinson projection (update `PJ_robin.c` from PROJ.4 SVN).
- Fix typo in `setup.py` (replace `"!= ['sdist','clean']"` with
  `"not in ['sdist','clean']"`).
- Clip coastlines for `"nplaea"`, `"npaeqd"`, `"splaea"`, `"spaeqd"` in
  stereographic coordinates to avoid South America disappearing in some
  south polar plots.
- Make sure `drawmeridians` can handle wrap-around (so that if projection is
  defined in -180 to 0 and user asks for meridians from 180 to 360 to be drawn,
  it should work). Only affects projections `"mill"`, `"gall"`, `"merc"` and
  `"cyl"`.

## 1.0.2 (git tag v1.0.2)

- Update included GEOS from 3.2.0 to 3.3.1 so it compiles with GCC 4.6.
- Add `colorbar` method that uses `axes_grid` toolkit to create colorbar axes.
- Replace `hasattr(arr,'mask')` with `numpy.ma.isMA(arr)`.
- Update docs and move to `matplotlib.github.com/basemap`.
- Add optional 1.25, 2.5 and 10 minute land/sea masks (derived from GSHHS
  coastline data). `"resolution"` and `"grid"` kwargs added to `drawlsmask` and
  maskoceans. `"resolution"` can be `"c"`, `"l"`, `"i"`, `"h"` or `"f"`,
  `"grid"` can be 1.25, 2.5, 5 or 10.
- Update `shapefile.py` from pyshp.googlecode.com to r72.
- Change default land-sea mask (now derived directly from GSHHS coastline data,
  default is 5 minutes use coastline resolution `"l"`). Default for plotting
  lakes is now True.
- Add `etopo` method (similar to `bluemarble`, but plots etopo relief image
  from www.ngdc.noaa.gov/mgg/global as a map background).
- Add `shadedrelief` method (similar to `bluemarble`, but plots shaded relief
  image from naturalearthdata.com as a map background).
- Replace `pyshapelib` with pure-Python `shapelib.py` from
  pyshp.googlecode.com. Allows full Python 3 compatibility.
- Fix `doc/conf.py` to use inheritance_diagram from Sphinx, not `matplotlib`.
- Fix `drawlsmask` so cylindrical projections work correctly when longitude
  range outside of -180 to 180.
- Python 3 compatibility.
- Add `lic_demo.py` to examples (line integral convolution, requires
  `scikit.vectorplot`).
- Remove deprecated `NetCDFFile`.
- Add `"zorder"` keyword to `drawmapscale`.
- Change default value for `"lakes"` kwarg in `drawlsmask` from False to True
  (**API change**).
- Change default value for `"inlands"` kwarg in `maskoceans` from False to
  True (**API change**).

## 1.0.1 (svn revision 8967)

- Regenerate C source with Cython 0.14.1.
- Add new `allsky` example from Tom Loredo.
- Add `"celestial"` keyword: if True, astronomical convention for longitude is
  used (negative for 'east', positive for 'west'); `celestial=True` implies
  `resolution=None` (no continents or coastlines).
- Improve placement of labels for parallels for pseudo-cylindrical projections
  like Mollweide and Hammer.
- Add support for Hammer projection (required adding inverse projection to
  PROJ.4 src in `src/PJ_hammer.c`).
- Update `src/pj_mutex.c` from PROJ.4 svn to fix a threading bug on Windows.
- If you try to transform NaNs to/from map projection coords, 1.e30 is returned
  (previously, this caused a segfault for some projections).
- Deprecate `NetCDFFile` function, will be removed in 1.0.2. Issue warning
  advising users to use `netcdf4-python` instead.
- Deleting an item from the dicts returned by `drawparallels` and
  `drawmeridians` removes the corresponding parallel or meridian (and
  associated labels) from the plot.
- Add a `remove` method to the tuples that are returned in the dicts returned
  by `drawparallels` and `drawmeridians`.

## 1.0.0 (svn revision 8531)

- Don't force `adjustable="box"` so `Basemap` is compatible with `AxesGrid`.
  Add `fcstmaps_axesgrid.py` example.
- Add support for plotting on unstructured grids using keyword `"tri"` in
  `pcolor`, `contour`, and `contourf` methods (which then forward to
  `tripcolor`, `tricontour`, and `tricontourf` axes methods).
  `examples/ploticos.py` added.
- Let continents that fill the whole map be filled.
- Add option for cubic spline interpolation in `interp` function (order=3)
  using `scipy.ndimage`.
- Add "near-sided perspective" projection for a satellite view at an
  arbitrary altitude.
- Patch from Stephane Raynaud to pass format string to `drawmapscale`, and
  allow `units="m"`.
- Update PROJ.4 source to version 4.7.0, `pyproj` to 1.8.6.
- Add `is_land` method to check whether a point is over land or water.
- GEOS-3.1.1 now required. 3.2.0 source included (instead of 2.2.3).
- `shiftgrid` no longer requires a cyclic point to be present (patch from Eric
  Bruning).
- Fix `date2index` bugs.
- Update `date2index` function with a bugfix from `netcdf4-python`.
- In `contourf` method, mask data outside map projection region (this prevents
  `contourf` from erroneously filling entire map).
- Add `nightshade` method to shade night regions on a map. `daynight.py`
  example added to illustrate usage.
- Add `lonmin` and `lonmax` instance variables.

## 0.99.4 (svn revision 7332)

- Replace `ax.frame` with `ax.spines` to maintain compatibility with
  `matplotlib` spines support.
- Add `"latmax"` kwarg to `drawparallels` and `drawmeridians` (patch from Chris
  Murphy).
- Add new example `plotmap_shaded.py` (shaded relief plot).
- Add new example `plothighsandlows.py`.
- Add `"fix_aspect"` kwarg to `Basemap.__init__`, when False `axes.set_aspect`
  is set to `"auto"` instead of default `"equal"`. Can be used to make plot
  fill whole plot region, even if the plot region doesn't match the aspect
  ratio of the map region.
- Add `date2index` function.
- Update `netcdftime` to 0.7.1.
- Add `maskoceans` function.
- Update `pupynere` to version 1.0.8 (supports writing large files).
- Add more informative error message in `readshapefile` when one of the
  shapefile components can't be found.

## 0.99.3 (svn revision 6780)

- If upper-right/lower-left corners nor width/height given for azimuthal
  equidistant (`"aeqd"`), the whole world is drawn in a circle (only works for
  perfect spheres, not ellipsoids).
- Make `setup.py` check for already installed `pyshapelib` (just like it does
  for `httplib2` and `pydap`).
- `Basemap` will now look for its data in `BASEMAPDATA`. If that environment
  variable is not set, it will fall back to its default location.
- If `readshapefile` is called with `drawbounds=True`, a `LineCollection`
  object is appended to the returned tuple.
- Make sure `drawmapscale` method returns a list of objects that can be
  iterated over to remove them from the plot.
- `fillcontinents` was returning just last `Polygon` instance. Now returns a
  list of all `Polygon` instances.
- Pass `bluemarble`/`warpimage` kwargs to `imshow` and return `Image` instance.

## 0.99.2 (svn revision 6541)

- Fix `drawlsmask` method so that it works for cylindrical projections with
  limits outside (-180, 180).
- Add `"scale"` keyword to `bluemarble` and `warpimage` methods to downsample
  image background.
- Make `"lat_ts"` default to 0 for Mercator.
- Now can specify just `lon_0` for all cylindrical projections (to produce
  global map centered on `lon_0`).
- Add `save_background.py` example, showing how to re-use a map background
  without redrawing coastlines.
- Add `embedding_map_in_wx.py` example (courtesy of Mauro Cavalcanti).
- Add masked array support to `shiftgrid` function (thanks to Jesper Larsen).
- Defer import of netcdf stuff till it is needed (in `NetCDFFile` function).
- Add McBryde-Thomas Flat Polar Quartic (`projection = "mbtfpq"`), Gall
  Stereographic Cylindrical (`projection = "gall"`) and van der Grinten
  (`projection = "vandg"`).
- Fix bugs in `warpimage` and `bluemarble` methods for several projections.
- Bugfix patch for `rotate_vector` from David Huard. David also contributed the
  beginnings of a test suite.
- Make sure scatter method sets pyplot color mappable.
- Add `cubed_sphere` example.
- Update `NetCDFFile` to use `pupynere` 1.0.2 (now can write as well as read!).
- Now works with GEOS version 3.
- Add `Basemap` instance variable `proj4string`.
- `testgdal.py` example now uses gdal to read topo data from a raster DEM file
  and ogr to read state boundaries from a shape file.
- `warpimage` method can now handle gray-scale images, and images specifed as
  URLs (for example, the Blue Marble images from
  http://earthobservatory.nasa.gov/Newsroom/BlueMarble/BlueMarble_monthlies.html).

## 0.99.1 (svn revision 5961)

- GEOS-2.2.3 patched for compatibility with GCC 4.3.
- Add `barbs` method to draw wind barbs on the map.
- Add `tissot` method for generating Tissot's indicatrix (see example
  `plot_tissot.py`).
- Fix processing of coastlines for gnomonic projection.
- Don't try to use PyNIO in `NetCDFFile` (it was causing too many suprises).
- Start of improved documentation using Sphinx/docutils. Can be viewed at
  http://matplotlib.sf.net/basemap/doc/html
- Change default behaviour of `num2date` and `date2num` to be the same as
  `matplotlib` counterparts.

## 0.99.0 (svn revision 5344)

- Fix to `warpimage` method for API change in `matplotlib` 0.98.0.
- Update `pyproj` to 1.8.6.
- Fix bug in `NetCDFFile` creating masked arrays when both `_FillValue` and
  `missing_value` exist.
- `drawparallels` and `drawmeridians` return a dictionary containing the
  `Line2D` and `Text` instances associated with each lat or lon.
- `drawcoastlines`, `drawcountries` and friends now have `PatchCollection`
  return values.
- Make sure `_nolabel_` set on coastlines, countries, states, rivers, parallels
  and meridians so they are not included in a legend.
- Add `drawmapscale` method to create a map scale bar similar to that available
  with the GMT's psbasemap.
- Now lives in `mpl_toolkits.basemap` (**API change**). Instead of
  `from matplotlib.toolkits.basemap import Basemap`, use
  `from mpl_toolkits.basemap import Basemap`. All examples changed. Uses
  `matplotlib mpl_toolkits` namespace package, so `basemap` can now be
  installed if `matplotlib` is installed as an egg. Python 2.3 support
  re-enabled.
- Change `_geos` to `_geoslib`, so as not to conflict with the Python module
  bundled with the GEOS library.
- Some fixes/enhancements for `"omerc"` projection (added `"no_rot"` flag).
- Add `warpimage` method for displaying an image background. Default is NASA's
  blue marble image, which is included.

## 0.9.9.1 (svn revision 4808)

- Require Python 2.4 (really only needed for building). Once namespace packages
  are re-enabled in matplotlib, Python 2.3 should work again.

## 0.9.9 (svn revision 4799)

- Updated PROJ.4 sources to version 4.6.0.
- Remove hidden dependency on `setuptools` (in `dap` module).
- Fix exception handling bug in code that looks for intersection between
  boundary feature and map projection region.
- `setup.py` now looks for GEOS library in a few standard places (/usr/local,
  /opt, $HOME, /sw) if the `GEOS_DIR` environment variable is not set. This is
  a workaround for a new Leopard 'feature' (sudo does not inherit enviroment
  variables).
- Add support for reading `Point` and `MultiPoint` shapes from ESRI shapefiles.
- Now automatically draws figure if running in interactive mode (so `draw()`
  does not need to be called explicitly in IPython).
- Add `num2date` and `date2num` functions, which use included `netcdftime`
  module.

## 0.9.8 (svn revision 4526)

- Fixes for filling continents in orthographic projection.
- Add `"maskandscale"` kwarg to `NetCDFFile` to optionally turn off automatic
  masking and rescaling of variable data.
- `NetCDFFile` will try to use PyNIO if installed and the file cannot be read
  with `pupynere`. This allows GRIB1, GRIB2, HDF4 and HDFEOS2 files to be read.
- `"fmt"` kwarg to `drawparallels` and `drawmeridians` can now be a custom
  string formatting function (example `customticks.py` demonstrates usage).
- Remove `"linestyle"` kwarg from `drawparallels` and `drawmeridians` (it never
  did anything anyway since it was overridden by the `"dashes"` kwarg).
- Modify `NetCDFFile` to use `dap` module to read remote datasets over http.
  Include `dap` and `httplib2` modules.
- Modify `NetCDFFile` to automatically apply `"scale_factor"` and
  `"add_offset"`, and return masked arrays masked where `data == missing_value`
  or `data == _FillValue`.
- Add `"fill_color"` option to `drawmapboundary`, to optionally fill the map
  projection background a certain color.
- Add `"sstanom"` colormap from
  http://www.ghrsst-pp.org/GHRSST-PP-Data-Tools.html

## 0.9.7 (svn revision 4422)

- Fix bug in `drawlsmask` for `"moll"`, `"robin"` and `"sinu"` projections.
- Add `"lake_color"` keyword to `fillcontinents`.
- Fix a bug in the `"tmerc"` projection.
- Add pure-Python `NetCDFFile` reader from Roberto De Almeida to `basemap`
  namespace (from `matplotlib.toolkits.basemap import NetCDFFile`).
- Add support for full-resolution boundaries (will be a separate download).
  Full-res files (totaling around 100 MB) available in SVN.
- High-resolution boundaries now included.
- Postpone processing of countries, states and river boundaries until a draw is
  requested. Only the coastlines are processed in `__init__`.
- Use a Cython interface to the GEOS library (http://geos.refractions.net,
  LGPL-2.1 license) to find geometries that are within map projection region.
  This speeds up instance creation for small map regions and high resolution
  coastlines. Boundary datasets now in binary format (I/O is faster). Requires
  GEOS version 2.2.3, source code included.
- Remove all numerix imports.
- Fix `rotate_vector` so it works in S. Hem and for non-orthogonal grids.
  Support for masked velocity vectors also added. (EF)
- Numpification. (EF)

## 0.9.6 (svn revision 3888)

- Fix `addcyclic` function so it handles masked arrays.
- Labelling of meridians and parallels now works with very small map regions
  (less than 0.2 degrees square).
- Subregions of the globe may be specified with
  `llcrnrlat, llcrnrlon, urcrnrlat, urcrnrlon` keywords for `"ortho"` and
  `"geos"` (illustrated by `examples/geos_demo_2.py`).
- Add `"labelstyle"` keyword to `drawparallels` and `drawmeridians`. If set to
  `"+/-"`, labels are given prefixed by `"+"` or `"-"`, instead of suffixed
  with `"N"`, `"S"`, `"E"` or `"W"`. Useful for astronomical plots, where there
  is no such thing as north, south, east or west.
- Add support for geostationary satellite projection (`projection = "geos"`),
  contributed by Scott Sinclair. Illustrated by `examples/geos_demo.py`.
- Add a bunch of extra colormaps (mostly from GMT), and a script to plot them
  (`examples/show_colormaps.py`). To import new colormaps, use
  `from matplotlib.toolkits.basemap import cm`.
- Orthographic projection only defined for perfect sphere, raise an error if
  user tries to use an ellipsoid.
- Print a warning in `contour` and `contourf` in situations that may
  result in a screwy looking plot (`x` not monotonically increasing,
  because the data wraps around the edge of the plot). The warning
  suggests using the `shiftgrid` function to recenter the data on the
  map projection region.
- Fix `setup.py` so it works properly with Python 2.3, and on Windows.
- Add `"zorder"` keyword to `drawmeridians`, `drawparallels`, `drawcoastlines`,
  `drawstates`, `drawcountries`, `drawrivers`, `readshapefile` and
  `fillcontinents`.
- `numpy` now required.
- Added `"srs"` (spatial reference system) instance variable.
- Update `pyproj` to version 1.8.3. Now uses `pyproj.Geod` for Great Circle
  calculations. PROJ.4 data files now included.
- Make sure axes ticks are always turned off, unless `noticks=False` is set
  when creating a `Basemap` instance.

## 0.9.5 (svn revision 3083) - 2007-03-14

- Fix examples to conform to 'one show() per script' rule.
- Intermediate coastlines now installed by default. `basemap-data` is no longer
  a separate package (couldn't figure out how to manage the egg). If the `"h"`
  res boundaries are needed, the data files must be manually put in place by
  the user. `BASEMAP_DATA_PATH` environment variable is no longer used.
- Reorganize data files so that `bdist_egg` includes the data. `setup-data.py`
  now used to install the high-res data files.
- Make sure PROJ.4 returns 1.e30 instead of `HUGE_VAL` for undefined
  transformations (`HUGE_VAL` is inf on most platforms, which gets
  embedded in the postscript, resulting in a un-renderable file).
- Use `typedef Py_ssize_t` if necessary (ADS).
- Rename `pyproj.so` to `_pyproj.so`, add a Python wrapper `pyproj.py` and
  move `pyproj` into `matplotlib.toolkits.basemap`.

## 0.9.4 - 2006-11-17

- Add "Tissot's indicatrix" example.
- `"extent"` keyword was erroneously being passed to `pcolor`.
- Update PROJ.4 source files to version 4.5.0.
- Update `pyproj` to version 1.8.0 (better error handling).

## 0.9.3 - 2006-10-17

- Update `pyproj.c` to be compatible with Python 2.5.
- Add new example `ccsm_popgrid.py` (contributed by Ivan Lima).

## 0.9.2 - 2006-08-31

- Fix several bugs in `drawlsmask` method.
- Remove buggy optimizations for cylindrical projections not crossing the
  Greenwich meridian.
- Can now specify map projection region in `Basemap.__init__` by setting
  width and height in projection coordinates (in meters) instead of
  specifying lat/lon of upper-right and lower-left corners (**API change**).

## 0.9.1 - 2006-07-27

- Make sure `llcrnrlat` and `llcrnrlon` are not at poles for Mercator.
- Use Eric Firing's new `quiver` in `Basemap.quiver` method.
- `interp` functions now work with masked arrays.
- Add some sanity checks for projection parameters.
- Change from classic to new-style classes.
- Remove deprecated `createfigure` method.
- Fix some creeping `numpy`'isms (which caused breakage when `numarray` or
  `Numeric` were used).

## 0.9.0 - 2006-06-09

- Update for new `matplotlib` aspect ratio handling. Now maps will always have
  the correct aspect ratio.
- If `"resolution"` keyword is set to None when a `Basemap` instance is
  created, no boundary data sets are needed (methods to draw boundaries, like
  `drawcoastlines`, will raise an exception).
- Rename `proj4` module to `pyproj` to avoid conflicts with `proj4` module
  from CDAT.
- Deprecate `createfigure` method, since maps will now automatically have the
  correct aspect ratio.
- Add new projections Xpstere, Xplaea, Xpaeqd (where X can be n or s). These
  are special-case, polar-centric versions of the stereographic, lambert
  azimuthal equal area and azimuthal equidistant projections that don't require
  you specify the lat/lon values of the lower-left and upper-right corners.
- Fix bugs in `plot`, `scatter` and `mapboundary` methods for Miller,
  cylindrical and Mercator projections.
- `"crude"` and `"low"` resolution boundary datasets now installed by default.
  `basemap_data` package now only needed to get `"intermediate"` and `"high"`
  resolution datasets.
- Move all packages under single `lib/` directory so that setuptools' "develop"
  command works properly.
- Add sinusoidal projection.
- Bilinear interpolation routines can return masked arrays with values outside
  range of data coordinates masked.
- New examples:
  - `warpimage.py`: warping an image to different map projections.
  - `polarmaps.py`: simplified polar projections.
  - `garp.py`: 'World According to Garp' maps.
- Add `pcolormesh` method.
- Add `drawlsmask` method for masking oceans and/or land areas.
- Add 5-minute land-sea mask dataset.

## 0.8.2 - 2006-02-22

- Minor bugfixes, mostly in examples.

## 0.8.1 - 2006-02-03

- Huge speedups for `numpy` (no significant differences for `Numeric` and
  `numarray`).

## 0.8.0 - 2006-01-14

- Add `numpy` compatibility.

## 0.7.2.1 - 2005-11-18

- There was a problem running examples that read pickle files. The pickle files
  were created with `numarray`, and the data would not be read correctly using
  `Numeric`. Fixed so that pickles are created with `Numeric` and `Numeric` is
  used to read them.

## 0.7.2 - 2005-10-18

- No longer requires `numarray` (`interp` function no longer uses
  `numarray.nd_image`). This means that `interp` does not accept `"mode"` and
  `"cval"` any longer (**API change**). `"order"` keyword must be 0 or 1.
- Modify to work with the new `ContourSet` returned by `contour` and
  `contourf`.
- Turn off axes frame by default for non-rectangular projections (`"ortho"`,
  `"robin"` and `"moll"`).
- Add `createfigure` method to create a figure with the same aspect ratio as
  the map using `pylab.figure`.
- Reset subplot.params defaults so that default axes rectangle will have both
  a width and height of 0.9 (this ensures that the figure size determines that
  aspect ratio of the plot).
- Make `readshapefile` method raise an exception if the vertices look like they
  are not in geographic (lat/lon) coordinates.

## 0.7.1 - 2005-09-21

- Fix several bugs in meridian/parallel labelling and cylindrical projections
  that crossed Greenwich were not being handled properly.
- Add `"fmt"` kwarg to `drawmeridians` and `drawparallels` (default is `"%g"`).
- Fix bug in `readshapefile` that prevented boundaries from being drawn for
  `"cyl"`, `"merc"` or `"miller"` projections when the map region did not cross
  the Greenwich meridian.
- Modify `imshow` method so `"origin"` keyword is accepted (it was always set
  to `"lower"` previously).
- Add `testgdal.py` example showing how to plot raster geospatial data with
  `gdal` module (gdal.maptools.org).
- Meridians and parallels labelled correctly when
  `rcParams["text.usetext"] = True`.

## 0.7.0 - 2005-09-14

- Optimizations to reduce the time it takes to create a `Basemap` class
  instance (now nearly 4 times faster when using `resolution="i"`).
- Add `"h"` (high) resolution boundary data.
- Add datasets for major rivers, "drawrivers" class method.
- Fix some errors in boundaries datasets.
- Boundary datasets now installed in a separate package.
- Should now handle Numeric to numarray conversions internally, so removed
  warning when `rcParams["numerix"] != "numarray"`.
- Change default `"area_thresh"` so it depends on coastline resolution
  (10000 for `"c"` declining to 10 for `"h"`).

## 0.6.2 - 2005-09-01

- Warning issued if numerix = 'Numeric' (a user reported crashes due to botched
  Numeric --> numarray conversions).
- Changes to PROJ.4 wrapper to make `Basemap` instances pickle-able.

## 0.6.1 - 2005-08-14

- Add `hurrtracks.py` example (plot hurricane tracks from shapefile).
- Now includes `pyshapelib` (extracted from Thuban source).

## 0.6.0 - 2005-08-11

- SF bug #1254163 (make sure lat/lon limits of projection are floats).
  `wiki_example.py` and `plotclimdiv.py` examples added.
- Add `readshapefile` method for reading and plotting data from ESRI shapefiles
  (requires `pyshapelib` from Thuban). `plotclimdiv.py` is an example that
  illustrates this.

## 0.5.2 - 2005-06-28

- Fix bug in meridian labelling when lon > 360 or lon < -180.
- Add `"ax"` keyword to `Basemap.__init__`. This will set default axis
  instance, which can be overridden by using `"ax"` keyword in method calls
  (**API change**).

## 0.5.1 - 2005-06-26

- Add `"ax"` keyword to most `Basemap` methods to allow use of a pre-existing
  `Axes` instance. Default is still to use the current instance.
- Full control of font properties for parallel and meridian labels (now uses
  unicode instead of mathtext for degree symbol). Replace `"font"` and
  `"fontsize"` keyword args for `drawparallels` and `drawmeridians` replaced by
  `**kwargs`, which is passed directly to `Axes.text` method (**API change**).

## 0.5.0 - 2005-06-02

- Add Orthographic, Mollweide and Robinson projections.
- Add `drawmapboundary` method to draw a line around the map projection
  region.
- Add `"suppress_ticks"` keyword to `Basemap.__init__` It's True by default,
  but can be set to False if you want to label ticks in native map projection
  coordinates.
- Add `rotate_vector` method to rotate vectors to map projection coordinates
  (without interpolation, as in `transform_vector` method).
- Modified `pcolor`, `contour`, `contourf` methods to use masked arrays.
- Now requires matplotlib v0.81.
- `drawparallels` and `drawmeridians` methods now take optional keyword
  arguments `"xoffset"` and `"yoffset"`, which control how far from the edge of
  the map labels are drawn.
- Make `llcrnrlon`, `llcrnrlat`, `urcrnrlon` and `urcrnrlat` optional
  keyword arguments in `Basemap.__init__` (**API change**).

## 0.4.3 - 2005-05-11

- Add Oblique Mercator.
- Great circle calculations now use Vincenty's equations for an ellipsoid.

## 0.4.2 - 2005-05-10

- Remove `"preserve_magnitude"` keyword from `transform_vector`. Now
  `transform_vector` does a simple rotation of the vector from geographic
  to map coordinates, preserving the vector magnitude (**API change**).
- Fix minor bugs in Miller and Mercator projections.
- Add Gnomonic, Cassini-Soldner and Polyconic projections (now 13 projections
  supported).

## 0.4.1 - 2005-05-09

- Fix Miller projection being erroneously referred to by the name `"miller"`
  instead of `"mill"`.
- Add the ability to specify the major and minor sphere radii by specifying the
  `"rsphere"` keyword in `__init__` to be a tuple instead of a scalar.

## 0.4.0 - 2005-05-05

- Add support for miller cylindrical, equidistant conic, and azimuthal
  equidistant projections. 
- Fix bugs in coastline drawing and continent filling methods.

## 0.3.3 - 2005-04-27

- Modify `fillcontinents` to not fill lakes (they are actually still filled,
  but with axis background color).

## 0.3.2 - 2005-04-20

- Code cleanups.
- Docstring typo fixes.
- Environment variable `BASEMAP_DATA_PATH` can now be used to point to data
  files in a non-standard install (i.e using `--prefix` or `--home`).

## 0.3.0 - 2005-04-15

- Add `transform_scalar` and `transform_vector` methods for interpolating
  scalar and vector fields to a map projection grid. 
- Add `shiftgrid` and `addcyclic` convenience functions.
- Add `quiver_demo.py` example illustrating how to plot wind vectors on a map.
- Update examples to use `transform_scalar` instead of calling `interp`
  directly.
- Change Mercator `x` coordinate units from degrees to meters.
- `setup.py` now installs data in Python version-numbered directory (so you
  can have separate copies for different Python versions).
- Fix aspect ratio of mercator plots.
- Add `set_axes_limits`, `plot`, `scatter`, `contourf`, `contour`, `pcolor` and
  `quiver` methods.
- `axes` instance no longer a method argument to any `Basemap` method,
  `gca` is called to obtain the current axes instance instead (**API change**).

## 0.2.1 - 2005-04-10

- Add `gcpoints` and `drawgreatcircles` methods.
- Add intermediate resolution coastline and political boundary databases.
- Fix bug in filling continents.
- Fix bug in drawing parallels/meridians.
- Add `ireland.py` example.

## 0.2.0 - 2005-04-04

- `drawparallels` and `drawmeridians` can now draw labels at edge of map.

## 0.1.2 - 2005-03-31

- Now can handle negative longitudes (patch from Michael Brady).
- `basemap.interp` can now handle irregular (but still rectilinear) lat/lon
  grids.

## 0.1.0 - 2005-02-03

- First release on SF.

## 2005-02-02

- Change `LineCollections.color` to `set_color` in a try/except block.
  (color was a typo in 0.71 and is deprecated).

## 2005-01-30

- No user visible changes. Uses new pyrex generated C extension interface
  to PROJ.4 which is twice as fast as the Thuban one.

## 2005-01-24

- Some code reorganisation.
- Fix bugs in S. Hem. projections.
- `Basemap` instance variables `xmin, xmax, ymin, ymax` renamed to
  `llcrnrx, urcrnrx, llcrnry, urcrnry`.
- Add `__call__` and makegrid methods to `Basemap` class.

## 2005-01-21

- PROJ.4 is now called via a C-library interface, instead of using the proj
  command-line tool via `os.popen`.

## 2005-01-19

- Fix glitches in drawing of parallels and meridians.


[#564]:
https://github.com/matplotlib/basemap/pull/564
[#563]:
https://github.com/matplotlib/basemap/pull/563
[#561]:
https://github.com/matplotlib/basemap/issues/561
[#560]:
https://github.com/matplotlib/basemap/pull/560
[#559]:
https://github.com/matplotlib/basemap/pull/559
[#555]:
https://github.com/matplotlib/basemap/issues/555
[#548]:
https://github.com/matplotlib/basemap/pull/548
[#547]:
https://github.com/matplotlib/basemap/issues/547
[#546]:
https://github.com/matplotlib/basemap/issues/546
[#541]:
https://github.com/matplotlib/basemap/pull/541
[#539]:
https://github.com/matplotlib/basemap/issues/539
[#538]:
https://github.com/matplotlib/basemap/pull/538
[#537]:
https://github.com/matplotlib/basemap/pull/537
[#536]:
https://github.com/matplotlib/basemap/issues/536
[#535]:
https://github.com/matplotlib/basemap/issues/535
[#533]:
https://github.com/matplotlib/basemap/issues/533
[#532]:
https://github.com/matplotlib/basemap/pull/532
[#531]:
https://github.com/matplotlib/basemap/issues/531
[#530]:
https://github.com/matplotlib/basemap/issues/530
[#527]:
https://github.com/matplotlib/basemap/issues/527
[#526]:
https://github.com/matplotlib/basemap/issues/526
[#525]:
https://github.com/matplotlib/basemap/issues/525
[#522]:
https://github.com/matplotlib/basemap/issues/522
[#521]:
https://github.com/matplotlib/basemap/issues/521
[#518]:
https://github.com/matplotlib/basemap/issues/518
[#513]:
https://github.com/matplotlib/basemap/issues/513
[#512]:
https://github.com/matplotlib/basemap/issues/512
[#510]:
https://github.com/matplotlib/basemap/issues/510
[#505]:
https://github.com/matplotlib/basemap/pull/505
[#491]:
https://github.com/matplotlib/basemap/issues/491
[#489]:
https://github.com/matplotlib/basemap/issues/489
[#488]:
https://github.com/matplotlib/basemap/issues/488
[#487]:
https://github.com/matplotlib/basemap/issues/487
[#476]:
https://github.com/matplotlib/basemap/pull/476
[#461]:
https://github.com/matplotlib/basemap/issues/461
[#456]:
https://github.com/matplotlib/basemap/issues/456
[#445]:
https://github.com/matplotlib/basemap/issues/445
[#436]:
https://github.com/matplotlib/basemap/issues/436
[#422]:
https://github.com/matplotlib/basemap/issues/422
[#405]:
https://github.com/matplotlib/basemap/issues/405
[#403]:
https://github.com/matplotlib/basemap/issues/403
[#383]:
https://github.com/matplotlib/basemap/issues/383
[#362]:
https://github.com/matplotlib/basemap/issues/362
[#228]:
https://github.com/matplotlib/basemap/issues/228
[#179]:
https://github.com/matplotlib/basemap/issues/179

[Unreleased]:
https://github.com/matplotlib/basemap/compare/v1.3.6...develop
[1.3.6]:
https://github.com/matplotlib/basemap/compare/v1.3.5...v1.3.6
[1.3.5]:
https://github.com/matplotlib/basemap/compare/v1.3.4...v1.3.5
[1.3.4]:
https://github.com/matplotlib/basemap/compare/v1.3.3...v1.3.4
[1.3.3]:
https://github.com/matplotlib/basemap/compare/v1.3.2...v1.3.3
[1.3.2]:
https://github.com/matplotlib/basemap/compare/v1.3.1...v1.3.2
[1.3.1]:
https://github.com/matplotlib/basemap/compare/v1.3.0...v1.3.1
[1.3.0]:
https://github.com/matplotlib/basemap/compare/v1.2.2rel...v1.3.0
[1.2.2]:
https://github.com/matplotlib/basemap/compare/v1.2.1rel...v1.2.2rel
[1.2.1]:
https://github.com/matplotlib/basemap/compare/v1.2.0rel...v1.2.1rel
[1.2.0]:
https://github.com/matplotlib/basemap/compare/v1.1.0...v1.2.0rel
[1.1.0]:
https://github.com/matplotlib/basemap/compare/v1.0.7rel...v1.1.0
[1.0.7]:
https://github.com/matplotlib/basemap/compare/v1.0.6rel...v1.0.7rel
[1.0.6]:
https://github.com/matplotlib/basemap/compare/v1.0.5rel...v1.0.6rel
[1.0.5]:
https://github.com/matplotlib/basemap/compare/v1.0.4rel...v1.0.5rel
[1.0.4]:
https://github.com/matplotlib/basemap/compare/v1.0.3rel...v1.0.4rel
[1.0.3]:
https://github.com/matplotlib/basemap/tree/v1.0.3rel

[CVE-2022-24303]:
https://nvd.nist.gov/vuln/detail/CVE-2022-24303
[CVE-2022-22817]:
https://nvd.nist.gov/vuln/detail/CVE-2022-22817
[CVE-2022-22816]:
https://nvd.nist.gov/vuln/detail/CVE-2022-22816
[CVE-2022-22815]:
https://nvd.nist.gov/vuln/detail/CVE-2022-22815
[CVE-2021-41496]:
https://nvd.nist.gov/vuln/detail/CVE-2021-41496
[CVE-2021-41495]:
https://nvd.nist.gov/vuln/detail/CVE-2021-41495
[CVE-2021-34552]:
https://nvd.nist.gov/vuln/detail/CVE-2021-34552
[CVE-2021-34141]:
https://nvd.nist.gov/vuln/detail/CVE-2021-34141
[CVE-2021-33430]:
https://nvd.nist.gov/vuln/detail/CVE-2021-33430
[CVE-2021-28678]:
https://nvd.nist.gov/vuln/detail/CVE-2021-28678
[CVE-2021-28677]:
https://nvd.nist.gov/vuln/detail/CVE-2021-28677
[CVE-2021-28676]:
https://nvd.nist.gov/vuln/detail/CVE-2021-28676
[CVE-2021-28675]:
https://nvd.nist.gov/vuln/detail/CVE-2021-28675
[CVE-2021-27923]:
https://nvd.nist.gov/vuln/detail/CVE-2021-27923
[CVE-2021-27922]:
https://nvd.nist.gov/vuln/detail/CVE-2021-27922
[CVE-2021-27921]:
https://nvd.nist.gov/vuln/detail/CVE-2021-27921
[CVE-2021-25293]:
https://nvd.nist.gov/vuln/detail/CVE-2021-25293
[CVE-2021-25292]:
https://nvd.nist.gov/vuln/detail/CVE-2021-25292
[CVE-2021-25291]:
https://nvd.nist.gov/vuln/detail/CVE-2021-25291
[CVE-2021-25290]:
https://nvd.nist.gov/vuln/detail/CVE-2021-25290
[CVE-2021-25288]:
https://nvd.nist.gov/vuln/detail/CVE-2021-25288
[CVE-2021-25287]:
https://nvd.nist.gov/vuln/detail/CVE-2021-25287
[CVE-2021-23437]:
https://nvd.nist.gov/vuln/detail/CVE-2021-23437
[CVE-2020-35655]:
https://nvd.nist.gov/vuln/detail/CVE-2020-35655
[CVE-2020-35654]:
https://nvd.nist.gov/vuln/detail/CVE-2020-35654
[CVE-2020-35653]:
https://nvd.nist.gov/vuln/detail/CVE-2020-35653
[CVE-2020-11538]:
https://nvd.nist.gov/vuln/detail/CVE-2020-11538
[CVE-2020-10994]:
https://nvd.nist.gov/vuln/detail/CVE-2020-10994
[CVE-2020-10379]:
https://nvd.nist.gov/vuln/detail/CVE-2020-10379
[CVE-2020-10378]:
https://nvd.nist.gov/vuln/detail/CVE-2020-10378
[CVE-2020-10177]:
https://nvd.nist.gov/vuln/detail/CVE-2020-10177
[CVE-2020-5313]:
https://nvd.nist.gov/vuln/detail/CVE-2020-5313
[CVE-2020-5312]:
https://nvd.nist.gov/vuln/detail/CVE-2020-5312
[CVE-2020-5310]:
https://nvd.nist.gov/vuln/detail/CVE-2020-5310
[CVE-2019-19911]:
https://nvd.nist.gov/vuln/detail/CVE-2019-19911
[CVE-2019-16865]:
https://nvd.nist.gov/vuln/detail/CVE-2019-16865
