from __future__ import (absolute_import, division, print_function)

"""
Module for plotting data on maps with matplotlib.

Contains the :class:`Basemap` class (which does most of the
heavy lifting), and the following functions:

:func:`interp`:  bilinear interpolation between rectilinear grids.

:func:`maskoceans`:  mask 'wet' points of an input array.

:func:`shiftgrid`:  shifts global lat/lon grids east or west.

:func:`addcyclic`: Add cyclic (wraparound) point in longitude.
"""
from distutils.version import LooseVersion

try:
    from urllib import urlretrieve
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlretrieve, urlopen

from matplotlib import __version__ as _matplotlib_version
try:
    from inspect import cleandoc as dedent
except ImportError:
    # Deprecated as of version 3.1. Not quite the same
    # as textwrap.dedent.
    from matplotlib.cbook import dedent
# check to make sure matplotlib is not too old.
_matplotlib_version = LooseVersion(_matplotlib_version)
_mpl_required_version = LooseVersion('0.98')
if _matplotlib_version < _mpl_required_version:
    msg = dedent("""
    your matplotlib is too old - basemap requires version %s or
    higher, you have version %s""" %
    (_mpl_required_version,_matplotlib_version))
    raise ImportError(msg)
from matplotlib import rcParams, is_interactive
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Ellipse, Circle, Polygon, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox
import pyproj
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.image import imread
import sys, os, math
from .proj import Proj
import numpy as np
import numpy.ma as ma
import _geoslib
import functools


__version__ = "1.3.8"

# basemap data files now installed in lib/matplotlib/toolkits/basemap/data
# check to see if environment variable BASEMAPDATA set to a directory,
# and if so look for the data there.
if 'BASEMAPDATA' in os.environ:
    basemap_datadir = os.environ['BASEMAPDATA']
    if not os.path.isdir(basemap_datadir):
        raise RuntimeError('Path in environment BASEMAPDATA not a directory')
else:
    from mpl_toolkits import basemap_data
    basemap_datadir = os.path.abspath(list(basemap_data.__path__)[0])

# module variable that sets the default value for the 'latlon' kwarg.
# can be set to True by user so plotting functions can take lons,lats
# in degrees by default, instead of x,y (map projection coords in meters).
latlon_default = False

# supported map projections.
_projnames = {'cyl'      : 'Cylindrical Equidistant',
             'merc'     : 'Mercator',
             'tmerc'    : 'Transverse Mercator',
             'omerc'    : 'Oblique Mercator',
             'mill'     : 'Miller Cylindrical',
             'gall'     : 'Gall Stereographic Cylindrical',
             'cea'      : 'Cylindrical Equal Area',
             'lcc'      : 'Lambert Conformal',
             'laea'     : 'Lambert Azimuthal Equal Area',
             'nplaea'   : 'North-Polar Lambert Azimuthal',
             'splaea'   : 'South-Polar Lambert Azimuthal',
             'eqdc'     : 'Equidistant Conic',
             'aeqd'     : 'Azimuthal Equidistant',
             'npaeqd'   : 'North-Polar Azimuthal Equidistant',
             'spaeqd'   : 'South-Polar Azimuthal Equidistant',
             'aea'      : 'Albers Equal Area',
             'stere'    : 'Stereographic',
             'npstere'  : 'North-Polar Stereographic',
             'spstere'  : 'South-Polar Stereographic',
             'cass'     : 'Cassini-Soldner',
             'poly'     : 'Polyconic',
             'ortho'    : 'Orthographic',
             'geos'     : 'Geostationary',
             'nsper'    : 'Near-Sided Perspective',
             'sinu'     : 'Sinusoidal',
             'moll'     : 'Mollweide',
             'hammer'   : 'Hammer',
             'robin'    : 'Robinson',
             'kav7'     : 'Kavrayskiy VII',
             'eck4'     : 'Eckert IV',
             'vandg'    : 'van der Grinten',
             'mbtfpq'   : 'McBryde-Thomas Flat-Polar Quartic',
             'gnom'     : 'Gnomonic',
             'rotpole'  : 'Rotated Pole',
             }
supported_projections = []
for _items in _projnames.items():
    supported_projections.append(" %-17s%-40s\n" % (_items))
supported_projections = ''.join(supported_projections)

_cylproj = ['cyl','merc','mill','gall','cea']
_pseudocyl = ['moll','robin','eck4','kav7','sinu','mbtfpq','vandg','hammer']
_dg2rad = math.radians(1.)
_rad2dg = math.degrees(1.)

# projection specific parameters.
projection_params = {'cyl'      : 'corners only (no width/height)',
             'merc'     : 'corners plus lat_ts (no width/height)',
             'tmerc'    : 'lon_0,lat_0,k_0',
             'omerc'    : 'lon_0,lat_0,lat_1,lat_2,lon_1,lon_2,no_rot,k_0',
             'mill'     : 'corners only (no width/height)',
             'gall'     : 'corners only (no width/height)',
             'cea'      : 'corners only plus lat_ts (no width/height)',
             'lcc'      : 'lon_0,lat_0,lat_1,lat_2,k_0',
             'laea'     : 'lon_0,lat_0',
             'nplaea'   : 'bounding_lat,lon_0,lat_0,no corners or width/height',
             'splaea'   : 'bounding_lat,lon_0,lat_0,no corners or width/height',
             'eqdc'     : 'lon_0,lat_0,lat_1,lat_2',
             'aeqd'     : 'lon_0,lat_0',
             'npaeqd'   : 'bounding_lat,lon_0,lat_0,no corners or width/height',
             'spaeqd'   : 'bounding_lat,lon_0,lat_0,no corners or width/height',
             'aea'      : 'lon_0,lat_0,lat_1',
             'stere'    : 'lon_0,lat_0,lat_ts,k_0',
             'npstere'  : 'bounding_lat,lon_0,lat_0,no corners or width/height',
             'spstere'  : 'bounding_lat,lon_0,lat_0,no corners or width/height',
             'cass'     : 'lon_0,lat_0',
             'poly'     : 'lon_0,lat_0',
             'ortho'    : 'lon_0,lat_0,llcrnrx,llcrnry,urcrnrx,urcrnry,no width/height',
             'geos'     : 'lon_0,satellite_height,llcrnrx,llcrnry,urcrnrx,urcrnry,no width/height',
             'nsper'    : 'lon_0,satellite_height,llcrnrx,llcrnry,urcrnrx,urcrnry,no width/height',
             'sinu'     : 'lon_0,lat_0,no corners or width/height',
             'moll'     : 'lon_0,lat_0,no corners or width/height',
             'hammer'   : 'lon_0,lat_0,no corners or width/height',
             'robin'    : 'lon_0,lat_0,no corners or width/height',
             'eck4'    : 'lon_0,lat_0,no corners or width/height',
             'kav7'    : 'lon_0,lat_0,no corners or width/height',
             'vandg'    : 'lon_0,lat_0,no corners or width/height',
             'mbtfpq'   : 'lon_0,lat_0,no corners or width/height',
             'gnom'     : 'lon_0,lat_0',
             'rotpole'  : 'lon_0,o_lat_p,o_lon_p,corner lat/lon or corner x,y (no width/height)'
             }

# create dictionary that maps epsg codes to Basemap kwargs.
epsgf = open(os.path.join(basemap_datadir, 'epsg'))
epsg_dict={}
for line in epsgf:
    if line.startswith("#"):
        continue
    l = line.split()
    code = l[0].strip("<>")
    parms = ' '.join(l[1:-1])
    _kw_args={}
    for s in l[1:-1]:
        try:
            k,v = s.split('=')
        except:
            pass
        k = k.strip("+")
        if k=='proj':
            if v == 'longlat': v = 'cyl'
            if v not in _projnames:
                continue
            k='projection'
        if k=='k':
            k='k_0'
        if k in ['projection','lat_1','lat_2','lon_0','lat_0',\
                 'a','b','k_0','lat_ts','ellps','datum']:
            if k not in ['projection','ellps','datum']:
                v = float(v)
            _kw_args[k]=v
    if 'projection' in _kw_args:
        if 'a' in _kw_args:
            if 'b' in _kw_args:
                _kw_args['rsphere']=(_kw_args['a'],_kw_args['b'])
                del _kw_args['b']
            else:
                _kw_args['rsphere']=_kw_args['a']
            del _kw_args['a']
        if 'datum' in _kw_args:
            if _kw_args['datum'] == 'NAD83':
                _kw_args['ellps'] = 'GRS80'
            elif _kw_args['datum'] == 'NAD27':
                _kw_args['ellps'] = 'clrk66'
            elif _kw_args['datum'] == 'WGS84':
                _kw_args['ellps'] = 'WGS84'
            del _kw_args['datum']
        # supported epsg projections.
        # omerc not supported yet, since we can't handle
        # alpha,gamma and lonc keywords.
        if _kw_args['projection'] != 'omerc':
            epsg_dict[code]=_kw_args
epsgf.close()

# The __init__ docstring is pulled out here because it is so long;
# Having it in the usual place makes it hard to get from the
# __init__ argument list to the code that uses the arguments.
_Basemap_init_doc = """

 Sets up a basemap with specified map projection.
 and creates the coastline data structures in map projection
 coordinates.

 Calling a Basemap class instance with the arguments lon, lat will
 convert lon/lat (in degrees) to x/y map projection coordinates
 (in meters). The inverse transformation is done if the optional keyword
 ``inverse`` is set to True.

 The desired projection is set with the projection keyword. Default is ``cyl``.
 Supported values for the projection keyword are:

 ==============   ====================================================
 Value            Description
 ==============   ====================================================
%(supported_projections)s
 ==============   ====================================================

 For most map projections, the map projection region can either be
 specified by setting these keywords:

 .. tabularcolumns:: |l|L|

 ==============   ====================================================
 Keyword          Description
 ==============   ====================================================
 llcrnrlon        longitude of lower left hand corner of the desired map
                  domain (degrees).
 llcrnrlat        latitude of lower left hand corner of the desired map
                  domain (degrees).
 urcrnrlon        longitude of upper right hand corner of the desired map
                  domain (degrees).
 urcrnrlat        latitude of upper right hand corner of the desired map
                  domain (degrees).
 ==============   ====================================================

 or these

 .. tabularcolumns:: |l|L|

 ==============   ====================================================
 Keyword          Description
 ==============   ====================================================
 width            width of desired map domain in projection coordinates
                  (meters).
 height           height of desired map domain in projection coordinates
                  (meters).
 lon_0            center of desired map domain (in degrees).
 lat_0            center of desired map domain (in degrees).
 ==============   ====================================================

 For ``sinu``, ``moll``, ``hammer``, ``npstere``, ``spstere``, ``nplaea``, ``splaea``,
 ``npaeqd``, ``spaeqd``, ``robin``, ``eck4``, ``kav7``, or ``mbtfpq``, the values of
 llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, width and height are ignored
 (because either they are computed internally, or entire globe is
 always plotted).

 For the cylindrical projections (``cyl``, ``merc``, ``mill``, ``cea``  and ``gall``),
 the default is to use
 llcrnrlon=-180,llcrnrlat=-90, urcrnrlon=180 and urcrnrlat=90). For all other
 projections except ``ortho``, ``geos`` and ``nsper``, either the lat/lon values of the
 corners or width and height must be specified by the user.

 For ``ortho``, ``geos`` and ``nsper``, the lat/lon values of the corners may be specified,
 or the x/y values of the corners (llcrnrx,llcrnry,urcrnrx,urcrnry) in the
 coordinate system of the global projection (with x=0,y=0 at the center
 of the global projection).  If the corners are not specified,
 the entire globe is plotted.

 For ``rotpole``, the lat/lon values of the corners on the unrotated sphere
 may be provided as llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat, or the lat/lon
 values of the corners on the rotated sphere can be given as
 llcrnrx,llcrnry,urcrnrx,urcrnry.

 Other keyword arguments:

 .. tabularcolumns:: |l|L|

 ==============   ====================================================
 Keyword          Description
 ==============   ====================================================
 resolution       resolution of boundary database to use. Can be ``c``
                  (crude), ``l`` (low), ``i`` (intermediate), ``h``
                  (high), ``f`` (full) or None.
                  If None, no boundary data will be read in (and
                  class methods such as drawcoastlines will raise an
                  if invoked).
                  Resolution drops off by roughly 80%% between datasets.
                  Higher res datasets are much slower to draw.
                  Default ``c``. Coastline data is from the GSHHS
                  (http://www.soest.hawaii.edu/wessel/gshhs/gshhs.html).
                  State, country and river datasets from the Generic
                  Mapping Tools (http://gmt.soest.hawaii.edu).
 area_thresh      coastline or lake with an area smaller than
                  area_thresh in km^2 will not be plotted.
                  Default 10000,1000,100,10,1 for resolution
                  ``c``, ``l``, ``i``, ``h``, ``f``.
 rsphere          radius of the sphere used to define map projection
                  (default 6370997 meters, close to the arithmetic mean
                  radius of the earth). If given as a sequence, the
                  first two elements are interpreted as the radii
                  of the major and minor axes of an ellipsoid.
                  Note: sometimes an ellipsoid is specified by the
                  major axis and an inverse flattening parameter (if).
                  The minor axis (b) can be computed from the major
                  axis (a) and the inverse flattening parameter using
                  the formula if = a/(a-b).
 ellps            string describing ellipsoid ('GRS80' or 'WGS84',
                  for example). If both rsphere and ellps are given,
                  rsphere is ignored. Default None. See pyproj.pj_ellps
                  for allowed values.
 suppress_ticks   suppress automatic drawing of axis ticks and labels
                  in map projection coordinates.  Default True,
                  so parallels and meridians can be labelled instead.
                  If parallel or meridian labelling is requested
                  (using drawparallels and drawmeridians methods),
                  automatic tick labelling will be supressed even if
                  suppress_ticks=False.  suppress_ticks=False
                  is useful if you want to use your own custom tick
                  formatter, or  if you want to let matplotlib label
                  the axes in meters using map projection
                  coordinates.
 fix_aspect       fix aspect ratio of plot to match aspect ratio
                  of map projection region (default True).
 anchor           determines how map is placed in axes rectangle
                  (passed to axes.set_aspect). Default is ``C``,
                  which means map is centered.
                  Allowed values are
                  ``C``, ``SW``, ``S``, ``SE``, ``E``, ``NE``,
                  ``N``, ``NW``, and ``W``.
 celestial        use astronomical conventions for longitude (i.e.
                  negative longitudes to the east of 0). Default False.
                  Implies resolution=None.
 ax               set default axes instance
                  (default None - matplotlib.pyplot.gca() may be used
                  to get the current axes instance).
                  If you do not want matplotlib.pyplot to be imported,
                  you can either set this to a pre-defined axes
                  instance, or use the ``ax`` keyword in each Basemap
                  method call that does drawing. In the first case,
                  all Basemap method calls will draw to the same axes
                  instance.  In the second case, you can draw to
                  different axes with the same Basemap instance.
                  You can also use the ``ax`` keyword in individual
                  method calls to selectively override the default
                  axes instance.
 ==============   ====================================================

 The following keywords are map projection parameters which all default to
 None.  Not all parameters are used by all projections, some are ignored.
 The module variable ``projection_params`` is a dictionary which
 lists which parameters apply to which projections.

 .. tabularcolumns:: |l|L|

 ================ ====================================================
 Keyword          Description
 ================ ====================================================
 lat_ts           latitude of true scale. Optional for stereographic,
                  cylindrical equal area and mercator projections.
                  default is lat_0 for stereographic projection.
                  default is 0 for mercator and cylindrical equal area
                  projections.
 lat_1            first standard parallel for lambert conformal,
                  albers equal area and equidistant conic.
                  Latitude of one of the two points on the projection
                  centerline for oblique mercator. If lat_1 is not given, but
                  lat_0 is, lat_1 is set to lat_0 for lambert
                  conformal, albers equal area and equidistant conic.
 lat_2            second standard parallel for lambert conformal,
                  albers equal area and equidistant conic.
                  Latitude of one of the two points on the projection
                  centerline for oblique mercator. If lat_2 is not
                  given it is set to lat_1 for lambert conformal,
                  albers equal area and equidistant conic.
 lon_1            Longitude of one of the two points on the projection
                  centerline for oblique mercator.
 lon_2            Longitude of one of the two points on the projection
                  centerline for oblique mercator.
 k_0              Scale factor at natural origin (used
                  by 'tmerc', 'omerc', 'stere' and 'lcc').
 no_rot           only used by oblique mercator.
                  If set to True, the map projection coordinates will
                  not be rotated to true North.  Default is False
                  (projection coordinates are automatically rotated).
 lat_0            central latitude (y-axis origin) - used by all
                  projections.
 lon_0            central meridian (x-axis origin) - used by all
                  projections.
 o_lat_p          latitude of rotated pole (only used by 'rotpole')
 o_lon_p          longitude of rotated pole (only used by 'rotpole')
 boundinglat      bounding latitude for pole-centered projections
                  (npstere,spstere,nplaea,splaea,npaeqd,spaeqd).
                  These projections are square regions centered
                  on the north or south pole.
                  The longitude lon_0 is at 6-o'clock, and the
                  latitude circle boundinglat is tangent to the edge
                  of the map at lon_0.
 round            cut off pole-centered projection at boundinglat
                  (so plot is a circle instead of a square). Only
                  relevant for npstere,spstere,nplaea,splaea,npaeqd
                  or spaeqd projections. Default False.
 satellite_height height of satellite (in m) above equator -
                  only relevant for geostationary
                  and near-sided perspective (``geos`` or ``nsper``)
                  projections. Default 35,786 km.
 ================ ====================================================

 Useful instance variables:

 .. tabularcolumns:: |l|L|

 ================ ====================================================
 Variable Name    Description
 ================ ====================================================
 projection       map projection. Print the module variable
                  ``supported_projections`` to see a list of allowed
                  values.
 epsg             EPSG code defining projection (see
                  http://spatialreference.org for a list of
                  EPSG codes and their definitions).
 aspect           map aspect ratio
                  (size of y dimension / size of x dimension).
 llcrnrlon        longitude of lower left hand corner of the
                  selected map domain.
 llcrnrlat        latitude of lower left hand corner of the
                  selected map domain.
 urcrnrlon        longitude of upper right hand corner of the
                  selected map domain.
 urcrnrlat        latitude of upper right hand corner of the
                  selected map domain.
 llcrnrx          x value of lower left hand corner of the
                  selected map domain in map projection coordinates.
 llcrnry          y value of lower left hand corner of the
                  selected map domain in map projection coordinates.
 urcrnrx          x value of upper right hand corner of the
                  selected map domain in map projection coordinates.
 urcrnry          y value of upper right hand corner of the
                  selected map domain in map projection coordinates.
 rmajor           equatorial radius of ellipsoid used (in meters).
 rminor           polar radius of ellipsoid used (in meters).
 resolution       resolution of boundary dataset being used (``c``
                  for crude, ``l`` for low, etc.).
                  If None, no boundary dataset is associated with the
                  Basemap instance.
 proj4string      the string describing the map projection that is
                  used by PROJ.4.
 ================ ====================================================

 **Converting from Geographic (lon/lat) to Map Projection (x/y) Coordinates**

 Calling a Basemap class instance with the arguments lon, lat will
 convert lon/lat (in degrees) to x/y map projection
 coordinates (in meters).  If optional keyword ``inverse`` is
 True (default is False), the inverse transformation from x/y
 to lon/lat is performed.

 For cylindrical equidistant projection (``cyl``), this
 does nothing (i.e. x,y == lon,lat).

 For non-cylindrical projections, the inverse transformation
 always returns longitudes between -180 and 180 degrees. For
 cylindrical projections (self.projection == ``cyl``, ``mill``,
 ``cea``, ``gall`` or ``merc``)
 the inverse transformation will return longitudes between
 self.llcrnrlon and self.llcrnrlat.

 Input arguments lon, lat can be either scalar floats, sequences
 or numpy arrays.

 **Example Usage:**

 >>> from mpl_toolkits.basemap import Basemap
 >>> import numpy as np
 >>> import matplotlib.pyplot as plt
 >>> # read in topo data (on a regular lat/lon grid)
 >>> etopo = np.loadtxt('etopo20data.gz')
 >>> lons  = np.loadtxt('etopo20lons.gz')
 >>> lats  = np.loadtxt('etopo20lats.gz')
 >>> # create Basemap instance for Robinson projection.
 >>> m = Basemap(projection='robin',lon_0=0.5*(lons[0]+lons[-1]))
 >>> # compute map projection coordinates for lat/lon grid.
 >>> x, y = m(*np.meshgrid(lons,lats))
 >>> # make filled contour plot.
 >>> cs = m.contourf(x,y,etopo,30,cmap=plt.cm.jet)
 >>> m.drawcoastlines() # draw coastlines
 >>> m.drawmapboundary() # draw a line around the map region
 >>> m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0]) # draw parallels
 >>> m.drawmeridians(np.arange(0.,420.,60.),labels=[0,0,0,1]) # draw meridians
 >>> plt.title('Robinson Projection') # add a title
 >>> plt.show()

 [this example (simpletest.py) plus many others can be found in the
 examples directory of source distribution.  The "OO" version of this
 example (which does not use matplotlib.pyplot) is called "simpletest_oo.py".]
""" % locals()

# unsupported projection error message.
_unsupported_projection = ["'%s' is an unsupported projection.\n"]
_unsupported_projection.append("The supported projections are:\n")
_unsupported_projection.append(supported_projections)
_unsupported_projection = ''.join(_unsupported_projection)

def _validated_ll(param, name, minval, maxval):
    param = float(param)
    if param > maxval or param < minval:
        raise ValueError('%s must be between %f and %f degrees' %
                                           (name, minval, maxval))
    return param


def _validated_or_none(param, name, minval, maxval):
    if param is None:
        return None
    return _validated_ll(param, name, minval, maxval)


def _insert_validated(d, param, name, minval, maxval):
    if param is not None:
        d[name] = _validated_ll(param, name, minval, maxval)

def _transform(plotfunc):
    # shift data and longitudes to map projection region, then compute
    # transformation to map projection coordinates.
    @functools.wraps(plotfunc)
    def with_transform(self,x,y,data,*args,**kwargs):
        # input coordinates are latitude/longitude, not map projection coords.
        if kwargs.pop('latlon', latlon_default):
            # shift data to map projection region for
            # cylindrical and pseudo-cylindrical projections.
            if self.projection in _cylproj or self.projection in _pseudocyl:
                x, data = self.shiftdata(x, data,
                                         fix_wrap_around=plotfunc.__name__ not in ["scatter"])
            # convert lat/lon coords to map projection coords.
            x, y = self(x,y)
        return plotfunc(self,x,y,data,*args,**kwargs)
    return with_transform

def _transform1d(plotfunc):
    # shift data and longitudes to map projection region, then compute
    # transformation to map projection coordinates.
    @functools.wraps(plotfunc)
    def with_transform(self,x,y,*args,**kwargs):
        x = np.asarray(x)
        # input coordinates are latitude/longitude, not map projection coords.
        if kwargs.pop('latlon', latlon_default):
            # shift data to map projection region for
            # cylindrical and pseudo-cylindrical projections.
            if self.projection in _cylproj or self.projection in _pseudocyl:
                if x.ndim == 1:
                    x = self.shiftdata(x, fix_wrap_around=plotfunc.__name__ not in ["scatter"])
                elif x.ndim == 0:
                    if x > 180:
                        x = x - 360.
            # convert lat/lon coords to map projection coords.
            x, y = self(x,y)
        return plotfunc(self,x,y,*args,**kwargs)
    return with_transform

def _transformuv(plotfunc):
    # shift data and longitudes to map projection region, then compute
    # transformation to map projection coordinates. Works when call
    # signature has two data arrays instead of one.
    @functools.wraps(plotfunc)
    def with_transform(self,x,y,u,v,*args,**kwargs):
        # input coordinates are latitude/longitude, not map projection coords.
        if kwargs.pop('latlon', latlon_default):
            # shift data to map projection region for
            # cylindrical and pseudo-cylindrical projections.
            if self.projection in _cylproj or self.projection in _pseudocyl:
                x1, u = self.shiftdata(x, u)
                x, v = self.shiftdata(x, v)
            # convert lat/lon coords to map projection coords.
            x, y = self(x,y)
        return plotfunc(self,x,y,u,v,*args,**kwargs)
    return with_transform

class Basemap(object):

    def __init__(self, llcrnrlon=None, llcrnrlat=None,
                       urcrnrlon=None, urcrnrlat=None,
                       llcrnrx=None, llcrnry=None,
                       urcrnrx=None, urcrnry=None,
                       width=None, height=None,
                       projection='cyl', resolution='c',
                       area_thresh=None, rsphere=6370997.0,
                       ellps=None, lat_ts=None,
                       lat_1=None, lat_2=None,
                       lat_0=None, lon_0=None,
                       lon_1=None, lon_2=None,
                       o_lon_p=None, o_lat_p=None,
                       k_0=None,
                       no_rot=False,
                       suppress_ticks=True,
                       satellite_height=35786000,
                       boundinglat=None,
                       fix_aspect=True,
                       anchor='C',
                       celestial=False,
                       round=False,
                       epsg=None,
                       ax=None):
        # docstring is added after __init__ method definition

        # set epsg code if given, set to 4326 for projection='cyl':
        if epsg is not None:
            self.epsg = epsg
        elif projection == 'cyl':
            self.epsg = 4326
        # replace kwarg values with those implied by epsg code,
        # if given.
        if hasattr(self,'epsg'):
            if str(self.epsg) not in epsg_dict:
                raise ValueError('%s is not a supported EPSG code' %
                        self.epsg)
            epsg_params = epsg_dict[str(self.epsg)]
            for k in epsg_params:
                if k == 'projection':
                    projection = epsg_params[k]
                elif k == 'rsphere':
                    rsphere = epsg_params[k]
                elif k == 'ellps':
                    ellps = epsg_params[k]
                elif k == 'lat_1':
                    lat_1 = epsg_params[k]
                elif k == 'lat_2':
                    lat_2 = epsg_params[k]
                elif k == 'lon_0':
                    lon_0 = epsg_params[k]
                elif k == 'lat_0':
                    lat_0 = epsg_params[k]
                elif k == 'lat_ts':
                    lat_ts = epsg_params[k]
                elif k == 'k_0':
                    k_0 = epsg_params[k]

        # fix aspect to ratio to match aspect ratio of map projection
        # region
        self.fix_aspect = fix_aspect
        # where to put plot in figure (default is 'C' or center)
        self.anchor = anchor
        # geographic or celestial coords?
        self.celestial = celestial
        # map projection.
        self.projection = projection
        # bounding lat (for pole-centered plots)
        self.boundinglat = boundinglat
        # is a round pole-centered plot desired?
        self.round = round
        # full disk projection?
        self._fulldisk = False # default value

        # set up projection parameter dict.
        projparams = {}
        projparams['proj'] = projection
        # if ellps keyword specified, it over-rides rsphere.
        if ellps is not None:
            try:
                elldict = pyproj.pj_ellps[ellps]
            except KeyError:
                raise ValueError(
                'illegal ellps definition, allowed values are %s' %
                pyproj.pj_ellps.keys())
            projparams['a'] = elldict['a']
            if 'b' in elldict:
                projparams['b'] = elldict['b']
            else:
                projparams['b'] = projparams['a']*(1.0-(1.0/elldict['rf']))
        else:
            try:
                if rsphere[0] > rsphere[1]:
                    projparams['a'] = rsphere[0]
                    projparams['b'] = rsphere[1]
                else:
                    projparams['a'] = rsphere[1]
                    projparams['b'] = rsphere[0]
            except:
                if projection == 'tmerc':
                # use bR_a instead of R because of obscure bug
                # in proj4 for tmerc projection.
                    projparams['bR_a'] = rsphere
                else:
                    projparams['R'] = rsphere
        # set units to meters.
        projparams['units']='m'
        # check for sane values of lon_0, lat_0, lat_ts, lat_1, lat_2
        lat_0 = _validated_or_none(lat_0, 'lat_0', -90, 90)
        lat_1 = _validated_or_none(lat_1, 'lat_1', -90, 90)
        lat_2 = _validated_or_none(lat_2, 'lat_2', -90, 90)
        lat_ts = _validated_or_none(lat_ts, 'lat_ts', -90, 90)
        lon_0 = _validated_or_none(lon_0, 'lon_0', -360, 720)
        lon_1 = _validated_or_none(lon_1, 'lon_1', -360, 720)
        lon_2 = _validated_or_none(lon_2, 'lon_2', -360, 720)
        llcrnrlon = _validated_or_none(llcrnrlon, 'llcrnrlon', -360, 720)
        urcrnrlon = _validated_or_none(urcrnrlon, 'urcrnrlon', -360, 720)
        llcrnrlat = _validated_or_none(llcrnrlat, 'llcrnrlat', -90, 90)
        urcrnrlat = _validated_or_none(urcrnrlat, 'urcrnrlat', -90, 90)

        _insert_validated(projparams, lat_0, 'lat_0', -90, 90)
        _insert_validated(projparams, lat_1, 'lat_1', -90, 90)
        _insert_validated(projparams, lat_2, 'lat_2', -90, 90)
        _insert_validated(projparams, lat_ts, 'lat_ts', -90, 90)
        _insert_validated(projparams, lon_0, 'lon_0', -360, 720)
        _insert_validated(projparams, lon_1, 'lon_1', -360, 720)
        _insert_validated(projparams, lon_2, 'lon_2', -360, 720)
        if projection in ['geos','nsper']:
            projparams['h'] = satellite_height
        # check for sane values of projection corners.
        using_corners = (None not in [llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat])
        if using_corners:
            self.llcrnrlon = _validated_ll(llcrnrlon, 'llcrnrlon', -360, 720)
            self.urcrnrlon = _validated_ll(urcrnrlon, 'urcrnrlon', -360, 720)
            self.llcrnrlat = _validated_ll(llcrnrlat, 'llcrnrlat', -90, 90)
            self.urcrnrlat = _validated_ll(urcrnrlat, 'urcrnrlat', -90, 90)

        # for each of the supported projections,
        # compute lat/lon of domain corners
        # and set values in projparams dict as needed.

        if projection in ['lcc', 'eqdc', 'aea']:
            if projection == 'lcc' and k_0 is not None:
                projparams['k_0']=k_0
            # if lat_0 is given, but not lat_1,
            # set lat_1=lat_0
            if lat_1 is None and lat_0 is not None:
                lat_1 = lat_0
                projparams['lat_1'] = lat_1
            if lat_1 is None or lon_0 is None:
                raise ValueError('must specify lat_1 or lat_0 and lon_0 for %s basemap (lat_2 is optional)' % _projnames[projection])
            if lat_2 is None:
                projparams['lat_2'] = lat_1
            if not using_corners:
                using_cornersxy = (None not in [llcrnrx,llcrnry,urcrnrx,urcrnry])
                if using_cornersxy:
                    llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat = _choosecornersllur(llcrnrx,llcrnry,urcrnrx,urcrnry,**projparams)
                    self.llcrnrlon = llcrnrlon; self.llcrnrlat = llcrnrlat
                    self.urcrnrlon = urcrnrlon; self.urcrnrlat = urcrnrlat
                else:
                    if width is None or height is None:
                        raise ValueError('must either specify lat/lon values of corners (llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat) in degrees or width and height in meters')
                    if lon_0 is None or lat_0 is None:
                        raise ValueError('must specify lon_0 and lat_0 when using width, height to specify projection region')
                    llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat = _choosecorners(width,height,**projparams)
                    self.llcrnrlon = llcrnrlon; self.llcrnrlat = llcrnrlat
                    self.urcrnrlon = urcrnrlon; self.urcrnrlat = urcrnrlat
        elif projection == 'stere':
            if k_0 is not None:
                projparams['k_0']=k_0
            if lat_0 is None or lon_0 is None:
                raise ValueError('must specify lat_0 and lon_0 for Stereographic basemap (lat_ts is optional)')
            if not using_corners:
                if width is None or height is None:
                    raise ValueError('must either specify lat/lon values of corners (llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat) in degrees or width and height in meters')
                if lon_0 is None or lat_0 is None:
                    raise ValueError('must specify lon_0 and lat_0 when using width, height to specify projection region')
                llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat = _choosecorners(width,height,**projparams)
                self.llcrnrlon = llcrnrlon; self.llcrnrlat = llcrnrlat
                self.urcrnrlon = urcrnrlon; self.urcrnrlat = urcrnrlat
        elif projection in ['spstere', 'npstere',
                            'splaea', 'nplaea',
                            'spaeqd', 'npaeqd']:
            if (projection == 'splaea' and boundinglat >= 0) or\
               (projection == 'nplaea' and boundinglat <= 0):
                msg='boundinglat cannot extend into opposite hemisphere'
                raise ValueError(msg)
            if boundinglat is None or lon_0 is None:
                raise ValueError('must specify boundinglat and lon_0 for %s basemap' % _projnames[projection])
            if projection[0] == 's':
                sgn = -1
            else:
                sgn = 1
            rootproj = projection[2:]
            projparams['proj'] = rootproj
            if rootproj == 'stere':
                projparams['lat_ts'] = sgn * 90.
            projparams['lat_0'] = sgn * 90.
            self.llcrnrlon = lon_0 - sgn*45.
            self.urcrnrlon = lon_0 + sgn*135.
            proj = pyproj.Proj(projparams)
            x,y = proj(lon_0,boundinglat)
            lon,self.llcrnrlat = proj(math.sqrt(2.)*y,0.,inverse=True)
            self.urcrnrlat = self.llcrnrlat
            if width is not None or height is not None:
                sys.stdout.write('warning: width and height keywords ignored for %s projection' % _projnames[projection])
        elif projection == 'laea':
            if lat_0 is None or lon_0 is None:
                raise ValueError('must specify lat_0 and lon_0 for Lambert Azimuthal basemap')
            if not using_corners:
                if width is None or height is None:
                    raise ValueError('must either specify lat/lon values of corners (llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat) in degrees or width and height in meters')
                if lon_0 is None or lat_0 is None:
                    raise ValueError('must specify lon_0 and lat_0 when using width, height to specify projection region')
                llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat = _choosecorners(width,height,**projparams)
                self.llcrnrlon = llcrnrlon; self.llcrnrlat = llcrnrlat
                self.urcrnrlon = urcrnrlon; self.urcrnrlat = urcrnrlat
        elif projection in ['tmerc','gnom','cass','poly'] :
            if projection == 'tmerc' and k_0 is not None:
                projparams['k_0']=k_0
            if projection == 'gnom' and 'R' not in projparams:
                raise ValueError('gnomonic projection only works for perfect spheres - not ellipsoids')
            if lat_0 is None or lon_0 is None:
                raise ValueError('must specify lat_0 and lon_0 for Transverse Mercator, Gnomonic, Cassini-Soldnerr and Polyconic basemap')
            if not using_corners:
                if width is None or height is None:
                    raise ValueError('must either specify lat/lon values of corners (llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat) in degrees or width and height in meters')
                if lon_0 is None or lat_0 is None:
                    raise ValueError('must specify lon_0 and lat_0 when using width, height to specify projection region')
                llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat = _choosecorners(width,height,**projparams)
                self.llcrnrlon = llcrnrlon; self.llcrnrlat = llcrnrlat
                self.urcrnrlon = urcrnrlon; self.urcrnrlat = urcrnrlat
        elif projection == 'ortho':
            if 'R' not in projparams:
                raise ValueError('orthographic projection only works for perfect spheres - not ellipsoids')
            if lat_0 is None or lon_0 is None:
                raise ValueError('must specify lat_0 and lon_0 for Orthographic basemap')
            if (lat_0 == 90 or lat_0 == -90) and\
               None in [llcrnrx,llcrnry,urcrnrx,urcrnry]:
                # for ortho plot centered on pole, set boundinglat to equator.
                # (so meridian labels can be drawn in this special case).
                self.boundinglat = 0
                self.round = True
            if width is not None or height is not None:
                sys.stdout.write('warning: width and height keywords ignored for %s projection' % _projnames[self.projection])
            if not using_corners:
                llcrnrlon = -180.
                llcrnrlat = -90.
                urcrnrlon = 180
                urcrnrlat = 90.
                self._fulldisk = True
            else:
                self._fulldisk = False
            self.llcrnrlon = llcrnrlon; self.llcrnrlat = llcrnrlat
            self.urcrnrlon = urcrnrlon; self.urcrnrlat = urcrnrlat
            # FIXME: won't work for points exactly on equator??
            if np.abs(lat_0) < 1.e-2: lat_0 = 1.e-2
            projparams['lat_0'] = lat_0
        elif projection == 'geos':
            if lat_0 is not None and lat_0 != 0:
                raise ValueError('lat_0 must be zero for Geostationary basemap')
            if lon_0 is None:
                raise ValueError('must specify lon_0 for Geostationary basemap')
            if width is not None or height is not None:
                sys.stdout.write('warning: width and height keywords ignored for %s projection' % _projnames[self.projection])
            if not using_corners:
                llcrnrlon = -180.
                llcrnrlat = -90.
                urcrnrlon = 180
                urcrnrlat = 90.
                self._fulldisk = True
            else:
                self._fulldisk = False
            self.llcrnrlon = llcrnrlon; self.llcrnrlat = llcrnrlat
            self.urcrnrlon = urcrnrlon; self.urcrnrlat = urcrnrlat
        elif projection == 'nsper':
            if 'R' not in projparams:
                raise ValueError('near-sided perspective projection only works for perfect spheres - not ellipsoids')
            if lat_0 is None or lon_0 is None:
                msg='must specify lon_0 and lat_0 for near-sided perspective Basemap'
                raise ValueError(msg)
            if width is not None or height is not None:
                sys.stdout.write('warning: width and height keywords ignored for %s projection' % _projnames[self.projection])
            if not using_corners:
                llcrnrlon = -180.
                llcrnrlat = -90.
                urcrnrlon = 180
                urcrnrlat = 90.
                self._fulldisk = True
            else:
                self._fulldisk = False
            self.llcrnrlon = llcrnrlon; self.llcrnrlat = llcrnrlat
            self.urcrnrlon = urcrnrlon; self.urcrnrlat = urcrnrlat
        elif projection in _pseudocyl:
            if lon_0 is None:
                raise ValueError('must specify lon_0 for %s projection' % _projnames[self.projection])
            if width is not None or height is not None:
                sys.stdout.write('warning: width and height keywords ignored for %s projection' % _projnames[self.projection])
            llcrnrlon = lon_0-180.
            llcrnrlat = -90.
            urcrnrlon = lon_0+180
            urcrnrlat = 90.
            self.llcrnrlon = llcrnrlon; self.llcrnrlat = llcrnrlat
            self.urcrnrlon = urcrnrlon; self.urcrnrlat = urcrnrlat
        elif projection == 'omerc':
            if k_0 is not None:
                projparams['k_0']=k_0
            if lat_1 is None or lon_1 is None or lat_2 is None or lon_2 is None:
                raise ValueError('must specify lat_1,lon_1 and lat_2,lon_2 for Oblique Mercator basemap')
            projparams['lat_1'] = lat_1
            projparams['lon_1'] = lon_1
            projparams['lat_2'] = lat_2
            projparams['lon_2'] = lon_2
            projparams['lat_0'] = lat_0
            if no_rot:
                projparams['no_rot']=''
            #if not using_corners:
            #    raise ValueError, 'cannot specify map region with width and height keywords for this projection, please specify lat/lon values of corners'
            if not using_corners:
                if width is None or height is None:
                    raise ValueError('must either specify lat/lon values of corners (llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat) in degrees or width and height in meters')
                if lon_0 is None or lat_0 is None:
                    raise ValueError('must specify lon_0 and lat_0 when using width, height to specify projection region')
                llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat = _choosecorners(width,height,**projparams)
                self.llcrnrlon = llcrnrlon; self.llcrnrlat = llcrnrlat
                self.urcrnrlon = urcrnrlon; self.urcrnrlat = urcrnrlat
        elif projection == 'aeqd':
            if lat_0 is None or lon_0 is None:
                raise ValueError('must specify lat_0 and lon_0 for Azimuthal Equidistant basemap')
            if not using_corners:
                if width is None or height is None:
                    self._fulldisk = True
                    llcrnrlon = -180.
                    llcrnrlat = -90.
                    urcrnrlon = 180
                    urcrnrlat = 90.
                else:
                    self._fulldisk = False
                if lon_0 is None or lat_0 is None:
                    raise ValueError('must specify lon_0 and lat_0 when using width, height to specify projection region')
                if not self._fulldisk:
                    llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat = _choosecorners(width,height,**projparams)
                self.llcrnrlon = llcrnrlon; self.llcrnrlat = llcrnrlat
                self.urcrnrlon = urcrnrlon; self.urcrnrlat = urcrnrlat
        elif projection in _cylproj:
            if projection == 'merc' or projection == 'cea':
                if lat_ts is None:
                    lat_ts = 0.
                    projparams['lat_ts'] = lat_ts
            if not using_corners:
                llcrnrlat = -90.
                urcrnrlat = 90.
                if lon_0 is not None:
                    llcrnrlon = lon_0-180.
                    urcrnrlon = lon_0+180.
                else:
                    llcrnrlon = -180.
                    urcrnrlon = 180
                if projection == 'merc':
                    # clip plot region to be within -89.99S to 89.99N
                    # (mercator is singular at poles)
                    if llcrnrlat < -89.99: llcrnrlat = -89.99
                    if llcrnrlat > 89.99: llcrnrlat = 89.99
                    if urcrnrlat < -89.99: urcrnrlat = -89.99
                    if urcrnrlat > 89.99: urcrnrlat = 89.99
                self.llcrnrlon = llcrnrlon; self.llcrnrlat = llcrnrlat
                self.urcrnrlon = urcrnrlon; self.urcrnrlat = urcrnrlat
            if width is not None or height is not None:
                sys.stdout.write('warning: width and height keywords ignored for %s projection' % _projnames[self.projection])
            if lon_0 is not None:
                projparams['lon_0'] = lon_0
            else:
                projparams['lon_0']=0.5*(llcrnrlon+urcrnrlon)
        elif projection == 'rotpole':
            if lon_0 is None or o_lon_p is None or o_lat_p is None:
                msg='must specify lon_0,o_lat_p,o_lon_p for rotated pole Basemap'
                raise ValueError(msg)
            if width is not None or height is not None:
                sys.stdout.write('warning: width and height keywords ignored for %s projection' % _projnames[self.projection])
            projparams['lon_0']=lon_0
            projparams['o_lon_p']=o_lon_p
            projparams['o_lat_p']=o_lat_p
            projparams['o_proj']='longlat'
            projparams['proj']='ob_tran'
            if not using_corners and None in [llcrnrx,llcrnry,urcrnrx,urcrnry]:
                raise ValueError('must specify lat/lon values of corners in degrees')
            if None not in [llcrnrx,llcrnry,urcrnrx,urcrnry]:
                p = pyproj.Proj(projparams)
                llcrnrx = _dg2rad*llcrnrx; llcrnry = _dg2rad*llcrnry
                urcrnrx = _dg2rad*urcrnrx; urcrnry = _dg2rad*urcrnry
                llcrnrlon, llcrnrlat = p(llcrnrx,llcrnry,inverse=True)
                urcrnrlon, urcrnrlat = p(urcrnrx,urcrnry,inverse=True)
                self.llcrnrlon = llcrnrlon; self.llcrnrlat = llcrnrlat
                self.urcrnrlon = urcrnrlon; self.urcrnrlat = urcrnrlat
        else:
            raise ValueError(_unsupported_projection % projection)

        # initialize proj4
        proj = Proj(projparams,self.llcrnrlon,self.llcrnrlat,self.urcrnrlon,self.urcrnrlat)

        # make sure axis ticks are suppressed.
        self.noticks = suppress_ticks
        # map boundary not yet drawn.
        self._mapboundarydrawn = False

        # make Proj instance a Basemap instance variable.
        self.projtran = proj
        # copy some Proj attributes.
        atts = ['rmajor','rminor','esq','flattening','ellipsoid','projparams']
        for att in atts:
            self.__dict__[att] = proj.__dict__[att]
        # these only exist for geostationary projection.
        if hasattr(proj,'_width'):
            self.__dict__['_width'] = proj.__dict__['_width']
        if hasattr(proj,'_height'):
            self.__dict__['_height'] = proj.__dict__['_height']
        # spatial reference string (useful for georeferencing output
        # images with gdal_translate).
        if hasattr(self,'_proj4'):
            #self.srs = proj._proj4.srs
            self.srs = proj._proj4.pjinitstring
        else:
            pjargs = []
            for key,value in self.projparams.items():
                # 'cyl' projection translates to 'eqc' in PROJ.4
                if projection == 'cyl' and key == 'proj':
                    value = 'eqc'
                # ignore x_0 and y_0 settings for 'cyl' projection
                # (they are not consistent with what PROJ.4 uses)
                elif projection == 'cyl' and key in ['x_0','y_0']:
                    continue
                pjargs.append('+'+key+"="+str(value)+' ')
            self.srs = ''.join(pjargs)
        self.proj4string = self.srs
        # set instance variables defining map region.
        self.xmin = proj.xmin
        self.xmax = proj.xmax
        self.ymin = proj.ymin
        self.ymax = proj.ymax
        if projection == 'cyl':
            self.aspect = (self.urcrnrlat-self.llcrnrlat)/(self.urcrnrlon-self.llcrnrlon)
        else:
            self.aspect = (proj.ymax-proj.ymin)/(proj.xmax-proj.xmin)
        if projection in ['geos','ortho','nsper'] and \
           None not in [llcrnrx,llcrnry,urcrnrx,urcrnry]:
            self.llcrnrx = llcrnrx+0.5*proj.xmax
            self.llcrnry = llcrnry+0.5*proj.ymax
            self.urcrnrx = urcrnrx+0.5*proj.xmax
            self.urcrnry = urcrnry+0.5*proj.ymax
            self._fulldisk = False
        else:
            self.llcrnrx = proj.llcrnrx
            self.llcrnry = proj.llcrnry
            self.urcrnrx = proj.urcrnrx
            self.urcrnry = proj.urcrnry

        if self.projection == 'rotpole':
            lon0,lat0 = self(0.5*(self.llcrnrx + self.urcrnrx),\
                             0.5*(self.llcrnry + self.urcrnry),\
                             inverse=True)
            self.projparams['lat_0']=lat0

        # if ax == None, pyplot.gca may be used.
        self.ax = ax
        self.lsmask = None
        # This will record hashs of Axes instances.
        self._initialized_axes = set()

        # set defaults for area_thresh.
        self.resolution = resolution
        # celestial=True implies resolution=None (no coastlines).
        if self.celestial:
            self.resolution=None
        if area_thresh is None and self.resolution is not None:
            if resolution == 'c':
                area_thresh = 10000.
            elif resolution == 'l':
                area_thresh = 1000.
            elif resolution == 'i':
                area_thresh = 100.
            elif resolution == 'h':
                area_thresh = 10.
            elif resolution == 'f':
                area_thresh = 1.
            else:
                raise ValueError("boundary resolution must be one of 'c','l','i','h' or 'f'")
        self.area_thresh = area_thresh
        # define map boundary polygon (in lat/lon coordinates)
        blons, blats, self._boundarypolyll, self._boundarypolyxy = self._getmapboundary()
        self.boundarylats = blats
        self.boundarylons = blons
        # set min/max lats for projection domain.
        if self.projection in _cylproj:
            self.latmin = self.llcrnrlat
            self.latmax = self.urcrnrlat
            self.lonmin = self.llcrnrlon
            self.lonmax = self.urcrnrlon
        elif self.projection in ['ortho','geos','nsper'] + _pseudocyl:
            self.latmin = -90.
            self.latmax = 90.
            self.lonmin = self.llcrnrlon
            self.lonmax = self.urcrnrlon
        else:
            lons, lats = self.makegrid(1001,1001)
            lats = ma.masked_where(lats > 1.e20,lats)
            lons = ma.masked_where(lons > 1.e20,lons)
            self.latmin = lats.min()
            self.latmax = lats.max()
            self.lonmin = lons.min()
            self.lonmax = lons.max()
            NPole = _geoslib.Point(self(0.,90.))
            SPole = _geoslib.Point(self(0.,-90.))
            if lat_0 is None:
                lon_0, lat_0 =\
                self(0.5*(self.xmin+self.xmax),
                     0.5*(self.ymin+self.ymax),inverse=True)
            Dateline = _geoslib.Point(self(180.,lat_0))
            Greenwich = _geoslib.Point(self(0.,lat_0))
            hasNP = NPole.within(self._boundarypolyxy)
            hasSP = SPole.within(self._boundarypolyxy)
            hasPole = hasNP or hasSP
            hasDateline = Dateline.within(self._boundarypolyxy)
            hasGreenwich = Greenwich.within(self._boundarypolyxy)
            # projection crosses dateline (and not Greenwich or pole).
            if not hasPole and hasDateline and not hasGreenwich:
                if self.lonmin < 0 and self.lonmax > 0.:
                    lons = np.where(lons < 0, lons+360, lons)
                    self.lonmin = lons.min()
                    self.lonmax = lons.max()
        # read in coastline polygons, only keeping those that
        # intersect map boundary polygon.
        if self.resolution is not None:
            self.coastsegs, self.coastpolygontypes =\
            self._readboundarydata('gshhs',as_polygons=True)
            # reformat for use in matplotlib.patches.Polygon.
            self.coastpolygons = []
            for seg in self.coastsegs:
                x, y = list(zip(*seg))
                self.coastpolygons.append((x,y))
            # replace coastsegs with line segments (instead of polygons)
            self.coastsegs, types =\
            self._readboundarydata('gshhs',as_polygons=False)
            self.coastsegs = [sg for sg in self.coastsegs if len(sg) > 0]
        # create geos Polygon structures for land areas.
        # currently only used in is_land method.
        self.landpolygons=[]
        self.lakepolygons=[]
        if self.resolution is not None and len(self.coastpolygons) > 0:
            #self.islandinlakepolygons=[]
            #self.lakeinislandinlakepolygons=[]
            x, y = list(zip(*self.coastpolygons))
            for x,y,typ in zip(x,y,self.coastpolygontypes):
                b = np.asarray([x,y]).T
                if typ == 1: self.landpolygons.append(_geoslib.Polygon(b))
                if typ == 2: self.lakepolygons.append(_geoslib.Polygon(b))
                #if typ == 3: self.islandinlakepolygons.append(_geoslib.Polygon(b))
                #if typ == 4: self.lakeinislandinlakepolygons.append(_geoslib.Polygon(b))

    # set __init__'s docstring
    __init__.__doc__ = _Basemap_init_doc

    def __call__(self,x,y,inverse=False):
        """
        Calling a Basemap class instance with the arguments lon, lat will
        convert lon/lat (in degrees) to x/y map projection
        coordinates (in meters).  If optional keyword ``inverse`` is
        True (default is False), the inverse transformation from x/y
        to lon/lat is performed.

        For cylindrical equidistant projection (``cyl``), this
        does nothing (i.e. x,y == lon,lat).

        For non-cylindrical projections, the inverse transformation
        always returns longitudes between -180 and 180 degrees. For
        cylindrical projections (self.projection == ``cyl``,
        ``cea``, ``mill``, ``gall`` or ``merc``)
        the inverse transformation will return longitudes between
        self.llcrnrlon and self.llcrnrlat.

        Input arguments lon, lat can be either scalar floats,
        sequences, or numpy arrays.
        """
        if self.celestial:
            # don't assume center of map is at greenwich
            # (only relevant for cyl or pseudo-cyl projections)
            if self.projection in _pseudocyl or self.projection in _cylproj:
                lon_0=self.projparams['lon_0']
            else:
                lon_0 = 0.
        if self.celestial and not inverse:
            try:
                x = 2.*lon_0-x
            except TypeError:
                x = [2*lon_0-xx for xx in x]
        if self.projection == 'rotpole' and inverse:
            try:
                x = _dg2rad*x
            except TypeError:
                x = [_dg2rad*xx for xx in x]
            try:
                y = _dg2rad*y
            except TypeError:
                y = [_dg2rad*yy for yy in y]
        xout,yout = self.projtran(x,y,inverse=inverse)
        if self.celestial and inverse:
            try:
                xout = -2.*lon_0-xout
            except:
                xout = [-2.*lon_0-xx for xx in xout]
        if self.projection == 'rotpole' and not inverse:
            try:
                xout = _rad2dg*xout
                xout = np.where(xout < 0., xout+360, xout)
            except TypeError:
                xout = [_rad2dg*xx for xx in xout]
                xout = [xx+360. if xx < 0 else xx for xx in xout]
            try:
                yout = _rad2dg*yout
            except TypeError:
                yout = [_rad2dg*yy for yy in yout]
        return xout,yout

    def makegrid(self,nx,ny,returnxy=False):
        """
        return arrays of shape (ny,nx) containing lon,lat coordinates of
        an equally spaced native projection grid.

        If ``returnxy = True``, the x,y values of the grid are returned also.
        """
        return self.projtran.makegrid(nx,ny,returnxy=returnxy)

    def _readboundarydata(self,name,as_polygons=False):
        """
        read boundary data, clip to map projection region.
        """
        msg = dedent("""
        Unable to open boundary dataset file. Only the 'crude', 'low' and
        'intermediate' resolution datasets are installed by default. If you
        are requesting a 'high' or 'full' resolution dataset, you need to
        install the `basemap-data-hires` package.""")
        # only gshhs coastlines can be polygons.
        if name != 'gshhs': as_polygons=False
        try:
            bdatfile = open(os.path.join(basemap_datadir,name+'_'+self.resolution+'.dat'),'rb')
            bdatmetafile = open(os.path.join(basemap_datadir,name+'meta_'+self.resolution+'.dat'),'r')
        except:
            raise IOError(msg)
        polygons = []
        polygon_types = []
        # coastlines are polygons, other boundaries are line segments.
        if name == 'gshhs':
            Shape = _geoslib.Polygon
        else:
            Shape = _geoslib.LineString
        # see if map projection region polygon contains a pole.
        NPole = _geoslib.Point(self(0.,90.))
        SPole = _geoslib.Point(self(0.,-90.))
        boundarypolyxy = self._boundarypolyxy
        boundarypolyll = self._boundarypolyll
        hasNP = NPole.within(boundarypolyxy)
        hasSP = SPole.within(boundarypolyxy)
        containsPole = hasNP or hasSP
        # these projections cannot cross pole.
        if containsPole and\
            self.projection in _cylproj + _pseudocyl + ['geos']:
            raise ValueError('%s projection cannot cross pole'%(self.projection))
        # make sure some projections have has containsPole=True
        # we will compute the intersections in stereographic
        # coordinates, then transform back. This is
        # because these projections are only defined on a hemisphere, and
        # some boundary features (like Eurasia) would be undefined otherwise.
        tostere =\
        ['omerc','ortho','gnom','nsper','nplaea','npaeqd','splaea','spaeqd']
        if self.projection in tostere and name == 'gshhs':
            containsPole = True
            lon_0=self.projparams['lon_0']
            lat_0=self.projparams['lat_0']
            re = self.projparams['R']
            # center of stereographic projection restricted to be
            # nearest one of 6 points on the sphere (every 90 deg lat/lon).
            lon0 = 90.*(np.around(lon_0/90.))
            lat0 = 90.*(np.around(lat_0/90.))
            if np.abs(int(lat0)) == 90: lon0=0.
            maptran = pyproj.Proj(proj='stere',lon_0=lon0,lat_0=lat0,R=re)
            # boundary polygon for ortho/gnom/nsper projection
            # in stereographic coordinates.
            b = self._boundarypolyll.boundary
            blons = b[:,0]; blats = b[:,1]
            b[:,0], b[:,1] = maptran(blons, blats)
            boundarypolyxy = _geoslib.Polygon(b)
        for line in bdatmetafile:
            linesplit = line.split()
            area = float(linesplit[1])
            south = float(linesplit[3])
            north = float(linesplit[4])
            crossdatelineE=False; crossdatelineW=False
            if name == 'gshhs':
                id = linesplit[7]
                if id.endswith('E'):
                    crossdatelineE = True
                elif id.endswith('W'):
                    crossdatelineW = True
            # make sure south/north limits of dateline crossing polygons
            # (Eurasia) are the same, since they will be merged into one.
            # (this avoids having one filtered out and not the other).
            if crossdatelineE:
                south_save=south
                north_save=north
            if crossdatelineW:
                south=south_save
                north=north_save
            if area < 0.: area = 1.e30
            useit = self.latmax>=south and self.latmin<=north and area>self.area_thresh
            if useit:
                typ = int(linesplit[0])
                npts = int(linesplit[2])
                offsetbytes = int(linesplit[5])
                bytecount = int(linesplit[6])
                bdatfile.seek(offsetbytes,0)
                # read in binary string convert into an npts by 2
                # numpy array (first column is lons, second is lats).
                polystring = bdatfile.read(bytecount)
                # binary data is little endian.
                b = np.array(np.frombuffer(polystring,dtype='<f4'),'f8')
                b.shape = (npts,2)
                b2 = b.copy()
                # merge polygons that cross dateline.
                poly = Shape(b)
                # hack to try to avoid having Antartica filled polygon
                # covering entire map (if skipAnart = False, this happens
                # for ortho lon_0=-120, lat_0=60, for example).
                skipAntart = self.projection in tostere and south < -89 and \
                 not hasSP
                if crossdatelineE and not skipAntart:
                    if not poly.is_valid(): poly=poly.fix()
                    polyE = poly
                    continue
                elif crossdatelineW and not skipAntart:
                    if not poly.is_valid(): poly=poly.fix()
                    b = poly.boundary
                    b[:,0] = b[:,0]+360.
                    poly = Shape(b)
                    poly = poly.union(polyE)
                    if not poly.is_valid(): poly=poly.fix()
                    b = poly.boundary
                    b2 = b.copy()
                    # fix Antartica.
                    if name == 'gshhs' and south < -89:
                        b = b[3:,:]
                        b2 = b.copy()
                        poly = Shape(b)
                # if map boundary polygon is a valid one in lat/lon
                # coordinates (i.e. it does not contain either pole),
                # the intersections of the boundary geometries
                # and the map projection region can be computed before
                # transforming the boundary geometry to map projection
                # coordinates (this saves time, especially for small map
                # regions and high-resolution boundary geometries).
                if not containsPole:
                    # close Antarctica.
                    if name == 'gshhs' and south < -89:
                        lons2 = b[:,0]
                        lats = b[:,1]
                        lons1 = lons2 - 360.
                        lons3 = lons2 + 360.
                        lons = lons1.tolist()+lons2.tolist()+lons3.tolist()
                        lats = lats.tolist()+lats.tolist()+lats.tolist()
                        lonstart,latstart = lons[0], lats[0]
                        lonend,latend = lons[-1], lats[-1]
                        lons.insert(0,lonstart)
                        lats.insert(0,-90.)
                        lons.append(lonend)
                        lats.append(-90.)
                        b = np.empty((len(lons),2),np.float64)
                        b[:,0] = lons; b[:,1] = lats
                        poly = Shape(b)
                        if not poly.is_valid(): poly=poly.fix()
                        # if polygon instersects map projection
                        # region, process it.
                        if poly.intersects(boundarypolyll):
                            if name != 'gshhs' or as_polygons:
                                geoms = poly.intersection(boundarypolyll)
                            else:
                                # convert polygons to line segments
                                poly = _geoslib.LineString(poly.boundary)
                                geoms = poly.intersection(boundarypolyll)
                            # iterate over geometries in intersection.
                            for psub in geoms:
                                b = psub.boundary
                                blons = b[:,0]; blats = b[:,1]
                                bx, by = self(blons, blats)
                                polygons.append(list(zip(bx,by)))
                                polygon_types.append(typ)
                    else:
                        # create duplicate polygons shifted by -360 and +360
                        # (so as to properly treat polygons that cross
                        # Greenwich meridian).
                        b2[:,0] = b[:,0]-360
                        poly1 = Shape(b2)
                        b2[:,0] = b[:,0]+360
                        poly2 = Shape(b2)
                        polys = [poly1,poly,poly2]
                        for poly in polys:
                            # try to fix "non-noded intersection" errors.
                            if not poly.is_valid(): poly=poly.fix()
                            # if polygon instersects map projection
                            # region, process it.
                            if poly.intersects(boundarypolyll):
                                if name != 'gshhs' or as_polygons:
                                    geoms = poly.intersection(boundarypolyll)
                                else:
                                    # convert polygons to line segments
                                    # note: use fix method here or Eurasia
                                    # line segments sometimes disappear.
                                    poly = _geoslib.LineString(poly.fix().boundary)
                                    geoms = poly.intersection(boundarypolyll)
                                # iterate over geometries in intersection.
                                for psub in geoms:
                                    b = psub.boundary
                                    blons = b[:,0]; blats = b[:,1]
                                    # transformation from lat/lon to
                                    # map projection coordinates.
                                    bx, by = self(blons, blats)
                                    if not as_polygons or len(bx) > 4:
                                        polygons.append(list(zip(bx,by)))
                                        polygon_types.append(typ)
                # if map boundary polygon is not valid in lat/lon
                # coordinates, compute intersection between map
                # projection region and boundary geometries in map
                # projection coordinates.
                else:
                    # transform coordinates from lat/lon
                    # to map projection coordinates.
                    # special case for ortho/gnom/nsper, compute coastline polygon
                    # vertices in stereographic coords.
                    if name == 'gshhs' and as_polygons and self.projection in tostere:
                        b[:,0], b[:,1] = maptran(b[:,0], b[:,1])
                    else:
                        b[:,0], b[:,1] = self(b[:,0], b[:,1])
                    goodmask = np.logical_and(b[:,0]<1.e20,b[:,1]<1.e20)
                    # if less than two points are valid in
                    # map proj coords, skip this geometry.
                    if np.sum(goodmask) <= 1: continue
                    if name != 'gshhs' or (name == 'gshhs' and not as_polygons):
                        # if not a polygon,
                        # just remove parts of geometry that are undefined
                        # in this map projection.
                        bx = np.compress(goodmask, b[:,0])
                        by = np.compress(goodmask, b[:,1])
                        # split coastline segments that jump across entire plot.
                        xd = (bx[1:]-bx[0:-1])**2
                        yd = (by[1:]-by[0:-1])**2
                        dist = np.sqrt(xd+yd)
                        split = dist > 0.1*(self.xmax-self.xmin)
                        if np.sum(split) and self.projection not in _cylproj:
                            ind = (np.compress(split,np.squeeze(split*np.indices(xd.shape)))+1).tolist()
                            iprev = 0
                            ind.append(len(xd))
                            for i in ind:
                                # don't add empty lists.
                                if len(list(range(iprev,i))):
                                    polygons.append(list(zip(bx[iprev:i],by[iprev:i])))
                                iprev = i
                        else:
                            polygons.append(list(zip(bx,by)))
                        polygon_types.append(typ)
                        continue
                    # create a GEOS geometry object.
                    if name == 'gshhs' and not as_polygons:
                        # convert polygons to line segments
                        poly = _geoslib.LineString(poly.boundary)
                    else:
                        # this is a workaround to avoid
                        # GEOS_ERROR: CGAlgorithmsDD::orientationIndex encountered NaN/Inf numbers
                        b[np.isposinf(b)] = 1e20
                        b[np.isneginf(b)] = -1e20
                        poly = Shape(b)
                    # this is a workaround to avoid
                    # "GEOS_ERROR: TopologyException:
                    # found non-noded intersection between ..."
                    if not poly.is_valid(): poly=poly.fix()
                    # if geometry instersects map projection
                    # region, and doesn't have any invalid points, process it.
                    if goodmask.all() and poly.intersects(boundarypolyxy):
                        # if geometry intersection calculation fails,
                        # just move on.
                        try:
                            geoms = poly.intersection(boundarypolyxy)
                        except:
                            continue
                        # iterate over geometries in intersection.
                        for psub in geoms:
                            b = psub.boundary
                            # if projection in ['ortho','gnom','nsper'],
                            # transform polygon from stereographic
                            # to ortho/gnom/nsper coordinates.
                            if self.projection in tostere:
                                # if coastline polygon covers more than 99%
                                # of map region for fulldisk projection,
                                # it's probably bogus, so skip it.
                                #areafrac = psub.area()/boundarypolyxy.area()
                                #if self.projection == ['ortho','nsper']:
                                #    if name == 'gshhs' and\
                                #       self._fulldisk and\
                                #       areafrac > 0.99: continue
                                # inverse transform from stereographic
                                # to lat/lon.
                                b[:,0], b[:,1] = maptran(b[:,0], b[:,1], inverse=True)
                                # orthographic/gnomonic/nsper.
                                b[:,0], b[:,1]= self(b[:,0], b[:,1])
                            if not as_polygons or len(b) > 4:
                                polygons.append(list(zip(b[:,0],b[:,1])))
                                polygon_types.append(typ)
        bdatfile.close()
        bdatmetafile.close()
        return polygons, polygon_types

    def _getmapboundary(self):
        """
        create map boundary polygon (in lat/lon and x/y coordinates)
        """
        nx = 100; ny = 100
        maptran = self
        if self.projection in ['ortho','geos','nsper']:
            # circular region.
            thetas = np.linspace(0.,2.*np.pi,2*nx*ny)[:-1]
            rminor = self._height
            rmajor = self._width
            x = rmajor*np.cos(thetas) + rmajor
            y = rminor*np.sin(thetas) + rminor
            b = np.empty((len(x),2),np.float64)
            b[:,0]=x; b[:,1]=y
            boundaryxy = _geoslib.Polygon(b)
            # compute proj instance for full disk, if necessary.
            if not self._fulldisk:
                projparms = self.projparams.copy()
                del projparms['x_0']
                del projparms['y_0']
                if self.projection == 'ortho':
                    llcrnrx = -self.rmajor
                    llcrnry = -self.rmajor
                    urcrnrx = -llcrnrx
                    urcrnry = -llcrnry
                else:
                    llcrnrx = -self._width
                    llcrnry = -self._height
                    urcrnrx = -llcrnrx
                    urcrnry = -llcrnry
                projparms['x_0']=-llcrnrx
                projparms['y_0']=-llcrnry
                maptran = pyproj.Proj(projparms)
        elif self.projection == 'aeqd' and self._fulldisk:
            # circular region.
            thetas = np.linspace(0.,2.*np.pi,2*nx*ny)[:-1]
            rminor = self._height
            rmajor = self._width
            x = rmajor*np.cos(thetas) + rmajor
            y = rminor*np.sin(thetas) + rminor
            b = np.empty((len(x),2),np.float64)
            b[:,0]=x; b[:,1]=y
            boundaryxy = _geoslib.Polygon(b)
        elif self.projection in _pseudocyl:
            nx = 10*nx; ny = 10*ny
            # quasi-elliptical region.
            lon_0 = self.projparams['lon_0']
            # left side
            lats1 = np.linspace(-89.9999,89.9999,ny).tolist()
            lons1 = len(lats1)*[lon_0-179.9]
            # top.
            lons2 = np.linspace(lon_0-179.9,lon_0+179.9,nx).tolist()
            lats2 = len(lons2)*[89.9999]
            # right side
            lats3 = np.linspace(89.9999,-89.9999,ny).tolist()
            lons3 = len(lats3)*[lon_0+179.9]
            # bottom.
            lons4 = np.linspace(lon_0+179.9,lon_0-179.9,nx).tolist()
            lats4 = len(lons4)*[-89.9999]
            lons = np.array(lons1+lons2+lons3+lons4,np.float64)
            lats = np.array(lats1+lats2+lats3+lats4,np.float64)
            x, y = maptran(lons,lats)
            b = np.empty((len(x),2),np.float64)
            b[:,0]=x; b[:,1]=y
            boundaryxy = _geoslib.Polygon(b)
        else: # all other projections are rectangular.
            nx = 100*nx; ny = 100*ny
            # left side (x = xmin, ymin <= y <= ymax)
            yy = np.linspace(self.ymin, self.ymax, ny)[:-1]
            x = len(yy)*[self.xmin]; y = yy.tolist()
            # top (y = ymax, xmin <= x <= xmax)
            xx = np.linspace(self.xmin, self.xmax, nx)[:-1]
            x = x + xx.tolist()
            y = y + len(xx)*[self.ymax]
            # right side (x = xmax, ymin <= y <= ymax)
            yy = np.linspace(self.ymax, self.ymin, ny)[:-1]
            x = x + len(yy)*[self.xmax]; y = y + yy.tolist()
            # bottom (y = ymin, xmin <= x <= xmax)
            xx = np.linspace(self.xmax, self.xmin, nx)[:-1]
            x = x + xx.tolist()
            y = y + len(xx)*[self.ymin]
            x = np.array(x,np.float64)
            y = np.array(y,np.float64)
            b = np.empty((4,2),np.float64)
            b[:,0]=[self.xmin,self.xmin,self.xmax,self.xmax]
            b[:,1]=[self.ymin,self.ymax,self.ymax,self.ymin]
            boundaryxy = _geoslib.Polygon(b)
        if self.projection in _cylproj:
            # make sure map boundary doesn't quite include pole.
            if self.urcrnrlat > 89.9999:
                urcrnrlat = 89.9999
            else:
                urcrnrlat = self.urcrnrlat
            if self.llcrnrlat < -89.9999:
                llcrnrlat = -89.9999
            else:
                llcrnrlat = self.llcrnrlat
            lons = [self.llcrnrlon, self.llcrnrlon, self.urcrnrlon, self.urcrnrlon]
            lats = [llcrnrlat, urcrnrlat, urcrnrlat, llcrnrlat]
            self.boundarylonmin = min(lons)
            self.boundarylonmax = max(lons)
            x, y = self(lons, lats)
            b = np.empty((len(x),2),np.float64)
            b[:,0]=x; b[:,1]=y
            boundaryxy = _geoslib.Polygon(b)
        else:
            if self.projection not in _pseudocyl:
                lons, lats = maptran(x,y,inverse=True)
                # fix lons so there are no jumps.
                n = 1
                lonprev = lons[0]
                for lon,lat in zip(lons[1:],lats[1:]):
                    if np.abs(lon-lonprev) > 90.:
                        if lonprev < 0:
                            lon = lon - 360.
                        else:
                            lon = lon + 360
                        lons[n] = lon
                    lonprev = lon
                    n = n + 1
                self.boundarylonmin = lons.min()
                self.boundarylonmax = lons.max()
                # for circular full disk projections where boundary is
                # a latitude circle, set boundarylonmax and boundarylonmin
                # to cover entire world (so parallels will be drawn).
                if self._fulldisk and \
                   np.abs(self.boundarylonmax-self.boundarylonmin) < 1.:
                   self.boundarylonmin = -180.
                   self.boundarylonmax = 180.
        b = np.empty((len(lons),2),np.float64)
        b[:,0] = lons; b[:,1] = lats
        boundaryll = _geoslib.Polygon(b)
        return lons, lats, boundaryll, boundaryxy


    def drawmapboundary(self,color='k',linewidth=1.0,fill_color=None,\
                        zorder=None,ax=None):
        """
        draw boundary around map projection region, optionally
        filling interior of region.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keyword          Description
        ==============   ====================================================
        linewidth        line width for boundary (default 1.)
        color            color of boundary line (default black)
        fill_color       fill the map region background with this
                         color (default is to fill with axis
                         background color). If set to the string
                         'none', no filling is done.
        zorder           sets the zorder for filling map background
                         (default 0).
        ax               axes instance to use
                         (default None, use default axes instance).
        ==============   ====================================================

        returns matplotlib.collections.PatchCollection representing map boundary.
        """
        # get current axes instance (if none specified).
        ax = ax or self._check_ax()
        # if no fill_color given, use axes background color.
        # if fill_color is string 'none', really don't fill.
        if fill_color is None:
            if _matplotlib_version >= '2.0':
                fill_color = ax.get_facecolor()
            else:
                fill_color = ax.get_axis_bgcolor()
        elif fill_color == 'none' or fill_color == 'None':
            fill_color = None
        limb = None
        if self.projection in ['ortho','geos','nsper'] or (self.projection=='aeqd' and\
           self._fulldisk):
            limb = Ellipse((self._width,self._height),2.*self._width,2.*self._height)
        if self.projection in ['ortho','geos','nsper','aeqd'] and self._fulldisk:
            # elliptical region.
            ax.set_frame_on(False)
        elif self.projection in _pseudocyl:  # elliptical region.
            ax.set_frame_on(False)
            nx = 100; ny = 100
            if self.projection == 'vandg':
                nx = 10*nx; ny = 10*ny
            # quasi-elliptical region.
            lon_0 = self.projparams['lon_0']
            # left side
            lats1 = np.linspace(-89.9999,89.99999,ny).tolist()
            lons1 = len(lats1)*[lon_0-179.9]
            # top.
            lons2 = np.linspace(lon_0-179.9999,lon_0+179.9999,nx).tolist()
            lats2 = len(lons2)*[89.9999]
            # right side
            lats3 = np.linspace(89.9999,-89.9999,ny).tolist()
            lons3 = len(lats3)*[lon_0+179.9999]
            # bottom.
            lons4 = np.linspace(lon_0+179.9999,lon_0-179.9999,nx).tolist()
            lats4 = len(lons4)*[-89.9999]
            lons = np.array(lons1+lons2+lons3+lons4,np.float64)
            lats = np.array(lats1+lats2+lats3+lats4,np.float64)
            x, y = self(lons,lats)
            xy = list(zip(x,y))
            limb = Polygon(xy)
        elif self.round:
            ax.set_frame_on(False)
            limb = Circle((0.5*(self.xmax+self.xmin),0.5*(self.ymax+self.ymin)),
                    radius=0.5*(self.xmax-self.xmin),fc='none')
        else: # all other projections are rectangular.
            ax.set_frame_on(True)
            for spine in ax.spines.values():
                spine.set_linewidth(linewidth)
                spine.set_edgecolor(color)
                if zorder is not None:
                    spine.set_zorder(zorder)
            if self.projection not in ['geos','ortho','nsper']:
                limb = ax.patch

        if limb is not None:
            if limb is not ax.patch:
                ax.add_patch(limb)
            self._mapboundarydrawn = limb
            if fill_color is None:
                limb.set_fill(False)
            else:
                limb.set_facecolor(fill_color)
                limb.set_zorder(0)
            limb.set_edgecolor(color)
            limb.set_linewidth(linewidth)
            if zorder is not None:
                limb.set_zorder(zorder)
            limb.set_clip_on(True)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        return limb

    def fillcontinents(self,color='0.8',lake_color=None,ax=None,zorder=None,alpha=None):
        """
        Fill continents.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keyword          Description
        ==============   ====================================================
        color            color to fill continents (default gray).
        lake_color       color to fill inland lakes (default axes background).
        ax               axes instance (overrides default axes instance).
        zorder           sets the zorder for the continent polygons (if not
                         specified, uses default zorder for a Polygon patch).
                         Set to zero if you want to paint over the filled
                         continents).
        alpha            sets alpha transparency for continent polygons
        ==============   ====================================================

        After filling continents, lakes are re-filled with
        axis background color.

        returns a list of matplotlib.patches.Polygon objects.
        """
        if self.resolution is None:
            raise AttributeError('there are no boundary datasets associated with this Basemap instance')
        # get current axes instance (if none specified).
        ax = ax or self._check_ax()
        # get axis background color.
        if _matplotlib_version >= '2.0':
            axisbgc = ax.get_facecolor()
        else:
            axisbgc = ax.get_axis_bgcolor()
        npoly = 0
        polys = []
        for x,y in self.coastpolygons:
            xa = np.array(x,np.float32)
            ya = np.array(y,np.float32)
        # check to see if all four corners of domain in polygon (if so,
        # don't draw since it will just fill in the whole map).
        # ** turn this off for now since it prevents continents that
        # fill the whole map from being filled **
            #delx = 10; dely = 10
            #if self.projection in ['cyl']:
            #    delx = 0.1
            #    dely = 0.1
            #test1 = np.fabs(xa-self.urcrnrx) < delx
            #test2 = np.fabs(xa-self.llcrnrx) < delx
            #test3 = np.fabs(ya-self.urcrnry) < dely
            #test4 = np.fabs(ya-self.llcrnry) < dely
            #hasp1 = np.sum(test1*test3)
            #hasp2 = np.sum(test2*test3)
            #hasp4 = np.sum(test2*test4)
            #hasp3 = np.sum(test1*test4)
            #if not hasp1 or not hasp2 or not hasp3 or not hasp4:
            if 1:
                xy = list(zip(xa.tolist(),ya.tolist()))
                if self.coastpolygontypes[npoly] not in [2,4]:
                    poly = Polygon(xy,facecolor=color,edgecolor=color,linewidth=0)
                else: # lakes filled with background color by default
                    if lake_color is None:
                        poly = Polygon(xy,facecolor=axisbgc,edgecolor=axisbgc,linewidth=0)
                    else:
                        poly = Polygon(xy,facecolor=lake_color,edgecolor=lake_color,linewidth=0)
                if zorder is not None:
                    poly.set_zorder(zorder)
                if alpha is not None:
                    poly.set_alpha(alpha)
                ax.add_patch(poly)
                polys.append(poly)
            npoly = npoly + 1
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip continent polygons to map limbs
        polys,c = self._cliplimb(ax,polys)
        return polys

    def _cliplimb(self,ax,coll):
        if not self._mapboundarydrawn:
            return coll, None
        c = self._mapboundarydrawn
        if c not in ax.patches:
            p = ax.add_patch(c)
            #p.set_clip_on(False)
        try:
            coll.set_clip_path(c)
        except:
            for item in coll:
                item.set_clip_path(c)
        return coll,c

    def drawcoastlines(self,linewidth=1.,linestyle='solid',color='k',antialiased=1,ax=None,zorder=None):
        """
        Draw coastlines.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keyword          Description
        ==============   ====================================================
        linewidth        coastline width (default 1.)
        linestyle        coastline linestyle (default solid)
        color            coastline color (default black)
        antialiased      antialiasing switch for coastlines (default True).
        ax               axes instance (overrides default axes instance)
        zorder           sets the zorder for the coastlines (if not specified,
                         uses default zorder for
                         matplotlib.patches.LineCollections).
        ==============   ====================================================

        returns a matplotlib.patches.LineCollection object.
        """
        if self.resolution is None:
            raise AttributeError('there are no boundary datasets associated with this Basemap instance')
        # get current axes instance (if none specified).
        ax = ax or self._check_ax()
        coastlines = LineCollection(self.coastsegs,antialiaseds=(antialiased,))
        coastlines.set_color(color)
        coastlines.set_linestyle(linestyle)
        coastlines.set_linewidth(linewidth)
        coastlines.set_label('_nolabel_')
        if zorder is not None:
            coastlines.set_zorder(zorder)
        ax.add_collection(coastlines)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip to map limbs
        coastlines,c = self._cliplimb(ax,coastlines)
        return coastlines

    def drawcountries(self,linewidth=0.5,linestyle='solid',color='k',antialiased=1,ax=None,zorder=None):
        """
        Draw country boundaries.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keyword          Description
        ==============   ====================================================
        linewidth        country boundary line width (default 0.5)
        linestyle        coastline linestyle (default solid)
        color            country boundary line color (default black)
        antialiased      antialiasing switch for country boundaries (default
                         True).
        ax               axes instance (overrides default axes instance)
        zorder           sets the zorder for the country boundaries (if not
                         specified uses default zorder for
                         matplotlib.patches.LineCollections).
        ==============   ====================================================

        returns a matplotlib.patches.LineCollection object.
        """
        if self.resolution is None:
            raise AttributeError('there are no boundary datasets associated with this Basemap instance')
        # read in country line segments, only keeping those that
        # intersect map boundary polygon.
        if not hasattr(self,'cntrysegs'):
            self.cntrysegs, types = self._readboundarydata('countries')
        # get current axes instance (if none specified).
        ax = ax or self._check_ax()
        countries = LineCollection(self.cntrysegs,antialiaseds=(antialiased,))
        countries.set_color(color)
        countries.set_linestyle(linestyle)
        countries.set_linewidth(linewidth)
        countries.set_label('_nolabel_')
        if zorder is not None:
            countries.set_zorder(zorder)
        ax.add_collection(countries)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip countries to map limbs
        countries,c = self._cliplimb(ax,countries)
        return countries

    def drawstates(self,linewidth=0.5,linestyle='solid',color='k',antialiased=1,ax=None,zorder=None):
        """
        Draw state boundaries in Americas.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keyword          Description
        ==============   ====================================================
        linewidth        state boundary line width (default 0.5)
        linestyle        coastline linestyle (default solid)
        color            state boundary line color (default black)
        antialiased      antialiasing switch for state boundaries
                         (default True).
        ax               axes instance (overrides default axes instance)
        zorder           sets the zorder for the state boundaries (if not
                         specified, uses default zorder for
                         matplotlib.patches.LineCollections).
        ==============   ====================================================

        returns a matplotlib.patches.LineCollection object.
        """
        if self.resolution is None:
            raise AttributeError('there are no boundary datasets associated with this Basemap instance')
        # read in state line segments, only keeping those that
        # intersect map boundary polygon.
        if not hasattr(self,'statesegs'):
            self.statesegs, types = self._readboundarydata('states')
        # get current axes instance (if none specified).
        ax = ax or self._check_ax()
        states = LineCollection(self.statesegs,antialiaseds=(antialiased,))
        states.set_color(color)
        states.set_linestyle(linestyle)
        states.set_linewidth(linewidth)
        states.set_label('_nolabel_')
        if zorder is not None:
            states.set_zorder(zorder)
        ax.add_collection(states)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip states to map limbs
        states,c = self._cliplimb(ax,states)
        return states

    def drawcounties(self,linewidth=0.1,linestyle='solid',color='k',antialiased=1,
                     facecolor='none',ax=None,zorder=None,drawbounds=False):
        """
        Draw county boundaries in US. The county boundary shapefile
        originates with the NOAA Coastal Geospatial Data Project
        (http://coastalgeospatial.noaa.gov/data_gis.html).

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keyword          Description
        ==============   ====================================================
        linewidth        county boundary line width (default 0.1)
        linestyle        coastline linestyle (default solid)
        color            county boundary line color (default black)
        facecolor        fill color of county (default is no fill)
        antialiased      antialiasing switch for county boundaries
                         (default True).
        ax               axes instance (overrides default axes instance)
        zorder           sets the zorder for the county boundaries (if not
                         specified, uses default zorder for
                         matplotlib.patches.LineCollections).
        ==============   ====================================================

        returns a matplotlib.patches.LineCollection object.
        """
        ax = ax or self._check_ax()
        gis_file = os.path.join(basemap_datadir,'UScounties')
        county_info = self.readshapefile(gis_file,'counties',\
                      default_encoding='latin-1',drawbounds=drawbounds)
        counties = [coords for coords in self.counties]
        counties = PolyCollection(counties)
        counties.set_linestyle(linestyle)
        counties.set_linewidth(linewidth)
        counties.set_edgecolor(color)
        counties.set_facecolor(facecolor)
        counties.set_label('counties')
        if zorder:
            counties.set_zorder(zorder)
        ax.add_collection(counties)
        return counties

    def drawrivers(self,linewidth=0.5,linestyle='solid',color='k',antialiased=1,ax=None,zorder=None):
        """
        Draw major rivers.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keyword          Description
        ==============   ====================================================
        linewidth        river boundary line width (default 0.5)
        linestyle        coastline linestyle (default solid)
        color            river boundary line color (default black)
        antialiased      antialiasing switch for river boundaries (default
                         True).
        ax               axes instance (overrides default axes instance)
        zorder           sets the zorder for the rivers (if not
                         specified uses default zorder for
                         matplotlib.patches.LineCollections).
        ==============   ====================================================

        returns a matplotlib.patches.LineCollection object.
        """
        if self.resolution is None:
            raise AttributeError('there are no boundary datasets associated with this Basemap instance')
        # read in river line segments, only keeping those that
        # intersect map boundary polygon.
        if not hasattr(self,'riversegs'):
            self.riversegs, types = self._readboundarydata('rivers')
        # get current axes instance (if none specified).
        ax = ax or self._check_ax()
        rivers = LineCollection(self.riversegs,antialiaseds=(antialiased,))
        rivers.set_color(color)
        rivers.set_linestyle(linestyle)
        rivers.set_linewidth(linewidth)
        rivers.set_label('_nolabel_')
        if zorder is not None:
            rivers.set_zorder(zorder)
        ax.add_collection(rivers)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip rivers to map limbs
        rivers,c = self._cliplimb(ax,rivers)
        return rivers

    def is_land(self,xpt,ypt):
        """
        Returns True if the given x,y point (in projection coordinates) is
        over land, False otherwise.  The definition of land is based upon
        the GSHHS coastline polygons associated with the class instance.
        Points over lakes inside land regions are not counted as land points.
        """
        if self.resolution is None: return None
        landpt = False
        for poly in self.landpolygons:
            landpt = _geoslib.Point((xpt,ypt)).within(poly)
            if landpt: break
        lakept = False
        for poly in self.lakepolygons:
            lakept = _geoslib.Point((xpt,ypt)).within(poly)
            if lakept: break
        return landpt and not lakept

    def readshapefile(self,shapefile,name,drawbounds=True,zorder=None,
                      linewidth=0.5,color='k',antialiased=1,ax=None,
                      default_encoding='utf-8'):
        """
        Read in shape file, optionally draw boundaries on map.

        .. note::
          - Assumes shapes are 2D
          - only works for Point, MultiPoint, Polyline and Polygon shapes.
          - vertices/points must be in geographic (lat/lon) coordinates.

        Mandatory Arguments:

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Argument         Description
        ==============   ====================================================
        shapefile        path to shapefile components.  Example:
                         shapefile='/home/jeff/esri/world_borders' assumes
                         that world_borders.shp, world_borders.shx and
                         world_borders.dbf live in /home/jeff/esri.
        name             name for Basemap attribute to hold the shapefile
                         vertices or points in map projection
                         coordinates. Class attribute name+'_info' is a list
                         of dictionaries, one for each shape, containing
                         attributes of each shape from dbf file, For
                         example, if name='counties', self.counties
                         will be a list of x,y vertices for each shape in
                         map projection  coordinates and self.counties_info
                         will be a list of dictionaries with shape
                         attributes.  Rings in individual Polygon
                         shapes are split out into separate polygons, and
                         additional keys 'RINGNUM' and 'SHAPENUM' are added
                         to the shape attribute dictionary.
        ==============   ====================================================

        The following optional keyword arguments are only relevant for Polyline
        and Polygon shape types, for Point and MultiPoint shapes they are
        ignored.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keyword          Description
        ==============   ====================================================
        drawbounds       draw boundaries of shapes (default True).
        zorder           shape boundary zorder (if not specified,
                         default for mathplotlib.lines.LineCollection
                         is used).
        linewidth        shape boundary line width (default 0.5)
        color            shape boundary line color (default black)
        antialiased      antialiasing switch for shape boundaries
                         (default True).
        ax               axes instance (overrides default axes instance)
        ==============   ====================================================

        A tuple (num_shapes, type, min, max) containing shape file info
        is returned.
        num_shapes is the number of shapes, type is the type code (one of
        the SHPT* constants defined in the shapelib module, see
        http://shapelib.maptools.org/shp_api.html) and min and
        max are 4-element lists with the minimum and maximum values of the
        vertices. If ``drawbounds=True`` a
        matplotlib.patches.LineCollection object is appended to the tuple.
        """
        import shapefile as shp
        from shapefile import Reader
        shp.default_encoding = default_encoding
        if not os.path.exists('%s.shp'%shapefile):
            raise IOError('cannot locate %s.shp'%shapefile)
        if not os.path.exists('%s.shx'%shapefile):
            raise IOError('cannot locate %s.shx'%shapefile)
        if not os.path.exists('%s.dbf'%shapefile):
            raise IOError('cannot locate %s.dbf'%shapefile)
        # open shapefile, read vertices for each object, convert
        # to map projection coordinates (only works for 2D shape types).
        try:
            shf = Reader(shapefile, encoding=default_encoding)
        except:
            raise IOError('error reading shapefile %s.shp' % shapefile)
        fields = shf.fields
        coords = []; attributes = []
        msg=dedent("""
        shapefile must have lat/lon vertices  - it looks like this one has vertices
        in map projection coordinates. You can convert the shapefile to geographic
        coordinates using the shpproj utility from the shapelib tools
        (http://shapelib.maptools.org/shapelib-tools.html)""")
        shptype = shf.shapes()[0].shapeType
        bbox = shf.bbox.tolist()
        info = (shf.numRecords,shptype,bbox[0:2]+[0.,0.],bbox[2:]+[0.,0.])
        npoly = 0
        for shprec in shf.shapeRecords():
            shp = shprec.shape; rec = shprec.record
            npoly = npoly + 1
            if shptype != shp.shapeType:
                raise ValueError('readshapefile can only handle a single shape type per file')
            if shptype not in [1,3,5,8]:
                raise ValueError('readshapefile can only handle 2D shape types')
            verts = shp.points
            if shptype in [1,8]: # a Point or MultiPoint shape.
                lons, lats = list(zip(*verts))
                if max(lons) > 721. or min(lons) < -721. or max(lats) > 90.01 or min(lats) < -90.01:
                    raise ValueError(msg)
                # if latitude is slightly greater than 90, truncate to 90
                lats = [max(min(lat, 90.0), -90.0) for lat in lats]
                if len(verts) > 1: # MultiPoint
                    x,y = self(lons, lats)
                    coords.append(list(zip(x,y)))
                else: # single Point
                    x,y = self(lons[0], lats[0])
                    coords.append((x,y))
                attdict={}
                for r,key in zip(rec,fields[1:]):
                    attdict[key[0]]=r
                attributes.append(attdict)
            else: # a Polyline or Polygon shape.
                parts = shp.parts.tolist()
                ringnum = 0
                for indx1,indx2 in zip(parts,parts[1:]+[len(verts)]):
                    ringnum = ringnum + 1
                    lons, lats = list(zip(*verts[indx1:indx2]))
                    if max(lons) > 721. or min(lons) < -721. or max(lats) > 90.01 or min(lats) < -90.01:
                        raise ValueError(msg)
                    # if latitude is slightly greater than 90, truncate to 90
                    lats = [max(min(lat, 90.0), -90.0) for lat in lats]
                    x, y = self(lons, lats)
                    coords.append(list(zip(x,y)))
                    attdict={}
                    for r,key in zip(rec,fields[1:]):
                        attdict[key[0]]=r
                    # add information about ring number to dictionary.
                    attdict['RINGNUM'] = ringnum
                    attdict['SHAPENUM'] = npoly
                    attributes.append(attdict)
        # draw shape boundaries for polylines, polygons  using LineCollection.
        if shptype not in [1,8] and drawbounds:
            # get current axes instance (if none specified).
            ax = ax or self._check_ax()
            # make LineCollections for each polygon.
            lines = LineCollection(coords,antialiaseds=(1,))
            lines.set_color(color)
            lines.set_linewidth(linewidth)
            lines.set_label('_nolabel_')
            if zorder is not None:
               lines.set_zorder(zorder)
            ax.add_collection(lines)
            # set axes limits to fit map region.
            self.set_axes_limits(ax=ax)
            # clip boundaries to map limbs
            lines,c = self._cliplimb(ax,lines)
            info = info + (lines,)
        self.__dict__[name]=coords
        self.__dict__[name+'_info']=attributes
        return info

    def drawparallels(self,circles,color='k',textcolor='k',linewidth=1.,zorder=None, \
                      dashes=[1,1],labels=[0,0,0,0],labelstyle=None, \
                      fmt='%g',xoffset=None,yoffset=None,ax=None,latmax=None,
                      **text_kwargs):
        """
        Draw and label parallels (latitude lines) for values (in degrees)
        given in the sequence ``circles``.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keyword          Description
        ==============   ====================================================
        color            color to draw parallels (default black).
        textcolor        color to draw labels (default black).
        linewidth        line width for parallels (default 1.)
        zorder           sets the zorder for parallels (if not specified,
                         uses default zorder for matplotlib.lines.Line2D
                         objects).
        dashes           dash pattern for parallels (default [1,1], i.e.
                         1 pixel on, 1 pixel off).
        labels           list of 4 values (default [0,0,0,0]) that control
                         whether parallels are labelled where they intersect
                         the left, right, top or bottom of the plot. For
                         example labels=[1,0,0,1] will cause parallels
                         to be labelled where they intersect the left and
                         and bottom of the plot, but not the right and top.
        labelstyle       if set to "+/-", north and south latitudes are
                         labelled with "+" and "-", otherwise they are
                         labelled with "N" and "S".
        fmt              a format string to format the parallel labels
                         (default '%g') **or** a function that takes a
                         latitude value in degrees as it's only argument
                         and returns a formatted string.
        xoffset          label offset from edge of map in x-direction
                         (default is 0.01 times width of map in map
                         projection coordinates).
        yoffset          label offset from edge of map in y-direction
                         (default is 0.01 times height of map in map
                         projection coordinates).
        ax               axes instance (overrides default axes instance)
        latmax           absolute value of latitude to which meridians are drawn
                         (default is 80).
        \**text_kwargs   additional keyword arguments controlling text
                         for labels that are passed on to
                         the text method of the axes instance (see
                         matplotlib.pyplot.text documentation).
        ==============   ====================================================

        returns a dictionary whose keys are the parallel values, and
        whose values are tuples containing lists of the
        matplotlib.lines.Line2D and matplotlib.text.Text instances
        associated with each parallel. Deleting an item from the
        dictionary removes the corresponding parallel from the plot.
        """
        text_kwargs['color']=textcolor # pass textcolor kwarg on to ax.text
        # if celestial=True, don't use "N" and "S" labels.
        if labelstyle is None and self.celestial:
            labelstyle="+/-"
        # get current axes instance (if none specified).
        ax = ax or self._check_ax()
        # don't draw meridians past latmax, always draw parallel at latmax.
        if latmax is None: latmax = 80.
        # offset for labels.
        if yoffset is None:
            yoffset = (self.urcrnry-self.llcrnry)/100.
            if self.aspect > 1:
                yoffset = self.aspect*yoffset
            else:
                yoffset = yoffset/self.aspect
        if xoffset is None:
            xoffset = (self.urcrnrx-self.llcrnrx)/100.

        if self.projection in _cylproj + _pseudocyl:
            lons = np.linspace(self.llcrnrlon, self.urcrnrlon, 10001)
        elif self.projection in ['tmerc']:
            lon_0 = self.projparams['lon_0']
            # tmerc only defined within +/- 90 degrees of lon_0
            lons = np.linspace(lon_0-90,lon_0+90,100001)
        else:
            lonmin = self.boundarylonmin; lonmax = self.boundarylonmax
            lons = np.linspace(lonmin, lonmax, 10001)
        # make sure latmax degree parallel is drawn if projection not merc or cyl or miller
        try:
            circlesl = list(circles)
        except:
            circlesl = circles
        if self.projection not in _cylproj + _pseudocyl:
            if max(circlesl) > 0 and latmax not in circlesl:
                circlesl.append(latmax)
            if min(circlesl) < 0 and -latmax not in circlesl:
                circlesl.append(-latmax)
        xdelta = 0.01*(self.xmax-self.xmin)
        ydelta = 0.01*(self.ymax-self.ymin)
        linecolls = {}
        for circ in circlesl:
            lats = circ*np.ones(len(lons),np.float32)
            x,y = self(lons,lats)
            # remove points outside domain.
            # leave a little slop around edges (3*xdelta)
            # don't really know why, but this appears to be needed to
            # or lines sometimes don't reach edge of plot.
            testx = np.logical_and(x>=self.xmin-3*xdelta,x<=self.xmax+3*xdelta)
            x = np.compress(testx, x)
            y = np.compress(testx, y)
            testy = np.logical_and(y>=self.ymin-3*ydelta,y<=self.ymax+3*ydelta)
            x = np.compress(testy, x)
            y = np.compress(testy, y)
            lines = []
            if len(x) > 1 and len(y) > 1:
                # split into separate line segments if necessary.
                # (not necessary for cylindrical or pseudocylindricl projections)
                xd = (x[1:]-x[0:-1])**2
                yd = (y[1:]-y[0:-1])**2
                dist = np.sqrt(xd+yd)
                if self.projection not in ['cyl','rotpole']:
                    split = dist > self.rmajor/10.
                else:
                    split = dist > 1.
                if np.sum(split) and self.projection not in _cylproj:
                    ind = (np.compress(split,np.squeeze(split*np.indices(xd.shape)))+1).tolist()
                    xl = []
                    yl = []
                    iprev = 0
                    ind.append(len(xd))
                    for i in ind:
                        xl.append(x[iprev:i])
                        yl.append(y[iprev:i])
                        iprev = i
                else:
                    xl = [x]
                    yl = [y]
                # draw each line segment.
                for x,y in zip(xl,yl):
                    # skip if only a point.
                    if len(x) > 1 and len(y) > 1:
                        l = Line2D(x,y,linewidth=linewidth)
                        l.set_color(color)
                        l.set_dashes(dashes)
                        l.set_label('_nolabel_')
                        if zorder is not None:
                            l.set_zorder(zorder)
                        ax.add_line(l)
                        lines.append(l)
            linecolls[circ] = (lines,[])
        # draw labels for parallels
        # parallels not labelled for fulldisk orthographic or geostationary
        if self.projection in ['ortho','geos','nsper','vandg','aeqd'] and max(labels):
            if self.projection == 'vandg' or self._fulldisk:
                sys.stdout.write('Warning: Cannot label parallels on %s basemap' % _projnames[self.projection])
                labels = [0,0,0,0]
        # search along edges of map to see if parallels intersect.
        # if so, find x,y location of intersection and draw a label there.
        dx = (self.xmax-self.xmin)/1000.
        dy = (self.ymax-self.ymin)/1000.
        if self.projection in _pseudocyl:
            lon_0 = self.projparams['lon_0']
        for dolab,side in zip(labels,['l','r','t','b']):
            if not dolab: continue
            # for cylindrical projections, don't draw parallels on top or bottom.
            if self.projection in _cylproj + _pseudocyl and side in ['t','b']: continue
            if side in ['l','r']:
                nmax = int((self.ymax-self.ymin)/dy+1)
                yy = np.linspace(self.llcrnry,self.urcrnry,nmax)
                if side == 'l':
                    if self.projection in _pseudocyl:
                        lats = np.linspace(-89.99,89,99,nmax)
                        if self.celestial:
                            lons = (self.projparams['lon_0']+180.)*np.ones(len(lats),lats.dtype)
                        else:
                            lons = (self.projparams['lon_0']-180.)*np.ones(len(lats),lats.dtype)
                        xx, yy = self(lons, lats)
                    else:
                        xx = self.llcrnrx*np.ones(yy.shape,yy.dtype)
                        lons,lats = self(xx,yy,inverse=True)
                        lons = lons.tolist(); lats = lats.tolist()
                else:
                    if self.projection in _pseudocyl:
                        lats = np.linspace(-89.99,89,99,nmax)
                        if self.celestial:
                           lons = (self.projparams['lon_0']-180.)*np.ones(len(lats),lats.dtype)
                        else:
                           lons = (self.projparams['lon_0']+180.)*np.ones(len(lats),lats.dtype)
                        xx, yy = self(lons, lats)
                    else:
                        xx = self.urcrnrx*np.ones(yy.shape,yy.dtype)
                        lons,lats = self(xx,yy,inverse=True)
                        lons = lons.tolist(); lats = lats.tolist()
                if max(lons) > 1.e20 or max(lats) > 1.e20:
                    raise ValueError('inverse transformation undefined - please adjust the map projection region')
                # adjust so 0 <= lons < 360
                lons = [(lon+360) % 360 for lon in lons]
            else:
                nmax = int((self.xmax-self.xmin)/dx+1)
                xx = np.linspace(self.llcrnrx,self.urcrnrx,nmax)
                if side == 'b':
                    lons,lats = self(xx,self.llcrnry*np.ones(xx.shape,np.float32),inverse=True)
                    lons = lons.tolist(); lats = lats.tolist()
                else:
                    lons,lats = self(xx,self.urcrnry*np.ones(xx.shape,np.float32),inverse=True)
                    lons = lons.tolist(); lats = lats.tolist()
                if max(lons) > 1.e20 or max(lats) > 1.e20:
                    raise ValueError('inverse transformation undefined - please adjust the map projection region')
                # adjust so 0 <= lons < 360
                lons = [(lon+360) % 360 for lon in lons]
            for lat in circles:
                # don't label parallels for round polar plots
                if self.round: continue
                # find index of parallel (there may be two, so
                # search from left and right).
                nl = _searchlist(lats,lat)
                nr = _searchlist(lats[::-1],lat)
                if nr != -1: nr = len(lons)-nr-1
                latlab = _setlatlab(fmt,lat,labelstyle)
                # parallels can intersect each map edge twice.
                for i,n in enumerate([nl,nr]):
                    # don't bother if close to the first label.
                    if i and abs(nr-nl) < 100: continue
                    if n >= 0:
                        t = None
                        if side == 'l':
                            if self.projection in _pseudocyl:
                                if self.celestial:
                                    xlab,ylab = self(lon_0+179.9,lat)
                                else:
                                    xlab,ylab = self(lon_0-179.9,lat)
                            else:
                                xlab = self.llcrnrx
                            xlab = xlab-xoffset
                            if self.projection in _pseudocyl:
                                if lat>0:
                                   t=ax.text(xlab,yy[n],latlab,horizontalalignment='right',verticalalignment='bottom',**text_kwargs)
                                elif lat<0:
                                   t=ax.text(xlab,yy[n],latlab,horizontalalignment='right',verticalalignment='top',**text_kwargs)
                                else:
                                   t=ax.text(xlab,yy[n],latlab,horizontalalignment='right',verticalalignment='center',**text_kwargs)
                            else:
                               t=ax.text(xlab,yy[n],latlab,horizontalalignment='right',verticalalignment='center',**text_kwargs)
                        elif side == 'r':
                            if self.projection in _pseudocyl:
                                if self.celestial:
                                   xlab,ylab = self(lon_0-179.9,lat)
                                else:
                                   xlab,ylab = self(lon_0+179.9,lat)
                            else:
                                xlab = self.urcrnrx
                            xlab = xlab+xoffset
                            if self.projection in _pseudocyl:
                                if lat>0:
                                   t=ax.text(xlab,yy[n],latlab,horizontalalignment='left',verticalalignment='bottom',**text_kwargs)
                                elif lat<0:
                                   t=ax.text(xlab,yy[n],latlab,horizontalalignment='left',verticalalignment='top',**text_kwargs)
                                else:
                                   t=ax.text(xlab,yy[n],latlab,horizontalalignment='left',verticalalignment='center',**text_kwargs)
                            else:
                               t=ax.text(xlab,yy[n],latlab,horizontalalignment='left',verticalalignment='center',**text_kwargs)
                        elif side == 'b':
                            t = ax.text(xx[n],self.llcrnry-yoffset,latlab,horizontalalignment='center',verticalalignment='top',**text_kwargs)
                        else:
                            t = ax.text(xx[n],self.urcrnry+yoffset,latlab,horizontalalignment='center',verticalalignment='bottom',**text_kwargs)
                        if t is not None: linecolls[lat][1].append(t)

        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        keys = list(linecolls.keys()); vals = list(linecolls.values())
        for k,v in zip(keys,vals):
            if v == ([], []):
                del linecolls[k]
            # add a remove method to each tuple.
            else:
                linecolls[k] = _tup(linecolls[k])
        # override __delitem__ in dict to call remove() on values.
        pardict = _dict(linecolls)
        # clip parallels for round polar plots (and delete labels).
        for lines, _ in pardict.values():
            self._cliplimb(ax, lines)
        return pardict

    def drawmeridians(self,meridians,color='k',textcolor='k',linewidth=1., zorder=None,\
                      dashes=[1,1],labels=[0,0,0,0],labelstyle=None,\
                      fmt='%g',xoffset=None,yoffset=None,ax=None,latmax=None,
                      **text_kwargs):
        """
        Draw and label meridians (longitude lines) for values (in degrees)
        given in the sequence ``meridians``.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keyword          Description
        ==============   ====================================================
        color            color to draw meridians (default black).
        textcolor        color to draw labels (default black).
        linewidth        line width for meridians (default 1.)
        zorder           sets the zorder for meridians (if not specified,
                         uses default zorder for matplotlib.lines.Line2D
                         objects).
        dashes           dash pattern for meridians (default [1,1], i.e.
                         1 pixel on, 1 pixel off).
        labels           list of 4 values (default [0,0,0,0]) that control
                         whether meridians are labelled where they intersect
                         the left, right, top or bottom of the plot. For
                         example labels=[1,0,0,1] will cause meridians
                         to be labelled where they intersect the left and
                         and bottom of the plot, but not the right and top.
        labelstyle       if set to "+/-", east and west longitudes are
                         labelled with "+" and "-", otherwise they are
                         labelled with "E" and "W".
        fmt              a format string to format the meridian labels
                         (default '%g') **or** a function that takes a
                         longitude value in degrees as it's only argument
                         and returns a formatted string.
        xoffset          label offset from edge of map in x-direction
                         (default is 0.01 times width of map in map
                         projection coordinates).
        yoffset          label offset from edge of map in y-direction
                         (default is 0.01 times height of map in map
                         projection coordinates).
        ax               axes instance (overrides default axes instance)
        latmax           absolute value of latitude to which meridians are drawn
                         (default is 80).
        \**text_kwargs   additional keyword arguments controlling text
                         for labels that are passed on to
                         the text method of the axes instance (see
                         matplotlib.pyplot.text documentation).
        ==============   ====================================================

        returns a dictionary whose keys are the meridian values, and
        whose values are tuples containing lists of the
        matplotlib.lines.Line2D and matplotlib.text.Text instances
        associated with each meridian. Deleting an item from the
        dictionary removes the correpsonding meridian from the plot.
        """
        text_kwargs['color']=textcolor # pass textcolor kwarg on to ax.text
        # for cylindrical projections, try to handle wraparound (i.e. if
        # projection is defined in -180 to 0 and user asks for meridians from
        # 180 to 360 to be drawn, it should work)
        if self.projection in _cylproj or self.projection in _pseudocyl:
            def addlon(meridians,madd):
                minside = (madd >= self.llcrnrlon and madd <= self.urcrnrlon)
                if minside and madd not in meridians: meridians.append(madd)
                return meridians
            merids = list(meridians)
            meridians = []
            for m in merids:
                meridians = addlon(meridians,m)
                meridians = addlon(meridians,m+360)
                meridians = addlon(meridians,m-360)
            meridians.sort()
        # if celestial=True, don't use "E" and "W" labels.
        if labelstyle is None and self.celestial:
            labelstyle="+/-"
        # get current axes instance (if none specified).
        ax = ax or self._check_ax()
        # don't draw meridians past latmax, always draw parallel at latmax.
        if latmax is None: latmax = 80. # unused w/ cyl, merc or miller proj.
        # offset for labels.
        if yoffset is None:
            yoffset = (self.urcrnry-self.llcrnry)/100.
            if self.aspect > 1:
                yoffset = self.aspect*yoffset
            else:
                yoffset = yoffset/self.aspect
        if xoffset is None:
            xoffset = (self.urcrnrx-self.llcrnrx)/100.

        lats = np.linspace(self.latmin,self.latmax,10001)
        if self.projection not in _cylproj + _pseudocyl:
            testlat = np.logical_and(lats>-latmax,lats<latmax)
            lats = np.compress(testlat,lats)

        xdelta = 0.01*(self.xmax-self.xmin)
        ydelta = 0.01*(self.ymax-self.ymin)
        linecolls = {}
        for merid in meridians:
            lons = merid*np.ones(len(lats),np.float32)
            x,y = self(lons,lats)
            # remove points outside domain.
            # leave a little slop around edges (3*xdelta)
            # don't really know why, but this appears to be needed to
            # or lines sometimes don't reach edge of plot.
            testx = np.logical_and(x>=self.xmin-3*xdelta,x<=self.xmax+3*xdelta)
            x = np.compress(testx, x)
            y = np.compress(testx, y)
            testy = np.logical_and(y>=self.ymin-3*ydelta,y<=self.ymax+3*ydelta)
            x = np.compress(testy, x)
            y = np.compress(testy, y)
            lines = []
            if len(x) > 1 and len(y) > 1:
                # split into separate line segments if necessary.
                # (not necessary for mercator or cylindrical or miller).
                xd = (x[1:]-x[0:-1])**2
                yd = (y[1:]-y[0:-1])**2
                dist = np.sqrt(xd+yd)
                if self.projection not in ['cyl','rotpole']:
                    split = dist > self.rmajor/10.
                else:
                    split = dist > 1.
                if np.sum(split) and self.projection not in _cylproj:
                    ind = (np.compress(split,np.squeeze(split*np.indices(xd.shape)))+1).tolist()
                    xl = []
                    yl = []
                    iprev = 0
                    ind.append(len(xd))
                    for i in ind:
                        xl.append(x[iprev:i])
                        yl.append(y[iprev:i])
                        iprev = i
                else:
                    xl = [x]
                    yl = [y]
                # draw each line segment.
                for x,y in zip(xl,yl):
                    # skip if only a point.
                    if len(x) > 1 and len(y) > 1:
                        l = Line2D(x,y,linewidth=linewidth)
                        l.set_color(color)
                        l.set_dashes(dashes)
                        l.set_label('_nolabel_')
                        if zorder is not None:
                            l.set_zorder(zorder)
                        ax.add_line(l)
                        lines.append(l)
            linecolls[merid] = (lines,[])
        # draw labels for meridians.
        # meridians not labelled for sinusoidal, hammer, mollweide,
        # VanDerGrinten or full-disk orthographic/geostationary.
        if self.projection in ['sinu','moll','hammer','vandg'] and max(labels):
            sys.stdout.write('Warning: Cannot label meridians on %s basemap' % _projnames[self.projection])
            labels = [0,0,0,0]
        if self.projection in ['ortho','geos','nsper','aeqd'] and max(labels):
            if self._fulldisk and self.boundinglat is None:
                sys.stdout.write(dedent(
                """'Warning: Cannot label meridians on full-disk
                Geostationary, Orthographic or Azimuthal equidistant basemap
                """))
                labels = [0,0,0,0]
        # search along edges of map to see if parallels intersect.
        # if so, find x,y location of intersection and draw a label there.
        dx = (self.xmax-self.xmin)/1000.
        dy = (self.ymax-self.ymin)/1000.
        if self.projection in _pseudocyl:
            lon_0 = self.projparams['lon_0']
            xmin,ymin = self(lon_0-179.9,-90)
            xmax,ymax = self(lon_0+179.9,90)
        for dolab,side in zip(labels,['l','r','t','b']):
            if not dolab or self.round: continue
            # for cylindrical projections, don't draw meridians on left or right.
            if self.projection in _cylproj + _pseudocyl and side in ['l','r']: continue
            if side in ['l','r']:
                nmax = int((self.ymax-self.ymin)/dy+1)
                yy = np.linspace(self.llcrnry,self.urcrnry,nmax)
                if side == 'l':
                    lons,lats = self(self.llcrnrx*np.ones(yy.shape,np.float32),yy,inverse=True)
                    lons = lons.tolist(); lats = lats.tolist()
                else:
                    lons,lats = self(self.urcrnrx*np.ones(yy.shape,np.float32),yy,inverse=True)
                    lons = lons.tolist(); lats = lats.tolist()
                if max(lons) > 1.e20 or max(lats) > 1.e20:
                    raise ValueError('inverse transformation undefined - please adjust the map projection region')
                # adjust so 0 <= lons < 360
                lons = [(lon+360) % 360 for lon in lons]
            else:
                nmax = int((self.xmax-self.xmin)/dx+1)
                if self.projection in _pseudocyl:
                    xx = np.linspace(xmin,xmax,nmax)
                else:
                    xx = np.linspace(self.llcrnrx,self.urcrnrx,nmax)
                if side == 'b':
                    lons,lats = self(xx,self.llcrnry*np.ones(xx.shape,np.float32),inverse=True)
                    lons = lons.tolist(); lats = lats.tolist()
                else:
                    lons,lats = self(xx,self.urcrnry*np.ones(xx.shape,np.float32),inverse=True)
                    lons = lons.tolist(); lats = lats.tolist()
                if max(lons) > 1.e20 or max(lats) > 1.e20:
                    raise ValueError('inverse transformation undefined - please adjust the map projection region')
                # adjust so 0 <= lons < 360
                lons = [(lon+360) % 360 for lon in lons]
            for lon in meridians:
                # adjust so 0 <= lon < 360
                lon2 = (lon+360) % 360
                # find index of meridian (there may be two, so
                # search from left and right).
                nl = _searchlist(lons,lon2)
                nr = _searchlist(lons[::-1],lon2)
                if nr != -1: nr = len(lons)-nr-1
                lonlab = _setlonlab(fmt,lon2,labelstyle)
                # meridians can intersect each map edge twice.
                for i,n in enumerate([nl,nr]):
                    lat = lats[n]/100.
                    # no meridians > latmax for projections other than merc,cyl,miller.
                    if self.projection not in _cylproj and lat > latmax: continue
                    # don't bother if close to the first label.
                    if i and abs(nr-nl) < 100: continue
                    if n >= 0:
                        t = None
                        if side == 'l':
                            t = ax.text(self.llcrnrx-xoffset,yy[n],lonlab,horizontalalignment='right',verticalalignment='center',**text_kwargs)
                        elif side == 'r':
                            t = ax.text(self.urcrnrx+xoffset,yy[n],lonlab,horizontalalignment='left',verticalalignment='center',**text_kwargs)
                        elif side == 'b':
                            t = ax.text(xx[n],self.llcrnry-yoffset,lonlab,horizontalalignment='center',verticalalignment='top',**text_kwargs)
                        else:
                            t = ax.text(xx[n],self.urcrnry+yoffset,lonlab,horizontalalignment='center',verticalalignment='bottom',**text_kwargs)

                        if t is not None: linecolls[lon][1].append(t)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # remove empty values from linecolls dictionary
        keys = list(linecolls.keys()); vals = list(linecolls.values())
        for k,v in zip(keys,vals):
            if v == ([], []):
                del linecolls[k]
            else:
            # add a remove method to each tuple.
                linecolls[k] = _tup(linecolls[k])
        # override __delitem__ in dict to call remove() on values.
        meridict = _dict(linecolls)
        # for round polar plots, clip meridian lines and label them.
        if self.round:
            # label desired?
            label = False
            for lab in labels:
                if lab: label = True
            for merid in meridict:
                if not label: continue
                # label
                lonlab = _setlonlab(fmt,merid,labelstyle)
                x,y = self(merid,self.boundinglat)
                r = np.sqrt((x-0.5*(self.xmin+self.xmax))**2+
                            (y-0.5*(self.ymin+self.ymax))**2)
                r = r + np.sqrt(xoffset**2+yoffset**2)
                if self.projection.startswith('np'):
                    pole = 1
                elif self.projection.startswith('sp'):
                    pole = -1
                elif self.projection == 'ortho' and self.round:
                    pole = 1
                if pole == 1:
                    theta = (np.pi/180.)*(merid-self.projparams['lon_0']-90)
                    if self.projection == 'ortho' and\
                       self.projparams['lat_0'] == -90:
                        theta = (np.pi/180.)*(-merid+self.projparams['lon_0']+90)
                    x = r*np.cos(theta)+0.5*(self.xmin+self.xmax)
                    y = r*np.sin(theta)+0.5*(self.ymin+self.ymax)
                    if x > 0.5*(self.xmin+self.xmax)+xoffset:
                        horizalign = 'left'
                    elif x < 0.5*(self.xmin+self.xmax)-xoffset:
                        horizalign = 'right'
                    else:
                        horizalign = 'center'
                    if y > 0.5*(self.ymin+self.ymax)+yoffset:
                        vertalign = 'bottom'
                    elif y < 0.5*(self.ymin+self.ymax)-yoffset:
                        vertalign = 'top'
                    else:
                        vertalign = 'center'
                    # labels [l,r,t,b]
                    if labels[0] and not labels[1] and x >= 0.5*(self.xmin+self.xmax)+xoffset: continue
                    if labels[1] and not labels[0] and x <= 0.5*(self.xmin+self.xmax)-xoffset: continue
                    if labels[2] and not labels[3] and y <= 0.5*(self.ymin+self.ymax)-yoffset: continue
                    if labels[3] and not labels[2]and y >= 0.5*(self.ymin+self.ymax)+yoffset: continue
                elif pole == -1:
                    theta = (np.pi/180.)*(-merid+self.projparams['lon_0']+90)
                    x = r*np.cos(theta)+0.5*(self.xmin+self.xmax)
                    y = r*np.sin(theta)+0.5*(self.ymin+self.ymax)
                    if x > 0.5*(self.xmin+self.xmax)-xoffset:
                        horizalign = 'right'
                    elif x < 0.5*(self.xmin+self.xmax)+xoffset:
                        horizalign = 'left'
                    else:
                        horizalign = 'center'
                    if y > 0.5*(self.ymin+self.ymax)-yoffset:
                        vertalign = 'top'
                    elif y < 0.5*(self.ymin+self.ymax)+yoffset:
                        vertalign = 'bottom'
                    else:
                        vertalign = 'center'
                    # labels [l,r,t,b]
                    if labels[0] and not labels[1] and x <=  0.5*(self.xmin+self.xmax)+xoffset: continue
                    if labels[1] and not labels[0] and x >=  0.5*(self.xmin+self.xmax)-xoffset: continue
                    if labels[2] and not labels[3] and y >=  0.5*(self.ymin+self.ymax)-yoffset: continue
                    if labels[3] and not labels[2] and y <=  0.5*(self.ymin+self.ymax)+yoffset: continue
                t=ax.text(x,y,lonlab,horizontalalignment=horizalign,verticalalignment=vertalign,**text_kwargs)
                meridict[merid][1].append(t)
        for lines, _ in meridict.values():
            self._cliplimb(ax, lines)
        return meridict

    def tissot(self,lon_0,lat_0,radius_deg,npts,ax=None,**kwargs):
        """
        Draw a polygon centered at ``lon_0,lat_0``.  The polygon
        approximates a circle on the surface of the earth with radius
        ``radius_deg`` degrees latitude along longitude ``lon_0``,
        made up of ``npts`` vertices.
        The polygon represents a Tissot's indicatrix
        (http://en.wikipedia.org/wiki/Tissot's_Indicatrix),
        which when drawn on a map shows the distortion
        inherent in the map projection.

        .. note::
         Cannot handle situations in which the polygon intersects
         the edge of the map projection domain, and then re-enters the domain.

        Extra keyword ``ax`` can be used to override the default axis instance.

        Other \**kwargs passed on to matplotlib.patches.Polygon.

        returns a matplotlib.patches.Polygon object."""
        ax = kwargs.pop('ax', None) or self._check_ax()
        g = pyproj.Geod(a=self.rmajor,b=self.rminor)
        az12,az21,dist = g.inv(lon_0,lat_0,lon_0,lat_0+radius_deg)
        seg = [self(lon_0,lat_0+radius_deg)]
        delaz = 360./npts
        az = az12
        for n in range(npts):
            az = az+delaz
            lon, lat, az21 = g.fwd(lon_0, lat_0, az, dist)
            x,y = self(lon,lat)
            # add segment if it is in the map projection region.
            if x < 1.e20 and y < 1.e20:
                seg.append((x,y))
        poly = Polygon(seg,**kwargs)
        ax.add_patch(poly)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip polygons to map limbs
        poly,c = self._cliplimb(ax,poly)
        return poly

    def gcpoints(self,lon1,lat1,lon2,lat2,npoints):
        """
        compute ``points`` points along a great circle with endpoints
        ``(lon1,lat1)`` and ``(lon2,lat2)``.

        Returns arrays x,y with map projection coordinates.
        """
        gc = pyproj.Geod(a=self.rmajor,b=self.rminor)
        lonlats = gc.npts(lon1,lat1,lon2,lat2,npoints-2)
        lons=[lon1];lats=[lat1]
        for lon,lat in lonlats:
            lons.append(lon); lats.append(lat)
        lons.append(lon2); lats.append(lat2)
        x, y = self(lons, lats)
        return x,y

    def drawgreatcircle(self,lon1,lat1,lon2,lat2,del_s=100.,**kwargs):
        """
        Draw a great circle on the map from the longitude-latitude
        pair ``lon1,lat1`` to ``lon2,lat2``

        .. tabularcolumns:: |l|L|

        ==============   =======================================================
        Keyword          Description
        ==============   =======================================================
        del_s            points on great circle computed every del_s kilometers
                         (default 100).
        \**kwargs        other keyword arguments are passed on to :meth:`plot`
                         method of Basemap instance.
        ==============   =======================================================

        Returns a list with a single ``matplotlib.lines.Line2D`` object like a
        call to ``pyplot.plot()``.
        """
        # use great circle formula for a perfect sphere.
        gc = pyproj.Geod(a=self.rmajor,b=self.rminor)
        az12,az21,dist = gc.inv(lon1,lat1,lon2,lat2)
        npoints = int((dist+0.5*1000.*del_s)/(1000.*del_s))
        lonlats = gc.npts(lon1,lat1,lon2,lat2,npoints)
        lons = [lon1]; lats = [lat1]
        for lon, lat in lonlats:
            lons.append(lon)
            lats.append(lat)
        lons.append(lon2); lats.append(lat2)
        x, y = self(lons, lats)

        # Correct wrap around effect of great circles

        # get points
        _p = self.plot(x,y,**kwargs)
        p = _p[0].get_path()

        # since we know the difference between any two points, we can use this to find wrap arounds on the plot
        max_dist = 1000*del_s*2

        # calculate distances and compare with max allowable distance
        dists = np.abs(np.diff(p.vertices[:,0]))
        cuts = np.where( dists > max_dist )[0]

        # if there are any cut points, cut them and begin again at the next point
        for i,k in enumerate(cuts):
            # vertex to cut at
            cut_point = cuts[i]

            # create new vertices with a nan inbetween and set those as the path's vertices
            verts = np.concatenate(
                                       [p.vertices[:cut_point, :],
                                        [[np.nan, np.nan]],
                                        p.vertices[cut_point+1:, :]]
                                       )
            p.codes = None
            p.vertices = verts

        return _p

    def transform_scalar(self,datin,lons,lats,nx,ny,returnxy=False,checkbounds=False,order=1,masked=False):
        """
        Interpolate a scalar field (``datin``) from a lat/lon grid with
        longitudes = ``lons`` and latitudes = ``lats`` to a ``ny`` by ``nx``
        map projection grid.  Typically used to transform data to
        map projection coordinates for plotting on a map with
        the :meth:`imshow`.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Argument         Description
        ==============   ====================================================
        datin            input data on a lat/lon grid.
        lons, lats       rank-1 arrays containing longitudes and latitudes
                         (in degrees) of input data in increasing order.
                         For non-cylindrical projections (those other than
                         ``cyl``, ``merc``, ``cea``, ``gall`` and ``mill``) lons
                         must fit within range -180 to 180.
        nx, ny           The size of the output regular grid in map
                         projection coordinates
        ==============   ====================================================

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keyword          Description
        ==============   ====================================================
        returnxy         If True, the x and y values of the map
                         projection grid are also returned (Default False).
        checkbounds      If True, values of lons and lats are checked to see
                         that they lie within the map projection region.
                         Default is False, and data outside map projection
                         region is clipped to values on boundary.
        masked           If True, interpolated data is returned as a masked
                         array with values outside map projection region
                         masked (Default False).
        order            0 for nearest-neighbor interpolation, 1 for
                         bilinear, 3 for cubic spline (Default 1).
                         Cubic spline interpolation requires scipy.ndimage.
        ==============   ====================================================

        Returns ``datout`` (data on map projection grid).
        If returnxy=True, returns ``data,x,y``.
        """
        # check that lons, lats increasing
        delon = lons[1:]-lons[0:-1]
        delat = lats[1:]-lats[0:-1]
        if min(delon) < 0. or min(delat) < 0.:
            raise ValueError('lons and lats must be increasing!')
        # check that lons in -180,180 for non-cylindrical projections.
        if self.projection not in _cylproj:
            lonsa = np.array(lons)
            count = np.sum(lonsa < -180.00001) + np.sum(lonsa > 180.00001)
            if count > 1:
                raise ValueError('grid must be shifted so that lons are monotonically increasing and fit in range -180,+180 (see shiftgrid function)')
            # allow for wraparound point to be outside.
            elif count == 1 and math.fabs(lons[-1]-lons[0]-360.) > 1.e-4:
                raise ValueError('grid must be shifted so that lons are monotonically increasing and fit in range -180,+180 (see shiftgrid function)')
        if returnxy:
            lonsout, latsout, x, y = self.makegrid(nx,ny,returnxy=True)
        else:
            lonsout, latsout = self.makegrid(nx,ny)
        datout = interp(datin,lons,lats,lonsout,latsout,checkbounds=checkbounds,order=order,masked=masked)
        if returnxy:
            return datout, x, y
        else:
            return datout

    def transform_vector(self,uin,vin,lons,lats,nx,ny,returnxy=False,checkbounds=False,order=1,masked=False):
        """
        Rotate and interpolate a vector field (``uin,vin``) from a
        lat/lon grid with longitudes = ``lons`` and latitudes = ``lats``
        to a ``ny`` by ``nx`` map projection grid.

        The input vector field is defined in spherical coordinates (it
        has eastward and northward components) while the output
        vector field is rotated to map projection coordinates (relative
        to x and y). The magnitude of the vector is preserved.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Arguments        Description
        ==============   ====================================================
        uin, vin         input vector field on a lat/lon grid.
        lons, lats       rank-1 arrays containing longitudes and latitudes
                         (in degrees) of input data in increasing order.
                         For non-cylindrical projections (those other than
                         ``cyl``, ``merc``, ``cea``, ``gall`` and ``mill``) lons
                         must fit within range -180 to 180.
        nx, ny           The size of the output regular grid in map
                         projection coordinates
        ==============   ====================================================

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keyword          Description
        ==============   ====================================================
        returnxy         If True, the x and y values of the map
                         projection grid are also returned (Default False).
        checkbounds      If True, values of lons and lats are checked to see
                         that they lie within the map projection region.
                         Default is False, and data outside map projection
                         region is clipped to values on boundary.
        masked           If True, interpolated data is returned as a masked
                         array with values outside map projection region
                         masked (Default False).
        order            0 for nearest-neighbor interpolation, 1 for
                         bilinear, 3 for cubic spline (Default 1).
                         Cubic spline interpolation requires scipy.ndimage.
        ==============   ====================================================

        Returns ``uout, vout`` (vector field on map projection grid).
        If returnxy=True, returns ``uout,vout,x,y``.
        """
        # check that lons, lats increasing
        delon = lons[1:]-lons[0:-1]
        delat = lats[1:]-lats[0:-1]
        if min(delon) < 0. or min(delat) < 0.:
            raise ValueError('lons and lats must be increasing!')
        # check that lons in -180,180 for non-cylindrical projections.
        if self.projection not in _cylproj:
            lonsa = np.array(lons)
            count = np.sum(lonsa < -180.00001) + np.sum(lonsa > 180.00001)
            if count > 1:
                raise ValueError('grid must be shifted so that lons are monotonically increasing and fit in range -180,+180 (see shiftgrid function)')
            # allow for wraparound point to be outside.
            elif count == 1 and math.fabs(lons[-1]-lons[0]-360.) > 1.e-4:
                raise ValueError('grid must be shifted so that lons are monotonically increasing and fit in range -180,+180 (see shiftgrid function)')
        lonsout, latsout, x, y = self.makegrid(nx,ny,returnxy=True)
        # interpolate to map projection coordinates.
        uin = interp(uin,lons,lats,lonsout,latsout,checkbounds=checkbounds,order=order,masked=masked)
        vin = interp(vin,lons,lats,lonsout,latsout,checkbounds=checkbounds,order=order,masked=masked)
        # rotate from geographic to map coordinates.
        return self.rotate_vector(uin,vin,lonsout,latsout,returnxy=returnxy)

    def rotate_vector(self,uin,vin,lons,lats,returnxy=False):
        """
        Rotate a vector field (``uin,vin``) on a rectilinear grid
        with longitudes = ``lons`` and latitudes = ``lats`` from
        geographical (lat/lon) into map projection (x/y) coordinates.

        Differs from transform_vector in that no interpolation is done.
        The vector is returned on the same grid, but rotated into
        x,y coordinates.

        The input vector field is defined in spherical coordinates (it
        has eastward and northward components) while the output
        vector field is rotated to map projection coordinates (relative
        to x and y). The magnitude of the vector is preserved.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Arguments        Description
        ==============   ====================================================
        uin, vin         input vector field on a lat/lon grid.
        lons, lats       Arrays containing longitudes and latitudes
                         (in degrees) of input data in increasing order.
                         For non-cylindrical projections (those other than
                         ``cyl``, ``merc``, ``cyl``, ``gall`` and ``mill``) lons
                         must fit within range -180 to 180.
        ==============   ====================================================

        Returns ``uout, vout`` (rotated vector field).
        If the optional keyword argument
        ``returnxy`` is True (default is False),
        returns ``uout,vout,x,y`` (where ``x,y`` are the map projection
        coordinates of the grid defined by ``lons,lats``).
        """
        # if lons,lats are 1d and uin,vin are 2d, and
        # lats describes 1st dim of uin,vin, and
        # lons describes 2nd dim of uin,vin, make lons,lats 2d
        # with meshgrid.
        if lons.ndim == lats.ndim == 1 and uin.ndim == vin.ndim == 2 and\
           uin.shape[1] == vin.shape[1] == lons.shape[0] and\
           uin.shape[0] == vin.shape[0] == lats.shape[0]:
            lons, lats = np.meshgrid(lons, lats)
        else:
            if not lons.shape == lats.shape == uin.shape == vin.shape:
                raise TypeError("shapes of lons,lats and uin,vin don't match")
        x, y = self(lons, lats)
        # rotate from geographic to map coordinates.
        if ma.isMaskedArray(uin):
            mask = ma.getmaskarray(uin)
            masked = True
            uin = uin.filled(1)
            vin = vin.filled(1)
        else:
            masked = False

        # Map the (lon, lat) vector in the complex plane.
        uvc = uin + 1j*vin
        uvmag = np.abs(uvc)
        theta = np.angle(uvc)

        # Define a displacement (dlon, dlat) that moves all
        # positions (lons, lats) a small distance in the
        # direction of the original vector.
        dc = 1E-5 * np.exp(theta*1j)
        dlat = dc.imag * np.cos(np.radians(lats))
        dlon = dc.real

        # Deal with displacements that overshoot the North or South Pole.
        farnorth = np.abs(lats+dlat) >= 90.0
        somenorth = farnorth.any()
        if somenorth:
            dlon[farnorth] *= -1.0
            dlat[farnorth] *= -1.0

        # Add displacement to original location and find the native coordinates.
        lon1 = lons + dlon
        lat1 = lats + dlat
        xn, yn = self(lon1, lat1)

        # Determine the angle of the displacement in the native coordinates.
        vecangle = np.arctan2(yn-y, xn-x)
        if somenorth:
            vecangle[farnorth] += np.pi

        # Compute the x-y components of the original vector.
        uvcout = uvmag * np.exp(1j*vecangle)
        uout = uvcout.real
        vout = uvcout.imag

        if masked:
            uout = ma.array(uout, mask=mask)
            vout = ma.array(vout, mask=mask)
        if returnxy:
            return uout,vout,x,y
        else:
            return uout,vout

    def set_axes_limits(self,ax=None):
        """
        Final step in Basemap method wrappers of Axes plotting methods:

        Set axis limits, fix aspect ratio for map domain using current
        or specified axes instance.  This is done only once per axes
        instance.

        In interactive mode, this method always calls draw_if_interactive
        before returning.

        """
        # get current axes instance (if none specified).
        ax = ax or self._check_ax()

        # If we have already set the axes limits, and if the user
        # has not defeated this by turning autoscaling back on,
        # then all we need to do is plot if interactive.
        if (hash(ax) in self._initialized_axes
                                 and not ax.get_autoscalex_on()
                                 and not ax.get_autoscaley_on()):
            if is_interactive():
                import matplotlib.pyplot as plt
                plt.draw_if_interactive()
            return

        self._initialized_axes.add(hash(ax))
        # Take control of axis scaling:
        ax.set_autoscale_on(False)
        # update data limits for map domain.
        corners = ((self.llcrnrx, self.llcrnry), (self.urcrnrx, self.urcrnry))
        ax.update_datalim(corners)
        ax.set_xlim((self.llcrnrx, self.urcrnrx))
        ax.set_ylim((self.llcrnry, self.urcrnry))
        # if map boundary not yet drawn for elliptical maps, draw it with default values.
        if not self._mapboundarydrawn or self._mapboundarydrawn not in ax.patches:
            # elliptical map, draw boundary manually.
            if ((self.projection in ['ortho', 'geos', 'nsper', 'aeqd'] and
                    self._fulldisk) or self.round or
                    self.projection in _pseudocyl):
                # first draw boundary, no fill
                limb1 = self.drawmapboundary(fill_color='none', ax=ax)
                # draw another filled patch, with no boundary.
                limb2 = self.drawmapboundary(linewidth=0, ax=ax)
                self._mapboundarydrawn = limb2
        # for elliptical map, always turn off axis_frame.
        if ((self.projection in ['ortho', 'geos', 'nsper', 'aeqd'] and
                self._fulldisk) or self.round or
                self.projection in _pseudocyl):
            # turn off axes frame.
            ax.set_frame_on(False)
        # make sure aspect ratio of map preserved.
        # plot is re-centered in bounding rectangle.
        # (anchor instance var determines where plot is placed)
        if self.fix_aspect:
            ax.set_aspect('equal',anchor=self.anchor)
        else:
            ax.set_aspect('auto',anchor=self.anchor)
        # make sure axis ticks are turned off.
        if self.noticks:
            ax.set_xticks([])
            ax.set_yticks([])
        # force draw if in interactive mode.
        if is_interactive():
            import matplotlib.pyplot as plt
            plt.draw_if_interactive()


    def _save_use_hold(self, ax, kwargs):
        h = kwargs.pop('hold', None)
        if hasattr(ax, '_hold'):
            self._tmp_hold = ax._hold
            if h is not None:
                ax._hold = h

    def _restore_hold(self, ax):
        if hasattr(ax, '_hold'):
            ax._hold = self._tmp_hold

    @_transform1d
    def scatter(self, *args, **kwargs):
        """
        Plot points with markers on the map
        (see matplotlib.pyplot.scatter documentation).

        If ``latlon`` keyword is set to True, x,y are intrepreted as
        longitude and latitude in degrees.  Data and longitudes are
        automatically shifted to match map projection region for cylindrical
        and pseudocylindrical projections, and x,y are transformed to map
        projection coordinates. If ``latlon`` is False (default), x and y
        are assumed to be map projection coordinates.

        Extra keyword ``ax`` can be used to override the default axes instance.

        Other \**kwargs passed on to matplotlib.pyplot.scatter.
        """
        ax, plt = self._ax_plt_from_kw(kwargs)
        self._save_use_hold(ax, kwargs)
        try:
            ret =  ax.scatter(*args, **kwargs)
        finally:
            self._restore_hold(ax)
        # reset current active image (only if pyplot is imported).
        if plt:
            plt.sci(ret)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip to map limbs
        ret,c = self._cliplimb(ax,ret)
        return ret

    @_transform1d
    def plot(self, *args, **kwargs):
        """
        Draw lines and/or markers on the map
        (see matplotlib.pyplot.plot documentation).

        If ``latlon`` keyword is set to True, x,y are intrepreted as
        longitude and latitude in degrees.  Data and longitudes are
        automatically shifted to match map projection region for cylindrical
        and pseudocylindrical projections, and x,y are transformed to map
        projection coordinates. If ``latlon`` is False (default), x and y
        are assumed to be map projection coordinates.

        Extra keyword ``ax`` can be used to override the default axis instance.

        Other \**kwargs passed on to matplotlib.pyplot.plot.
        """
        ax = kwargs.pop('ax', None) or self._check_ax()
        self._save_use_hold(ax, kwargs)
        try:
            ret =  ax.plot(*args, **kwargs)
        finally:
            self._restore_hold(ax)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip to map limbs
        ret,c = self._cliplimb(ax,ret)
        return ret

    def imshow(self, *args, **kwargs):
        """
        Display an image over the map
        (see matplotlib.pyplot.imshow documentation).

        ``extent`` and ``origin`` keywords set automatically so image
        will be drawn over map region.

        Extra keyword ``ax`` can be used to override the default axis instance.

        Other \**kwargs passed on to matplotlib.pyplot.plot.

        returns an matplotlib.image.AxesImage instance.
        """
        ax, plt = self._ax_plt_from_kw(kwargs)
        kwargs['extent']=(self.llcrnrx,self.urcrnrx,self.llcrnry,self.urcrnry)
        # use origin='lower', unless overridden.
        if 'origin' not in kwargs:
            kwargs['origin']='lower'
        self._save_use_hold(ax, kwargs)
        try:
            ret =  ax.imshow(*args, **kwargs)
        finally:
            self._restore_hold(ax)
        # reset current active image (only if pyplot is imported).
        if plt:
            plt.sci(ret)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip image to map limbs
        ret,c = self._cliplimb(ax,ret)
        return ret

    @_transform
    def pcolor(self,x,y,data,**kwargs):
        """
        Make a pseudo-color plot over the map
        (see matplotlib.pyplot.pcolor documentation).

        If ``latlon`` keyword is set to True, x,y are intrepreted as
        longitude and latitude in degrees.  Data and longitudes are
        automatically shifted to match map projection region for cylindrical
        and pseudocylindrical projections, and x,y are transformed to map
        projection coordinates. If ``latlon`` is False (default), x and y
        are assumed to be map projection coordinates.

        If x or y are outside projection limb (i.e. they have values > 1.e20)
        they will be convert to masked arrays with those values masked.
        As a result, those values will not be plotted.

        If ``tri`` is set to ``True``, an unstructured grid is assumed
        (x,y,data must be 1-d) and matplotlib.pyplot.tripcolor is used.

        Extra keyword ``ax`` can be used to override the default axis instance.

        Other \**kwargs passed on to matplotlib.pyplot.pcolor (or tripcolor if
        ``tri=True``).

        Note: (taken from matplotlib.pyplot.pcolor documentation)
        Ideally the dimensions of x and y should be one greater than those of data;
        if the dimensions are the same, then the last row and column of data will be ignored.
        """
        ax, plt = self._ax_plt_from_kw(kwargs)
        self._save_use_hold(ax, kwargs)
        try:
            if kwargs.pop('tri', False):
                try:
                    import matplotlib.tri as tri
                except:
                    msg='need matplotlib > 0.99.1 to plot on unstructured grids'
                    raise ImportError(msg)
                # for unstructured grids, toss out points outside
                # projection limb (don't use those points in triangulation).
                if ma.isMA(data):
                    data = data.filled(fill_value=1.e30)
                    masked=True
                else:
                    masked=False
                mask = np.logical_or(x<1.e20,y<1.e20)
                x = np.compress(mask,x)
                y = np.compress(mask,y)
                data = np.compress(mask,data)
                if masked:
                    triang = tri.Triangulation(x, y)
                    z = data[triang.triangles]
                    mask = (z > 1.e20).sum(axis=-1)
                    triang.set_mask(mask)
                    ret = ax.tripcolor(triang,data,**kwargs)
                else:
                    ret = ax.tripcolor(x,y,data,**kwargs)
            else:
                # make x,y masked arrays
                # (masked where data is outside of projection limb)
                x = ma.masked_values(np.where(x > 1.e20,1.e20,x), 1.e20)
                y = ma.masked_values(np.where(y > 1.e20,1.e20,y), 1.e20)
                ret = ax.pcolor(x,y,data,**kwargs)
        finally:
            self._restore_hold(ax)
        # reset current active image (only if pyplot is imported).
        if plt:
            plt.sci(ret)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip to map limbs
        ret,c = self._cliplimb(ax,ret)
        if self.round:
            # for some reason, frame gets turned on.
            ax.set_frame_on(False)
        return ret

    @_transform
    def pcolormesh(self,x,y,data,**kwargs):
        """
        Make a pseudo-color plot over the map
        (see matplotlib.pyplot.pcolormesh documentation).

        If ``latlon`` keyword is set to True, x,y are intrepreted as
        longitude and latitude in degrees.  Data and longitudes are
        automatically shifted to match map projection region for cylindrical
        and pseudocylindrical projections, and x,y are transformed to map
        projection coordinates. If ``latlon`` is False (default), x and y
        are assumed to be map projection coordinates.

        Extra keyword ``ax`` can be used to override the default axis instance.

        Other \**kwargs passed on to matplotlib.pyplot.pcolormesh.

        Note: (taken from matplotlib.pyplot.pcolor documentation)
        Ideally the dimensions of x and y should be one greater than those of data;
        if the dimensions are the same, then the last row and column of data will be ignored.
        """
        ax, plt = self._ax_plt_from_kw(kwargs)
        # fix for invalid grid points
        if ((np.any(x > 1e20) or np.any(y > 1e20)) and
            x.ndim == 2 and y.ndim == 2):
            if x.shape != y.shape:
                raise ValueError('pcolormesh: x and y need same dimension')
            nx,ny = x.shape
            if nx < data.shape[0] or ny < data.shape[1]:
                raise ValueError('pcolormesh: data dimension needs to be at least that of x and y.')
            mask = (
                (x[:-1,:-1] > 1e20) |
                (x[1:,:-1] > 1e20) |
                (x[:-1,1:] > 1e20) |
                (x[1:,1:] > 1e20) |
                (y[:-1,:-1] > 1e20) |
                (y[1:,:-1] > 1e20) |
                (y[:-1,1:] > 1e20) |
                (y[1:,1:] > 1e20)
                )
            # we do not want to overwrite original array
            data = data[:nx-1,:ny-1].copy()
            data[mask] = np.nan
        self._save_use_hold(ax, kwargs)
        try:
            ret =  ax.pcolormesh(x,y,data,**kwargs)
        finally:
            self._restore_hold(ax)
        # reset current active image (only if pyplot is imported).
        if plt:
            plt.sci(ret)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip to map limbs
        ret,c = self._cliplimb(ax,ret)
        if self.round:
            # for some reason, frame gets turned on.
            ax.set_frame_on(False)
        return ret

    def hexbin(self,x,y,**kwargs):
        """
        Make a hexagonal binning plot of x versus y, where x, y are 1-D
        sequences of the same length, N. If C is None (the default), this is a
        histogram of the number of occurences of the observations at
        (x[i],y[i]).

        If C is specified, it specifies values at the coordinate (x[i],y[i]).
        These values are accumulated for each hexagonal bin and then reduced
        according to reduce_C_function, which defaults to the numpy mean function
        (np.mean). (If C is specified, it must also be a 1-D sequence of the
        same length as x and y.)

        x, y and/or C may be masked arrays, in which case only unmasked points
        will be plotted.

        (see matplotlib.pyplot.hexbin documentation).

        Extra keyword ``ax`` can be used to override the default axis instance.

        Other \**kwargs passed on to matplotlib.pyplot.hexbin
        """
        ax, plt = self._ax_plt_from_kw(kwargs)
        self._save_use_hold(ax, kwargs)
        try:
            # make x,y masked arrays
            # (masked where data is outside of projection limb)
            x = ma.masked_values(np.where(x > 1.e20,1.e20,x), 1.e20)
            y = ma.masked_values(np.where(y > 1.e20,1.e20,y), 1.e20)
            ret = ax.hexbin(x,y,**kwargs)
        finally:
            self._restore_hold(ax)
        # reset current active image (only if pyplot is imported).
        if plt:
            plt.sci(ret)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip to map limbs
        ret,c = self._cliplimb(ax,ret)
        return ret

    @_transform
    def contour(self,x,y,data,*args,**kwargs):
        """
        Make a contour plot over the map
        (see matplotlib.pyplot.contour documentation).

        If ``latlon`` keyword is set to True, x,y are intrepreted as
        longitude and latitude in degrees.  Data and longitudes are
        automatically shifted to match map projection region for cylindrical
        and pseudocylindrical projections, and x,y are transformed to map
        projection coordinates. If ``latlon`` is False (default), x and y
        are assumed to be map projection coordinates.

        Extra keyword ``ax`` can be used to override the default axis instance.

        If ``tri`` is set to ``True``, an unstructured grid is assumed
        (x,y,data must be 1-d) and matplotlib.pyplot.tricontour is used.

        Other \*args and \**kwargs passed on to matplotlib.pyplot.contour
        (or tricontour if ``tri=True``).
        """
        ax, plt = self._ax_plt_from_kw(kwargs)
        self._save_use_hold(ax, kwargs)
        try:
            if kwargs.pop('tri', False):
                try:
                    import matplotlib.tri as tri
                except:
                    msg='need matplotlib > 0.99.1 to plot on unstructured grids'
                    raise ImportError(msg)
                # for unstructured grids, toss out points outside
                # projection limb (don't use those points in triangulation).
                if ma.isMA(data):
                    data = data.filled(fill_value=1.e30)
                    masked=True
                else:
                    masked=False
                mask = np.logical_or(x<self.xmin,y<self.xmin) +\
                       np.logical_or(x>self.xmax,y>self.xmax)
                x = np.compress(mask,x)
                y = np.compress(mask,y)
                data = np.compress(mask,data)
                if masked:
                    triang = tri.Triangulation(x, y)
                    z = data[triang.triangles]
                    mask = (z > 1.e20).sum(axis=-1)
                    triang.set_mask(mask)
                    CS = ax.tricontour(triang,data,*args,**kwargs)
                else:
                    CS = ax.tricontour(x,y,data,*args,**kwargs)
            else:
                # make sure x is monotonically increasing - if not,
                # print warning suggesting that the data be shifted in longitude
                # with the shiftgrid function.
                # only do this check for global projections.
                if self.projection in _cylproj + _pseudocyl:
                    xx = x[x.shape[0]//2,:]
                    condition = (xx >= self.xmin) & (xx <= self.xmax)
                    xl = xx.compress(condition).tolist()
                    xs = xl[:]
                    xs.sort()
                    if xl != xs:
                        sys.stdout.write(dedent("""
                             WARNING: x coordinate not montonically increasing - contour plot
                             may not be what you expect.  If it looks odd, your can either
                             adjust the map projection region to be consistent with your data, or
                             (if your data is on a global lat/lon grid) use the shiftdata
                             method to adjust the data to be consistent with the map projection
                             region (see examples/shiftdata.py)."""))
                # mask for points more than one grid length outside projection limb.
                xx = ma.masked_where(x > 1.e20, x)
                yy = ma.masked_where(y > 1.e20, y)
                epsx = np.abs(xx[:,1:]-xx[:,0:-1]).max()
                epsy = np.abs(yy[1:,:]-yy[0:-1,:]).max()
                xymask = \
                np.logical_or(np.greater(x,self.xmax+epsx),np.greater(y,self.ymax+epsy))
                xymask = xymask + \
                np.logical_or(np.less(x,self.xmin-epsx),np.less(y,self.ymin-epsy))
                data = ma.asarray(data)
                # combine with data mask.
                mask = np.logical_or(ma.getmaskarray(data),xymask)
                data = ma.masked_array(data,mask=mask)
                CS = ax.contour(x,y,data,*args,**kwargs)
        finally:
            self._restore_hold(ax)
        # reset current active image (only if pyplot is imported).
        if plt and CS.get_array() is not None:
            plt.sci(CS)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip to map limbs
        CS.collections,c = self._cliplimb(ax,CS.collections)
        return CS

    @_transform
    def contourf(self,x,y,data,*args,**kwargs):
        """
        Make a filled contour plot over the map
        (see matplotlib.pyplot.contourf documentation).

        If ``latlon`` keyword is set to True, x,y are intrepreted as
        longitude and latitude in degrees.  Data and longitudes are
        automatically shifted to match map projection region for cylindrical
        and pseudocylindrical projections, and x,y are transformed to map
        projection coordinates. If ``latlon`` is False (default), x and y
        are assumed to be map projection coordinates.

        If x or y are outside projection limb (i.e. they have values > 1.e20),
        the corresponing data elements will be masked.

        Extra keyword 'ax' can be used to override the default axis instance.

        If ``tri`` is set to ``True``, an unstructured grid is assumed
        (x,y,data must be 1-d) and matplotlib.pyplot.tricontourf is used.

        Other \*args and \**kwargs passed on to matplotlib.pyplot.contourf
        (or tricontourf if ``tri=True``).
        """
        ax, plt = self._ax_plt_from_kw(kwargs)
        self._save_use_hold(ax, kwargs)
        try:
            if kwargs.get('tri', False):
                try:
                    import matplotlib.tri as tri
                except:
                    msg='need matplotlib > 0.99.1 to plot on unstructured grids'
                    raise ImportError(msg)
                # for unstructured grids, toss out points outside
                # projection limb (don't use those points in triangulation).
                if ma.isMA(data):
                    data = data.filled(fill_value=1.e30)
                    masked=True
                else:
                    masked=False
                mask = np.logical_or(x<1.e20,y<1.e20)
                x = np.compress(mask,x)
                y = np.compress(mask,y)
                data = np.compress(mask,data)
                if masked:
                    triang = tri.Triangulation(x, y)
                    z = data[triang.triangles]
                    mask = (z > 1.e20).sum(axis=-1)
                    triang.set_mask(mask)
                    CS = ax.tricontourf(triang,data,*args,**kwargs)
                else:
                    CS = ax.tricontourf(x,y,data,*args,**kwargs)
            else:
                # make sure x is monotonically increasing - if not,
                # print warning suggesting that the data be shifted in longitude
                # with the shiftgrid function.
                # only do this check for global projections.
                if self.projection in _cylproj + _pseudocyl:
                    xx = x[x.shape[0]//2,:]
                    condition = (xx >= self.xmin) & (xx <= self.xmax)
                    xl = xx.compress(condition).tolist()
                    xs = xl[:]
                    xs.sort()
                    if xl != xs:
                        sys.stdout.write(dedent("""
                             WARNING: x coordinate not montonically increasing - contour plot
                             may not be what you expect.  If it looks odd, your can either
                             adjust the map projection region to be consistent with your data, or
                             (if your data is on a global lat/lon grid) use the shiftgrid
                             function to adjust the data to be consistent with the map projection
                             region (see examples/contour_demo.py)."""))
                # mask for points more than one grid length outside projection limb.
                xx = ma.masked_where(x > 1.e20, x)
                yy = ma.masked_where(y > 1.e20, y)
                if self.projection != 'omerc':
                    epsx = np.abs(xx[:,1:]-xx[:,0:-1]).max()
                    epsy = np.abs(yy[1:,:]-yy[0:-1,:]).max()
                else: # doesn't work for omerc (FIXME)
                    epsx = 0.; epsy = 0
                xymask = \
                np.logical_or(np.greater(x,self.xmax+epsx),np.greater(y,self.ymax+epsy))
                xymask = xymask + \
                np.logical_or(np.less(x,self.xmin-epsx),np.less(y,self.ymin-epsy))
                data = ma.asarray(data)
                # combine with data mask.
                mask = np.logical_or(ma.getmaskarray(data),xymask)
                data = ma.masked_array(data,mask=mask)
                CS = ax.contourf(x,y,data,*args,**kwargs)
        finally:
            self._restore_hold(ax)
        # reset current active image (only if pyplot is imported).
        if plt and CS.get_array() is not None:
            plt.sci(CS)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip to map limbs
        CS.collections,c = self._cliplimb(ax,CS.collections)
        return CS

    @_transformuv
    def quiver(self, x, y, u, v, *args, **kwargs):
        """
        Make a vector plot (u, v) with arrows on the map.

        Arguments may be 1-D or 2-D arrays or sequences
        (see matplotlib.pyplot.quiver documentation for details).

        If ``latlon`` keyword is set to True, x,y are intrepreted as
        longitude and latitude in degrees.  Data and longitudes are
        automatically shifted to match map projection region for cylindrical
        and pseudocylindrical projections, and x,y are transformed to map
        projection coordinates. If ``latlon`` is False (default), x and y
        are assumed to be map projection coordinates.

        Extra keyword ``ax`` can be used to override the default axis instance.

        Other \*args and \**kwargs passed on to matplotlib.pyplot.quiver.
        """
        ax, plt = self._ax_plt_from_kw(kwargs)
        self._save_use_hold(ax, kwargs)
        try:
            ret =  ax.quiver(x,y,u,v,*args,**kwargs)
        finally:
            self._restore_hold(ax)
        if plt is not None and ret.get_array() is not None:
            plt.sci(ret)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip to map limbs
        ret,c = self._cliplimb(ax,ret)
        return ret

    @_transformuv
    def streamplot(self, x, y, u, v, *args, **kwargs):
        """
        Draws streamlines of a vector flow.
        (see matplotlib.pyplot.streamplot documentation).

        If ``latlon`` keyword is set to True, x,y are intrepreted as
        longitude and latitude in degrees.  Data and longitudes are
        automatically shifted to match map projection region for cylindrical
        and pseudocylindrical projections, and x,y are transformed to map
        projection coordinates. If ``latlon`` is False (default), x and y
        are assumed to be map projection coordinates.

        Extra keyword ``ax`` can be used to override the default axis instance.

        Other \*args and \**kwargs passed on to matplotlib.pyplot.streamplot.
        """
        if _matplotlib_version < '1.2':
            msg = dedent("""
            streamplot method requires matplotlib 1.2 or higher,
            you have %s""" % _matplotlib_version)
            raise NotImplementedError(msg)
        ax, plt = self._ax_plt_from_kw(kwargs)
        self._save_use_hold(ax, kwargs)
        try:
            ret =  ax.streamplot(x,y,u,v,*args,**kwargs)
        finally:
            self._restore_hold(ax)
        if plt is not None and ret.lines.get_array() is not None:
            plt.sci(ret.lines)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip to map limbs
        ret.lines,c = self._cliplimb(ax,ret.lines)
        ret.arrows,c = self._cliplimb(ax,ret.arrows)
        # streamplot arrows not returned in matplotlib 1.1.1, so clip all
        # FancyArrow patches attached to axes instance.
        if c is not None:
            for p in ax.patches:
                if isinstance(p,FancyArrowPatch): p.set_clip_path(c)
        return ret

    @_transformuv
    def barbs(self, x, y, u, v, *args, **kwargs):
        """
        Make a wind barb plot (u, v) with on the map.
        (see matplotlib.pyplot.barbs documentation).

        If ``latlon`` keyword is set to True, x,y are intrepreted as
        longitude and latitude in degrees.  Data and longitudes are
        automatically shifted to match map projection region for cylindrical
        and pseudocylindrical projections, and x,y are transformed to map
        projection coordinates. If ``latlon`` is False (default), x and y
        are assumed to be map projection coordinates.

        Extra keyword ``ax`` can be used to override the default axis instance.

        Other \*args and \**kwargs passed on to matplotlib.pyplot.barbs

        Returns two matplotlib.axes.Barbs instances, one for the Northern
        Hemisphere and one for the Southern Hemisphere.
        """
        if _matplotlib_version < '0.98.3':
            msg = dedent("""
            barb method requires matplotlib 0.98.3 or higher,
            you have %s""" % _matplotlib_version)
            raise NotImplementedError(msg)
        ax, plt = self._ax_plt_from_kw(kwargs)
        lons, lats = self(x, y, inverse=True)
        unh = ma.masked_where(lats <= 0, u)
        vnh = ma.masked_where(lats <= 0, v)
        ush = ma.masked_where(lats > 0, u)
        vsh = ma.masked_where(lats > 0, v)
        self._save_use_hold(ax, kwargs)
        try:
            retnh =  ax.barbs(x,y,unh,vnh,*args,**kwargs)
            kwargs['flip_barb']=True
            retsh =  ax.barbs(x,y,ush,vsh,*args,**kwargs)
        finally:
            self._restore_hold(ax)
        # Because there are two collections returned in general,
        # we can't set the current image...
        #if plt is not None and ret.get_array() is not None:
        #    plt.sci(retnh)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # clip to map limbs
        retnh,c = self._cliplimb(ax,retnh)
        retsh,c = self._cliplimb(ax,retsh)

        return retnh,retsh

    def drawlsmask(self,land_color="0.8",ocean_color="w",lsmask=None,
                   lsmask_lons=None,lsmask_lats=None,lakes=True,resolution='l',grid=5,**kwargs):
        """
        Draw land-sea mask image.

        .. note::
         The land-sea mask image cannot be overlaid on top
         of other images, due to limitations in matplotlib image handling
         (you can't specify the zorder of an image).

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keywords         Description
        ==============   ====================================================
        land_color       desired land color (color name or rgba tuple).
                         Default gray ("0.8").
        ocean_color      desired water color (color name or rgba tuple).
                         Default white.
        lsmask           An array of 0's for ocean pixels, 1's for
                         land pixels and 2's for lake/pond pixels.
                         Default is None
                         (default 5-minute resolution land-sea mask is used).
        lakes            Plot lakes and ponds (Default True)
        lsmask_lons      1d array of longitudes for lsmask (ignored
                         if lsmask is None). Longitudes must be ordered
                         from -180 W eastward.
        lsmask_lats      1d array of latitudes for lsmask (ignored
                         if lsmask is None). Latitudes must be ordered
                         from -90 S northward.
        resolution       gshhs coastline resolution used to define land/sea
                         mask (default 'l', available 'c','l','i','h' or 'f')
        grid             land/sea mask grid spacing in minutes (Default 5;
                         10, 2.5 and 1.25 are also available).
        \**kwargs        extra keyword arguments passed on to
                         :meth:`imshow`
        ==============   ====================================================

        If any of the lsmask, lsmask_lons or lsmask_lats keywords are not
        set, the built in GSHHS land-sea mask datasets are used.

        Extra keyword ``ax`` can be used to override the default axis instance.

        returns a matplotlib.image.AxesImage instance.
        """
        # convert land and water colors to integer rgba tuples with
        # values between 0 and 255.
        from matplotlib.colors import ColorConverter
        c = ColorConverter()
        # if conversion fails, assume it's because the color
        # given is already an rgba tuple with values between 0 and 255.
        try:
            cl = c.to_rgba(land_color)
            rgba_land = tuple([int(255*x) for x in cl])
        except:
            rgba_land = land_color
        try:
            co = c.to_rgba(ocean_color)
            rgba_ocean = tuple([int(255*x) for x in co])
        except:
            rgba_ocean = ocean_color
        # look for axes instance (as keyword, an instance variable
        # or from plt.gca().
        ax = kwargs.pop('ax', None) or self._check_ax()
        # Clear saved lsmask if new lsmask is passed
        if lsmask is not None or lsmask_lons is not None \
                or lsmask_lats is not None:
            # Make sure passed lsmask is not the same as cached mask
            if lsmask is not self.lsmask:
                self.lsmask = None
        # if lsmask,lsmask_lons,lsmask_lats keywords not given,
        # read default land-sea mask in from file.
        if lsmask is None or lsmask_lons is None or lsmask_lats is None:
            # if lsmask instance variable already set, data already
            # read in.
            if self.lsmask is None:
                # read in land/sea mask.
                lsmask_lons, lsmask_lats, lsmask =\
                _readlsmask(lakes=lakes,resolution=resolution,grid=grid)
                # instance variable lsmask is set on first invocation,
                # it contains the land-sea mask interpolated to the native
                # projection grid.  Further calls to drawlsmask will not
                # redo the interpolation (unless a new land-sea mask is passed
                # in via the lsmask, lsmask_lons, lsmask_lats keywords).

                # is it a cylindrical projection whose limits lie
                # outside the limits of the image?
                cylproj =  self.projection in _cylproj and \
                          (self.urcrnrlon > lsmask_lons[-1] or \
                           self.llcrnrlon < lsmask_lons[0])
                if cylproj:
                    # stack grids side-by-side (in longitiudinal direction), so
                    # any range of longitudes may be plotted on a world map.
                    # in versions of NumPy later than 1.10.0, concatenate will
                    # not stack these arrays as expected. If axis 1 is outside
                    # the dimensions of the array, concatenate will now raise
                    # an IndexError. Using hstack instead.
                    lsmask_lons = \
                            np.hstack((lsmask_lons,lsmask_lons[1:] + 360))
                    lsmask = \
                            np.hstack((lsmask,lsmask[:,1:]))
        else:
            if lakes: lsmask = np.where(lsmask==2,np.array(0,np.uint8),lsmask)

        # transform mask to nx x ny regularly spaced native projection grid
        # nx and ny chosen to have roughly the same horizontal
        # resolution as mask.
        if self.lsmask is None:
            nlons = len(lsmask_lons)
            nlats = len(lsmask_lats)
            if self.projection == 'cyl':
                dx = lsmask_lons[1]-lsmask_lons[0]
            else:
                dx = (np.pi/180.)*(lsmask_lons[1]-lsmask_lons[0])*self.rmajor
            nx = int((self.xmax-self.xmin)/dx)+1; ny = int((self.ymax-self.ymin)/dx)+1
        # interpolate rgba values from proj='cyl' (geographic coords)
        # to a rectangular map projection grid.
            mask,x,y = self.transform_scalar(lsmask,lsmask_lons,\
                       lsmask_lats,nx,ny,returnxy=True,order=0,masked=255)
            lsmask_lats.dtype
            # for these projections, points outside the projection
            # limb have to be set to transparent manually.
            if self.projection in _pseudocyl:
                lons, lats = self(x, y, inverse=True)
                lon_0 = self.projparams['lon_0']
                lats = lats[:,nx//2]
                lons1 = (lon_0+180.)*np.ones(lons.shape[0],np.float64)
                lons2 = (lon_0-180.)*np.ones(lons.shape[0],np.float64)
                xmax,ytmp = self(lons1,lats)
                xmin,ytmp = self(lons2,lats)
                for j in range(lats.shape[0]):
                    xx = x[j,:]
                    mask[j,:]=np.where(np.logical_or(xx<xmin[j],xx>xmax[j]),\
                                        255,mask[j,:])
            self.lsmask = mask
        ny, nx = self.lsmask.shape
        rgba = np.ones((ny,nx,4),np.uint8)
        rgba_land = np.array(rgba_land,np.uint8)
        rgba_ocean = np.array(rgba_ocean,np.uint8)
        for k in range(4):
            rgba[:,:,k] = np.where(self.lsmask,rgba_land[k],rgba_ocean[k])
        # make points outside projection limb transparent.
        rgba[:,:,3] = np.where(self.lsmask==255,0,rgba[:,:,3])
        # plot mask as rgba image.
        im = self.imshow(rgba,interpolation='nearest',ax=ax,**kwargs)
        # clip to map limbs.
        im,c = self._cliplimb(ax,im)
        return im

    def bluemarble(self,ax=None,scale=None,**kwargs):
        """
        display blue marble image (from http://visibleearth.nasa.gov)
        as map background.
        Default image size is 5400x2700, which can be quite slow and
        use quite a bit of memory.  The ``scale`` keyword can be used
        to downsample the image (``scale=0.5`` downsamples to 2700x1350).

        \**kwargs passed on to :meth:`imshow`.

        returns a matplotlib.image.AxesImage instance.
        """
        if ax is not None:
            return self.warpimage(image='bluemarble',ax=ax,scale=scale,**kwargs)
        else:
            return self.warpimage(image='bluemarble',scale=scale,**kwargs)

    def shadedrelief(self,ax=None,scale=None,**kwargs):
        """
        display shaded relief image (from http://www.shadedrelief.com)
        as map background.
        Default image size is 10800x5400, which can be quite slow and
        use quite a bit of memory.  The ``scale`` keyword can be used
        to downsample the image (``scale=0.5`` downsamples to 5400x2700).

        \**kwargs passed on to :meth:`imshow`.

        returns a matplotlib.image.AxesImage instance.
        """
        if ax is not None:
            return self.warpimage(image='shadedrelief',ax=ax,scale=scale,**kwargs)
        else:
            return self.warpimage(image='shadedrelief',scale=scale,**kwargs)

    def etopo(self,ax=None,scale=None,**kwargs):
        """
        display etopo relief image (from
        http://www.ngdc.noaa.gov/mgg/global/global.html)
        as map background.
        Default image size is 5400x2700, which can be quite slow and
        use quite a bit of memory.  The ``scale`` keyword can be used
        to downsample the image (``scale=0.5`` downsamples to 5400x2700).

        \**kwargs passed on to :meth:`imshow`.

        returns a matplotlib.image.AxesImage instance.
        """
        if ax is not None:
            return self.warpimage(image='etopo',ax=ax,scale=scale,**kwargs)
        else:
            return self.warpimage(image='etopo',scale=scale,**kwargs)

    def warpimage(self,image="bluemarble",scale=None,**kwargs):
        """
        Display an image (filename given by ``image`` keyword) as a map background.
        If image is a URL (starts with 'http'), it is downloaded to a temp
        file using urllib.urlretrieve.

        Default (if ``image`` not specified) is to display
        'blue marble next generation' image from http://visibleearth.nasa.gov/.

        Specified image must have pixels covering the whole globe in a regular
        lat/lon grid, starting and -180W and the South Pole.
        Works with the global images from
        http://earthobservatory.nasa.gov/Features/BlueMarble/BlueMarble_monthlies.php.

        The ``scale`` keyword can be used to downsample (rescale) the image.
        Values less than 1.0 will speed things up at the expense of image
        resolution.

        Extra keyword ``ax`` can be used to override the default axis instance.

        \**kwargs passed on to :meth:`imshow`.

        returns a matplotlib.image.AxesImage instance.
        """

        # fix PIL import on some versions of OSX and scipy
        try:
            from PIL import Image
        except ImportError:
            try:
                import Image
            except ImportError:
                msg = ('warpimage method requires PIL '
                       '(http://pillow.readthedocs.io)')
                raise ImportError(msg)

        from matplotlib.image import pil_to_array
        if self.celestial:
            msg='warpimage does not work in celestial coordinates'
            raise ValueError(msg)
        ax = kwargs.pop('ax', None) or self._check_ax()
        # default image file is blue marble next generation
        # from NASA (http://visibleearth.nasa.gov).
        if image == "bluemarble":
            file = os.path.join(basemap_datadir,'bmng.jpg')
        # display shaded relief image (from
        # http://www.shadedreliefdata.com)
        elif image == "shadedrelief":
            file = os.path.join(basemap_datadir,'shadedrelief.jpg')
        # display etopo image (from
        # http://www.ngdc.noaa.gov/mgg/image/globalimages.html)
        elif image == "etopo":
            file = os.path.join(basemap_datadir,'etopo1.jpg')
        else:
            file = image
        # if image is same as previous invocation, used cached data.
        # if not, regenerate rgba data.
        if not hasattr(self,'_bm_file') or self._bm_file != file:
            newfile = True
        else:
            newfile = False
        if file.startswith('http'):
            self._bm_file, headers = urlretrieve(file)
        else:
            self._bm_file = file
        # bmproj is True if map projection region is same as
        # image region.
        bmproj = self.projection == 'cyl' and \
                 self.llcrnrlon == -180 and self.urcrnrlon == 180 and \
                 self.llcrnrlat == -90 and self.urcrnrlat == 90
        # read in jpeg image to rgba array of normalized floats.
        if not hasattr(self,'_bm_rgba') or newfile:
            pilImage = Image.open(self._bm_file)
            if scale is not None:
                w, h = pilImage.size
                width = int(np.round(w*scale))
                height = int(np.round(h*scale))
                pilImage = pilImage.resize((width,height),Image.LANCZOS)
            if _matplotlib_version >= '1.2':
                # orientation of arrays returned by pil_to_array
                # changed (https://github.com/matplotlib/matplotlib/pull/616)
                self._bm_rgba = pil_to_array(pilImage)[::-1,:]
            else:
                self._bm_rgba = pil_to_array(pilImage)
            # define lat/lon grid that image spans.
            nlons = self._bm_rgba.shape[1]; nlats = self._bm_rgba.shape[0]
            delta = 360./float(nlons)
            self._bm_lons = np.arange(-180.+0.5*delta,180.,delta)
            self._bm_lats = np.arange(-90.+0.5*delta,90.,delta)
            # is it a cylindrical projection whose limits lie
            # outside the limits of the image?
            cylproj =  self.projection in _cylproj and \
                      (self.urcrnrlon > self._bm_lons[-1] or \
                       self.llcrnrlon < self._bm_lons[0])
            # if pil_to_array returns a 2D array, it's a grayscale image.
            # create an RGB image, with R==G==B.
            if self._bm_rgba.ndim == 2:
                tmp = np.empty(self._bm_rgba.shape+(3,),np.uint8)
                for k in range(3):
                    tmp[:,:,k] = self._bm_rgba
                self._bm_rgba = tmp
            if cylproj and not bmproj:
                # stack grids side-by-side (in longitiudinal direction), so
                # any range of longitudes may be plotted on a world map.
                self._bm_lons = \
                np.concatenate((self._bm_lons,self._bm_lons+360),0)
                self._bm_rgba = \
                np.concatenate((self._bm_rgba,self._bm_rgba),1)
            # convert to normalized floats.
            self._bm_rgba = self._bm_rgba.astype(np.float32)/255.
        if not bmproj: # interpolation necessary.
            if newfile or not hasattr(self,'_bm_rgba_warped'):
                # transform to nx x ny regularly spaced native
                # projection grid.
                # nx and ny chosen to have roughly the
                # same horizontal res as original image.
                if self.projection != 'cyl':
                    dx = 2.*np.pi*self.rmajor/float(nlons)
                    nx = int((self.xmax-self.xmin)/dx)+1
                    ny = int((self.ymax-self.ymin)/dx)+1
                else:
                    dx = 360./float(nlons)
                    nx = int((self.urcrnrlon-self.llcrnrlon)/dx)+1
                    ny = int((self.urcrnrlat-self.llcrnrlat)/dx)+1
                self._bm_rgba_warped = np.ones((ny,nx,4),np.float64)
                # interpolate rgba values from geographic coords (proj='cyl')
                # to map projection coords.
                # if masked=True, values outside of
                # projection limb will be masked.
                for k in range(self._bm_rgba.shape[2]):
                    self._bm_rgba_warped[:,:,k],x,y = \
                    self.transform_scalar(self._bm_rgba[:,:,k],\
                    self._bm_lons,self._bm_lats,nx,ny,returnxy=True)
                # for ortho,geos mask pixels outside projection limb.
                if self.projection in ['geos','ortho','nsper'] or \
                   (self.projection == 'aeqd' and self._fulldisk):
                    lonsr,latsr = self(x,y,inverse=True)
                    mask = ma.zeros((ny,nx,4),np.int8)
                    mask[:,:,0] = np.logical_or(lonsr>1.e20,latsr>1.e30)
                    for k in range(1,4):
                        mask[:,:,k] = mask[:,:,0]
                    self._bm_rgba_warped = \
                    ma.masked_array(self._bm_rgba_warped,mask=mask)
                    # make points outside projection limb transparent.
                    self._bm_rgba_warped = self._bm_rgba_warped.filled(0.)
                # treat pseudo-cyl projections such as mollweide, robinson and sinusoidal.
                elif self.projection in _pseudocyl and \
                     self.projection != 'hammer':
                    lonsr,latsr = self(x,y,inverse=True)
                    mask = ma.zeros((ny,nx,4),np.int8)
                    lon_0 = self.projparams['lon_0']
                    lonright = lon_0+180.
                    lonleft = lon_0-180.
                    x1 = np.array(ny*[0.5*(self.xmax + self.xmin)],np.float)
                    y1 = np.linspace(self.ymin, self.ymax, ny)
                    lons1, lats1 = self(x1,y1,inverse=True)
                    lats1 = np.where(lats1 < -89.999999, -89.999999, lats1)
                    lats1 = np.where(lats1 > 89.999999, 89.999999, lats1)
                    for j,lat in enumerate(lats1):
                        xmax,ymax = self(lonright,lat)
                        xmin,ymin = self(lonleft,lat)
                        mask[j,:,0] = np.logical_or(x[j,:]>xmax,x[j,:]<xmin)
                    for k in range(1,4):
                        mask[:,:,k] = mask[:,:,0]
                    self._bm_rgba_warped = \
                    ma.masked_array(self._bm_rgba_warped,mask=mask)
                    # make points outside projection limb transparent.
                    # FIXME: Probably not needed anymore
                    self._bm_rgba_warped = self._bm_rgba_warped.filled(0.)
            # plot warped rgba image.
            im = self.imshow(self._bm_rgba_warped,ax=ax,**kwargs)
            # for hammer projection, use clip path defined by
            # projection limb (patch created in drawmapboundary).
            # FIXME: Is this now redundant?
            if self.projection == 'hammer':
                if not self._mapboundarydrawn:
                    self.drawmapboundary(color='none',linewidth=None)
                im.set_clip_path(self._mapboundarydrawn)
        else:
            # bmproj True, no interpolation necessary.
            im = self.imshow(self._bm_rgba,ax=ax,**kwargs)
        # clip to map limbs
        im,c = self._cliplimb(ax,im)
        return im

    def arcgisimage(self,server='http://server.arcgisonline.com/ArcGIS',\
                 service='World_Imagery',xpixels=400,ypixels=None,\
                 dpi=96,verbose=False,**kwargs):
        """
        Retrieve an image using the ArcGIS Server REST API and display it on
        the map. In order to use this method, the Basemap instance must be
        created using the ``epsg`` keyword to define the map projection, unless
        the ``cyl`` projection is used (in which case the epsg code 4326 is
        assumed).

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keywords         Description
        ==============   ====================================================
        server           web map server URL (default
                         http://server.arcgisonline.com/ArcGIS).
        service          service (image type) hosted on server (default
                         'World_Imagery', which is NASA 'Blue Marble'
                         image).
        xpixels          requested number of image pixels in x-direction
                         (default 400).
        ypixels          requested number of image pixels in y-direction.
                         Default (None) is to infer the number from
                         from xpixels and the aspect ratio of the
                         map projection region.
        dpi              The device resolution of the exported image (dots per
                         inch, default 96).
        verbose          if True, print URL used to retrieve image (default
                         False).
        ==============   ====================================================

        Extra keyword ``ax`` can be used to override the default axis instance.

        returns a matplotlib.image.AxesImage instance.
        """


        # fix PIL import on some versions of OSX and scipy
        try:
            from PIL import Image
        except ImportError:
            try:
                import Image
            except ImportError:
                msg = ('arcgisimage method requires PIL '
                       '(http://pillow.readthedocs.io)')
                raise ImportError(msg)



        if not hasattr(self,'epsg'):
            msg = dedent("""
            Basemap instance must be creating using an EPSG code
            (http://spatialreference.org) in order to use the wmsmap method""")
            raise ValueError(msg)
        ax = kwargs.pop('ax', None) or self._check_ax()
        # find the x,y values at the corner points.
        p = pyproj.Proj(init="epsg:%s" % self.epsg, preserve_units=True)
        xmin,ymin = p(self.llcrnrlon,self.llcrnrlat)
        xmax,ymax = p(self.urcrnrlon,self.urcrnrlat)
        if self.projection in _cylproj:
            Dateline =\
            _geoslib.Point(self(180.,0.5*(self.llcrnrlat+self.urcrnrlat)))
            hasDateline = Dateline.within(self._boundarypolyxy)
            if hasDateline:
                msg=dedent("""
                arcgisimage cannot handle images that cross
                the dateline for cylindrical projections.""")
                raise ValueError(msg)
        # ypixels not given, find by scaling xpixels by the map aspect ratio.
        if ypixels is None:
            ypixels = int(self.aspect*xpixels)
        # construct a URL using the ArcGIS Server REST API.
        basemap_url = \
"%s/rest/services/%s/MapServer/export?\
bbox=%s,%s,%s,%s&\
bboxSR=%s&\
imageSR=%s&\
size=%s,%s&\
dpi=%s&\
format=png32&\
transparent=true&\
f=image" %\
(server,service,xmin,ymin,xmax,ymax,self.epsg,self.epsg,xpixels,ypixels,dpi)
        # print URL?
        if verbose: print(basemap_url)
        # return AxesImage instance.
        img = Image.open(urlopen(basemap_url))
        return self.imshow(img, ax=ax, origin='upper')

    def wmsimage(self,server,\
                 xpixels=400,ypixels=None,\
                 format='png',alpha=None,verbose=False,**kwargs):
        """
        Retrieve an image using from a WMS server using the
        Open Geospatial Consortium (OGC) standard interface
        and display on the map. Requires OWSLib
        (http://pypi.python.org/pypi/OWSLib).
        In order to use this method, the Basemap instance must be
        created using the ``epsg`` keyword to define the map projection, unless
        the ``cyl`` projection is used (in which case the epsg code 4326 is
        assumed).

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keywords         Description
        ==============   ====================================================
        server           WMS server URL.
        xpixels          requested number of image pixels in x-direction
                         (default 400).
        ypixels          requested number of image pixels in y-direction.
                         Default (None) is to infer the number from
                         from xpixels and the aspect ratio of the
                         map projection region.
        format           desired image format (default 'png')
        alpha            The alpha blending value,
                         between 0 (transparent) and 1 (opaque) (default None)
        verbose          if True, print WMS server info (default
                         False).
        \**kwargs        extra keyword arguments passed on to
                         OWSLib.wms.WebMapService.getmap.
        ==============   ====================================================

        Extra keyword ``ax`` can be used to override the default axis instance.

        returns a matplotlib.image.AxesImage instance.
        """
        try:
            from owslib.wms import WebMapService
        except ImportError:
            raise ImportError('OWSLib required to use wmsimage method')
        import io
        ax = kwargs.pop('ax', None) or self._check_ax()
        if not hasattr(self,'epsg'):
            msg = dedent("""
            Basemap instance must be creating using an EPSG code
            (http://spatialreference.org) in order to use the wmsmap method""")
            raise ValueError(msg)
        if 'layers' not in kwargs:
            raise ValueError('no layers specified')
        # find the x,y values at the corner points.
        p = pyproj.Proj(init="epsg:%s" % self.epsg, preserve_units=True)
        xmin,ymin = p(self.llcrnrlon,self.llcrnrlat)
        xmax,ymax = p(self.urcrnrlon,self.urcrnrlat)
        if self.projection in _cylproj:
            Dateline =\
            _geoslib.Point(self(180.,0.5*(self.llcrnrlat+self.urcrnrlat)))
            hasDateline = Dateline.within(self._boundarypolyxy)
            if hasDateline:
                msg=dedent("""
                wmsimage cannot handle images that cross
                the dateline for cylindrical projections.""")
                raise ValueError(msg)
        if self.projection == 'cyl':
            xmin = (180./np.pi)*xmin; xmax = (180./np.pi)*xmax
            ymin = (180./np.pi)*ymin; ymax = (180./np.pi)*ymax
        # ypixels not given, find by scaling xpixels by the map aspect ratio.
        if ypixels is None:
            ypixels = int(self.aspect*xpixels)
        if verbose: print(server)
        wms = WebMapService(server)
        if verbose:
            print('id: %s, version: %s' %
            (wms.identification.type,wms.identification.version))
            print('title: %s, abstract: %s' %
            (wms.identification.title,wms.identification.abstract))
            print('available layers:')
            print(list(wms.contents))
            print('projection options:')
            print(wms[kwargs['layers'][0]].crsOptions)
        # remove keys from kwargs that are over-ridden
        for k in ['format','bbox','service','size','srs']:
            if 'format' in kwargs: del kwargs['format']
        img = wms.getmap(service='wms',bbox=(xmin,ymin,xmax,ymax),
                         size=(xpixels,ypixels),format='image/%s'%format,
                         srs='EPSG:%s' % self.epsg, **kwargs)
        # return AxesImage instance.
        # this works for png and jpeg.
        return self.imshow(imread(io.BytesIO(img.read()),
                           format=format),origin='upper',alpha=alpha,ax=ax)

    def drawmapscale(self,lon,lat,lon0,lat0,length,barstyle='simple',\
                     units='km',fontsize=9,yoffset=None,labelstyle='simple',\
                     fontcolor='k',fillcolor1='w',fillcolor2='k',ax=None,\
                     format='%d',zorder=None,linecolor=None,linewidth=None):
        """
        Draw a map scale at ``lon,lat`` of length ``length``
        representing distance in the map
        projection coordinates at ``lon0,lat0``.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keywords         Description
        ==============   ====================================================
        units            the units of the length argument (Default km).
        barstyle         ``simple`` or ``fancy`` (roughly corresponding
                         to the styles provided by Generic Mapping Tools).
                         Default ``simple``.
        fontsize         for map scale annotations, default 9.
        fontcolor            for map scale annotations, default black.
        labelstyle       ``simple`` (default) or ``fancy``.  For
                         ``fancy`` the map scale factor (ratio betwee
                         the actual distance and map projection distance
                         at lon0,lat0) and the value of lon0,lat0 are also
                         displayed on the top of the scale bar. For
                         ``simple``, just the units are display on top
                         and the distance below the scale bar.
                         If equal to False, plot an empty label.
        format           a string formatter to format numeric values
        yoffset          yoffset controls how tall the scale bar is,
                         and how far the annotations are offset from the
                         scale bar.  Default is 0.02 times the height of
                         the map (0.02*(self.ymax-self.ymin)).
        fillcolor1(2)    colors of the alternating filled regions
                         (default white and black).  Only relevant for
                         'fancy' barstyle.
        zorder           sets the zorder for the map scale.
        linecolor        sets the color of the scale, by default, fontcolor
                         is used
        linewidth        linewidth for scale and ticks
        ==============   ====================================================

        Extra keyword ``ax`` can be used to override the default axis instance.
        """
        # get current axes instance (if none specified).
        ax = ax or self._check_ax()
        # not valid for cylindrical projection
        if self.projection == 'cyl':
            raise ValueError("cannot draw map scale for projection='cyl'")
        # convert length to meters
        lenlab = length
        if units == 'km':
            length = length*1000
        elif units == 'mi':
            length = length*1609.344
        elif units == 'nmi':
            length = length*1852
        elif units == 'ft':
            length = length*0.3048
        elif units != 'm':
            msg = "units must be 'm' (meters), 'km' (kilometers), "\
            "'mi' (miles), 'nmi' (nautical miles), or 'ft' (feet)"
            raise KeyError(msg)
        # reference point and center of scale.
        x0,y0 = self(lon0,lat0)
        xc,yc = self(lon,lat)
        # make sure lon_0 between -180 and 180
        lon_0 = ((lon0+360) % 360) - 360
        if lat0>0:
            if lon>0:
                lonlatstr = u'%g\N{DEGREE SIGN}N, %g\N{DEGREE SIGN}E' % (lat0,lon_0)
            elif lon<0:
                lonlatstr = u'%g\N{DEGREE SIGN}N, %g\N{DEGREE SIGN}W' % (lat0,lon_0)
            else:
                lonlatstr = u'%g\N{DEGREE SIGN}, %g\N{DEGREE SIGN}W' % (lat0,lon_0)
        else:
            if lon>0:
                lonlatstr = u'%g\N{DEGREE SIGN}S, %g\N{DEGREE SIGN}E' % (lat0,lon_0)
            elif lon<0:
                lonlatstr = u'%g\N{DEGREE SIGN}S, %g\N{DEGREE SIGN}W' % (lat0,lon_0)
            else:
                lonlatstr = u'%g\N{DEGREE SIGN}S, %g\N{DEGREE SIGN}' % (lat0,lon_0)
        # left edge of scale
        lon1,lat1 = self(x0-length/2,y0,inverse=True)
        x1,y1 = self(lon1,lat1)
        # right edge of scale
        lon4,lat4 = self(x0+length/2,y0,inverse=True)
        x4,y4 = self(lon4,lat4)
        x1 = x1-x0+xc; y1 = y1-y0+yc
        x4 = x4-x0+xc; y4 = y4-y0+yc
        if x1 > 1.e20 or x4 > 1.e20 or y1 > 1.e20 or y4 > 1.e20:
            raise ValueError("scale bar positioned outside projection limb")
        # scale factor for true distance
        gc = pyproj.Geod(a=self.rmajor,b=self.rminor)
        az12,az21,dist = gc.inv(lon1,lat1,lon4,lat4)
        scalefact = dist/length
        # label to put on top of scale bar.
        if labelstyle=='simple':
            labelstr = units
        elif labelstyle == 'fancy':
            labelstr = units+" (scale factor %4.2f at %s)"%(scalefact,lonlatstr)
        elif labelstyle == False:
            labelstr = ''
        else:
            raise KeyError("labelstyle must be 'simple' or 'fancy'")
        # default y offset is 2 percent of map height.
        if yoffset is None: yoffset = 0.02*(self.ymax-self.ymin)
        rets = [] # will hold all plot objects generated.
        # set linecolor
        if linecolor is None:
            linecolor = fontcolor
        # 'fancy' style
        if barstyle == 'fancy':
            #we need 5 sets of x coordinates (in map units)
            #quarter scale
            lon2,lat2 = self(x0-length/4,y0,inverse=True)
            x2,y2 = self(lon2,lat2)
            x2 = x2-x0+xc; y2 = y2-y0+yc
            #three quarter scale
            lon3,lat3 = self(x0+length/4,y0,inverse=True)
            x3,y3 = self(lon3,lat3)
            x3 = x3-x0+xc; y3 = y3-y0+yc
            #plot top line
            ytop = yc+yoffset/2
            ybottom = yc-yoffset/2
            ytick = ybottom - yoffset/2
            ytext = ytick - yoffset/2
            rets.append(self.plot([x1,x4],[ytop,ytop],color=linecolor, linewidth=linewidth)[0])
            #plot bottom line
            rets.append(self.plot([x1,x4],[ybottom,ybottom],color=linecolor, linewidth=linewidth)[0])
            #plot left edge
            rets.append(self.plot([x1,x1],[ybottom,ytop],color=linecolor, linewidth=linewidth)[0])
            #plot right edge
            rets.append(self.plot([x4,x4],[ybottom,ytop],color=linecolor, linewidth=linewidth)[0])
            #make a filled black box from left edge to 1/4 way across
            rets.append(ax.fill([x1,x2,x2,x1,x1],[ytop,ytop,ybottom,ybottom,ytop],\
                        ec=fontcolor,fc=fillcolor1)[0])
            #make a filled white box from 1/4 way across to 1/2 way across
            rets.append(ax.fill([x2,xc,xc,x2,x2],[ytop,ytop,ybottom,ybottom,ytop],\
                        ec=fontcolor,fc=fillcolor2)[0])
            #make a filled white box from 1/2 way across to 3/4 way across
            rets.append(ax.fill([xc,x3,x3,xc,xc],[ytop,ytop,ybottom,ybottom,ytop],\
                        ec=fontcolor,fc=fillcolor1)[0])
            #make a filled white box from 3/4 way across to end
            rets.append(ax.fill([x3,x4,x4,x3,x3],[ytop,ytop,ybottom,ybottom,ytop],\
                        ec=fontcolor,fc=fillcolor2)[0])
            #plot 3 tick marks at left edge, center, and right edge
            rets.append(self.plot([x1,x1],[ytick,ybottom],color=linecolor, linewidth=linewidth)[0])
            rets.append(self.plot([xc,xc],[ytick,ybottom],color=linecolor, linewidth=linewidth)[0])
            rets.append(self.plot([x4,x4],[ytick,ybottom],color=linecolor, linewidth=linewidth)[0])
            #label 3 tick marks
            rets.append(ax.text(x1,ytext,format % (0),\
            horizontalalignment='center',\
            verticalalignment='top',\
            fontsize=fontsize,color=fontcolor))
            rets.append(ax.text(xc,ytext,format % (0.5*lenlab),\
            horizontalalignment='center',\
            verticalalignment='top',\
            fontsize=fontsize,color=fontcolor))
            rets.append(ax.text(x4,ytext,format % (lenlab),\
            horizontalalignment='center',\
            verticalalignment='top',\
            fontsize=fontsize,color=fontcolor))
            #put units, scale factor on top
            rets.append(ax.text(xc,ytop+yoffset/2,labelstr,\
            horizontalalignment='center',\
            verticalalignment='bottom',\
            fontsize=fontsize,color=fontcolor))
        # 'simple' style
        elif barstyle == 'simple':
            rets.append(self.plot([x1,x4],[yc,yc],color=linecolor, linewidth=linewidth)[0])
            rets.append(self.plot([x1,x1],[yc-yoffset,yc+yoffset],color=linecolor, linewidth=linewidth)[0])
            rets.append(self.plot([x4,x4],[yc-yoffset,yc+yoffset],color=linecolor, linewidth=linewidth)[0])
            rets.append(ax.text(xc,yc-yoffset,format % lenlab,\
            verticalalignment='top',horizontalalignment='center',\
            fontsize=fontsize,color=fontcolor))
            #put units, scale factor on top
            rets.append(ax.text(xc,yc+yoffset,labelstr,\
            horizontalalignment='center',\
            verticalalignment='bottom',\
            fontsize=fontsize,color=fontcolor))
        else:
            raise KeyError("barstyle must be 'simple' or 'fancy'")
        if zorder is not None:
            for ret in rets:
                try:
                    ret.set_zorder(zorder)
                except:
                    pass
        return rets

    def colorbar(self,mappable=None,location='right',size="5%",pad='2%',fig=None,ax=None,**kwargs):
        """
        Add colorbar to axes associated with a map.
        The colorbar axes instance is created using the axes_grid toolkit.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keywords         Description
        ==============   ====================================================
        mappable         the Image, ContourSet, etc. to which the colorbar
                         applies.  Default None, matplotlib.pyplot.gci() is
                         used to retrieve the current image mappable.
        location         where to put colorbar ('top','bottom','left','right')
                         Default 'right'.
        size             width of colorbar axes (string 'N%', where N is
                         an integer describing the fractional width of the parent
                         axes). Default '5%'.
        pad              Padding between parent axes and colorbar axes in
                         same units as size. Default '2%'.
        fig              Figure instance the map axes instance is associated
                         with. Default None, and matplotlib.pyplot.gcf() is used
                         to retrieve the current active figure instance.
        ax               The axes instance which the colorbar will be
                         associated with.  Default None, searches for self.ax,
                         and if None uses matplotlib.pyplot.gca().
        \**kwargs        extra keyword arguments passed on to
                         colorbar method of the figure instance.
        ==============   ====================================================

        Returns a matplotlib colorbar instance.
        """
        # get current axes instance (if none specified).
        ax = ax or self._check_ax()
        # get current figure instance (if none specified).
        if fig is None or mappable is None:
            import matplotlib.pyplot as plt
        if fig is None:
            fig = plt.gcf()
        # get current mappable if none specified.
        if mappable is None:
            mappable = plt.gci()
        # create colorbar axes uses axes_grid toolkit.
        divider = make_axes_locatable(ax)
        if location in ['left','right']:
            orientation = 'vertical'
        elif location in ['top','bottom']:
            orientation = 'horizontal'
        else:
            raise ValueError('location must be top,bottom,left or right')
        cax = divider.append_axes(location, size=size, pad=pad)
        # create colorbar.
        cb = fig.colorbar(mappable,orientation=orientation,cax=cax,**kwargs)
        fig.sca(ax) # reset parent axes as current axes.
        return cb

    def nightshade(self,date,color="k",delta=0.25,alpha=0.5,ax=None,zorder=2):
        """
        Shade the regions of the map that are in darkness at the time
        specifed by ``date``.  ``date`` is a datetime instance,
        assumed to be UTC.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keywords         Description
        ==============   ====================================================
        color            color to shade night regions (default black).
        delta            day/night terminator is computed with a
                         a resolution of ``delta`` degrees (default 0.25).
        alpha            alpha transparency for shading (default 0.5, so
                         map background shows through).
        zorder           zorder for shading (default 2).
        ==============   ====================================================

        Extra keyword ``ax`` can be used to override the default axis instance.

        returns a matplotlib.contour.ContourSet instance.
        """
        from .solar import daynight_grid
        # make sure date is UTC, or naive with repect to time zones
        if date.utcoffset():
            raise ValueError('datetime instance must be UTC, not {0}'.format(date.tzname()))
        # get current axes instance (if none specified).
        ax = ax or self._check_ax()
        # create grid of day=0, night=1
        lons,lats,daynight = daynight_grid(date,delta,self.lonmin,self.lonmax)
        x,y = self(lons,lats)
        # contour the day-night grid, coloring the night area
        # with the specified color and transparency.
        CS = self.contourf(x,y,daynight,1,colors=[color],alpha=alpha,ax=ax)
        # set zorder on ContourSet collections show night shading
        # is on top.
        for c in CS.collections:
            c.set_zorder(zorder)
        # clip to map limbs
        CS.collections,c = self._cliplimb(ax,CS.collections)
        return CS

    def _check_ax(self):
        """
        Returns the axis on which to draw.
        Returns self.ax, or if self.ax=None returns plt.gca().
        """
        if self.ax is None:
            try:
                ax = plt.gca()
            except:
                import matplotlib.pyplot as plt
                ax = plt.gca()
            # associate an axes instance with this Basemap instance
            # the first time this method is called.
            #self.ax = ax
        else:
            ax = self.ax
        return ax

    def _ax_plt_from_kw(self, kw):
        """
        Return (ax, plt), where ax is the current axes, and plt is
        None or a reference to the pyplot module.

        plt will be None if ax was popped from kw or taken from self.ax;
        otherwise, pyplot was used and is returned.
        """
        plt = None
        _ax = kw.pop('ax', None)
        if _ax is None:
            _ax = self.ax
            if _ax is None:
                import matplotlib.pyplot as plt
                _ax = plt.gca()
        return _ax, plt

    def shiftdata(self,lonsin,datain=None,lon_0=None,fix_wrap_around=True):
        """
        Shift longitudes (and optionally data) so that they match map projection region.
        Only valid for cylindrical/pseudo-cylindrical global projections and data
        on regular lat/lon grids. longitudes and data can be 1-d or 2-d, if 2-d
        it is assumed longitudes are 2nd (rightmost) dimension.

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Arguments        Description
        ==============   ====================================================
        lonsin           original 1-d or 2-d longitudes.
        ==============   ====================================================

        .. tabularcolumns:: |l|L|

        ==============   ====================================================
        Keywords         Description
        ==============   ====================================================
        datain           original 1-d or 2-d data. Default None.
        lon_0            center of map projection region. Defaut None,
                         given by current map projection.
        fix_wrap_around  if True reindex (if required) longitudes (and data) to
                         avoid jumps caused by remapping of longitudes of
                         points from outside of the [lon_0-180, lon_0+180]
                         interval back into the interval.
                         If False do not reindex longitudes and data, but do
                         make sure that longitudes are in the
                         [lon_0-180, lon_0+180] range.
        ==============   ====================================================

        if datain given, returns ``dataout,lonsout`` (data and longitudes shifted to fit in interval
        [lon_0-180,lon_0+180]), otherwise just returns longitudes.  If
        transformed longitudes lie outside map projection region, data is
        masked and longitudes are set to 1.e30.
        """
        if lon_0 is None and 'lon_0' not in self.projparams:
            msg='lon_0 keyword must be provided'
            raise ValueError(msg)
        lonsin = np.asarray(lonsin)
        if lonsin.ndim not in [1,2]:
            raise ValueError('1-d or 2-d longitudes required')
        if datain is not None:
            # if it's a masked array, leave it alone.
            if not ma.isMA(datain): datain = np.asarray(datain)
            if datain.ndim not in [1,2]:
                raise ValueError('1-d or 2-d data required')
        if lon_0 is None:
            lon_0 = self.projparams['lon_0']
        # 2-d data.
        if lonsin.ndim == 2:
            nlats = lonsin.shape[0]
            nlons = lonsin.shape[1]
            lonsin1 = lonsin[0,:]
            lonsin1 = np.where(lonsin1 > lon_0+180, lonsin1-360 ,lonsin1)
            lonsin1 = np.where(lonsin1 < lon_0-180, lonsin1+360 ,lonsin1)
            if nlons > 1:
                londiff = np.abs(lonsin1[0:-1]-lonsin1[1:])
                londiff_sort = np.sort(londiff)
                thresh = 360.-londiff_sort[-2] if nlons > 2 else 360.-londiff_sort[-1]
                itemindex = nlons-np.where(londiff>=thresh)[0]
            else:
                lonsin[0, :] = lonsin1
                itemindex = 0

            # if no shift necessary, itemindex will be
            # empty, so don't do anything
            if fix_wrap_around and itemindex:
                # check to see if cyclic (wraparound) point included
                # if so, remove it.
                if np.abs(lonsin1[0]-lonsin1[-1]) < 1.e-4:
                    hascyclic = True
                    lonsin_save = lonsin.copy()
                    lonsin = lonsin[:,1:]
                    if datain is not None:
                       datain_save = datain.copy()
                       datain = datain[:,1:]
                else:
                    hascyclic = False
                lonsin = np.where(lonsin > lon_0+180, lonsin-360 ,lonsin)
                lonsin = np.where(lonsin < lon_0-180, lonsin+360 ,lonsin)
                lonsin = np.roll(lonsin,itemindex-1,axis=1)
                if datain is not None:
                    # np.roll works on ndarrays and on masked arrays
                    datain = np.roll(datain,itemindex-1,axis=1)
                # add cyclic point back at beginning.
                if hascyclic:
                    lonsin_save[:,1:] = lonsin
                    lonsin_save[:,0] = lonsin[:,-1]-360.
                    lonsin = lonsin_save
                    if datain is not None:
                        datain_save[:,1:] = datain
                        datain_save[:,0] = datain[:,-1]
                        datain = datain_save

        # 1-d data.
        elif lonsin.ndim == 1:
            nlons = len(lonsin)
            lonsin = np.where(lonsin > lon_0+180, lonsin-360 ,lonsin)
            lonsin = np.where(lonsin < lon_0-180, lonsin+360 ,lonsin)

            if nlons > 1:
                londiff = np.abs(lonsin[0:-1]-lonsin[1:])
                londiff_sort = np.sort(londiff)
                thresh = 360.-londiff_sort[-2] if nlons > 2 else 360.0 - londiff_sort[-1]
                itemindex = len(lonsin)-np.where(londiff>=thresh)[0]
            else:
                itemindex = 0

            if fix_wrap_around and itemindex:
                # check to see if cyclic (wraparound) point included
                # if so, remove it.
                if np.abs(lonsin[0]-lonsin[-1]) < 1.e-4:
                    hascyclic = True
                    lonsin_save = lonsin.copy()
                    lonsin = lonsin[1:]
                    if datain is not None:
                        datain_save = datain.copy()
                        datain = datain[1:]
                else:
                    hascyclic = False
                lonsin = np.roll(lonsin,itemindex-1)
                if datain is not None:
                    datain = np.roll(datain,itemindex-1)
                # add cyclic point back at beginning.
                if hascyclic:
                    lonsin_save[1:] = lonsin
                    lonsin_save[0] = lonsin[-1]-360.
                    lonsin = lonsin_save
                    if datain is not None:
                        datain_save[1:] = datain
                        datain_save[0] = datain[-1]
                        datain = datain_save

        # mask points outside
        # map region so they don't wrap back in the domain.
        mask = np.logical_or(lonsin<lon_0-180,lonsin>lon_0+180)
        lonsin = np.where(mask,1.e30,lonsin)
        if datain is not None and mask.any():
            datain = ma.masked_where(mask, datain)

        if datain is not None:
            return lonsin, datain
        else:
            return lonsin

### End of Basemap class

def _searchlist(a,x):
    """
    like bisect, but works for lists that are not sorted,
    and are not in increasing order.
    returns -1 if x does not fall between any two elements"""
    # make sure x is a float (and not an array scalar)
    x = float(x)
    itemprev = a[0]
    nslot = -1
    eps = 180.
    for n,item in enumerate(a[1:]):
        if item < itemprev:
            if itemprev-item>eps:
                if ((x>itemprev and x<=360.) or (x<item and x>=0.)):
                    nslot = n+1
                    break
            elif x <= itemprev and x > item and itemprev:
                nslot = n+1
                break
        else:
            if item-itemprev>eps:
                if ((x<itemprev and x>=0.) or (x>item and x<=360.)):
                    nslot = n+1
                    break
            elif x >= itemprev and x < item:
                nslot = n+1
                break
        itemprev = item
    return nslot

def interp(datain,xin,yin,xout,yout,checkbounds=False,masked=False,order=1):
    """
    Interpolate data (``datain``) on a rectilinear grid (with x = ``xin``
    y = ``yin``) to a grid with x = ``xout``, y= ``yout``.

    .. tabularcolumns:: |l|L|

    ==============   ====================================================
    Arguments        Description
    ==============   ====================================================
    datain           a rank-2 array with 1st dimension corresponding to
                     y, 2nd dimension x.
    xin, yin         rank-1 arrays containing x and y of
                     datain grid in increasing order.
    xout, yout       rank-2 arrays containing x and y of desired output grid.
    ==============   ====================================================

    .. tabularcolumns:: |l|L|

    ==============   ====================================================
    Keywords         Description
    ==============   ====================================================
    checkbounds      If True, values of xout and yout are checked to see
                     that they lie within the range specified by xin
                     and xin.
                     If False, and xout,yout are outside xin,yin,
                     interpolated values will be clipped to values on
                     boundary of input grid (xin,yin)
                     Default is False.
    masked           If True, points outside the range of xin and yin
                     are masked (in a masked array).
                     If masked is set to a number, then
                     points outside the range of xin and yin will be
                     set to that number. Default False.
    order            0 for nearest-neighbor interpolation, 1 for
                     bilinear interpolation, 3 for cublic spline
                     (default 1). order=3 requires scipy.ndimage.
    ==============   ====================================================

    .. note::
     If datain is a masked array and order=1 (bilinear interpolation) is
     used, elements of dataout will be masked if any of the four surrounding
     points in datain are masked.  To avoid this, do the interpolation in two
     passes, first with order=1 (producing dataout1), then with order=0
     (producing dataout2).  Then replace all the masked values in dataout1
     with the corresponding elements in dataout2 (using numpy.where).
     This effectively uses nearest neighbor interpolation if any of the
     four surrounding points in datain are masked, and bilinear interpolation
     otherwise.

    Returns ``dataout``, the interpolated data on the grid ``xout, yout``.
    """
    # xin and yin must be monotonically increasing.
    if xin[-1]-xin[0] < 0 or yin[-1]-yin[0] < 0:
        raise ValueError('xin and yin must be increasing!')
    if xout.shape != yout.shape:
        raise ValueError('xout and yout must have same shape!')
    # check that xout,yout are
    # within region defined by xin,yin.
    if checkbounds:
        if xout.min() < xin.min() or \
           xout.max() > xin.max() or \
           yout.min() < yin.min() or \
           yout.max() > yin.max():
            raise ValueError('yout or xout outside range of yin or xin')
    # compute grid coordinates of output grid.
    delx = xin[1:]-xin[0:-1]
    dely = yin[1:]-yin[0:-1]
    if max(delx)-min(delx) < 1.e-4 and max(dely)-min(dely) < 1.e-4:
        # regular input grid.
        xcoords = (len(xin)-1)*(xout-xin[0])/(xin[-1]-xin[0])
        ycoords = (len(yin)-1)*(yout-yin[0])/(yin[-1]-yin[0])
    else:
        # irregular (but still rectilinear) input grid.
        xoutflat = xout.flatten(); youtflat = yout.flatten()
        ix = (np.searchsorted(xin,xoutflat)-1).tolist()
        iy = (np.searchsorted(yin,youtflat)-1).tolist()
        xoutflat = xoutflat.tolist(); xin = xin.tolist()
        youtflat = youtflat.tolist(); yin = yin.tolist()
        xcoords = []; ycoords = []
        for n,i in enumerate(ix):
            if i < 0:
                xcoords.append(-1) # outside of range on xin (lower end)
            elif i >= len(xin)-1:
                xcoords.append(len(xin)) # outside range on upper end.
            else:
                xcoords.append(float(i)+(xoutflat[n]-xin[i])/(xin[i+1]-xin[i]))
        for m,j in enumerate(iy):
            if j < 0:
                ycoords.append(-1) # outside of range of yin (on lower end)
            elif j >= len(yin)-1:
                ycoords.append(len(yin)) # outside range on upper end
            else:
                ycoords.append(float(j)+(youtflat[m]-yin[j])/(yin[j+1]-yin[j]))
        xcoords = np.reshape(xcoords,xout.shape)
        ycoords = np.reshape(ycoords,yout.shape)
    # data outside range xin,yin will be clipped to
    # values on boundary.
    if masked:
        xmask = np.logical_or(np.less(xcoords,0),np.greater(xcoords,len(xin)-1))
        ymask = np.logical_or(np.less(ycoords,0),np.greater(ycoords,len(yin)-1))
        xymask = np.logical_or(xmask,ymask)
    xcoords = np.clip(xcoords,0,len(xin)-1)
    ycoords = np.clip(ycoords,0,len(yin)-1)
    # interpolate to output grid using bilinear interpolation.
    if order == 1:
        xi = xcoords.astype(np.int32)
        yi = ycoords.astype(np.int32)
        xip1 = xi+1
        yip1 = yi+1
        xip1 = np.clip(xip1,0,len(xin)-1)
        yip1 = np.clip(yip1,0,len(yin)-1)
        delx = xcoords-xi.astype(np.float32)
        dely = ycoords-yi.astype(np.float32)
        dataout = (1.-delx)*(1.-dely)*datain[yi,xi] + \
                  delx*dely*datain[yip1,xip1] + \
                  (1.-delx)*dely*datain[yip1,xi] + \
                  delx*(1.-dely)*datain[yi,xip1]
    elif order == 0:
        xcoordsi = np.around(xcoords).astype(np.int32)
        ycoordsi = np.around(ycoords).astype(np.int32)
        dataout = datain[ycoordsi,xcoordsi]
    elif order == 3:
        try:
            from scipy.ndimage import map_coordinates
        except ImportError:
            raise ValueError('scipy.ndimage must be installed if order=3')
        coords = [ycoords,xcoords]
        dataout = map_coordinates(datain,coords,order=3,mode='nearest')
    else:
        raise ValueError('order keyword must be 0, 1 or 3')
    if masked:
        newmask = ma.mask_or(ma.getmask(dataout), xymask)
        dataout = ma.masked_array(dataout, mask=newmask)
        if not isinstance(masked, bool):
            dataout = dataout.filled(masked)
    return dataout

def shiftgrid(lon0,datain,lonsin,start=True,cyclic=360.0):
    """
    Shift global lat/lon grid east or west.

    .. tabularcolumns:: |l|L|

    ==============   ====================================================
    Arguments        Description
    ==============   ====================================================
    lon0             starting longitude for shifted grid
                     (ending longitude if start=False). lon0 must be on
                     input grid (within the range of lonsin).
    datain           original data with longitude the right-most
                     dimension.
    lonsin           original longitudes.
    ==============   ====================================================

    .. tabularcolumns:: |l|L|

    ==============   ====================================================
    Keywords         Description
    ==============   ====================================================
    start            if True, lon0 represents the starting longitude
                     of the new grid. if False, lon0 is the ending
                     longitude. Default True.
    cyclic           width of periodic domain (default 360)
    ==============   ====================================================

    returns ``dataout,lonsout`` (data and longitudes on shifted grid).
    """
    if np.fabs(lonsin[-1]-lonsin[0]-cyclic) > 1.e-4:
        # Use all data instead of raise ValueError, 'cyclic point not included'
        start_idx = 0
    else:
        # If cyclic, remove the duplicate point
        start_idx = 1
    if lon0 < lonsin[0] or lon0 > lonsin[-1]:
        raise ValueError('lon0 outside of range of lonsin')
    i0 = np.argmin(np.fabs(lonsin-lon0))
    i0_shift = len(lonsin)-i0
    if ma.isMA(datain):
        dataout  = ma.zeros(datain.shape,datain.dtype)
    else:
        dataout  = np.zeros(datain.shape,datain.dtype)
    if ma.isMA(lonsin):
        lonsout = ma.zeros(lonsin.shape,lonsin.dtype)
    else:
        lonsout = np.zeros(lonsin.shape,lonsin.dtype)
    if start:
        lonsout[0:i0_shift] = lonsin[i0:]
    else:
        lonsout[0:i0_shift] = lonsin[i0:]-cyclic
    dataout[...,0:i0_shift] = datain[...,i0:]
    if start:
        lonsout[i0_shift:] = lonsin[start_idx:i0+start_idx]+cyclic
    else:
        lonsout[i0_shift:] = lonsin[start_idx:i0+start_idx]
    dataout[...,i0_shift:] = datain[...,start_idx:i0+start_idx]
    return dataout,lonsout

def addcyclic(*arr,**kwargs):
    """
    Adds cyclic (wraparound) points in longitude to one or several arrays,
    the last array being longitudes in degrees. e.g.

   ``data1out, data2out, lonsout = addcyclic(data1,data2,lons)``

    ==============   ====================================================
    Keywords         Description
    ==============   ====================================================
    axis             the dimension representing longitude (default -1,
                     or right-most)
    cyclic           width of periodic domain (default 360)
    ==============   ====================================================
    """
    # get (default) keyword arguments
    axis = kwargs.get('axis',-1)
    cyclic = kwargs.get('cyclic',360)
    # define functions
    def _addcyclic(a):
        """addcyclic function for a single data array"""
        npsel = np.ma if np.ma.is_masked(a) else np
        slicer = [slice(None)] * np.ndim(a)
        try:
            slicer[axis] = slice(0, 1)
        except IndexError:
            raise ValueError('The specified axis does not correspond to an '
                    'array dimension.')
        return npsel.concatenate((a,a[tuple(slicer)]),axis=axis)
    def _addcyclic_lon(a):
        """addcyclic function for a single longitude array"""
        # select the right numpy functions
        npsel = np.ma if np.ma.is_masked(a) else np
        # get cyclic longitudes
        clon = (np.take(a,[0],axis=axis)
                + cyclic * np.sign(np.diff(np.take(a,[0,-1],axis=axis),axis=axis)))
        # ensure the values do not exceed cyclic
        clonmod = npsel.where(clon<=cyclic,clon,np.mod(clon,cyclic))
        return npsel.concatenate((a,clonmod),axis=axis)
    # process array(s)
    if len(arr) == 1:
        return _addcyclic_lon(arr[-1])
    else:
        return list(map(_addcyclic,arr[:-1])) + [_addcyclic_lon(arr[-1])]

def _choosecorners(width,height,**kwargs):
    """
    private function to determine lat/lon values of projection region corners,
    given width and height of projection region in meters.
    """
    p = pyproj.Proj(kwargs)
    urcrnrlon, urcrnrlat = p(0.5*width,0.5*height, inverse=True)
    llcrnrlon, llcrnrlat = p(-0.5*width,-0.5*height, inverse=True)
    corners = llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat
    # test for invalid projection points on output
    if llcrnrlon > 1.e20 or urcrnrlon > 1.e20:
        raise ValueError('width and/or height too large for this projection, try smaller values')
    else:
        return corners

def _choosecornersllur(llcrnrx, llcrnry, urcrnrx, urcrnry,**kwargs):
    """
    private function to determine lat/lon values of projection region corners,
    given width and height of projection region in meters.
    """
    p = pyproj.Proj(kwargs)
    urcrnrlon, urcrnrlat = p(urcrnrx, urcrnry, inverse=True)
    llcrnrlon, llcrnrlat = p(llcrnrx, llcrnry, inverse=True)
    corners = llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat
    # test for invalid projection points on output
    if llcrnrlon > 1.e20 or urcrnrlon > 1.e20:
        raise ValueError('width and/or height too large for this projection, try smaller values')
    else:
        return corners

def maskoceans(lonsin,latsin,datain,inlands=True,resolution='l',grid=5):
    """
    mask data (``datain``), defined on a grid with latitudes ``latsin``
    longitudes ``lonsin`` so that points over water will not be plotted.

    .. tabularcolumns:: |l|L|

    ==============   ====================================================
    Arguments        Description
    ==============   ====================================================
    lonsin, latsin   rank-2 arrays containing longitudes and latitudes of
                     grid.
    datain           rank-2 input array on grid defined by ``lonsin`` and
                     ``latsin``.
    inlands          if False, masked only ocean points and not inland
                     lakes (Default True).
    resolution       gshhs coastline resolution used to define land/sea
                     mask (default 'l', available 'c','l','i','h' or 'f')
    grid             land/sea mask grid spacing in minutes (Default 5;
                     10, 2.5 and 1.25 are also available).
    ==============   ====================================================

    returns a masked array the same shape as datain with "wet" points masked.
    """
    # read in land/sea mask.
    lsmask_lons, lsmask_lats, lsmask =\
    _readlsmask(lakes=inlands,resolution=resolution,grid=grid)
    # nearest-neighbor interpolation to output grid.
    lsmasko = interp(lsmask,lsmask_lons,lsmask_lats,lonsin,latsin,masked=True,order=0)
    # mask input data.
    mask = lsmasko == 0
    return ma.masked_array(datain,mask=mask)

def _readlsmask(lakes=True,resolution='l',grid=5):
    # read in land/sea mask.
    if grid == 10:
        nlons = 2160
    elif grid == 5:
        nlons = 4320
    elif grid == 2.5:
        nlons = 8640
    elif grid == 1.25:
        nlons = 17280
    else:
        raise ValueError('grid for land/sea mask must be 10,5,2.5 or 1.25')
    nlats = nlons//2
    import gzip
    lsmaskf =\
    gzip.open(os.path.join(basemap_datadir,'lsmask_%smin_%s.bin' %\
        (grid,resolution)), 'rb')
    lsmask =\
    np.reshape(np.frombuffer(lsmaskf.read(),dtype=np.uint8),(nlats,nlons))
    if lakes:
        lsmask =\
        np.where(lsmask==2,np.array(0,dtype=np.uint8),lsmask)
    lsmaskf.close()
    delta = 360./nlons
    lsmask_lons = np.linspace(-180+0.5*delta,180-0.5*delta,nlons).astype(np.float32)
    lsmask_lats = np.linspace(-90+0.5*delta,90-0.5*delta,nlats).astype(np.float32)
    return lsmask_lons, lsmask_lats, lsmask

class _tup(tuple):
    # tuple with an added remove method.
    # used for objects returned by drawparallels and drawmeridians.
    def remove(self):
        for item in self:
            for x in item:
                x.remove()
class _dict(dict):
    # override __delitem__ to first call remove method on values.
    def __delitem__(self,key):
        self[key].remove()
        super(_dict, self).__delitem__(key)

def _setlonlab(fmt,lon,labelstyle):
    # set lon label string (called by Basemap.drawmeridians)
    try: # fmt is a function that returns a formatted string
        lonlab = fmt(lon)
    except: # fmt is a format string.
        if lon>180:
            if rcParams['text.usetex']:
                if labelstyle=='+/-':
                    lonlabstr = r'${\/-%s\/^{\circ}}$'%fmt
                else:
                    lonlabstr = r'${%s\/^{\circ}\/W}$'%fmt
            else:
                if labelstyle=='+/-':
                    lonlabstr = u'-%s\N{DEGREE SIGN}'%fmt
                else:
                    lonlabstr = u'%s\N{DEGREE SIGN}W'%fmt
            lonlab = lonlabstr%np.fabs(lon-360)
        elif lon<180 and lon != 0:
            if rcParams['text.usetex']:
                if labelstyle=='+/-':
                    lonlabstr = r'${\/+%s\/^{\circ}}$'%fmt
                else:
                    lonlabstr = r'${%s\/^{\circ}\/E}$'%fmt
            else:
                if labelstyle=='+/-':
                    lonlabstr = u'+%s\N{DEGREE SIGN}'%fmt
                else:
                    lonlabstr = u'%s\N{DEGREE SIGN}E'%fmt
            lonlab = lonlabstr%lon
        else:
            if rcParams['text.usetex']:
                lonlabstr = r'${%s\/^{\circ}}$'%fmt
            else:
                lonlabstr = u'%s\N{DEGREE SIGN}'%fmt
            lonlab = lonlabstr%lon
    return lonlab

def _setlatlab(fmt,lat,labelstyle):
    # set lat label string (called by Basemap.drawparallels)
    try: # fmt is a function that returns a formatted string
           latlab = fmt(lat)
    except: # fmt is a format string.
        if lat<0:
            if rcParams['text.usetex']:
                if labelstyle=='+/-':
                    latlabstr = r'${\/-%s\/^{\circ}}$'%fmt
                else:
                    latlabstr = r'${%s\/^{\circ}\/S}$'%fmt
            else:
                if labelstyle=='+/-':
                    latlabstr = u'-%s\N{DEGREE SIGN}'%fmt
                else:
                    latlabstr = u'%s\N{DEGREE SIGN}S'%fmt
            latlab = latlabstr%np.fabs(lat)
        elif lat>0:
            if rcParams['text.usetex']:
                if labelstyle=='+/-':
                    latlabstr = r'${\/+%s\/^{\circ}}$'%fmt
                else:
                    latlabstr = r'${%s\/^{\circ}\/N}$'%fmt
            else:
                if labelstyle=='+/-':
                    latlabstr = u'+%s\N{DEGREE SIGN}'%fmt
                else:
                    latlabstr = u'%s\N{DEGREE SIGN}N'%fmt
            latlab = latlabstr%lat
        else:
            if rcParams['text.usetex']:
                latlabstr = r'${%s\/^{\circ}}$'%fmt
            else:
                latlabstr = u'%s\N{DEGREE SIGN}'%fmt
            latlab = latlabstr%lat
    return latlab
