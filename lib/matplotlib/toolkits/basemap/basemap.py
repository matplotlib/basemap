from matplotlib import rcParams
from matplotlib import __version__ as matplotlib_version
# check to make sure matplotlib is not too old.
if matplotlib_version < '0.84':
    raise ImportError, 'your matplotlib is too old - basemap requires at least 0.84'
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
import sys, os, math, popen2
from proj import Proj
from greatcircle import GreatCircle, vinc_dist, vinc_pt
import matplotlib.numerix as NX
from matplotlib.numerix import ma
from matplotlib.mlab import linspace
from matplotlib.numerix.mlab import squeeze
from matplotlib.cbook import popd
from shapelib import ShapeFile
import dbflib

# look in sys.prefix for directory containing basemap files if
# BASEMAP_DATA_PATH env var not set.
_datadir = os.environ.get('BASEMAP_DATA_PATH')
if not _datadir:
   _datadir = os.path.join(sys.prefix,'share/basemap')

__version__ = '0.7.2'
__revision__ = '20051002'

# make sure subplots have same width and height
# so that figure size sets plot aspect ratio.
# (these are the same as matplotlib 0.84 defaults
#  except for subplot.left, which was changed from 0.125).
# print a warning when default is changed.
if rcParams['figure.subplot.left']!=0.1:
    print 'warning: figure.subplot.left rc value being reset to 0.1 in basemap'
    print 'use rcdefaults() to get the original value back'
    rcParams['figure.subplot.left']=0.1
if rcParams['figure.subplot.right']!=0.9 :
    print 'warning: figure.subplot.right rc value being reset to 0.9 in basemap'
    print 'use rcdefaults() to get the original value back'
    rcParams['figure.subplot.right']=0.9
if rcParams['figure.subplot.bottom']!=0.1:
    print 'warning: figure.subplot.bottom rc value being reset to 0.1 in basemap'
    print 'use rcdefaults() to get the original value back'
    rcParams['figure.subplot.bottom']=0.1
if rcParams['figure.subplot.top']!=0.9:
    print 'warning: figure.subplot.top rc value being reset to 0.9 in basemap'
    print 'use rcdefaults() to get the original value back'
    rcParams['figure.subplot.top']=0.9

class Basemap:

    """
 Set up a basemap with one of 17 supported map projections
 (cylindrical equidistant, mercator, polyconic, oblique mercator,
 transverse mercator, miller cylindrical, lambert conformal conic,
 azimuthal equidistant, equidistant conic, lambert azimuthal equal area,
 albers equal area conic, gnomonic, orthographic, mollweide,
 robinson, cassini-soldner or stereographic).
 Doesn't actually draw anything, but sets up the map projection class and
 creates the coastline, lake river and political boundary data
 structures in native map projection coordinates.
 Uses a pyrex interface to C-code from proj.4 (http://proj.maptools.org).

 Useful instance variables:

 projection - map projection ('cyl','merc','mill','lcc','eqdc','aea','aeqd',
  'laea', 'tmerc', 'omerc', 'cass', 'gnom', 'poly', 'ortho', 'robin',
  'moll' or 'stere')
 aspect - map aspect ratio (size of y dimension / size of x dimension).
 llcrnrlon - longitude of lower left hand corner of the desired map domain.
 llcrnrlon - latitude of lower left hand corner of the desired map domain.
 urcrnrlon - longitude of upper right hand corner of the desired map domain.
 urcrnrlon - latitude of upper right hand corner of the desired map domain.
 llcrnrx,llcrnry,urcrnrx,urcrnry - corners of map domain in projection coordinates.
 rmajor,rminor - equatorial and polar radii of ellipsoid used (in meters).

 Example Usage:

>>> from matplotlib.toolkits.basemap import Basemap
>>> import cPickle
>>> from pylab import *
>>> # read in topo data from pickle (on a regular lat/lon grid)
>>> topodict = cPickle.load(open('etopo20.pickle','rb'))
>>> etopo = topodict['data']; lons = topodict['lons']; lats = topodict['lats']
>>> # create Basemap instance for Robinson projection.
>>> m = Basemap(projection='robin',lon_0=0.5*(lons[0]+lons[-1]))
>>> # compute native map projection coordinates for lat/lon grid.
>>> x, y = m(*meshgrid(lons,lats))
>>> # create figure with same aspect ratio as map.
>>> fig = m.createfigure()
>>> # make filled contour plot.
>>> cs = m.contourf(x,y,etopo,30,cmap=cm.jet)
>>> m.drawcoastlines() # draw coastlines
>>> m.drawmapboundary() # draw a line around the map region
>>> m.drawparallels(arange(-90.,120.,30.),labels=[1,0,0,0]) # draw parallels
>>> m.drawmeridians(arange(0.,420.,60.),labels=[0,0,0,1]) # draw meridians
>>> title('Robinson Projection') # add a title
>>> show()

 [this example (simpletest.py) plus many others can be found in the
  examples directory of source distribution.  The "OO" version of this
  example (which does not use pylab) is called "simpletest_oo.py".]

 Contact: Jeff Whitaker <jeffrey.s.whitaker@noaa.gov>
    """

    def __init__(self,llcrnrlon=-180.,llcrnrlat=-90.,urcrnrlon=180.,urcrnrlat=90.,\
       projection='cyl',resolution='c',area_thresh=None,rsphere=6370997,\
       lat_ts=None,lat_1=None,lat_2=None,lat_0=None,lon_0=None,\
       lon_1=None,lon_2=None,suppress_ticks=True,ax=None):
        """
 create a Basemap instance.

 arguments:

 projection - map projection.  'cyl' - cylindrical equidistant, 'merc' -
  mercator, 'lcc' - lambert conformal conic, 'stere' - stereographic,
  'aea' - albers equal area conic, 'tmerc' - transverse mercator,
  'aeqd' - azimuthal equidistant, 'mill' - miller cylindrical,
  'eqdc' - equidistant conic, 'laea' - lambert azimuthal equal area,
  'cass' - cassini-soldner (transverse cylindrical equidistant),
  'poly' - polyconic, 'omerc' - oblique mercator, 'ortho' - orthographic,
  'moll' - mollweide, 'robin' - robinson,
  and 'gnom' - gnomonic are currently available.  Default 'cyl'.

 llcrnrlon - longitude of lower left hand corner of the desired map domain
  (Default -180).
 llcrnrlat - latitude of lower left hand corner of the desired map domain
  (Default -90).
 urcrnrlon - longitude of upper right hand corner of the desired map domain
  (Default 180).
 urcrnrlat - latitude of upper right hand corner of the desired map domain
  (Default 90).

 If the orthographic, mollweide or robinson projection is chosen
 the values of llcrnrlon,llcrnrlat,urcrnrlon and urcrnrlat are ignored,
 and the entire projection domain will be always be plotted.

 resolution - resolution of boundary database to use. Can be 'c' (crude),
  'l' (low), 'i' (intermediate) or 'h' (high).
  Resolution drops off by roughly 80%
  between datasets.  Higher res datasets are much slower to draw.
  Default 'c'. Coastline data is from the GSHHS
  (http://www.soest.hawaii.edu/wessel/gshhs/gshhs.html).
  State, country and river datasets from the Generic Mapping
  Tools ((http://gmt.soest.hawaii.edu).

 area_thresh - coastline or lake with an area smaller than area_thresh
  in km^2 will not be plotted.  Default 10000,1000,100,10 for resolution
  'c','l','i','h'.

 rsphere - radius of the sphere used to define map projection (default
  6370997 meters, close to the arithmetic mean radius of the earth). If
  given as a sequence, the first two elements are interpreted as
  the the radii of the major and minor axes of an ellipsoid. Note: sometimes
  an ellipsoid is specified by the major axis and an 'inverse flattening
  parameter' (if).  The minor axis (b) can be computed from the major axis (a)
  and the inverse flattening parameter using the formula if = a/(a-b).

 suppress_ticks - suppress automatic drawing of axis ticks and labels
 in map projection coordinates.  Default False, so parallels and meridians
 can be labelled instead. If parallel or meridian labelling is requested
 (using drawparallels and drawmeridians methods), automatic tick labelling
 will be supressed even is suppress_ticks=False.  Typically, you will
 only want to override the default if you want to label the axes in meters
 using native map projection coordinates.

 ax - set default axes instance (default None - pylab.gca() may be used
 to get the current axes instance).  If you don't want pylab to be imported,
 you can either set this to a pre-defined axes instance, or use the 'ax'
 keyword in each Basemap method call that does drawing. In the first case,
 all Basemap method calls will draw to the same axes instance.  In the
 second case, you can draw to different axes with the same Basemap instance.
 You can also use the 'ax' keyword in individual method calls to
 selectively override the default axes instance.

 The following parameters are map projection parameters which all default to
 None.  Not all parameters are used by all projections, some are ignored.

 lat_ts - latitude of natural origin (used for mercator and stereographic
  projections).
 lat_1 - first standard parallel for lambert conformal, albers
  equal area projection and equidistant conic projections. Latitude of one
  of the two points on the projection centerline for oblique mercator.
 lat_2 - second standard parallel for lambert conformal, albers
  equal area projection and equidistant conic projections. Latitude of one
  of the two points on the projection centerline for oblique mercator.
 lon_1 - Longitude of one of the two points on the projection centerline
  for oblique mercator.
 lon_2 - Longitude of one of the two points on the projection centerline
  for oblique mercator.
 lat_0 - central latitude (y-axis origin) - used by stereographic, polyconic,
  transverse mercator, miller cylindrical, cassini-soldner, oblique mercator,
  gnomonic, equidistant conic, orthographic and lambert azimuthal projections).
 lon_0 - central meridian (x-axis origin) - used by stereographic, polyconic,
  transverse mercator, miller cylindrical, cassini-soldner, mollweide, robinson,
  gnomonic, equidistant conic, orthographic and lambert azimuthal projections).
        """

        self.projection = projection
        # make sure lat/lon limits are converted to floats.
        self.llcrnrlon = float(llcrnrlon)
        self.llcrnrlat = float(llcrnrlat)
        self.urcrnrlon = float(urcrnrlon)
        self.urcrnrlat = float(urcrnrlat)
        # set min/max lats for projection domain.
        # only boundary segments within that range will be processed.
        self.latmax = None
        if projection in ['mill','cyl','merc']:
            # if self.crossgreenwich is False (projection
            # doesn't cross greenwich meridian) some time can
            # be saved in processing boundaries.
            if (self.llcrnrlon+360.)%360 >= (self.urcrnrlon+360.)%360:
                self.crossgreenwich = True
            else:
                self.crossgreenwich = False
                # if not crossgreenwich, adjust lon bounds
                # to be in range 0 to 360.
                self.urcrnrlon = (self.urcrnrlon+360.)%360
                self.llcrnrlon = (self.llcrnrlon+360.)%360
            self.latmin = self.llcrnrlat
            self.latmax = self.urcrnrlat
        elif projection in ['ortho','moll','robin']:
            self.latmin = -90.
            self.latmax = 90.
            self.crossgreenwich = True

        # set up projection parameter dict.
        # Create Proj class instance.
        projparams = {}
        projparams['proj'] = projection
        try:
            if rsphere[0] > rsphere[1]:
                projparams['a'] = rsphere[0]
                projparams['b'] = rsphere[1]
            else:
                projparams['a'] = rsphere[1]
                projparams['b'] = rsphere[0]
        except:
            projparams['R'] = rsphere

        if projection == 'lcc':
            if lat_1 is None or lon_0 is None:
                raise ValueError, 'must specify lat_1 and lon_0 for Lambert Conformal basemap'
            projparams['lat_1'] = lat_1
            if lat_2 != None:
                 projparams['lat_2'] = lat_2
            projparams['lon_0'] = lon_0
            proj = Proj(projparams,self.llcrnrlon,self.llcrnrlat,self.urcrnrlon,self.urcrnrlat)
        elif projection == 'eqdc':
            if lat_1 is None or lat_2 is None or lon_0 is None:
                raise ValueError, 'must specify lat_1, lat_2 and lon_0 for Equidistant Conic basemap'
            projparams['lat_1'] = lat_1
            projparams['lat_2'] = lat_2
            projparams['lon_0'] = lon_0
            proj = Proj(projparams,self.llcrnrlon,self.llcrnrlat,self.urcrnrlon,self.urcrnrlat)
        elif projection == 'aea':
            if lat_1 is None or lat_2 is None or lon_0 is None:
                raise ValueError, 'must specify lat_1, lat_2 and lon_0 for Albers Equal Area basemap'
            projparams['lat_1'] = lat_1
            projparams['lat_2'] = lat_2
            projparams['lon_0'] = lon_0
            proj = Proj(projparams,self.llcrnrlon,self.llcrnrlat,self.urcrnrlon,self.urcrnrlat)
        elif projection == 'stere':
            if lat_ts is None or lat_0 is None or lon_0 is None:
                raise ValueError, 'must specify lat_ts,lat_0 and lon_0 for Stereographic basemap'
            projparams['lat_ts'] = lat_ts
            projparams['lat_0'] = lat_0
            projparams['lon_0'] = lon_0
            proj = Proj(projparams,self.llcrnrlon,self.llcrnrlat,self.urcrnrlon,self.urcrnrlat)
        elif projection == 'laea':
            if lat_0 is None or lon_0 is None:
                raise ValueError, 'must specify lat_0 and lon_0 for Lambert Azimuthal basemap'
            projparams['lat_0'] = lat_0
            projparams['lon_0'] = lon_0
            proj = Proj(projparams,self.llcrnrlon,self.llcrnrlat,self.urcrnrlon,self.urcrnrlat)
        elif projection == 'merc':
            if lat_ts is None:
                raise ValueError, 'must specify lat_ts for Mercator basemap'
            projparams['lat_ts'] = lat_ts
            proj = Proj(projparams,self.llcrnrlon,self.llcrnrlat,self.urcrnrlon,self.urcrnrlat)
        elif projection in ['tmerc','gnom','cass','poly','ortho'] :
            if lat_0 is None or lon_0 is None:
                raise ValueError, 'must specify lat_0 and lon_0 for Transverse Mercator, Gnomonic, Cassini-Soldner, Orthographic or Polyconic basemap'
            projparams['lat_0'] = lat_0
            projparams['lon_0'] = lon_0
            proj = Proj(projparams,self.llcrnrlon,self.llcrnrlat,self.urcrnrlon,self.urcrnrlat)
        elif projection == 'moll' or projection == 'robin':
            if lon_0 is None:
                raise ValueError, 'must specify lon_0 for Robinson or Mollweide basemap'
            projparams['lon_0'] = lon_0
            proj = Proj(projparams,self.llcrnrlon,self.llcrnrlat,self.urcrnrlon,self.urcrnrlat)
        elif projection == 'omerc':
            if lat_1 is None or lon_1 is None or lat_2 is None or lon_2 is None:
                raise ValueError, 'must specify lat_1,lon_1 and lat_2,lon_2 for Oblique Mercator basemap'
            projparams['lat_1'] = lat_1
            projparams['lon_1'] = lon_1
            projparams['lat_2'] = lat_2
            projparams['lon_2'] = lon_2
            proj = Proj(projparams,self.llcrnrlon,self.llcrnrlat,self.urcrnrlon,self.urcrnrlat)
        elif projection == 'aeqd':
            if lat_0 is None or lon_0 is None:
                raise ValueError, 'must specify lat_0 and lon_0 for Azimuthal Equidistant basemap'
            projparams['lat_0'] = lat_0
            projparams['lon_0'] = lon_0
            proj = Proj(projparams,self.llcrnrlon,self.llcrnrlat,self.urcrnrlon,self.urcrnrlat)
        elif projection == 'mill':
            if lat_0 is not None:
                projparams['lat_0'] = lat_0
            if lon_0 is not None:
                projparams['lon_0'] = lon_0
            proj = Proj(projparams,self.llcrnrlon,self.llcrnrlat,self.urcrnrlon,self.urcrnrlat)
        elif projection == 'cyl':
            proj = Proj(projparams,self.llcrnrlon,self.llcrnrlat,self.urcrnrlon,self.urcrnrlat)
        else:
            raise ValueError, 'unsupported projection'

        # make sure axis ticks are suppressed.
        self.noticks = suppress_ticks

        # make Proj instance a Basemap instance variable.
        self.projtran = proj
        # copy some Proj attributes.
        atts = ['rmajor','rminor','esq','flattening','ellipsoid','projparams']
        for att in atts:
            self.__dict__[att] = proj.__dict__[att]
        # set instance variables defining map region.
        self.xmin = proj.xmin
        self.xmax = proj.xmax
        self.ymin = proj.ymin
        self.ymax = proj.ymax
        if projection == 'cyl':
            self.aspect = (self.urcrnrlat-self.llcrnrlat)/(self.urcrnrlon-self.llcrnrlon)
        else:
            self.aspect = (proj.ymax-proj.ymin)/(proj.xmax-proj.xmin)
        self.llcrnrx = proj.llcrnrx
        self.llcrnry = proj.llcrnry
        self.urcrnrx = proj.urcrnrx
        self.urcrnry = proj.urcrnry

        # find lat/lon bounds if it hasn't been done already.
        if self.latmax is None:
            lons, lats = self.makegrid(101,101)
            self.latmin = min(NX.ravel(lats))
            self.latmax = max(NX.ravel(lats))
            # too hard to figure out, so just assume it does.
            self.crossgreenwich = True

        # set defaults for area_thresh.
        if area_thresh is None:
            if resolution == 'c':
                area_thresh = 10000.
            elif resolution == 'l':
                area_thresh = 1000.
            elif resolution == 'i':
                area_thresh = 100.
            elif resolution == 'h':
                area_thresh = 10.
            else:
                raise ValueError, "bounday resolution must be one of 'c','l','i' or 'h'"
        # if ax == None, pylab.gca may be used.
        self.ax = ax
        # read in coastline data (only those polygons whose area > area_thresh).
        coastlons = []; coastlats = []; coastsegind = []; coastsegtype = []
        a = not self.crossgreenwich
        for line in open(os.path.join(_datadir,'gshhs_'+resolution+'.txt')):
            linesplit = line.split()
            if line.startswith('P'):
                area = float(linesplit[5])
                west,east,south,north = float(linesplit[6]),float(linesplit[7]),float(linesplit[8]),float(linesplit[9])
                type = int(linesplit[3])
                useit = self.latmax>=south and self.latmin<=north and area>area_thresh
                if useit:
                    coastsegind.append(len(coastlons))
                    coastsegtype.append(type)
                continue
            # lon/lat
            if useit:
                lon, lat = [float(val) for val in linesplit]
                coastlons.append(lon)
                coastlats.append(lat)
        coastsegtype.append(type)
        coastsegind.append(len(coastlons))

        # read in country boundary data.
        cntrylons = []; cntrylats = []; cntrysegind = []
        for line in open(os.path.join(_datadir,'countries_'+resolution+'.txt')):
            linesplit = line.split()
            if line.startswith('>'):
                west,east,south,north = float(linesplit[7]),float(linesplit[8]),float(linesplit[9]),float(linesplit[10])
                useit = self.latmax>=south and self.latmin<=north
                if useit: cntrysegind.append(len(cntrylons))
                continue
            # lon/lat
            if useit:
                lon, lat = [float(val) for val in linesplit]
                cntrylons.append(lon)
                cntrylats.append(lat)
        cntrysegind.append(len(cntrylons))

        # read in state boundaries (Americas only).
        statelons = []; statelats = []; statesegind = []
        for line in open(os.path.join(_datadir,'states_'+resolution+'.txt')):
            linesplit = line.split()
            if line.startswith('>'):
                west,east,south,north = float(linesplit[7]),float(linesplit[8]),float(linesplit[9]),float(linesplit[10])
                useit = self.latmax>=south and self.latmin<=north
                if useit: statesegind.append(len(statelons))
                continue
            # lon/lat
            if useit:
                lon, lat = [float(val) for val in linesplit]
                statelons.append(lon)
                statelats.append(lat)
        statesegind.append(len(statelons))

        # read in major rivers.
        riverlons = []; riverlats = []; riversegind = []
        for line in open(os.path.join(_datadir,'rivers_'+resolution+'.txt')):
            linesplit = line.split()
            if line.startswith('>'):
                west,east,south,north = float(linesplit[7]),float(linesplit[8]),float(linesplit[9]),float(linesplit[10])
                useit = self.latmax>=south and self.latmin<=north
                if useit: riversegind.append(len(riverlons))
                continue
            # lon/lat
            if useit:
                lon, lat = [float(val) for val in linesplit]
                riverlons.append(lon)
                riverlats.append(lat)
        riversegind.append(len(riverlons))

        # if projection straddles Greenwich,
        # extend longitudes around the earth a second time
        # so valid longitudes can range from -360 to 720.
        if self.crossgreenwich:
            coastlons2 = [lon+360. for lon in coastlons]
            cntrylons2 = [lon+360. for lon in cntrylons]
            statelons2 = [lon+360. for lon in statelons]
            riverlons2 = [lon+360. for lon in riverlons]
            coastlons3 = [lon-360. for lon in coastlons]
            cntrylons3 = [lon-360. for lon in cntrylons]
            statelons3 = [lon-360. for lon in statelons]
            riverlons3 = [lon-360. for lon in riverlons]

        # transform coastline polygons to native map coordinates.
        xc,yc = proj(NX.array(coastlons),NX.array(coastlats))
        if self.crossgreenwich:
           xc2,yc2 = proj(NX.array(coastlons2),NX.array(coastlats))
           xc3,yc3 = proj(NX.array(coastlons3),NX.array(coastlats))
        if self.crossgreenwich and projection == 'merc' or projection == 'mill':
            yc2 = yc
            yc3 = yc

        # set up segments in form needed for LineCollection,
        # ignoring 'inf' values that are off the map.
        segments = [zip(xc[i0:i1],yc[i0:i1]) for i0,i1 in zip(coastsegind[:-1],coastsegind[1:])]
        segmentsll = [zip(coastlons[i0:i1],coastlats[i0:i1]) for i0,i1 in zip(coastsegind[:-1],coastsegind[1:])]
        segtypes = [i for i in coastsegtype[:-1]]
        if self.crossgreenwich:
            segments2 = [zip(xc2[i0:i1],yc2[i0:i1]) for i0,i1 in zip(coastsegind[:-1],coastsegind[1:]) if max(xc2[i0:i1]) < 1.e20 and max(yc2[i0:i1]) < 1.e20]
            segmentsll2 = [zip(coastlons2[i0:i1],coastlats[i0:i1]) for i0,i1 in zip(coastsegind[:-1],coastsegind[1:]) if max(xc2[i0:i1]) < 1.e20 and max(yc2[i0:i1]) < 1.e20]
            segtypes2 = [i for i0,i1,i in zip(coastsegind[:-1],coastsegind[1:],coastsegtype[:-1]) if max(xc2[i0:i1]) < 1.e20 and max(yc2[i0:i1]) < 1.e20]
            segments3 = [zip(xc3[i0:i1],yc3[i0:i1]) for i0,i1 in zip(coastsegind[:-1],coastsegind[1:]) if max(xc3[i0:i1]) < 1.e20 and max(yc3[i0:i1]) < 1.e20]
            segmentsll3 = [zip(coastlons3[i0:i1],coastlats[i0:i1]) for i0,i1 in zip(coastsegind[:-1],coastsegind[1:]) if max(xc3[i0:i1]) < 1.e20 and max(yc3[i0:i1]) < 1.e20]
            segtypes3 = [i for i0,i1,i in zip(coastsegind[:-1],coastsegind[1:],coastsegtype[:-1]) if max(xc3[i0:i1]) < 1.e20 and max(yc3[i0:i1]) < 1.e20]
            self.coastsegs = segments+segments2+segments3
            self.coastsegsll = segmentsll+segmentsll2+segmentsll3
            self.coastsegtypes = segtypes+segtypes2+segtypes3
        else:
            self.coastsegs = segments
            self.coastsegsll = segmentsll
            self.coastsegtypes = segtypes

        # same as above for country segments.
        xc,yc = proj(NX.array(cntrylons,'f'),NX.array(cntrylats))
        if self.crossgreenwich:
            xc2,yc2 = proj(NX.array(cntrylons2,'f'),NX.array(cntrylats))
            xc3,yc3 = proj(NX.array(cntrylons3),NX.array(cntrylats))
        if self.crossgreenwich and projection == 'merc' or projection == 'mill':
            yc2=yc
            yc3=yc
        segments = [zip(xc[i0:i1],yc[i0:i1]) for i0,i1 in zip(cntrysegind[:-1],cntrysegind[1:])]
        if self.crossgreenwich:
            segments2 = [zip(xc2[i0:i1],yc2[i0:i1]) for i0,i1 in zip(cntrysegind[:-1],cntrysegind[1:]) if max(xc2[i0:i1]) < 1.e20 and max(yc2[i0:i1]) < 1.e20]
            segments3 = [zip(xc3[i0:i1],yc3[i0:i1]) for i0,i1 in zip(cntrysegind[:-1],cntrysegind[1:]) if max(xc3[i0:i1]) < 1.e20 and max(yc3[i0:i1]) < 1.e20]
            self.cntrysegs = segments+segments2+segments3
        else:
            self.cntrysegs = segments

        # same as above for state segments.
        xc,yc = proj(NX.array(statelons,'f'),NX.array(statelats))
        if self.crossgreenwich:
            xc2,yc2 = proj(NX.array(statelons2,'f'),NX.array(statelats))
            xc3,yc3 = proj(NX.array(statelons3),NX.array(statelats))
        if self.crossgreenwich and projection == 'merc' or projection == 'mill':
            yc2=yc
            yc3=yc
        segments = [zip(xc[i0:i1],yc[i0:i1]) for i0,i1 in zip(statesegind[:-1],statesegind[1:])]
        if self.crossgreenwich:
            segments2 = [zip(xc2[i0:i1],yc2[i0:i1]) for i0,i1 in zip(statesegind[:-1],statesegind[1:]) if max(xc2[i0:i1]) < 1.e20 and max(yc2[i0:i1]) < 1.e20]
            segments3 = [zip(xc3[i0:i1],yc3[i0:i1]) for i0,i1 in zip(statesegind[:-1],statesegind[1:]) if max(xc3[i0:i1]) < 1.e20 and max(yc3[i0:i1]) < 1.e20]
            self.statesegs = segments+segments2+segments3
        else:
            self.statesegs = segments

        # same as above for river segments.
        xc,yc = proj(NX.array(riverlons,'f'),NX.array(riverlats))
        if self.crossgreenwich:
            xc2,yc2 = proj(NX.array(riverlons2,'f'),NX.array(riverlats))
            xc3,yc3 = proj(NX.array(riverlons3),NX.array(riverlats))
        if self.crossgreenwich and projection == 'merc' or projection == 'mill':
            yc2=yc
            yc3=yc
        segments = [zip(xc[i0:i1],yc[i0:i1]) for i0,i1 in zip(riversegind[:-1],riversegind[1:])]
        if self.crossgreenwich:
            segments2 = [zip(xc2[i0:i1],yc2[i0:i1]) for i0,i1 in zip(riversegind[:-1],riversegind[1:]) if max(xc2[i0:i1]) < 1.e20 and max(yc2[i0:i1]) < 1.e20]
            segments3 = [zip(xc3[i0:i1],yc3[i0:i1]) for i0,i1 in zip(riversegind[:-1],riversegind[1:]) if max(xc3[i0:i1]) < 1.e20 and max(yc3[i0:i1]) < 1.e20]
            self.riversegs = segments+segments2+segments3
        else:
            self.riversegs = segments

        # store coast polygons for filling.
        self.coastpolygons = []
        coastpolygonsll = []
        self.coastpolygontypes = []
        if projection in ['merc','mill']:
            xsp,ysp = proj(0.,-89.9) # s. pole coordinates.
            xa,ya = proj(0.,-68.) # edge of antarctica.
            x0,y0 = proj(0.,0.)
            xm360,ym360 = proj(-360.,0.)
            x360,y360 = proj(360.,0.)
            x720,y720 = proj(720.,0.)
        for seg,segtype,segll in zip(self.coastsegs,self.coastsegtypes,self.coastsegsll):
            x = [lon for lon,lat in seg]
            y = [lat for lon,lat in seg]
            lons = [lon for lon,lat in segll]
            lats = [lat for lon,lat in segll]
            # the antarctic polygon is a nuisance, since it
            # spans all longitudes, it's not closed and will not be filled
            # without some projection dependant tweaking.
            if projection == 'cyl':
                if x[-1] == 0.000 and y[-1] < -68.: # close antarctica
                    x.append(0.)
                    y.append(-90.0000)
                    x.insert(0,360.)
                    y.insert(0,-90)
                    lons.append(0.)
                    lats.append(-90.)
                    lons.insert(0,360.)
                    lats.insert(0,-90.)
                if x[-1] == 360.000 and y[-1] < -68.:
                    x.append(360.)
                    y.append(-90)
                    x.insert(0,720.)
                    y.insert(0,-90)
                    lons.append(360.)
                    lats.append(-90.)
                    lons.insert(0,720.)
                    lats.insert(0,-90.)
                if x[-1] == -360.000 and y[-1] < -68.:
                    x.append(-360.)
                    y.append(-90)
                    x.insert(0,0.)
                    y.insert(0,-90)
                    lons.append(-360.)
                    lats.append(-90.)
                    lons.insert(0,0.)
                    lats.insert(0,-90.)
            elif projection in ['merc','mill']:
                if math.fabs(x[-1]-x0) < 1. and y[-1] < ya: # close antarctica
                    x.append(x0)
                    y.append(ysp)
                    x.insert(0,x360)
                    y.insert(0,ysp)
                    lons.append(0.)
                    lats.append(-90.)
                    lons.insert(0,360.)
                    lats.insert(0,-90.)
                if math.fabs(x[-1]-x360) < 1. and y[-1] < ya:
                    x.append(x360)
                    y.append(ysp)
                    x.insert(0,x720)
                    y.insert(0,ysp)
                    lons.append(360.)
                    lats.append(-90.)
                    lons.insert(0,720.)
                    lats.insert(0,-90.)
                if math.fabs(x[-1]-xm360) < 1. and y[-1] < ya:
                    x.append(xm360)
                    y.append(ysp)
                    x.insert(0,x0)
                    y.insert(0,ysp)
                    lons.append(-360.)
                    lats.append(-90.)
                    lons.insert(0,0.)
                    lats.insert(0,-90.)
            self.coastpolygons.append((x,y))
            coastpolygonsll.append((lons,lats))
            self.coastpolygontypes.append(segtype)

        # remove those segments/polygons that don't intersect map region.
        coastsegs = []
        coastsegtypes = []
        for seg,segtype in zip(self.coastsegs,self.coastsegtypes):
            if self._insidemap_seg(seg):
                coastsegs.append(seg)
                coastsegtypes.append(segtype)
        self.coastsegs = coastsegs
        self.coastsegtypes = coastsegtypes
        polygons = []
        polygonsll = []
        polygontypes = []
        for poly,polytype,polyll in zip(self.coastpolygons,self.coastpolygontypes,coastpolygonsll):
            if self._insidemap_poly(poly,polyll):
                polygons.append(poly)
                polygontypes.append(polytype)
                polygonsll.append(polyll)
        self.coastpolygons = polygons
        coastpolygonsll = polygonsll
        self.coastpolygontypes = polygontypes
        states = []; rivers = []; countries = []
        for seg in self.cntrysegs:
            if self._insidemap_seg(seg):
                countries.append(seg)
        for seg in self.statesegs:
            if self._insidemap_seg(seg):
                states.append(seg)
        for seg in self.riversegs:
            if self._insidemap_seg(seg):
                rivers.append(seg)
        self.statesegs = states
        self.riversegs = rivers
        self.cntryegs = countries

        # split up segments that go outside projection limb
        coastsegs = []
        coastsegtypes = []
        for seg,segtype in zip(self.coastsegs,self.coastsegtypes):
            xx = NX.array([x for x,y in seg],'f')
            yy = NX.array([y for x,y in seg],'f')
            i1,i2 = self._splitseg(xx,yy)
            if i1 and i2:
                for i,j in zip(i1,i2):
                    segment = zip(xx[i:j],yy[i:j])
                    coastsegs.append(segment)
                    coastsegtypes.append(segtype)
            else:
                coastsegs.append(seg)
                coastsegs.append(segtype)
        self.coastsegs = coastsegs
        self.coastsegtypes = coastsegtypes
        states = []
        for seg in self.statesegs:
            xx = NX.array([x for x,y in seg],'f')
            yy = NX.array([y for x,y in seg],'f')
            i1,i2 = self._splitseg(xx,yy)
            if i1 and i2:
                for i,j in zip(i1,i2):
                    segment = zip(xx[i:j],yy[i:j])
                    states.append(segment)
            else:
                states.append(seg)
        self.statesegs = states
        countries = []
        for seg in self.cntrysegs:
            xx = NX.array([x for x,y in seg],'f')
            yy = NX.array([y for x,y in seg],'f')
            i1,i2 = self._splitseg(xx,yy)
            if i1 and i2:
                for i,j in zip(i1,i2):
                    segment = zip(xx[i:j],yy[i:j])
                    countries.append(segment)
            else:
                countries.append(seg)
        self.cntrysegs = countries
        rivers = []
        for seg in self.riversegs:
            xx = NX.array([x for x,y in seg],'f')
            yy = NX.array([y for x,y in seg],'f')
            i1,i2 = self._splitseg(xx,yy)
            if i1 and i2:
                for i,j in zip(i1,i2):
                    segment = zip(xx[i:j],yy[i:j])
                    rivers.append(segment)
            else:
                rivers.append(seg)
        self.riversegs = rivers

        # split coastline segments that jump across entire plot.
        coastsegs = []
        coastsegtypes = []
        for seg,segtype in zip(self.coastsegs,self.coastsegtypes):
            xx = NX.array([x for x,y in seg],'f')
            yy = NX.array([y for x,y in seg],'f')
            xd = (xx[1:]-xx[0:-1])**2
            yd = (yy[1:]-yy[0:-1])**2
            dist = NX.sqrt(xd+yd)
            split = dist > 5000000.
            if NX.sum(split) and self.projection not in ['merc','cyl','mill']:
               ind = (NX.compress(split,squeeze(split*NX.indices(xd.shape)))+1).tolist()
               iprev = 0
               ind.append(len(xd))
               for i in ind:
                   coastsegs.append(zip(xx[iprev:i],yy[iprev:i]))
                   coastsegtypes.append(segtype)
                   iprev = i
            else:
                coastsegs.append(seg)
                coastsegtypes.append(segtype)
        self.coastsegs = coastsegs
        self.coastsegtypes = coastsegtypes

        # special treatment of coastline polygons for
        # orthographic, mollweide and robinson.
        # (polygon clipping along projection limb)
        if self.projection == 'ortho':
            lat_0 = math.radians(self.projparams['lat_0'])
            lon_0 = math.radians(self.projparams['lon_0'])
            rad = (2.*self.rmajor + self.rminor)/3.
            dtheta = 0.01
            coastpolygons = []
            coastpolygontypes = []
            for poly,polytype,polyll in zip(self.coastpolygons,self.coastpolygontypes,coastpolygonsll):
                x = poly[0]
                y = poly[1]
                lons = polyll[0]
                lats = polyll[1]
                mask = NX.logical_or(NX.greater(x,1.e20),NX.greater(y,1.e20))
                # replace values in polygons that are over the horizon.
                xsave = False
                ysave = False
                if NX.sum(mask):
                    i1,i2 = self._splitseg(x,y,mask=mask)
                    for i,j in zip(i1,i2):
                        if i and j != len(x):
                            dist,az1,alpha21=vinc_dist(self.flattening,rad,lat_0,lon_0,math.radians(lats[i]),math.radians(lons[i]))
                            lat1,lon1,az=vinc_pt(self.flattening,rad,lat_0,lon_0,az1,0.5*math.pi*rad)
                            dist,az2,alpha21=vinc_dist(self.flattening,rad,lat_0,lon_0,math.radians(lats[j]),math.radians(lons[j]))
                            lat2,lon2,az=vinc_pt(self.flattening,rad,lat_0,lon_0,az2,0.5*math.pi*rad)
                            gc = GreatCircle(self.rmajor,self.rminor,math.degrees(lon2),math.degrees(lat2),math.degrees(lon1),math.degrees(lat1))
                            npoints = int(gc.gcarclen/dtheta)+1
                            if npoints < 2: npoints=2
                            lonstmp, latstmp = gc.points(npoints)
                            xx, yy = self(lonstmp, latstmp)
                            xnew = x[i:j]
                            ynew = y[i:j]
                            xnew = x[i:j] + xx
                            ynew = y[i:j] + yy
                            coastpolygons.append((xnew,ynew))
                            coastpolygontypes.append(polytype)
                        elif i == 0:
                            xsave = x[0:j]
                            ysave = y[0:j]
                            lats_save = lats[0:j]
                            lons_save = lons[0:j]
                        elif j == len(x):
                            xnew = x[i:j] + xsave
                            ynew = y[i:j] + ysave
                            lonsnew = lons[i:j] + lons_save
                            latsnew = lats[i:j] + lats_save
                            dist,az1,alpha21=vinc_dist(self.flattening,rad,lat_0,lon_0,math.radians(latsnew[0]),math.radians(lonsnew[0]))
                            lat1,lon1,az=vinc_pt(self.flattening,rad,lat_0,lon_0,az1,0.5*math.pi*rad)
                            dist,az2,alpha21=vinc_dist(self.flattening,rad,lat_0,lon_0,math.radians(latsnew[-1]),math.radians(lonsnew[-1]))
                            lat2,lon2,az=vinc_pt(self.flattening,rad,lat_0,lon_0,az2,0.5*math.pi*rad)
                            gc = GreatCircle(self.rmajor,self.rminor,math.degrees(lon2),math.degrees(lat2),math.degrees(lon1),math.degrees(lat1))
                            npoints = int(gc.gcarclen/dtheta)+1
                            if npoints < 2: npoints=2
                            lonstmp, latstmp = gc.points(npoints)
                            xx, yy = self(lonstmp, latstmp)
                            xnew = xnew + xx
                            ynew = ynew + yy
                            coastpolygons.append((xnew,ynew))
                            coastpolygontypes.append(polytype)
                else:
                    coastpolygons.append(poly)
                    coastpolygontypes.append(polytype)
            self.coastpolygons = coastpolygons
            self.coastpolygontypes = coastpolygontypes
        elif self.projection in ['moll','robin']:
            lon_0 = self.projparams['lon_0']
            coastpolygons=[]
            for poly,polytype,polyll in zip(self.coastpolygons,self.coastpolygontypes,coastpolygonsll):
                x = poly[0]
                y = poly[1]
                lons = polyll[0]
                lats = polyll[1]
                xn=[]
                yn=[]
                # antarctic segment goes from 360 back to 0
                # reorder to go from lon_0-180 to lon_0+180.
                if lats[-1] < -68.0:
                    lons.reverse()
                    lats.reverse()
                    xx,yy = self(lons,lats)
                    xx = NX.array(xx); yy = NX.array(yy)
                    xdist = NX.fabs(xx[1:]-xx[0:-1])
                    if max(xdist) > 1000000:
                        nmin = NX.argmax(xdist)+1
                        xnew = NX.zeros(len(xx),'d')
                        ynew = NX.zeros(len(xx),'d')
                        lonsnew = len(xx)*[0]
                        latsnew = len(xx)*[0]
                        xnew[0:len(xx)-nmin] = xx[nmin:]
                        ynew[0:len(xx)-nmin] = yy[nmin:]
                        xnew[len(xx)-nmin:] = xx[0:nmin]
                        ynew[len(xx)-nmin:] = yy[0:nmin]
                        lonsnew[0:len(xx)-nmin] = lons[nmin:]
                        latsnew[0:len(xx)-nmin] = lats[nmin:]
                        lonsnew[len(xx)-nmin:] = lons[0:nmin]
                        latsnew[len(xx)-nmin:] = lats[0:nmin]
                        x = xnew.tolist(); y = ynew.tolist()
                        lons = lonsnew; lats = latsnew
                    else:
                        x.reverse()
                        y.reverse()
                    # close polygon (add lines along left and right edges down to S pole)
                    for phi in NX.arange(-89.999,lats[0],0.1):
                        xx,yy = self(lon_0-179.99,phi)
                        xn.append(xx); yn.append(yy)
                    xn = xn+x
                    yn = yn+y
                    for phi in NX.arange(lats[-1],-89.999,-0.1):
                        xx,yy = self(lon_0+179.99,phi)
                        xn.append(xx); yn.append(yy)
                # move points outside map to edge of map
                # along constant latitude.
                else:
                    for x,y,lon,lat in zip(x,y,lons,lats):
                        if lon > lon_0+180 or lon < lon_0-180:
                            if lon >= lon_0+180: lon=lon_0+180.
                            if lon <= lon_0-180: lon=lon_0-180.
                            xx,yy = self(lon,lat)
                            xn.append(xx); yn.append(yy)
                        else:
                            xn.append(x); yn.append(y)
                coastpolygons.append((xn,yn))
            self.coastpolygons = coastpolygons

    def createfigure(self,*args,**kwargs):
        """
 create a figure with same aspect ratio as map using pylab interface.
 returns a handle to the figure.
 see pylab.figure docs for details.
 Note:  Don't use this if you don't want pylab to be imported.
        See simpletest_oo.py example for how to create a figure
        with the right aspect ration using the non-pylab matplotlib
        OO interface.
        """
        import pylab as P
        if not kwargs.has_key('figsize'):
            figsize = rcParams['figure.figsize']
            xsize = figsize[0]; ysize = xsize
        else:
            figsize = kwargs['figsize']
            xsize = figsize[0]; ysize = figsize[1]
        if self.aspect <= 1:
            ysize = ysize*self.aspect
        else:
            xsize = xsize/self.aspect
        kwargs['figsize'] = (xsize,ysize)
        return P.figure(*args,**kwargs)

    def _splitseg(self,xx,yy,mask=None):
        """split segment up around missing values (outside projection limb)"""
        if mask is None:
            mask = NX.logical_or(NX.greater_equal(xx,1.e20),NX.greater_equal(yy,1.e20))
        i1=[]; i2=[]
        mprev = 1
        for i,m in enumerate(mask):
            if not m and mprev:
                i1.append(i)
            if m and not mprev:
                i2.append(i)
            mprev = m
        if not mprev: i2.append(len(mask))
        if len(i1) != len(i2):
            raise ValueError,'error in splitting boundary segments'
        return i1,i2

    def _insidemap_seg(self,seg):
        """returns True if any point in segment is inside map region"""
        xx = [x for x,y in seg]
        yy = [y for x,y in seg]
        isin = False
        for x,y in zip(xx,yy):
            if x >= self.xmin and x <= self.xmax and y >= self.ymin and y <= self.ymax:
                isin = True
                break
        return isin

    def _insidemap_poly(self,poly,polyll):
        """returns True if any point in polygon is inside map region"""
        isin = False
        xx = poly[0]; yy = poly[1]
        if self.projection in ['moll','robin']:
            lon_0 = self.projparams['lon_0']
            lons = polyll[0]
            for lon in lons:
                if lon < lon_0+180 and lon > lon_0-180:
                    isin = True
                    break
        else:
            for x,y in zip(xx,yy):
                if x >= self.xmin and x <= self.xmax and y >= self.ymin and y <= self.ymax:
                    isin = True
                    break
        return isin

    def __call__(self,x,y,inverse=False):
        """
 Calling a Basemap class instance with the arguments lon, lat will
 convert lon/lat (in degrees) to x/y native map projection
 coordinates (in meters).  If optional keyword 'inverse' is
 True (default is False), the inverse transformation from x/y
 to lon/lat is performed.

 For cylindrical equidistant projection ('cyl'), this
 does nothing (i.e. x,y == lon,lat).

 lon,lat can be either scalar floats or N arrays.
        """
        if not self.crossgreenwich and self.projection in ['merc','cyl','mill'] and not inverse:
            # for cylindrical projections that don't cross greenwich,
            # adjust lons to be in range 0 to 360.
            try:
                x = (360.+x)%360.
            except:
                x = [(360.+xx)%360. for xx in x]
        return self.projtran(x,y,inverse=inverse)

    def makegrid(self,nx,ny,returnxy=False):
        """
 return arrays of shape (ny,nx) containing lon,lat coordinates of
 an equally spaced native projection grid.
 if returnxy = True, the x,y values of the grid are returned also.
        """
        return self.projtran.makegrid(nx,ny,returnxy=returnxy)

    def drawmapboundary(self,color='k',linewidth=1.0,ax=None):
        """
 draw boundary around map projection region. If ax=None (default),
 default axis instance is used, otherwise specified axis instance is used.
        """
        # get current axes instance (if none specified).
        if ax is None and self.ax is None:
            try:
                ax = pylab.gca()
            except:
                import pylab
                ax = pylab.gca()
        elif ax is None and self.ax is not None:
            ax = self.ax
        x = []
        y = []
        dtheta = 0.01
        if self.projection == 'ortho': # circular region.
            r = (2.*self.rmajor+self.rminor)/3.
            r = r-1. # subtract 1 m to make sure it fits in plot region.
            for az in NX.arange(0.,2.*math.pi+dtheta,dtheta):
                x.append(r*math.cos(az)+0.5*self.xmax)
                y.append(r*math.sin(az)+0.5*self.ymax)
        elif self.projection in ['moll','robin']:  # elliptical region.
            # left side
            lats = NX.arange(-89.9,89.9+dtheta,dtheta).tolist()
            lons = len(lats)*[self.projparams['lon_0']-179.9]
            x,y = self(lons,lats)
            # top.
            lons = NX.arange(self.projparams['lon_0']-179.9,self.projparams['lon_0']+179+dtheta,dtheta).tolist()
            lats = len(lons)*[89.9]
            xx,yy = self(lons,lats)
            x = x+xx; y = y+yy
            # right side
            lats = NX.arange(89.9,-89.9-dtheta,-dtheta).tolist()
            lons = len(lats)*[self.projparams['lon_0']+179.9]
            xx,yy = self(lons,lats)
            x = x+xx; y = y+yy
            # bottom.
            lons = NX.arange(self.projparams['lon_0']+179.9,self.projparams['lon_0']-180-dtheta,-dtheta).tolist()
            lats = len(lons)*[-89.9]
            xx,yy = self(lons,lats)
            x = x+xx; y = y+yy
        else: # all other projections are rectangular.
            x = [self.llcrnrx+1.,self.llcrnrx+1.,self.urcrnrx-1.,self.urcrnrx-1.,self.llcrnrx+1.]
            y = [self.llcrnry+1.,self.urcrnry-1.,self.urcrnry-1.,self.llcrnry+1.,self.llcrnry+1.]
        ax.plot(x,y,color=color,linewidth=linewidth)
        # make sure axis ticks are turned off.
        if self.noticks:
            ax.set_xticks([])
            ax.set_yticks([])
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)

    def fillcontinents(self,color=0.8,ax=None):
        """
 Fill continents.

 color - color to fill continents (default gray).
 ax - axes instance (overrides default axes instance)

 After filling continents, lakes are re-filled with axis background color.
        """
        # get current axes instance (if none specified).
        if ax is None and self.ax is None:
            try:
                ax = pylab.gca()
            except:
                import pylab
                ax = pylab.gca()
        elif ax is None and self.ax is not None:
            ax = self.ax
        # get axis background color.
        axisbgc = ax.get_axis_bgcolor()
        np = 0
        for x,y in self.coastpolygons:
            xa = NX.array(x,'f')
            ya = NX.array(y,'f')
        # check to see if all four corners of domain in polygon (if so,
        # don't draw since it will just fill in the whole map).
            delx = 10; dely = 10
            if self.projection in ['cyl']:
                delx = 0.1
                dely = 0.1
            test1 = NX.fabs(xa-self.urcrnrx) < delx
            test2 = NX.fabs(xa-self.llcrnrx) < delx
            test3 = NX.fabs(ya-self.urcrnry) < dely
            test4 = NX.fabs(ya-self.llcrnry) < dely
            hasp1 = sum(test1*test3)
            hasp2 = sum(test2*test3)
            hasp4 = sum(test2*test4)
            hasp3 = sum(test1*test4)
            if not hasp1 or not hasp2 or not hasp3 or not hasp4:
                xy = zip(xa.tolist(),ya.tolist())
                if self.coastpolygontypes[np] != 2:
                    poly = Polygon(xy,facecolor=color,edgecolor=color,linewidth=0)
                else: # lakes filled with background color.
                    poly = Polygon(xy,facecolor=axisbgc,edgecolor=axisbgc,linewidth=0)
                ax.add_patch(poly)
            np = np + 1
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)

    def drawcoastlines(self,linewidth=1.,color='k',antialiased=1,ax=None):
        """
 Draw coastlines.

 linewidth - coastline width (default 1.)
 color - coastline color (default black)
 antialiased - antialiasing switch for coastlines (default True).
 ax - axes instance (overrides default axes instance)
        """
        # get current axes instance (if none specified).
        if ax is None and self.ax is None:
            try:
                ax = pylab.gca()
            except:
                import pylab
                ax = pylab.gca()
        elif ax is None and self.ax is not None:
            ax = self.ax
        coastlines = LineCollection(self.coastsegs,antialiaseds=(antialiased,))
        coastlines.set_color(color)
        coastlines.set_linewidth(linewidth)
        ax.add_collection(coastlines)
        # make sure axis ticks are turned off
        if self.noticks == True:
            ax.set_xticks([])
            ax.set_yticks([])
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)

    def drawcountries(self,linewidth=0.5,color='k',antialiased=1,ax=None):
        """
 Draw country boundaries.

 linewidth - country boundary line width (default 0.5)
 color - country boundary line color (default black)
 antialiased - antialiasing switch for country boundaries (default True).
 ax - axes instance (overrides default axes instance)
        """
        # get current axes instance (if none specified).
        if ax is None and self.ax is None:
            try:
                ax = pylab.gca()
            except:
                import pylab
                ax = pylab.gca()
        elif ax is None and self.ax is not None:
            ax = self.ax
        coastlines = LineCollection(self.cntrysegs,antialiaseds=(antialiased,))
        coastlines.set_color(color)
        coastlines.set_linewidth(linewidth)
        ax.add_collection(coastlines)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)

    def drawstates(self,linewidth=0.5,color='k',antialiased=1,ax=None):
        """
 Draw state boundaries in Americas.

 linewidth - state boundary line width (default 0.5)
 color - state boundary line color (default black)
 antialiased - antialiasing switch for state boundaries (default True).
 ax - axes instance (overrides default axes instance)
        """
        # get current axes instance (if none specified).
        if ax is None and self.ax is None:
            try:
                ax = pylab.gca()
            except:
                import pylab
                ax = pylab.gca()
        elif ax is None and self.ax is not None:
            ax = self.ax
        coastlines = LineCollection(self.statesegs,antialiaseds=(antialiased,))
        coastlines.set_color(color)
        coastlines.set_linewidth(linewidth)
        ax.add_collection(coastlines)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)

    def drawrivers(self,linewidth=0.5,color='k',antialiased=1,ax=None):
        """
 Draw major rivers.

 linewidth - river boundary line width (default 0.5)
 color - river boundary line color (default black)
 antialiased - antialiasing switch for river boundaries (default True).
 ax - axes instance (overrides default axes instance)
        """
        # get current axes instance (if none specified).
        if ax is None and self.ax is None:
            try:
                ax = pylab.gca()
            except:
                import pylab
                ax = pylab.gca()
        elif ax is None and self.ax is not None:
            ax = self.ax
        coastlines = LineCollection(self.riversegs,antialiaseds=(antialiased,))
        coastlines.set_color(color)
        coastlines.set_linewidth(linewidth)
        ax.add_collection(coastlines)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # make sure axis ticks are turned off.
        if self.noticks:
            ax.set_xticks([])
            ax.set_yticks([])

    def readshapefile(self,shapefile,name,drawbounds=True,
                      linewidth=0.5,color='k',antialiased=1,ax=None):
        """
 read in shape file, draw boundaries on map.
 Requires pyshapelib module from Thuban (http://thuban.intevation.org).

 Restrictions:
  - Assumes shapes are 2D
  - vertices must be in geographic (lat/lon) coordinates.

 shapefile - path to shapefile components.  Example:
  shapefile='/home/jeff/esri/world_borders' assumes that
  world_borders.shp, world_borders.shx and world_borders.dbf
  live in /home/jeff/esri.
 name - name for Basemap attribute to hold the shapefile
  vertices in native map projection coordinates.
  Class attribute name+'_info' is a list of dictionaries, one
  for each shape, containing attributes of each shape from dbf file.
  For example, if name='counties', self.counties
  will be a list of vertices for each shape in map projection
  coordinates and self.counties_info will be a list of dictionaries
  with shape attributes. Rings in individual shapes are split out
  into separate polygons.  Additional keys
  'RINGNUM' and 'SHAPENUM' are added to shape attribute dictionary.
 drawbounds - draw boundaries of shapes (default True)
 linewidth - shape boundary line width (default 0.5)
 color - shape boundary line color (default black)
 antialiased - antialiasing switch for shape boundaries (default True).
 ax - axes instance (overrides default axes instance)

 returns a tuple (num_shapes, type, min, max) containing shape file info.
 num_shapes is the number of shapes, type is the type code (one of
 the SHPT* constants defined in the shapelib module, see
 http://shapelib.maptools.org/shp_api.html) and min and
 max are 4-element lists with the min. and max. values of the
 vertices.
        """
        # open shapefile, read vertices for each object, convert
        # to map projection coordinates (only works for 2D shape types).
        try:
            shp = ShapeFile(shapefile)
        except:
            raise IOError, 'error reading shapefile %s.shp' % shapefile
        try:
            dbf = dbflib.open(shapefile)
        except:
            raise IOError, 'error reading dbffile %s.dbf' % shapefile
        info = shp.info()
        if info[1] not in [1,3,5,8]:
            raise ValueError, 'readshapefile can only handle 2D shape types'
        shpsegs = []
        shpinfo = []
        for npoly in range(shp.info()[0]):
            shp_object = shp.read_object(npoly)
            verts = shp_object.vertices()
            rings = len(verts)
            for ring in range(rings):
                lons, lats = zip(*verts[ring])
                if max(lons) > 721. or min(lons) < -721. or max(lats) > 91. or min(lats) < -91:
                    raise ValueError,'shapefile must have lat/lon vertices  - it looks like this one has vertices in map projection coordinates'
                if not self.crossgreenwich and self.projection in ['merc','cyl','mill']:
                    # for cylindrical projections that don't cross greenwich,
                    # adjust lons to be in range 0 to 360.
                    lons = [(lon + 360.)%360 for lon in lons]
                x, y = self(lons, lats)
                shpsegs.append(zip(x,y))
                if ring == 0:
                    shapedict = dbf.read_record(npoly)
                # add information about ring number to dictionary.
                shapedict['RINGNUM'] = ring+1
                shapedict['SHAPENUM'] = npoly+1
                shpinfo.append(shapedict)
        # draw shape boundaries using LineCollection.
        if drawbounds:
            # get current axes instance (if none specified).
            if ax is None and self.ax is None:
                try:
                    ax = pylab.gca()
                except:
                    import pylab
                    ax = pylab.gca()
            elif ax is None and self.ax is not None:
                ax = self.ax
            # make LineCollections for each polygon.
            lines = LineCollection(shpsegs,antialiaseds=(1,))
            lines.set_color(color)
            lines.set_linewidth(linewidth)
            ax.add_collection(lines)
            # make sure axis ticks are turned off
            if self.noticks == True:
                ax.set_xticks([])
                ax.set_yticks([])
            # set axes limits to fit map region.
            self.set_axes_limits(ax=ax)
        # save segments/polygons and shape attribute dicts as class attributes.
        self.__dict__[name]=shpsegs
        self.__dict__[name+'_info']=shpinfo
        shp.close()
        dbf.close()
        return info

    def drawparallels(self,circles,color='k',linewidth=1., \
                      linestyle='--',dashes=[1,1],labels=[0,0,0,0], \
                      fmt='%g',xoffset=None,yoffset=None,ax=None,**kwargs):
        """
 draw parallels (latitude lines).

 circles - list containing latitude values to draw (in degrees).
 color - color to draw parallels (default black).
 linewidth - line width for parallels (default 1.)
 linestyle - line style for parallels (default '--', i.e. dashed).
 dashes - dash pattern for parallels (default [1,1], i.e. 1 pixel on,
  1 pixel off).
 labels - list of 4 values (default [0,0,0,0]) that control whether
  parallels are labelled where they intersect the left, right, top or
  bottom of the plot. For example labels=[1,0,0,1] will cause parallels
  to be labelled where they intersect the left and bottom of the plot,
  but not the right and top.
 fmt is a format string to format the parallel labels (default '%g').
 xoffset - label offset from edge of map in x-direction
  (default is 0.01 times width of map in map projection coordinates).
 yoffset - label offset from edge of map in y-direction
  (default is 0.01 times height of map in map projection coordinates).
 ax - axes instance (overrides default axes instance)

 additional keyword arguments control text properties for labels (see
  NX.text documentation)
        """
        # get current axes instance (if none specified).
        if ax is None and self.ax is None:
            try:
                ax = pylab.gca()
            except:
                import pylab
                ax = pylab.gca()
        elif ax is None and self.ax is not None:
            ax = self.ax
        # don't draw meridians past latmax, always draw parallel at latmax.
        latmax = 80.
        # offset for labels.
        if yoffset is None:
            yoffset = (self.urcrnry-self.llcrnry)/100.
            if self.aspect > 1:
                yoffset = self.aspect*yoffset
            else:
                yoffset = yoffset/self.aspect
        if xoffset is None:
            xoffset = (self.urcrnrx-self.llcrnrx)/100.

        if self.projection in ['merc','cyl','mill','moll','robin']:
            lons = NX.arange(self.llcrnrlon,self.urcrnrlon+0.1,0.1).astype('f')
        else:
            lons = NX.arange(0,360.1,0.1).astype('f')
        # make sure latmax degree parallel is drawn if projection not merc or cyl or miller
        try:
            circlesl = circles.tolist()
        except:
            circlesl = circles
        if self.projection not in ['merc','cyl','mill','moll','robin']:
            if max(circlesl) > 0 and latmax not in circlesl:
                circlesl.append(latmax)
            if min(circlesl) < 0 and -latmax not in circlesl:
                circlesl.append(-latmax)
        xdelta = 0.01*(self.xmax-self.xmin)
        ydelta = 0.01*(self.ymax-self.ymin)
        for circ in circlesl:
            lats = circ*NX.ones(len(lons),'f')
            x,y = self(lons,lats)
            # remove points outside domain.
            testx = NX.logical_and(x>=self.xmin-xdelta,x<=self.xmax+xdelta)
            x = NX.compress(testx, x)
            y = NX.compress(testx, y)
            testy = NX.logical_and(y>=self.ymin-ydelta,y<=self.ymax+ydelta)
            x = NX.compress(testy, x)
            y = NX.compress(testy, y)
            if len(x) > 1 and len(y) > 1:
                # split into separate line segments if necessary.
                # (not necessary for mercator or cylindrical or miller).
                xd = (x[1:]-x[0:-1])**2
                yd = (y[1:]-y[0:-1])**2
                dist = NX.sqrt(xd+yd)
                split = dist > 500000.
                if NX.sum(split) and self.projection not in ['merc','cyl','mill','moll','robin']:
                   ind = (NX.compress(split,squeeze(split*NX.indices(xd.shape)))+1).tolist()
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
                        l = Line2D(x,y,linewidth=linewidth,linestyle=linestyle)
                        l.set_color(color)
                        l.set_dashes(dashes)
                        ax.add_line(l)
        # draw labels for parallels
        # parallels not labelled for orthographic, robinson or mollweide.
        if self.projection in ['ortho'] and max(labels):
            print 'Warning: Cannot label parallels on Orthographic basemap'
            labels = [0,0,0,0]
        # search along edges of map to see if parallels intersect.
        # if so, find x,y location of intersection and draw a label there.
        if self.projection == 'cyl':
            dx = 0.01; dy = 0.01
        else:
            dx = 1000; dy = 1000
        if self.projection == 'moll' or self.projection == 'robin':
            lon_0 = self.projparams['lon_0']
        for dolab,side in zip(labels,['l','r','t','b']):
            if not dolab: continue
            # for cylindrical projections, don't draw parallels on top or bottom.
            if self.projection in ['cyl','merc','mill','moll','robin'] and side in ['t','b']: continue
            if side in ['l','r']:
                nmax = int((self.ymax-self.ymin)/dy+1)
                yy = linspace(self.llcrnry,self.urcrnry,nmax)
                # mollweide inverse transform undefined at South Pole
                if self.projection == 'moll' and yy[0] < 1.e-4:
                    yy[0] = 1.e-4
                if side == 'l':
                    lons,lats = self(self.llcrnrx*NX.ones(yy.shape,'f'),yy,inverse=True)
                else:
                    lons,lats = self(self.urcrnrx*NX.ones(yy.shape,'f'),yy,inverse=True)
                if max(lons) > 1.e20 or max(lats) > 1.e20:
                    raise ValueError,'inverse transformation undefined - please adjust the map projection region'
                # adjust so 0 <= lons < 360
                lons = [(lon+360) % 360 for lon in lons]
            else:
                nmax = int((self.xmax-self.xmin)/dx+1)
                xx = linspace(self.llcrnrx,self.urcrnrx,nmax)
                if side == 'b':
                    lons,lats = self(xx,self.llcrnry*NX.ones(xx.shape,'f'),inverse=True)
                else:
                    lons,lats = self(xx,self.urcrnry*NX.ones(xx.shape,'f'),inverse=True)
                if max(lons) > 1.e20 or max(lats) > 1.e20:
                    raise ValueError,'inverse transformation undefined - please adjust the map projection region'
                # adjust so 0 <= lons < 360
                lons = [(lon+360) % 360 for lon in lons]
            for lat in circles:
                # find index of parallel (there may be two, so
                # search from left and right).
                nl = _searchlist(lats,lat)
                nr = _searchlist(lats[::-1],lat)
                if nr != -1: nr = len(lons)-nr-1
                if lat<0:
                    if rcParams['text.usetex']:
        	        latlabstr = r'${%s\/^{\circ}\/S}$'%fmt
                    else:
                        latlabstr = u'%s\N{DEGREE SIGN}S'%fmt
                    latlab = latlabstr%NX.fabs(lat)
                elif lat>0:
                    if rcParams['text.usetex']:
        	        latlabstr = r'${%s\/^{\circ}\/N}$'%fmt
                    else:
                        latlabstr = u'%s\N{DEGREE SIGN}N'%fmt
                    latlab = latlabstr%lat
                else:
                    if rcParams['text.usetex']:
        	        latlabstr = r'${%s\/^{\circ}}$'%fmt
                    else:
                        latlabstr = u'%s\N{DEGREE SIGN}'%fmt
                    latlab = latlabstr%lat
                # parallels can intersect each map edge twice.
                for i,n in enumerate([nl,nr]):
                    # don't bother if close to the first label.
                    if i and abs(nr-nl) < 100: continue
                    if n >= 0:
                        if side == 'l':
                            if self.projection in ['moll','robin']:
                                xlab,ylab = self(lon_0-179.9,lat)
                            else:
                                xlab = self.llcrnrx
                            xlab = xlab-xoffset
                            ax.text(xlab,yy[n],latlab,horizontalalignment='right',verticalalignment='center',**kwargs)
                        elif side == 'r':
                            if self.projection in ['moll','robin']:
                                xlab,ylab = self(lon_0+179.9,lat)
                            else:
                                xlab = self.urcrnrx
                            xlab = xlab+xoffset
                            ax.text(xlab,yy[n],latlab,horizontalalignment='left',verticalalignment='center',**kwargs)
                        elif side == 'b':
                            ax.text(xx[n],self.llcrnry-yoffset,latlab,horizontalalignment='center',verticalalignment='top',**kwargs)
                        else:
                            ax.text(xx[n],self.urcrnry+yoffset,latlab,horizontalalignment='center',verticalalignment='bottom',**kwargs)

        # make sure axis ticks are turned off is parallels labelled.
        if self.noticks or max(labels):
            ax.set_xticks([])
            ax.set_yticks([])
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)

    def drawmeridians(self,meridians,color='k',linewidth=1., \
                      linestyle='--',dashes=[1,1],labels=[0,0,0,0],\
                      fmt='%g',xoffset=None,yoffset=None,ax=None,**kwargs):
        """
 draw meridians (longitude lines).

 meridians - list containing longitude values to draw (in degrees).
 color - color to draw meridians (default black).
 linewidth - line width for meridians (default 1.)
 linestyle - line style for meridians (default '--', i.e. dashed).
 dashes - dash pattern for meridians (default [1,1], i.e. 1 pixel on,
  1 pixel off).
 labels - list of 4 values (default [0,0,0,0]) that control whether
  meridians are labelled where they intersect the left, right, top or
  bottom of the plot. For example labels=[1,0,0,1] will cause meridians
  to be labelled where they intersect the left and bottom of the plot,
  but not the right and top.
 fmt is a format string to format the meridian labels (default '%g').
 xoffset - label offset from edge of map in x-direction
  (default is 0.01 times width of map in map projection coordinates).
 yoffset - label offset from edge of map in y-direction
  (default is 0.01 times height of map in map projection coordinates).
 ax - axes instance (overrides default axes instance)

 additional keyword arguments control text properties for labels (see
  NX.text documentation)
        """
        # get current axes instance (if none specified).
        if ax is None and self.ax is None:
            try:
                ax = pylab.gca()
            except:
                import pylab
                ax = pylab.gca()
        elif ax is None and self.ax is not None:
            ax = self.ax
        # don't draw meridians past latmax, always draw parallel at latmax.
        latmax = 80. # not used for cyl, merc or miller projections.
        # offset for labels.
        if yoffset is None:
            yoffset = (self.urcrnry-self.llcrnry)/100.
            if self.aspect > 1:
                yoffset = self.aspect*yoffset
            else:
                yoffset = yoffset/self.aspect
        if xoffset is None:
            xoffset = (self.urcrnrx-self.llcrnrx)/100.

        if self.projection not in ['merc','cyl','mill','moll','robin']:
            lats = NX.arange(-latmax,latmax+0.1,0.1).astype('f')
        else:
            lats = NX.arange(-90,90.1,0.1).astype('f')
        xdelta = 0.1*(self.xmax-self.xmin)
        ydelta = 0.1*(self.ymax-self.ymin)
        for merid in meridians:
            if not self.crossgreenwich and self.projection in ['merc','cyl','mill']:
                # for cylindrical projections that don't cross greenwich,
                # adjust meridians to be in range 0 to 360.
                merid = (merid + 360.)%360
            lons = merid*NX.ones(len(lats),'f')
            x,y = self(lons,lats)
            # remove points outside domain.
            testx = NX.logical_and(x>=self.xmin-xdelta,x<=self.xmax+xdelta)
            x = NX.compress(testx, x)
            y = NX.compress(testx, y)
            testy = NX.logical_and(y>=self.ymin-ydelta,y<=self.ymax+ydelta)
            x = NX.compress(testy, x)
            y = NX.compress(testy, y)
            if len(x) > 1 and len(y) > 1:
                # split into separate line segments if necessary.
                # (not necessary for mercator or cylindrical or miller).
                xd = (x[1:]-x[0:-1])**2
                yd = (y[1:]-y[0:-1])**2
                dist = NX.sqrt(xd+yd)
                split = dist > 500000.
                if NX.sum(split) and self.projection not in ['merc','cyl','mill','moll','robin']:
                   ind = (NX.compress(split,squeeze(split*NX.indices(xd.shape)))+1).tolist()
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
                        l = Line2D(x,y,linewidth=linewidth,linestyle=linestyle)
                        l.set_color(color)
                        l.set_dashes(dashes)
                        ax.add_line(l)
        # draw labels for meridians.
        # meridians not labelled for orthographic, robinson or mollweide
        if self.projection in ['ortho','moll'] and max(labels):
            print 'Warning: Cannot label meridians on Mollweide or Orthographic basemap'
            labels = [0,0,0,0]
        # search along edges of map to see if parallels intersect.
        # if so, find x,y location of intersection and draw a label there.
        if self.projection == 'cyl':
            dx = 0.01; dy = 0.01
        else:
            dx = 1000; dy = 1000
        if self.projection == 'moll' or self.projection == 'robin':
            lon_0 = self.projparams['lon_0']
            xmin,ymin = self(lon_0-179.9,-90)
            xmax,ymax = self(lon_0+179.9,90)
        for dolab,side in zip(labels,['l','r','t','b']):
            if not dolab: continue
            # for cylindrical projections, don't draw meridians on left or right.
            if self.projection in ['cyl','merc','mill','robin','moll'] and side in ['l','r']: continue
            if side in ['l','r']:
                nmax = int((self.ymax-self.ymin)/dy+1)
                yy = linspace(self.llcrnry,self.urcrnry,nmax)
                if side == 'l':
                    lons,lats = self(self.llcrnrx*NX.ones(yy.shape,'f'),yy,inverse=True)
                else:
                    lons,lats = self(self.urcrnrx*NX.ones(yy.shape,'f'),yy,inverse=True)
                if max(lons) > 1.e20 or max(lats) > 1.e20:
                    raise ValueError,'inverse transformation undefined - please adjust the map projection region'
                # adjust so 0 <= lons < 360
                lons = [(lon+360) % 360 for lon in lons]
            else:
                nmax = int((self.xmax-self.xmin)/dx+1)
                xx = linspace(self.llcrnrx,self.urcrnrx,nmax)
                if side == 'b':
                    lons,lats = self(xx,self.llcrnry*NX.ones(xx.shape,'f'),inverse=True)
                else:
                    lons,lats = self(xx,self.urcrnry*NX.ones(xx.shape,'f'),inverse=True)
                if max(lons) > 1.e20 or max(lats) > 1.e20:
                    raise ValueError,'inverse transformation undefined - please adjust the map projection region'
                # adjust so 0 <= lons < 360
                lons = [(lon+360) % 360 for lon in lons]
            for lon in meridians:
                # adjust so 0 <= lon < 360
                lon = (lon+360) % 360
                # find index of meridian (there may be two, so
                # search from left and right).
                nl = _searchlist(lons,lon)
                nr = _searchlist(lons[::-1],lon)
                if nr != -1: nr = len(lons)-nr-1
                if lon>180:
                    if rcParams['text.usetex']:
        	        lonlabstr = r'${%s\/^{\circ}\/W}$'%fmt
                    else:
                        lonlabstr = u'%s\N{DEGREE SIGN}W'%fmt
                    lonlab = lonlabstr%NX.fabs(lon-360)
                elif lon<180 and lon != 0:
                    if rcParams['text.usetex']:
        	        lonlabstr = r'${%s\/^{\circ}\/E}$'%fmt
                    else:
                        lonlabstr = u'%s\N{DEGREE SIGN}E'%fmt
                    lonlab = lonlabstr%lon
                else:
                    if rcParams['text.usetex']:
        	        lonlabstr = r'${%s\/^{\circ}}$'%fmt
                    else:
                        lonlabstr = u'%s\N{DEGREE SIGN}'%fmt
                    lonlab = lonlabstr%lon
                # meridians can intersect each map edge twice.
                for i,n in enumerate([nl,nr]):
                    lat = lats[n]/100.
                    # no meridians > latmax for projections other than merc,cyl,miller.
                    if self.projection not in ['merc','cyl','mill'] and lat > latmax: continue
                    # don't bother if close to the first label.
                    if i and abs(nr-nl) < 100: continue
                    if n >= 0:
                        if side == 'l':
                            ax.text(self.llcrnrx-xoffset,yy[n],lonlab,horizontalalignment='right',verticalalignment='center',**kwargs)
                        elif side == 'r':
                            ax.text(self.urcrnrx+xoffset,yy[n],lonlab,horizontalalignment='left',verticalalignment='center',**kwargs)
                        elif side == 'b':
                            if self.projection != 'robin' or (xx[n] > xmin and xx[n] < xmax):
                                ax.text(xx[n],self.llcrnry-yoffset,lonlab,horizontalalignment='center',verticalalignment='top',**kwargs)
                        else:
                            if self.projection != 'robin' or (xx[n] > xmin and xx[n] < xmax):
                                ax.text(xx[n],self.urcrnry+yoffset,lonlab,horizontalalignment='center',verticalalignment='bottom',**kwargs)

        # make sure axis ticks are turned off if meridians labelled.
        if self.noticks or max(labels):
            ax.set_xticks([])
            ax.set_yticks([])
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)

    def gcpoints(self,lon1,lat1,lon2,lat2,npoints):
        """
 compute npoints points along a great circle with endpoints
 (lon1,lat1) and (lon2,lat2).  Returns arrays x,y
 with map projection coordinates.
        """
        gc = GreatCircle(self.rmajor,self.rminor,lon1,lat1,lon2,lat2)
        lons, lats = gc.points(npoints)
        x, y = self(lons, lats)
        return x,y

    def drawgreatcircle(self,lon1,lat1,lon2,lat2,dtheta=0.02,**kwargs):
        """
 draw a great circle on the map.

 lon1,lat1 - longitude,latitude of one endpoint of the great circle.
 lon2,lat2 - longitude,latitude of the other endpoint of the great circle.
 dtheta - points on great circle computed every dtheta radians (default 0.02).

 Other keyword arguments (**kwargs) control plotting of great circle line,
 see NX.plot documentation for details.

 Note:  cannot handle situations in which the great circle intersects
 the edge of the map projection domain, and then re-enters the domain.
        """
        # use great circle formula for a perfect sphere.
        gc = GreatCircle(self.rmajor,self.rminor,lon1,lat1,lon2,lat2)
        if gc.antipodal:
            raise ValueError,'cannot draw great circle whose endpoints are antipodal'
        # points have spacing of dtheta radians.
        npoints = int(gc.gcarclen/dtheta)+1
        lons, lats = gc.points(npoints)
        x, y = self(lons, lats)
        self.plot(x,y,**kwargs)

    def transform_scalar(self,datin,lons,lats,nx,ny,returnxy=False,checkbounds=False,order=1):
        """
 interpolate a scalar field (datin) from a lat/lon grid with longitudes =
 lons and latitudes = lats to a (ny,nx) native map projection grid.

 lons, lats must be rank-1 arrays containing longitudes and latitudes
 (in degrees) of datin grid in increasing order
 (i.e. from Greenwich meridian eastward, and South Pole northward).

 if returnxy=True, the x and y values of the native map projection grid
 are also returned.

 If checkbounds=True, values of xout and yout are checked to see that
 they lie within the range specified by xin and xin.  Default is False.
 If checkbounds=False, and xout,yout are outside xin,yin, interpolated
 values will be clipped to values on boundary of input grid (xin,yin).

 The order keyword can be 0 for nearest-neighbor interpolation,
 or 1 for bilinear interpolation (default 1).

 data on a lat/lon grid must be transformed to map projection coordinates
 before it can be plotted on the map with imshow.
        """
        if returnxy:
            lonsout, latsout, x, y = self.makegrid(nx,ny,returnxy=True)
            datout = interp(datin,lons,lats,lonsout,latsout,checkbounds=checkbounds,order=order)
            return datout, x, y
        else:
            lonsout, latsout = self.makegrid(nx,ny)
            datout = interp(datin,lons,lats,lonsout,latsout,checkbounds=checkbounds,order=order)
            return datout

    def transform_vector(self,uin,vin,lons,lats,nx,ny,returnxy=False,checkbounds=False,order=1):
        """
 rotate and interpolate a vector field (uin,vin) from a lat/lon grid
 with longitudes = lons and latitudes = lats to a
 (ny,nx) native map projection grid.

 lons, lats must be rank-1 arrays containing longitudes and latitudes
 (in degrees) of datin grid in increasing order
 (i.e. from Greenwich meridian eastward, and South Pole northward).

 The input vector field is defined in spherical coordinates (it
 has eastward and northward components) while the output
 vector field is rotated to map projection coordinates (relative
 to x and y). The magnitude of the vector is preserved.

 if returnxy=True, the x and y values of the native map projection grid
 are also returned (default False).

 If checkbounds=True, values of xout and yout are checked to see that
 they lie within the range specified by xin and xin.  Default is False.
 If checkbounds=False, and xout,yout are outside xin,yin, interpolated
 values will be clipped to values on boundary of input grid (xin,yin).

 The order keyword can be 0 for nearest-neighbor interpolation,
 or 1 for bilinear interpolation (default 1).
        """
        lonsout, latsout, x, y = self.makegrid(nx,ny,returnxy=True)
        # interpolate to map projection coordinates.
        uin = interp(uin,lons,lats,lonsout,latsout,checkbounds=checkbounds,order=order)
        vin = interp(vin,lons,lats,lonsout,latsout,checkbounds=checkbounds,order=order)
        # rotate from geographic to map coordinates.
        delta = 0.1 # incement in latitude used to estimate derivatives.
        xn,yn = self(lonsout,NX.where(latsout+delta<90.,latsout+delta,latsout-delta))
        dxdlat = NX.where(latsout+delta<90.,(xn-x)/(latsout+delta),(x-xn)/(latsout+delta))
        dydlat = NX.where(latsout+delta<90.,(yn-y)/(latsout+delta),(y-yn)/(latsout+delta))
        # northangle is the angle between true north and the y axis.
        northangle = NX.arctan2(dxdlat,dydlat)
        uout = uin*NX.cos(northangle) + vin*NX.sin(northangle)
        vout = vin*NX.cos(northangle) - uin*NX.sin(northangle)
        if returnxy:
            return uout,vout,x,y
        else:
            return uout,vout

    def rotate_vector(self,uin,vin,lons,lats,returnxy=False):
        """
 rotate a vector field (uin,vin) on a rectilinear lat/lon grid
 with longitudes = lons and latitudes = lats from geographical into map
 projection coordinates.

 Differs from transform_vector in that no interpolation is done,
 the vector is returned on the same lat/lon grid, but rotated into
 x,y coordinates.

 lons, lats must be rank-2 arrays containing longitudes and latitudes
 (in degrees) of grid.

 if returnxy=True, the x and y values of the lat/lon grid
 are also returned (default False).

 The input vector field is defined in spherical coordinates (it
 has eastward and northward components) while the output
 vector field is rotated to map projection coordinates (relative
 to x and y). The magnitude of the vector is preserved.
        """
        x, y = self(lons, lats)
        # rotate from geographic to map coordinates.
        delta = 0.1 # incement in latitude used to estimate derivatives.
        xn,yn = self(lons,NX.where(lats+delta<90.,lats+delta,lats-delta))
        dxdlat = NX.where(lats+delta<90.,(xn-x)/(lats+delta),(x-xn)/(lats+delta))
        dydlat = NX.where(lats+delta<90.,(yn-y)/(lats+delta),(y-yn)/(lats+delta))
        # northangle is the angle between true north and the y axis.
        northangle = NX.arctan2(dxdlat,dydlat)
        uout = uin*NX.cos(northangle) + vin*NX.sin(northangle)
        vout = vin*NX.cos(northangle) - uin*NX.sin(northangle)
        if returnxy:
            return uout,vout,x,y
        else:
            return uout,vout

    def set_axes_limits(self,ax=None):
        """
 Set axis limits for map domain using current or specified axes instance.
        """
        # get current axes instance (if none specified).
        if ax is None and self.ax is None:
            try:
                ax = pylab.gca()
            except:
                import pylab
                ax = pylab.gca()
        elif ax is None and self.ax is not None:
            ax = self.ax
        corners = ((self.llcrnrx,self.llcrnry), (self.urcrnrx,self.urcrnry))
        ax.update_datalim( corners )
        ax.set_xlim((self.llcrnrx, self.urcrnrx))
        ax.set_ylim((self.llcrnry, self.urcrnry))
        # turn off axes frame for non-rectangular projections.
        if self.projection in ['ortho','moll','robin']:
            ax.set_frame_on(False)

    def scatter(self, *args, **kwargs):
        """
 Plot points with markers on the map (see NX.scatter documentation).
 extra keyword 'ax' can be used to override the default axes instance.
        """
        if not self.crossgreenwich and self.projection in ['merc','cyl','mill']:
            # for cylindrical projections that don't cross greenwich,
            # adjust lons to be in range 0 to 360.
            args = list(args)
            args[0] = [(lon + 360.)%360 for lon in args[0]]
            args = tuple(args)
        if not kwargs.has_key('ax') and self.ax is None:
            try:
                ax = pylab.gca()
            except:
                import pylab
                ax = pylab.gca()
        elif not kwargs.has_key('ax') and self.ax is not None:
            ax = self.ax
        else:
            ax = popd(kwargs,'ax')
        # allow callers to override the hold state by passing hold=True|False
        b = ax.ishold()
        h = popd(kwargs, 'hold', None)
        if h is not None:
            ax.hold(h)
        try:
            ret =  ax.scatter(*args, **kwargs)
            try:
                pylab.draw_if_interactive()
            except:
                pass
        except:
            ax.hold(b)
            raise
        ax.hold(b)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # make sure axis ticks are turned off.
        if self.noticks:
            ax.set_xticks([])
            ax.set_yticks([])
        return ret

    def plot(self, *args, **kwargs):
        """
 Draw lines and/or markers on the map (see NX.plot documentation).
 extra keyword 'ax' can be used to override the default axis instance.
        """
        if not self.crossgreenwich and self.projection in ['merc','cyl','mill']:
            # for cylindrical projections that don't cross greenwich,
            # adjust lons to be in range 0 to 360.
            args = list(args)
            args[0] = [(lon + 360.)%360 for lon in args[0]]
            args = tuple(args)
        if not kwargs.has_key('ax') and self.ax is None:
            try:
                ax = pylab.gca()
            except:
                import pylab
                ax = pylab.gca()
        elif not kwargs.has_key('ax') and self.ax is not None:
            ax = self.ax
        else:
            ax = popd(kwargs,'ax')
        # allow callers to override the hold state by passing hold=True|False
        b = ax.ishold()
        h = popd(kwargs, 'hold', None)
        if h is not None:
            ax.hold(h)
        try:
            ret =  ax.plot(*args, **kwargs)
            try:
                pylab.draw_if_interactive()
            except:
                pass
        except:
            ax.hold(b)
            raise
        ax.hold(b)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # make sure axis ticks are turned off.
        if self.noticks:
            ax.set_xticks([])
            ax.set_yticks([])
        return ret

    def imshow(self, *args, **kwargs):
        """
 Display an image over the map (see NX.imshow documentation).
 extent and origin keywords set automatically so image will be drawn
 over map region.
 extra keyword 'ax' can be used to override the default axis instance.
        """
        if not kwargs.has_key('ax') and self.ax is None:
            try:
                ax = pylab.gca()
            except:
                import pylab
                ax = pylab.gca()
        elif not kwargs.has_key('ax') and self.ax is not None:
            ax = self.ax
        else:
            ax = popd(kwargs,'ax')
        kwargs['extent']=(self.llcrnrx,self.urcrnrx,self.llcrnry,self.urcrnry)
        # use origin='lower', unless overridden.
        if not kwargs.has_key('origin'):
            kwargs['origin']='lower'
        # allow callers to override the hold state by passing hold=True|False
        b = ax.ishold()
        h = popd(kwargs, 'hold', None)
        if h is not None:
            ax.hold(h)
        try:
            ret =  ax.imshow(*args, **kwargs)
            try:
                pylab.draw_if_interactive()
            except:
                pass
        except:
            ax.hold(b)
            raise
        ax.hold(b)
        # reset current active image (only if pylab is imported).
        try:
            pylab.gci._current = ret
        except:
            pass
        return ret

    def pcolor(self,x,y,data,**kwargs):
        """
 Make a pseudo-color plot over the map.
 see NX.pcolor documentation for definition of **kwargs
 extra keyword 'ax' can be used to override the default axis instance.
        """
        if not kwargs.has_key('ax') and self.ax is None:
            try:
                ax = pylab.gca()
            except:
                import pylab
                ax = pylab.gca()
        elif not kwargs.has_key('ax') and self.ax is not None:
            ax = self.ax
        else:
            ax = popd(kwargs,'ax')
        kwargs['extent']=(self.llcrnrx,self.urcrnrx,self.llcrnry,self.urcrnry)
        # make x,y masked arrays
        # (masked where data is outside of projection limb)
        x = ma.masked_values(NX.where(x > 1.e20,1.e20,x), 1.e20)
        y = ma.masked_values(NX.where(y > 1.e20,1.e20,y), 1.e20)
        # allow callers to override the hold state by passing hold=True|False
        b = ax.ishold()
        h = popd(kwargs, 'hold', None)
        if h is not None:
            ax.hold(h)
        try:
            ret =  ax.pcolor(x,y,data,**kwargs)
            try:
                pylab.draw_if_interactive()
            except:
                pass
        except:
            ax.hold(b)
            raise
        ax.hold(b)
        # reset current active image (only if pylab is imported).
        try:
            pylab.gci._current = ret
        except:
            pass
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # make sure axis ticks are turned off.
        if self.noticks:
            ax.set_xticks([])
            ax.set_yticks([])
        return ret

    def contour(self,x,y,data,*args,**kwargs):
        """
 Make a contour plot over the map (see NX.contour documentation).
 extra keyword 'ax' can be used to override the default axis instance.
        """
        if not kwargs.has_key('ax') and self.ax is None:
            try:
                ax = pylab.gca()
            except:
                import pylab
                ax = pylab.gca()
        elif not kwargs.has_key('ax') and self.ax is not None:
            ax = self.ax
        else:
            ax = popd(kwargs,'ax')
        # mask for points outside projection limb.
        xymask = NX.logical_or(NX.greater(x,1.e20),NX.greater(y,1.e20))
        data = ma.asarray(data)
        # combine with data mask.
        mask = NX.logical_or(ma.getmaskarray(data),xymask)
        data = ma.masked_array(data,mask=mask)
        # allow callers to override the hold state by passing hold=True|False
        b = ax.ishold()
        h = popd(kwargs, 'hold', None)
        if h is not None:
            ax.hold(h)
        try:
            CS = ax.contour(x,y,data,*args,**kwargs)
            try:
                pylab.draw_if_interactive()
            except:
                pass
        except:
            ax.hold(b)
            raise
        ax.hold(b)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # make sure axis ticks are turned off.
        if self.noticks:
            ax.set_xticks([])
            ax.set_yticks([])
        # reset current active image (only if pylab is imported).
        try:
            if CS._A is not None: pylab.gci._current = CS
        except:
            pass
        return CS

    def contourf(self,x,y,data,*args,**kwargs):
        """
 Make a filled contour plot over the map (see NX.contourf documentation).
 extra keyword 'ax' can be used to override the default axis instance.
        """
        if not kwargs.has_key('ax') and self.ax is None:
            try:
                ax = pylab.gca()
            except:
                import pylab
                ax = pylab.gca()
        elif not kwargs.has_key('ax') and self.ax is not None:
            ax = self.ax
        else:
            ax = popd(kwargs,'ax')
        # mask for points outside projection limb.
        xymask = NX.logical_or(NX.greater(x,1.e20),NX.greater(y,1.e20))
        data = ma.asarray(data)
        # combine with data mask.
        mask = NX.logical_or(ma.getmaskarray(data),xymask)
        data = ma.masked_array(data,mask=mask)
        # allow callers to override the hold state by passing hold=True|False
        b = ax.ishold()
        h = popd(kwargs, 'hold', None)
        if h is not None:
            ax.hold(h)
        try:
            CS = ax.contourf(x,y,data,*args,**kwargs)
            try:
                pylab.draw_if_interactive()
            except:
                pass
        except:
            ax.hold(b)
            raise
        ax.hold(b)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # make sure axis ticks are turned off.
        if self.noticks:
            ax.set_xticks([])
            ax.set_yticks([])
        # reset current active image (only if pylab is imported).
        try:
            if CS._A is not None: pylab.gci._current = CS
        except:
            pass
        return CS

    def quiver(self, x, y, u, v, scale=None, **kwargs):
        """
 Make a vector plot (u, v) with arrows on the map projection grid (x,y)
 If scale is specified, it is used to scale the vectors. If scale=None
 (default) arrows are scaled to longest one is equal to the maximum
 distance between grid points.

 Extra keyword arguments (**kwargs) passed to quiver Axes method (see
 NX.quiver documentation for details).
 extra keyword 'ax' can be used to override the default axis instance.
        """
        if not kwargs.has_key('ax') and self.ax is None:
            try:
                ax = pylab.gca()
            except:
                import pylab
                ax = pylab.gca()
        elif not kwargs.has_key('ax') and self.ax is not None:
            ax = self.ax
        else:
            ax = popd(kwargs,'ax')
        ny = x.shape[0]; nx = x.shape[1]
        if scale is None:
            scale = max([(self.xmax-self.xmin)/(nx-1),(self.ymax-self.ymin)/(ny-1)])
        else:
            scale = scale
        # allow callers to override the hold state by passing hold=True|False
        b = ax.ishold()
        h = popd(kwargs, 'hold', None)
        if h is not None:
            ax.hold(h)
        try:
            ret =  ax.quiver(x,y,u,v,scale,**kwargs)
            try:
                pylab.draw_if_interactive()
            except:
                pass
        except:
            ax.hold(b)
            raise
        ax.hold(b)
        # set axes limits to fit map region.
        self.set_axes_limits(ax=ax)
        # make sure axis ticks are turned off.
        if self.noticks:
            ax.set_xticks([])
            ax.set_yticks([])
        return ret

def _searchlist(a,x):
    """
 like bisect, but works for lists that are not sorted,
 and are not in increasing order.
 returns -1 if x does not fall between any two elements"""
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

def interp(datain,xin,yin,xout,yout,checkbounds=False,order=1):
    """
 dataout = interp(datain,xin,yin,xout,yout,order=1)

 interpolate data (datain) on a rectilinear grid (with x=xin
 y=yin) to a grid with x=xout, y=yout.

 datain is a rank-2 array with 1st dimension corresponding to y,
 2nd dimension x.

 xin, yin are rank-1 arrays containing x and y of
 datain grid in increasing order.

 xout, yout are rank-2 arrays containing x and y of desired output grid.

 If checkbounds=True, values of xout and yout are checked to see that
 they lie within the range specified by xin and xin.  Default is False.
 If checkbounds=False, and xout,yout are outside xin,yin, interpolated
 values will be clipped to values on boundary of input grid (xin,yin).

 The order keyword can be 0 for nearest-neighbor interpolation,
 or 1 for bilinear interpolation (default 1).
    """
    # xin and yin must be monotonically increasing.
    if xin[-1]-xin[0] < 0 or yin[-1]-yin[0] < 0:
        raise ValueError, 'xin and yin must be increasing!'
    if xout.shape != yout.shape:
        raise ValueError, 'xout and yout must have same shape!'
    # check that xout,yout are
    # within region defined by xin,yin.
    if checkbounds:
        if min(NX.ravel(xout)) < min(xin) or \
           max(NX.ravel(xout)) > max(xin) or \
           min(NX.ravel(yout)) < min(yin) or \
           max(NX.ravel(yout)) > max(yin):
            raise ValueError, 'yout or xout outside range of yin or xin'
    # compute grid coordinates of output grid.
    xoutflat = NX.ravel(xout)
    youtflat = NX.ravel(yout)
    datainflat = NX.ravel(datain)
    if max(xin)-min(xin) < 1.e-4 and max(yin)-min(yin) < 1.e-4:
        # regular input grid.
        xcoords = (len(xin)-1)*(xoutflat-xin[0])/(xin[-1]-xin[0])
        ycoords = (len(yin)-1)*(youtflat-yin[0])/(yin[-1]-yin[0])
    else:
        # irregular (but still rectilinear) input grid.
        ix = NX.searchsorted(xin,xoutflat)-1
        iy = NX.searchsorted(yin,youtflat)-1
        xcoords = NX.zeros(ix.shape,'f')
        ycoords = NX.zeros(iy.shape,'f')
        for n,i in enumerate(ix):
            if i < 0:
                xcoords[n] = -1 # outside of range on xin (lower end)
            elif i >= len(xin)-1:
                xcoords[n] = len(xin) # outside range on upper end.
            else:
                xcoords[n] = float(i)+(xoutflat[n]-xin[i])/(xin[i+1]-xin[i])
        for m,j in enumerate(iy):
            if j < 0:
                ycoords[m] = -1 # outside of range of yin (on lower end)
            elif j >= len(yin)-1:
                ycoords[m] = len(yin) # outside range on upper end
            else:
                ycoords[m] = float(j)+(youtflat[m]-yin[j])/(yin[j+1]-yin[j])
    # data outside range xin,yin will be clipped to
    # values on boundary.
    xcoords = NX.clip(xcoords,0,len(xin)-1)
    ycoords = NX.clip(ycoords,0,len(yin)-1)
    # interpolate to output grid using bilinear interpolation.
    if order == 1:
        xi = xcoords.astype('i')
        yi = ycoords.astype('i')
        xip1 = xi+1
        yip1 = yi+1
        xip1 = NX.clip(xip1,0,len(xin)-1)
        yip1 = NX.clip(yip1,0,len(yin)-1)
        delx = xcoords-xi.astype('f')
        dely = ycoords-yi.astype('f')
        coords = yi*datain.shape[1]+xi
        data11 = NX.take(datainflat,coords)
        coords = yip1*datain.shape[1]+xip1
        data22 = NX.take(datainflat,coords)
        coords = yi*datain.shape[1]+xip1
        data12 = NX.take(datainflat,coords)
        coords = yip1*datain.shape[1]+xi
        data21 = NX.take(datainflat,coords)
        dataout = (1.-delx)*(1.-dely)*data11 + \
                  delx*dely*data22 + \
                  (1.-delx)*dely*data21 + \
                  delx*(1.-dely)*data12
        dataout = NX.reshape(dataout,xout.shape).astype(datain.typecode())
    elif order == 0:
        xcoordsi = NX.around(xcoords).astype('i')
        ycoordsi = NX.around(ycoords).astype('i')
        coords = ycoordsi*datain.shape[1]+xcoordsi
        dataout = NX.take(datainflat,coords)
        dataout = NX.reshape(dataout,xout.shape).astype(datain.typecode())
    else:
        raise ValueError,'order keyword must be 0 or 1'
    return dataout

def shiftgrid(lon0,datain,lonsin,start=True):
    """
 shift global lat/lon grid east or west.
 assumes wraparound (or cyclic point) is included.

 lon0:  starting longitude for shifted grid
        (ending longitude if start=False). lon0 must be on
        input grid (with the range of lonsin).
 datain:  original data.
 lonsin:  original longitudes.
 start[True]: if True, lon0 represents he starting longitude
 of the new grid. if False, lon0 is the ending longitude.

 returns dataout,lonsout (data and longitudes on shifted grid).
    """
    if NX.fabs(lonsin[-1]-lonsin[0]-360.) > 1.e-4:
        raise ValueError, 'cyclic point not included'
    if lon0 < lonsin[0] or lon0 > lonsin[-1]:
        raise ValueError, 'lon0 outside of range of lonsin'
    i0 = NX.argsort(NX.fabs(lonsin-lon0))[0]
    dataout = NX.zeros(datain.shape,datain.typecode())
    lonsout = NX.zeros(lonsin.shape,lonsin.typecode())
    if start:
        lonsout[0:len(lonsin)-i0] = lonsin[i0:]
    else:
        lonsout[0:len(lonsin)-i0] = lonsin[i0:]-360.
    dataout[:,0:len(lonsin)-i0] = datain[:,i0:]
    if start:
        lonsout[len(lonsin)-i0:] = lonsin[1:i0+1]+360.
    else:
        lonsout[len(lonsin)-i0:] = lonsin[1:i0+1]
    dataout[:,len(lonsin)-i0:] = datain[:,1:i0+1]
    return dataout,lonsout

def addcyclic(arrin,lonsin):
   """
 Add cyclic (wraparound) point in longitude.
   """
   nlats = arrin.shape[0]
   nlons = arrin.shape[1]
   arrout  = NX.zeros((nlats,nlons+1),arrin.typecode())
   arrout[:,0:nlons] = arrin[:,:]
   arrout[:,nlons] = arrin[:,0]
   lonsout = NX.zeros(nlons+1,lonsin.typecode())
   lonsout[0:nlons] = lonsin[:]
   lonsout[nlons]  = lonsin[-1] + lonsin[1]-lonsin[0]
   return arrout,lonsout
