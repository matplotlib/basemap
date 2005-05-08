from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from numarray import nd_image
import sys, os, pylab, math
from proj import Proj
from greatcircle import GreatCircle

# look in sys.prefix for directory containing basemap files if
# BASEMAP_DATA_PATH env var not set.
_datadir = os.environ.get('BASEMAP_DATA_PATH')
if not _datadir:
   _datadir = os.path.join(sys.prefix,'share/basemap-py'+repr(sys.version_info[0])+repr(sys.version_info[1])) 

class Basemap:

    """
 Set up a basemap with one of 10 supported map projections
 (cylindrical equidistant, mercator,
 transverse mercator, miller cylindrical, lambert conformal conic,
 azimuthal equidistant, equidistant conic, lambert azimuthal equal area,
 albers equal area conic or stereographic).
 Doesn't actually draw anything, but sets up the map projection class and
 creates the coastline and political boundary polygons in native map 
 projection coordinates.  Requires matplotlib and numarray.
 Uses a pyrex interface to C-code from proj.4 (http://proj.maptools.org).
 
 Useful instance variables:
 
 projection - map projection ('cyl','merc','mill','lcc','eqdc','aea','aeqd',
  'laea', 'tmerc' or 'stere')
 aspect - map aspect ratio (size of y dimension / size of x dimension).
 llcrnrlon - longitude of lower left hand corner of the desired map domain.
 llcrnrlon - latitude of lower left hand corner of the desired map domain.      
 urcrnrlon - longitude of upper right hand corner of the desired map domain.
 urcrnrlon - latitude of upper right hand corner of the desired map domain.
 llcrnrx,llcrnry,urcrnrx,urcrnry - corners of map domain in projection coordinates.

 Example Usage:
 (this example plus others can be found in the examples directory of
  the source distribution)

>>> from matplotlib.toolkits.basemap import Basemap
>>> import cPickle
>>> from import pylab *
>>> # read in topo data from pickle (on a regular lat/lon grid)
>>> topodict = cPickle.load(open('etopo20.pickle','rb'))
>>> etopo = topodict['data']; lons = topodict['lons']; lats = topodict['lats']
>>> # setup map projection (global cylindrical equidistant is default)
>>> m = Basemap(lons[0],lats[0],lons[-1],lats[-1])
>>> # setup figure with same aspect ratio as map.
>>> xsize = rcParams['figure.figsize'][0]
>>> fig=figure(figsize=(xsize,m.aspect*xsize))
>>> fig.add_axes([0.1,0.1,0.8,0.8])
>>> im = m.imshow(etopo) # plot image over map.
>>> # draw coastlines and fill continents.
>>> m.drawcoastlines()
>>> m.fillcontinents()
>>> # draw parallels, label on bottom.
>>> circles = arange(-90.,120.,30.)
>>> m.drawparallels(circles,labels=[1,0,0,0])
>>> # draw meridians, label on left.
>>> meridians = arange(0.,390.,60.)
>>> m.drawmeridians(meridians,labels=[0,0,0,1])
>>> title('Cylindrical Equidistant')
>>> show()

 Version: 0.4.1 (20050509)
 Contact: Jeff Whitaker <jeffrey.s.whitaker@noaa.gov>
    """

    def __init__(self,llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat,\
        resolution='c',area_thresh=10000.,projection='cyl',rsphere=6370997,\
        lat_ts=None,lat_1=None,lat_2=None,lat_0=None,lon_0=None):
        """
 create a Basemap instance.
 
 mandatory input arguments:
 
 llcrnrlon - longitude of lower left hand corner of the desired map domain.
 llcrnrlon - latitude of lower left hand corner of the desired map domain.      
 urcrnrlon - longitude of upper right hand corner of the desired map domain.
 urcrnrlon - latitude of upper right hand corner of the desired map domain.
 
 optional keyword parameters:
 
 resolution - resolution of coastline database to use. Can be 'c' (crude), 
  'l' (low), or 'i' (intermediate). Resolution drops off by roughly 80%
  between datasets.  Higher res datasets are much slower to draw.
  Default 'c'. Coastline data is from the GSHHS
  (http://www.soest.hawaii.edu/wessel/gshhs/gshhs.html).
 area_thresh - coastline with an area smaller than area_thresh in km^2
  will not be plotted.  Default 10,000.
 projection - map projection.  'cyl' - cylindrical equidistant, 'merc' -
  mercator, 'lcc' - lambert conformal conic, 'stere' - stereographic,
  'aea' - albers equal area conic, 'tmerc' - transverse mercator,  
  'aeqd' - azimuthal equidistant, 'mill' - miller cylindrical,
  'eqdc' - equidistant conic, and 'laea' - lambert azimuthal equal area
  are currently available.  Default 'cyl'.
 rsphere - radius of the sphere used to define map projection (default
  6370997 meters, close to the arithmetic mean radius of the earth). If
  given as a sequence, the first two elements are interpreted as
  the the radii of the major and minor axes of an ellipsoid. Note: sometimes
  an ellipsoid is specified by the major axis and an 'inverse flattening
  parameter' (if).  The minor axis (b) can be computed from the major axis (a) 
  and the inverse flattening parameter using the formula if = a/(a-b).

 The following parameters are map projection parameters which all default to 
 None.  Not all parameters are used by all projections, some are ignored.
 
 lat_ts - latitude of natural origin (used for mercator and stereographic
  projections).
 lat_1 - first standard parallel for lambert conformal, albers
  equal area projection and equidistant conic projections.
 lat_2 - second standard parallel for lambert conformal, albers
  equal area projection and equidistant conic projections.
 lat_0 - central latitude (y-axis origin) - used by stereographic 
  transverse mercator, miller cylindrical and lambert azimuthal projections).
 lon_0 - central meridian (x-axis origin - used by lambert conformal,
  lambert azimuthal, equidistant conic, transverse mercator, oblique mercator,
  miller cylindrical and stereographic projections).
        """     

        # read in coastline data.
        coastlons = []; coastlats = []; coastsegind = []; coastsegarea = []; coastsegtype = []
        i = 0  # the current ind
        for line in open(os.path.join(_datadir,'gshhs_'+resolution+'.txt')):
            linesplit = line.split()
            if line.startswith('P'):
                coastsegind.append(i)
                coastsegtype.append(int(linesplit[3]))
                coastsegarea.append(float(linesplit[5]))
                continue
            # lon/lat
            lon, lat = [float(val) for val in linesplit]
            coastlons.append(lon)
            coastlats.append(lat)
            i += 1

        # read in country boundary data.
        cntrylons = []; cntrylats = []; cntrysegind = []
        i = 0  # the current ind
        for line in open(os.path.join(_datadir,'countries_'+resolution+'.txt')):
            linesplit = line.split()
            if line.startswith('>'):
                cntrysegind.append(i)
                continue
            # lon/lat
            lon, lat = [float(val) for val in linesplit]
            cntrylons.append(lon)
            cntrylats.append(lat)
            i += 1

        # read in state boundaries (Americas only).
        statelons = []; statelats = []; statesegind = []
        i = 0  # the current ind
        for line in open(os.path.join(_datadir,'states_'+resolution+'.txt')):
            linesplit = line.split()
            if line.startswith('>'):
                statesegind.append(i)
                continue
            # lon/lat
            lon, lat = [float(val) for val in linesplit]
            statelons.append(lon)
            statelats.append(lat)
            i += 1

        # extend longitudes around the earth a second time
        # (in case projection region straddles Greenwich meridian).
        # also include negative longitudes, so valid longitudes
        # can range from -360 to 720.
        coastlons2 = [lon+360. for lon in coastlons]
        cntrylons2 = [lon+360. for lon in cntrylons]
        statelons2 = [lon+360. for lon in statelons]
        coastlons3 = [lon-360. for lon in coastlons]
        cntrylons3 = [lon-360. for lon in cntrylons]
        statelons3 = [lon-360. for lon in statelons]

        # set up projections using Proj class.
        self.projection = projection
        self.llcrnrlon = llcrnrlon
        self.llcrnrlat = llcrnrlat
        self.urcrnrlon = urcrnrlon
        self.urcrnrlat = urcrnrlat
        projparams = {}
        projparams['proj'] = projection
        # rsphere instance variable is the
        # arithmetic mean of polar and equatorial radii.
        try:
            projparams['a'] = rsphere[0]
            projparams['b'] = rsphere[1]
            self.rsphere = 0.5*(rsphere[0]+rsphere[1])
        except:
            projparams['R'] = rsphere
            self.rsphere = rsphere

        if projection == 'lcc':
            if lat_1 is None or lon_0 is None:
                raise ValueError, 'must specify lat_1 and lon_0 for Lambert Conformal basemap'
            projparams['lat_1'] = lat_1
            if lat_2 != None:
                 projparams['lat_2'] = lat_2
            projparams['lon_0'] = lon_0
            proj = Proj(projparams,llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat)
        elif projection == 'eqdc':
            if lat_1 is None or lat_2 is None or lon_0 is None:
                raise ValueError, 'must specify lat_1, lat_2 and lon_0 for Equidistant Conic basemap'
            projparams['lat_1'] = lat_1
            projparams['lat_2'] = lat_2
            projparams['lon_0'] = lon_0
            proj = Proj(projparams,llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat)
        elif projection == 'aea':
            if lat_1 is None or lat_2 is None or lon_0 is None:
                raise ValueError, 'must specify lat_1, lat_2 and lon_0 for Albers Equal Area basemap'
            projparams['lat_1'] = lat_1
            projparams['lat_2'] = lat_2
            projparams['lon_0'] = lon_0
            proj = Proj(projparams,llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat)
        elif projection == 'stere':
            if lat_ts is None or lat_0 is None or lon_0 is None:
                raise ValueError, 'must specify lat_ts,lat_0 and lon_0 for Stereographic basemap'
            projparams['lat_ts'] = lat_ts
            projparams['lat_0'] = lat_0
            projparams['lon_0'] = lon_0
            proj = Proj(projparams,llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat)
        elif projection == 'laea':
            if lat_0 is None or lon_0 is None:
                raise ValueError, 'must specify lat_0 and lon_0 for Lambert Azimuthal basemap'
            projparams['lat_0'] = lat_0
            projparams['lon_0'] = lon_0
            proj = Proj(projparams,llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat)
        elif projection == 'merc':
            if lat_ts is None:
                raise ValueError, 'must specify lat_ts for Mercator basemap'
            projparams['lat_ts'] = lat_ts
            proj = Proj(projparams,llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat)
        elif projection == 'tmerc':
            if lat_0 is None or lon_0 is None:
                raise ValueError, 'must specify lat_0 and lon_0 for Transverse Mercator basemap'
            projparams['lat_0'] = lat_0
            projparams['lon_0'] = lon_0
            proj = Proj(projparams,llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat)
        elif projection == 'aeqd':
            if lat_0 is None or lon_0 is None:
                raise ValueError, 'must specify lat_0 and lon_0 for Azimuthal Equidistant basemap'
            projparams['lat_0'] = lat_0
            projparams['lon_0'] = lon_0
            proj = Proj(projparams,llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat)
        elif projection == 'mill':
            if lat_0 is not None:
                projparams['lat_0'] = lat_0
            if lon_0 is not None:
                projparams['lon_0'] = lon_0
            proj = Proj(projparams,llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat)
        elif projection == 'cyl':
            proj = Proj(projparams,llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat)
        else:
            raise ValueError, 'unsupported projection'

        # make Proj instance a Basemap instance variable.
        self.projtran = proj
        # set instance variables defining map region.
        self.xmin = proj.xmin
        self.xmax = proj.xmax
        self.ymin = proj.ymin
        self.ymax = proj.ymax
        if projection == 'cyl':
            self.aspect = (urcrnrlat-llcrnrlat)/(urcrnrlon-llcrnrlon)
        else:
            self.aspect = (proj.ymax-proj.ymin)/(proj.xmax-proj.xmin)
        self.llcrnrx = proj.llcrnrx
        self.llcrnry = proj.llcrnry
        self.urcrnrx = proj.urcrnrx
        self.urcrnry = proj.urcrnry

        # transform coastline polygons to native map coordinates.
        xc,yc = proj(pylab.array(coastlons,'f'),pylab.array(coastlats,'f'))
        xc2,yc2 = proj(pylab.array(coastlons2,'f'),pylab.array(coastlats,'f'))
        xc3,yc3 = proj(pylab.array(coastlons3,'f'),pylab.array(coastlats,'f'))
        if projection == 'merc' or projection == 'mill': 
            yc2 = yc
            yc3 = yc

        # set up segments in form needed for LineCollection,
        # ignoring 'inf' values that are off the map, and skipping
        # polygons that have an area > area_thresh..
        segments = [zip(xc[i0:i1],yc[i0:i1]) for a,i0,i1 in zip(coastsegarea[:-1],coastsegind[:-1],coastsegind[1:]) if a > area_thresh]
        segtypes = [i for a,i in zip(coastsegarea[:-1],coastsegtype[:-1]) if a > area_thresh]
        segments2 = [zip(xc2[i0:i1],yc2[i0:i1]) for a,i0,i1 in zip(coastsegarea[:-1],coastsegind[:-1],coastsegind[1:]) if a > area_thresh and max(xc2[i0:i1]) < 1.e20 and max(yc2[i0:i1]) < 1.e20]
        segtypes2 = [i for a,i0,i1,i in zip(coastsegarea[:-1],coastsegind[:-1],coastsegind[1:],coastsegtype[:-1]) if a > area_thresh and max(xc2[i0:i1]) < 1.e20 and max(yc2[i0:i1]) < 1.e20]
        segments3 = [zip(xc3[i0:i1],yc3[i0:i1]) for a,i0,i1 in zip(coastsegarea[:-1],coastsegind[:-1],coastsegind[1:]) if a > area_thresh and max(xc3[i0:i1]) < 1.e20 and max(yc3[i0:i1]) < 1.e20]
        segtypes3 = [i for a,i0,i1,i in zip(coastsegarea[:-1],coastsegind[:-1],coastsegind[1:],coastsegtype[:-1]) if a > area_thresh and max(xc3[i0:i1]) < 1.e20 and max(yc3[i0:i1]) < 1.e20]
        self.coastsegs = segments+segments2+segments3
        self.coastsegtypes = segtypes+segtypes2+segtypes3

        # same as above for country polygons.
        xc,yc = proj(pylab.array(cntrylons,'f'),pylab.array(cntrylats,'f'))
        xc2,yc2 = proj(pylab.array(cntrylons2,'f'),pylab.array(cntrylats,'f'))
        xc3,yc3 = proj(pylab.array(cntrylons3,'f'),pylab.array(cntrylats,'f'))
        if projection == 'merc' or projection == 'mill': 
            yc2=yc
            yc3=yc
        segments = [zip(xc[i0:i1],yc[i0:i1]) for i0,i1 in zip(cntrysegind[:-1],cntrysegind[1:])]
        segments2 = [zip(xc2[i0:i1],yc2[i0:i1]) for i0,i1 in zip(cntrysegind[:-1],cntrysegind[1:]) if max(xc2[i0:i1]) < 1.e20 and max(yc2[i0:i1]) < 1.e20]
        segments3 = [zip(xc3[i0:i1],yc3[i0:i1]) for i0,i1 in zip(cntrysegind[:-1],cntrysegind[1:]) if max(xc3[i0:i1]) < 1.e20 and max(yc3[i0:i1]) < 1.e20]
        self.cntrysegs = segments+segments2+segments3

        # same as above for state polygons.
        xc,yc = proj(pylab.array(statelons,'f'),pylab.array(statelats,'f'))
        xc2,yc2 = proj(pylab.array(statelons2,'f'),pylab.array(statelats,'f'))
        xc3,yc3 = proj(pylab.array(statelons3,'f'),pylab.array(statelats,'f'))
        if projection == 'merc' or projection == 'mill': 
            yc2=yc
            yc3=yc
        segments = [zip(xc[i0:i1],yc[i0:i1]) for i0,i1 in zip(statesegind[:-1],statesegind[1:])]
        segments2 = [zip(xc2[i0:i1],yc2[i0:i1]) for i0,i1 in zip(statesegind[:-1],statesegind[1:]) if max(xc2[i0:i1]) < 1.e20 and max(yc2[i0:i1]) < 1.e20]
        segments3 = [zip(xc3[i0:i1],yc3[i0:i1]) for i0,i1 in zip(statesegind[:-1],statesegind[1:]) if max(xc3[i0:i1]) < 1.e20 and max(yc3[i0:i1]) < 1.e20]
        self.statesegs = segments+segments2+segments3

        # store coast polygons for filling.
	self.coastpolygons = []
        self.coastpolygontypes = []
        if projection == 'merc' or projection ==  'mill': 
            xsp,ysp = proj(0.,-89.9) # s. pole coordinates.
            xa,ya = proj(0.,-68.0) # edge of antarctica.
	    x0,y0 = proj(0.,0.)
	    xm360,ym360 = proj(-360.,0.)
	    x360,y360 = proj(360.,0.)
	    x720,y720 = proj(720.,0.)
        for seg,segtype in zip(self.coastsegs,self.coastsegtypes):
            x = [lon for lon,lat in seg]
            y = [lat for lon,lat in seg]
            # the antarctic polygon is a nuisance, since it
            # spans all longitudes, it's not closed and will not be filled
            # without some projection dependant tweaking.
            if projection == 'cyl':
                if x[-1] == 0.000 and y[-1] < -68.: # close antarctica
                    x.append(0.)
                    y.append(-90.0000)
                    x.insert(0,360.)
                    y.insert(0,-90)
                if x[-1] == 360.000 and y[-1] < -68.: 
                    x.append(360.)
                    y.append(-90)
                    x.insert(0,720.)
                    y.insert(0,-90)
                if x[-1] == -360.000 and y[-1] < -68.: 
                    x.append(-360.)
                    y.append(-90)
                    x.insert(0,0.)
                    y.insert(0,-90)
            elif projection == 'merc' or projection == 'mill':
                if math.fabs(x[-1]-x0) < 1. and y[-1] < ya: # close antarctica
                    x.append(x0)
                    y.append(ysp)
                    x.insert(0,x360)
                    y.insert(0,ysp)
                if math.fabs(x[-1]-x360) < 1. and y[-1] < ya: 
                    x.append(x360)
                    y.append(ysp)
                    x.insert(0,x720)
                    y.insert(0,ysp)
                if math.fabs(x[-1]-xm360) < 1. and y[-1] < ya: 
                    x.append(xm360)
                    y.append(ysp)
                    x.insert(0,x0)
                    y.insert(0,ysp)
            self.coastpolygons.append((x,y))
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
        polygontypes = []
        for poly,polytype in zip(self.coastpolygons,self.coastpolygontypes):
            if self._insidemap_poly(poly):
                polygons.append(poly)
                polygontypes.append(polytype)
        self.coastpolygons = polygons
        self.coastpolygontypes = polygontypes
        states = []
        for seg in self.statesegs:
            if self._insidemap_seg(seg):
                states.append(seg)
        self.statesegs = states
        countries = []
        for seg in self.cntrysegs:
            if self._insidemap_seg(seg):
                countries.append(seg)
        self.cntrysegs = countries

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

    def _insidemap_poly(self,poly):
        """returns True if any point in polygon is inside map region"""
        isin = False
        for x,y in zip(poly[0],poly[1]):
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
        return self.projtran(x,y,inverse=inverse)
 
    def makegrid(self,nx,ny,returnxy=False):
        """
 return arrays of shape (ny,nx) containing lon,lat coordinates of
 an equally spaced native projection grid.
 if returnxy = True, the x,y values of the grid are returned also.
        """
        return self.projtran.makegrid(nx,ny,returnxy=returnxy)

    def fillcontinents(self,color=0.8):
        """
 Fill continents.

 color - color to fill continents (default gray).
        """
        # get current axes instance.
        ax = pylab.gca()
        # get axis background color.
        axisbgc = ax.get_axis_bgcolor()
        np = 0
        for x,y in self.coastpolygons:
            xa = pylab.array(x,'f')
            ya = pylab.array(y,'f')
        # check to see if all four corners of domain in polygon (if so,
        # don't draw since it will just fill in the whole map).
            delx = 10; dely = 10
            if self.projection in ['cyl']:
                delx = 0.1
                dely = 0.1
            test1 = pylab.fabs(xa-self.urcrnrx) < delx
            test2 = pylab.fabs(xa-self.llcrnrx) < delx
            test3 = pylab.fabs(ya-self.urcrnry) < dely
            test4 = pylab.fabs(ya-self.llcrnry) < dely
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
        self.set_axes_limits()

    def drawcoastlines(self,linewidth=1.,color='k',antialiased=1):
        """
 Draw coastlines.

 linewidth - coastline width (default 1.)
 color - coastline color (default black)
 antialiased - antialiasing switch for coastlines (default True).
        """
        # get current axes instance.
        ax = pylab.gca()
        coastlines = LineCollection(self.coastsegs,antialiaseds=(antialiased,))
        try:
            coastlines.set_color(color)
        except: # this was a typo that existed in matplotlib-0.71 and earlier
            coastlines.color(color)
        coastlines.set_linewidth(linewidth)
        ax.add_collection(coastlines)
        # set axes limits to fit map region.
        self.set_axes_limits()

    def drawcountries(self,linewidth=0.5,color='k',antialiased=1):
        """
 Draw country boundaries.

 linewidth - country boundary line width (default 0.5)
 color - country boundary line color (default black)
 antialiased - antialiasing switch for country boundaries (default True).
        """
        # get current axes instance.
        ax = pylab.gca()
        coastlines = LineCollection(self.cntrysegs,antialiaseds=(antialiased,))
        try:
            coastlines.set_color(color)
        except: # this was a typo that existed in matplotlib-0.71 and earlier
            coastlines.color(color)
        coastlines.set_linewidth(linewidth)
        ax.add_collection(coastlines)
        # set axes limits to fit map region.
        self.set_axes_limits()

    def drawstates(self,linewidth=0.5,color='k',antialiased=1):
        """
 Draw state boundaries in Americas.

 linewidth - state boundary line width (default 0.5)
 color - state boundary line color (default black)
 antialiased - antialiasing switch for state boundaries (default True).
        """
        # get current axes instance.
        ax = pylab.gca()
        coastlines = LineCollection(self.statesegs,antialiaseds=(antialiased,))
        try:
            coastlines.set_color(color)
        except: # this was a typo that existed in matplotlib-0.71 and earlier
            coastlines.color(color)
        coastlines.set_linewidth(linewidth)
        ax.add_collection(coastlines)
        # set axes limits to fit map region.
        self.set_axes_limits()

    def drawparallels(self,circles,color='k',linewidth=1., \
                      linestyle='--',dashes=[1,1],labels=[0,0,0,0], \
                      font='rm',fontsize=12):
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
  but not the right and top. Labels are drawn using mathtext.
 font - mathtext font used for labels ('rm','tt','it' or 'cal', default 'rm').
 fontsize - font size in points for labels (default 12).
        """
        # get current axes instance.
        ax = pylab.gca()
        # don't draw meridians past latmax, always draw parallel at latmax.
        latmax = 80.
        # offset for labels.
	yoffset = (self.urcrnry-self.llcrnry)/100./self.aspect
	xoffset = (self.urcrnrx-self.llcrnrx)/100.

        if self.projection in ['merc','cyl','mill']:
            lons = pylab.arange(self.llcrnrlon,self.urcrnrlon+0.1,0.1).astype('f')
        else:
            lons = pylab.arange(0,360.1,0.1).astype('f')
        # make sure latmax degree parallel is drawn if projection not merc or cyl or miller
        try:
            circlesl = circles.tolist()
        except:
            circlesl = circles
        if self.projection not in ['merc','cyl','mill']:
            if max(circlesl) > 0 and latmax not in circlesl: 
                circlesl.append(latmax)
            if min(circlesl) < 0 and -latmax not in circlesl: 
                circlesl.append(-latmax)
        xdelta = 0.1*(self.xmax-self.xmin)
        ydelta = 0.1*(self.ymax-self.ymin)
        for circ in circlesl:
            lats = circ*pylab.ones(len(lons),'f')
            x,y = self(lons,lats)
            # remove points outside domain.
            testx = pylab.logical_and(x>=self.xmin-xdelta,x<=self.xmax+xdelta)
            x = pylab.compress(testx, x)
            y = pylab.compress(testx, y)
            testy = pylab.logical_and(y>=self.ymin-ydelta,y<=self.ymax+ydelta)
            x = pylab.compress(testy, x)
            y = pylab.compress(testy, y)
            if len(x) > 1 and len(y) > 1:
                # split into separate line segments if necessary.
                # (not necessary for mercator or cylindrical or miller).
                xd = (x[1:]-x[0:-1])**2
                yd = (y[1:]-y[0:-1])**2
                dist = pylab.sqrt(xd+yd)
                split = dist > 500000.
                if pylab.asum(split) and self.projection not in ['merc','cyl','mill']:
                   ind = (pylab.compress(split,pylab.squeeze(split*pylab.indices(xd.shape)))+1).tolist()
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
        # search along edges of map to see if parallels intersect.
        # if so, find x,y location of intersection and draw a label there.
        if self.projection == 'cyl':
            dx = 0.01; dy = 0.01
        else:
            dx = 1000; dy = 1000
        for dolab,side in zip(labels,['l','r','t','b']):
            if not dolab: continue
            # for cyl, merc or miller, don't draw parallels on top or bottom.
            if self.projection in ['cyl','merc','mill'] and side in ['t','b']: continue
            if side in ['l','r']:
	        nmax = int((self.ymax-self.ymin)/dy+1)
                if self.urcrnry < self.llcrnry:
	            yy = self.llcrnry-dy*pylab.arange(nmax)
                else:
	            yy = self.llcrnry+dy*pylab.arange(nmax)
                if side == 'l':
	            lons,lats = self(self.llcrnrx*pylab.ones(yy.shape,'f'),yy,inverse=True)
                else:
	            lons,lats = self(self.urcrnrx*pylab.ones(yy.shape,'f'),yy,inverse=True)
                lons = pylab.where(lons < 0, lons+360, lons)
                lons = [int(lon*10) for lon in lons.tolist()]
                lats = [int(lat*10) for lat in lats.tolist()]
            else:
	        nmax = int((self.xmax-self.xmin)/dx+1)
                if self.urcrnrx < self.llcrnrx:
	            xx = self.llcrnrx-dx*pylab.arange(nmax)
                else:
	            xx = self.llcrnrx+dx*pylab.arange(nmax)
                if side == 'b':
	            lons,lats = self(xx,self.llcrnry*pylab.ones(xx.shape,'f'),inverse=True)
                else:
	            lons,lats = self(xx,self.urcrnry*pylab.ones(xx.shape,'f'),inverse=True)
                lons = pylab.where(lons < 0, lons+360, lons)
                lons = [int(lon*10) for lon in lons.tolist()]
                lats = [int(lat*10) for lat in lats.tolist()]
	    for lat in circles:
                # find index of parallel (there may be two, so
                # search from left and right).
                try:
                    nl = lats.index(int(lat*10))
                except:
                    nl = -1
                try:
                    nr = len(lats)-lats[::-1].index(int(lat*10))-1
                except:
                    nr = -1
                if lat<0:
        	    latlab = r'$\%s{%g\/^{\circ}\/S}$'%(font,pylab.fabs(lat))
        	elif lat>0:
        	    latlab = r'$\%s{%g\/^{\circ}\/N}$'%(font,lat)
        	else:
        	    latlab = r'$\%s{%g\/^{\circ}}$'%(font,lat)
                # parallels can intersect each map edge twice.
                for i,n in enumerate([nl,nr]):
                    # don't bother if close to the first label.
                    if i and abs(nr-nl) < 100: continue
                    if n >= 0:
                        if side == 'l':
        	            pylab.text(self.llcrnrx-xoffset,yy[n],latlab,horizontalalignment='right',verticalalignment='center',fontsize=fontsize)
                        elif side == 'r':
        	            pylab.text(self.urcrnrx+xoffset,yy[n],latlab,horizontalalignment='left',verticalalignment='center',fontsize=fontsize)
                        elif side == 'b':
        	            pylab.text(xx[n],self.llcrnry-yoffset,latlab,horizontalalignment='center',verticalalignment='top',fontsize=fontsize)
                        else:
        	            pylab.text(xx[n],self.urcrnry+yoffset,latlab,horizontalalignment='center',verticalalignment='bottom',fontsize=fontsize)

        # make sure axis ticks are turned off
        ax.set_xticks([]) 
        ax.set_yticks([])
        # set axes limits to fit map region.
        self.set_axes_limits()

    def drawmeridians(self,meridians,color='k',linewidth=1., \
                      linestyle='--',dashes=[1,1],labels=[0,0,0,0],\
                      font='rm',fontsize=12):
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
  but not the right and top. Labels are drawn using mathtext.
 font - mathtext font used for labels ('rm','tt','it' or 'cal', default 'rm').
 fontsize - font size in points for labels (default 12).
        """
        # get current axes instance.
        ax = pylab.gca()
        # don't draw meridians past latmax, always draw parallel at latmax.
        latmax = 80. # not used for cyl, merc or miller projections.
        # offset for labels.
	yoffset = (self.urcrnry-self.llcrnry)/100./self.aspect
	xoffset = (self.urcrnrx-self.llcrnrx)/100.

        if self.projection not in ['merc','cyl','mill']:
            lats = pylab.arange(-latmax,latmax+0.1,0.1).astype('f')
        else:
            lats = pylab.arange(-90,90.1,0.1).astype('f')
        xdelta = 0.1*(self.xmax-self.xmin)
        ydelta = 0.1*(self.ymax-self.ymin)
        for merid in meridians:
            lons = merid*pylab.ones(len(lats),'f')
            x,y = self(lons,lats)
            # remove points outside domain.
            testx = pylab.logical_and(x>=self.xmin-xdelta,x<=self.xmax+xdelta)
            x = pylab.compress(testx, x)
            y = pylab.compress(testx, y)
            testy = pylab.logical_and(y>=self.ymin-ydelta,y<=self.ymax+ydelta)
            x = pylab.compress(testy, x)
            y = pylab.compress(testy, y)
            if len(x) > 1 and len(y) > 1:
                # split into separate line segments if necessary.
                # (not necessary for mercator or cylindrical or miller).
                xd = (x[1:]-x[0:-1])**2
                yd = (y[1:]-y[0:-1])**2
                dist = pylab.sqrt(xd+yd)
                split = dist > 500000.
                if pylab.asum(split) and self.projection not in ['merc','cyl','mill']:
                   ind = (pylab.compress(split,pylab.squeeze(split*pylab.indices(xd.shape)))+1).tolist()
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
        # search along edges of map to see if parallels intersect.
        # if so, find x,y location of intersection and draw a label there.
        if self.projection == 'cyl':
            dx = 0.01; dy = 0.01
        else:
            dx = 1000; dy = 1000
        for dolab,side in zip(labels,['l','r','t','b']):
            if not dolab: continue
            # for cyl, merc or miller, don't draw meridians on left or right.
            if self.projection in ['cyl','merc','mill'] and side in ['l','r']: continue
            if side in ['l','r']:
	        nmax = int((self.ymax-self.ymin)/dy+1)
                if self.urcrnry < self.llcrnry:
	            yy = self.llcrnry-dy*pylab.arange(nmax)
                else:
	            yy = self.llcrnry+dy*pylab.arange(nmax)
                if side == 'l':
	            lons,lats = self(self.llcrnrx*pylab.ones(yy.shape,'f'),yy,inverse=True)
                else:
	            lons,lats = self(self.urcrnrx*pylab.ones(yy.shape,'f'),yy,inverse=True)
                lons = pylab.where(lons < 0, lons+360, lons)
                lons = [int(lon*10) for lon in lons.tolist()]
                lats = [int(lat*10) for lat in lats.tolist()]
            else:
	        nmax = int((self.xmax-self.xmin)/dx+1)
                if self.urcrnrx < self.llcrnrx:
	            xx = self.llcrnrx-dx*pylab.arange(nmax)
                else:
	            xx = self.llcrnrx+dx*pylab.arange(nmax)
                if side == 'b':
	            lons,lats = self(xx,self.llcrnry*pylab.ones(xx.shape,'f'),inverse=True)
                else:
	            lons,lats = self(xx,self.urcrnry*pylab.ones(xx.shape,'f'),inverse=True)
                lons = pylab.where(lons < 0, lons+360, lons)
                lons = [int(lon*10) for lon in lons.tolist()]
                lats = [int(lat*10) for lat in lats.tolist()]
            for lon in meridians:
                if lon<0: lon=lon+360.
                # find index of meridian (there may be two, so
                # search from left and right).
                try:
                    nl = lons.index(int(lon*10))
                except:
                    nl = -1
                try:
                    nr = len(lons)-lons[::-1].index(int(lon*10))-1
                except:
                    nr = -1
        	if lon>180:
        	    lonlab = r'$\%s{%g\/^{\circ}\/W}$'%(font,pylab.fabs(lon-360))
        	elif lon<180 and lon != 0:
        	    lonlab = r'$\%s{%g\/^{\circ}\/E}$'%(font,lon)
        	else:
        	    lonlab = r'$\%s{%g\/^{\circ}}$'%(font,lon)
                # meridians can intersect each map edge twice.
                for i,n in enumerate([nl,nr]):
                    lat = lats[n]/10.
                    # no meridians > latmax for projections other than merc,cyl,miller.
                    if self.projection not in ['merc','cyl','mill'] and lat > latmax: continue
                    # don't bother if close to the first label.
                    if i and abs(nr-nl) < 100: continue
                    if n >= 0:
                        if side == 'l':
        	            pylab.text(self.llcrnrx-xoffset,yy[n],lonlab,horizontalalignment='right',verticalalignment='center',fontsize=fontsize)
                        elif side == 'r':
        	            pylab.text(self.urcrnrx+xoffset,yy[n],lonlab,horizontalalignment='left',verticalalignment='center',fontsize=fontsize)
                        elif side == 'b':
        	            pylab.text(xx[n],self.llcrnry-yoffset,lonlab,horizontalalignment='center',verticalalignment='top',fontsize=fontsize)
                        else:
        	            pylab.text(xx[n],self.urcrnry+yoffset,lonlab,horizontalalignment='center',verticalalignment='bottom',fontsize=fontsize)

        # make sure axis ticks are turned off
        ax.set_xticks([]) 
        ax.set_yticks([])
        # set axes limits to fit map region.
        self.set_axes_limits()

    def gcpoints(self,lon1,lat1,lon2,lat2,npoints):
        """
 compute npoints points along a great circle with endpoints
 (lon1,lat1) and (lon2,lat2).  Returns numarrays x,y
 with map projection coordinates.
        """
        gc = GreatCircle(lon1,lat1,lon2,lat2)
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
 see pylab plot documentation for details.

 Note:  cannot handle situations in which the great circle intersects
 the edge of the map projection domain, and then re-enters the domain.
 Assumes a perfect sphere, doesn't take into account eccentricity of the
 ellipsoid (this implies an O(10km) error for a great circle extending
 all the way around the earth).
        """
        # get current axes instance.
        ax = pylab.gca()
        gc = GreatCircle(lon1,lat1,lon2,lat2)
        if gc.antipodal:
            raise ValueError,'cannot draw great circle whose endpoints are antipodal'
        # points have spacing of dtheta radians.
        npoints = int(gc.distance/dtheta)+1
        lons, lats = gc.points(npoints)
        x, y = self(lons, lats)
        self.plot(x,y,**kwargs)

    def transform_scalar(self,datin,lons,lats,nx,ny,returnxy=False,**kwargs):
        """
 transform a scalar field (datin) from a lat/lon grid with longitudes
 lons and latitudes lats to a (ny,nx) native map projection grid.

 lons, lats must be rank-1 arrays containing longitudes and latitudes
 (in degrees) of datin grid in increasing order
 (i.e. from Greenwich meridian eastward, and South Pole northward).

 if returnxy=True, the x and y values of the native map projection grid
 are also returned.

 See interp documentation for meaning of extra keyword arguments (**kwargs).
 
 data on a lat/lon grid must be transformed to map projection coordinates
 before it can be plotted on the map (with the contour, contourf,
 imshow or pcolor class methods).
        """
        if returnxy:
            lonsout, latsout, x, y = self.makegrid(nx,ny,returnxy=True)
            datout = interp(datin,lons,lats,lonsout,latsout,**kwargs)
            return datout, x, y
        else:
            lonsout, latsout = self.makegrid(nx,ny)
            datout = interp(datin,lons,lats,lonsout,latsout,**kwargs)
            return datout

    def transform_vector(self,uin,vin,lons,lats,nx,ny,returnxy=False,**kwargs):
        """
 transform a vector field (uin,vin) from a lat/lon grid with longitudes
 lons and latitudes lats to a (ny,nx) native map projection grid.

 lons, lats must be rank-1 arrays containing longitudes and latitudes
 (in degrees) of datin grid in increasing order
 (i.e. from Greenwich meridian eastward, and South Pole northward).

 The input vector field is defined in spherical coordinates (it
 has eastward and northward components) while the output
 vector field is rotated to map projection coordinates (relative
 to x and y). The magnitude of the vector is preserved.

 if returnxy=True, the x and y values of the native map projection grid
 are also returned (default False).

 See interp documentation for meaning of extra keyword arguments (**kwargs).

 vectors on a lat/lon grid must be transformed to map projection coordinates
 before they be plotted on the map (with the quiver class method).
        """
        lonsout, latsout, x, y = self.makegrid(nx,ny,returnxy=True)
        # interpolate to map projection coordinates.
        uin = interp(uin,lons,lats,lonsout,latsout,**kwargs)
        vin = interp(vin,lons,lats,lonsout,latsout,**kwargs)
        # rotate from geographic to map coordinates.
        delta = 0.1 # incement in latitude used to estimate derivatives.
        xn,yn = self(lonsout,pylab.where(latsout+delta<90.,latsout+delta,latsout-delta))
        dxdlat = pylab.where(latsout+delta<90.,(xn-x)/(latsout+delta),(x-xn)/(latsout+delta))
        dydlat = pylab.where(latsout+delta<90.,(yn-y)/(latsout+delta),(y-yn)/(latsout+delta))
        # northangle is the angle between true north and the y axis.
        northangle = pylab.arctan2(dxdlat,dydlat)
        uout = uin*pylab.cos(northangle) + vin*pylab.sin(northangle)
        vout = vin*pylab.cos(northangle) - uin*pylab.sin(northangle)
        if returnxy:
            return uout,vout,x,y
        else:
            return uout,vout

    def set_axes_limits(self):
        """
 Set axis limits for map domain using current axes instance.
        """
        # get current axes instance.
        ax = pylab.gca()
        corners = ((self.llcrnrx,self.llcrnry), (self.urcrnrx,self.urcrnry))
        ax.update_datalim( corners )                                          
        ax.set_xlim((self.llcrnrx, self.urcrnrx))
        ax.set_ylim((self.llcrnry, self.urcrnry))

    def scatter(self, *args, **kwargs):
        """
 Plot points with markers on the map (see pylab scatter documentation).
        """
        pylab.scatter(*args, **kwargs)
        self.set_axes_limits()

    def plot(self, *args, **kwargs):
        """
 Draw lines and/or markers on the map (see pylab plot documentation).
        """
        pylab.plot(*args, **kwargs)
        self.set_axes_limits()

    def imshow(self, *args, **kwargs):
        """
 Display an image over the map (see pylab imshow documentation).
 extent and origin keywords set automatically so image will be drawn
 over map region.
        """
        kwargs['extent']=(self.llcrnrx,self.urcrnrx,self.llcrnry,self.urcrnry)
        kwargs['origin']='lower'
        return pylab.imshow(*args,  **kwargs)

    def pcolor(self, *args, **kwargs):
        """
 Make a pseudo-color plot over the map (see pylab pcolor documentation).
        """
        pylab.pcolor(*args, **kwargs)
        self.set_axes_limits()

    def contour(self, *args, **kwargs):
        """
 Make a contour plot over the map (see pylab contour documentation).
        """
        levels, colls = pylab.contour(*args, **kwargs)
        self.set_axes_limits()
        return levels,colls

    def contourf(self, *args, **kwargs):
        """
 Make a filled contour plot over the map (see pylab documentation).
        """
        levels, colls = pylab.contourf(*args, **kwargs)
        self.set_axes_limits()
        return levels,colls

    def quiver(self, x, y, u, v, scale=None, **kwargs):
        """
 Make a vector plot (u, v) with arrows on the map projection grid (x,y)
 If scale is specified, it is used to scale the vectors. If scale=None 
 (default) arrows are scaled to longest one is equal to the maximum
 distance between grid points.   

 Extra keyword arguments (**kwargs) passed to pylab.quiver (see pylab 
 quiver documentation for details).
        """
        ny = x.shape[0]; nx = x.shape[1]
        if scale is None:
            scale = max([(self.xmax-self.xmin)/(nx-1),(self.ymax-self.ymin)/(ny-1)])
        else:
            scale = scale
        pylab.quiver(x,y,u,v,scale, **kwargs)
        self.set_axes_limits()

def interp(datain,lonsin,latsin,lonsout,latsout,checkbounds=False,mode='nearest',cval=0.0,order=3):
    """
 dataout = interp(datain,lonsin,latsin,lonsout,latsout,mode='constant',cval=0.0,order=3)

 interpolate data (datain) on a rectilinear lat/lon grid (with lons=lonsin
 lats=latsin) to a grid with lons=lonsout, lats=latsout.

 datain is a rank-2 array with 1st dimension corresponding to longitude,
 2nd dimension latitude.

 lonsin, latsin are rank-1 arrays containing longitudes and latitudes
 of datain grid in increasing order (i.e. from Greenwich meridian eastward, and
 South Pole northward)

 lonsout, latsout are rank-2 arrays containing lons and lats of desired
 output grid (typically a native map projection grid).

 If checkbounds=True, values of lonsout and latsout are checked to see that
 they lie within the range specified by lonsin and latsin.  Default is
 False, and values outside the borders are handled in the manner described
 by the 'mode' parameter (default mode='nearest', which means the nearest
 boundary value is used). See section 20.2 of the numarray docs for 
 information on the 'mode' keyword.

 See numarray.nd_image.map_coordinates documentation for information on
 the other optional keyword parameters.  The order keyword can be 0 
 for nearest neighbor interpolation (nd_image only allows 1-6) - if
 order=0 bounds checking is done even if checkbounds=False.
    """
    # lonsin and latsin must be monotonically increasing.
    if lonsin[-1]-lonsin[0] < 0 or latsin[-1]-latsin[0] < 0:
        raise ValueError, 'lonsin and latsin must be increasing!'
    # optionally, check that lonsout,latsout are 
    # within region defined by lonsin,latsin.
    # (this check is always done if nearest neighbor 
    # interpolation (order=0) requested).
    if checkbounds or order == 0:
        if min(pylab.ravel(lonsout)) < min(lonsin) or \
           max(pylab.ravel(lonsout)) > max(lonsin) or \
           min(pylab.ravel(latsout)) < min(latsin) or \
           max(pylab.ravel(latsout)) > max(latsin):
            raise ValueError, 'latsout or lonsout outside range of latsin or lonsin'
    # compute grid coordinates of output grid.
    delon = lonsin[1:]-lonsin[0:-1]
    delat = latsin[1:]-latsin[0:-1]
    if max(delat)-min(delat) < 1.e-4 and max(delon)-min(delon) < 1.e-4:
        # regular input grid.
        xcoords = (len(lonsin)-1)*(lonsout-lonsin[0])/(lonsin[-1]-lonsin[0])
        ycoords = (len(latsin)-1)*(latsout-latsin[0])/(latsin[-1]-latsin[0])
    else:
        # irregular (but still rectilinear) input grid.
        lonsoutflat = pylab.ravel(lonsout)
        latsoutflat = pylab.ravel(latsout)
        ix = pylab.searchsorted(lonsin,lonsoutflat)-1
        iy = pylab.searchsorted(latsin,latsoutflat)-1
        xcoords = pylab.zeros(ix.shape,'f')
        ycoords = pylab.zeros(iy.shape,'f')
        for n,i in enumerate(ix):
            if i < 0:
                xcoords[n] = -1 # outside of range on lonsin (lower end)
            elif i >= len(lonsin)-1:
                xcoords[n] = len(lonsin) # outside range on upper end.
            else:
                xcoords[n] = float(i)+(lonsoutflat[n]-lonsin[i])/(lonsin[i+1]-lonsin[i])
        xcoords = pylab.reshape(xcoords,lonsout.shape)
        for m,j in enumerate(iy):
            if j < 0:
                ycoords[m] = -1 # outside of range of latsin (on lower end)
            elif j >= len(latsin)-1:
                ycoords[m] = len(latsin) # outside range on upper end
            else:
                ycoords[m] = float(j)+(latsoutflat[m]-latsin[j])/(latsin[j+1]-latsin[j])
        ycoords = pylab.reshape(ycoords,latsout.shape)
    coords = [ycoords,xcoords]
    # interpolate to output grid using numarray.nd_image spline filter.
    if order:
        return nd_image.map_coordinates(datain,coords,mode=mode,cval=cval,order=order)
    else:
        # nearest neighbor interpolation if order=0.
        # uses index arrays, so first convert to numarray.
        datatmp = pylab.array(datain,datain.typecode())
        xi = pylab.around(xcoords).astype('i')
        yi = pylab.around(ycoords).astype('i')
        return datatmp[yi,xi]

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
    if pylab.fabs(lonsin[-1]-lonsin[0]-360.) > 1.e-4:
        raise ValueError, 'cyclic point not included'
    if lon0 < lonsin[0] or lon0 > lonsin[-1]:
        raise ValueError, 'lon0 outside of range of lonsin'
    i0 = pylab.argsort(pylab.fabs(lonsin-lon0))[0]
    dataout = pylab.zeros(datain.shape,datain.typecode())
    lonsout = pylab.zeros(lonsin.shape,lonsin.typecode())
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
   arrout  = pylab.zeros((nlats,nlons+1),arrin.typecode())
   arrout[:,0:nlons] = arrin[:,:]
   arrout[:,nlons] = arrin[:,0]
   lonsout = pylab.zeros(nlons+1,lonsin.typecode())
   lonsout[0:nlons] = lonsin[:]
   lonsout[nlons]  = lonsin[-1] + lonsin[1]-lonsin[0]
   return arrout,lonsout
