from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
import numarray as N
from numarray import nd_image
import sys, os, MLab
from proj import Proj

_datadir = os.path.join(sys.prefix,'share/basemap')

class Basemap:

    """
 Set up a basemap with a given map projection (cylindrical equidistant,
 mercator, lambert conformal conic, lambert azimuthal equal area,
 albers equal area conic and stereographic are currently available).
 Doesn't actually draw anything, but sets up the map projection class and
 creates the coastline and political boundary polygons in native map 
 projection coordinates.  Requires matplotlib and numarray.
 Uses a pyrex interface to C-code from proj.4 (http://proj.maptools.org).
 
 Useful instance variables:
 
 projection - map projection ('cyl','merc','lcc','aea','laea' or 'stere')
 aspect - map aspect ratio (size of y dimension / size of x dimension).
 llcrnrlon - longitude of lower left hand corner of the desired map domain.
 llcrnrlon - latitude of lower left hand corner of the desired map domain.      
 urcrnrlon - longitude of upper right hand corner of the desired map domain.
 urcrnrlon - latitude of upper right hand corner of the desired map domain.
 llcrnrx,llcrnry,urcrnrx,urcrnry - corners of map domain in projection coordinates.

 Example Usage:
 (this example plus others can be run by running test.py in the examples dir)

>>> from matplotlib.toolkits.basemap import Basemap
>>> import cPickle
>>> from pylab import *
>>> # read in topo data from pickle (on a regular lat/lon grid)
>>> topodict = cPickle.load(open('etopo20.pickle','rb'))
>>> etopo = topodict['data']; lons = topodict['lons']; lats = topodict['lats']
>>> m = Basemap(lons[0],lats[0],lons[-1],lats[-1])
>>> xsize = rcParams['figure.figsize'][0]
>>> fig=figure(figsize=(xsize,m.aspect*xsize))
>>> fig.add_axes([0.1,0.1,0.8,0.8])
>>> im = imshow(etopo,extent=(m.llcrnrx,m.urcrnrx,m.llcrnry,m.urcrnry),origin='lower')
>>> ax = gca() # get current axis instance
>>> # draw coastlines and fill continents.
>>> m.drawcoastlines(ax)
>>> m.fillcontinents(ax)
>>> # draw parallels
>>> circles = arange(-90.,120.,30.)
>>> m.drawparallels(ax,circles)
>>> # draw meridians
>>> meridians = arange(0.,390.,60.)
>>> m.drawmeridians(ax,meridians)
>>> ax.set_xticks([]) # no ticks
>>> ax.set_yticks([])
>>> title('Cylindrical Equidistant')
>>> show()

 Version: 0.1 (20050203)
 Contact: Jeff Whitaker <jeffrey.s.whitaker@noaa.gov>
    """

    def __init__(self,llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat,\
        resolution='c',area_thresh=10000.,projection='cyl',rsphere=6371009.,\
        lat_ts=None,lat_1=None,lat_2=None,lat_0=None,lon_0=None):
        """
 create a Basemap instance.
 
 mandatory input arguments:
 
 llcrnrlon - longitude of lower left hand corner of the desired map domain.
 llcrnrlon - latitude of lower left hand corner of the desired map domain.      
 urcrnrlon - longitude of upper right hand corner of the desired map domain.
 urcrnrlon - latitude of upper right hand corner of the desired map domain.
 
 optional keyword parameters:
 
 resolution - resolution of coastline database to use. Can be 'c' (crude, 
  roughly 25 km resolution) or 'l' (low, roughly 5 km resolution). Default 'c'.
  Coastline data is from the GSHHS 
  (http://www.soest.hawaii.edu/wessel/gshhs/gshhs.html).
 area_thresh - coastline with an area smaller than area_thresh in km^2
  will not be plotted.  Default 10,000.
 projection - map projection.  'cyl' - cylindrical equidistant, 'merc' -
  mercator, 'lcc' - lambert conformal conic, 'stere' - stereographic,
  'aea' - albers equal area conic, and 
  'laea' - lambert azimuthal equal area currently available.  Default 'cyl'.
 rsphere - radius of the sphere used to define map projection (default
  6371009 meters, close to the arithmetic mean radius of the earth).
  
 The following parameters are map projection parameters which all default to 
 None.  Not all parameters are used by all projections, some are ignored.
 
 lat_ts - latitude of natural origin (used for mercator and stereographic
  projections).
 lat_1 - first standard parallel for lambert conformal and albers
  equal area projections.
 lat_2 - second standard parallel for lambert conformal and albers
  equal area projections.
 lat_0 - central latitude (y-axis origin) - used by stereographic and
  lambert azimuthal projections).
 lon_0 - central meridian (x-axis origin - used by lambert conformal
  and lambert azimuthal and stereographic projections).
        """     

        # read in coastline data.
        coastlons = []; coastlats = []; coastsegind = []; coastsegarea = []
        i = 0  # the current ind
        for line in open(os.path.join(_datadir,'gshhs_'+resolution+'.txt')):
            linesplit = line.split()
            if line.startswith('P'):
                coastsegind.append(i)
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
            if None or lat_0 is None or lon_0 is None:
                raise ValueError, 'must specify lat_0 and lon_0 for Lambert Azimuthal basemap'
            projparams['lat_0'] = lat_0
            projparams['lon_0'] = lon_0
            proj = Proj(projparams,llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat)
        elif projection == 'merc':
            if lat_ts is None:
                raise ValueError, 'must specify lat_ts for Mercator basemap'
            projparams['lat_ts'] = lat_ts
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
        if projection == 'cyl' or projection == 'merc':
            self.aspect = (urcrnrlat-llcrnrlat)/(urcrnrlon-llcrnrlon)
        else:
            self.aspect = (proj.ymax-proj.ymin)/(proj.xmax-proj.xmin)
        self.llcrnrx = proj.llcrnrx
        self.llcrnry = proj.llcrnry
        self.urcrnrx = proj.urcrnrx
        self.urcrnry = proj.urcrnry

        # transform coastline polygons to native map coordinates.
        xc,yc = proj(N.array(coastlons,'f'),N.array(coastlats,'f'))
        xc2,yc2 = proj(N.array(coastlons2,'f'),N.array(coastlats,'f'))
        xc3,yc3 = proj(N.array(coastlons3,'f'),N.array(coastlats,'f'))
        if projection == 'merc': yc2=yc
        # set up segments in form needed for LineCollection,
        # ignoring 'inf' values that are off the map, and skipping
        # polygons that have an area > area_thresh..
        segments = [zip(xc[i0:i1],yc[i0:i1]) for a,i0,i1 in zip(coastsegarea[:-1],coastsegind[:-1],coastsegind[1:]) if a > area_thresh]
        segments2 = [zip(xc2[i0:i1],yc2[i0:i1]) for a,i0,i1 in zip(coastsegarea[:-1],coastsegind[:-1],coastsegind[1:]) if a > area_thresh and max(xc2[i0:i1]) < 1.e20 and max(yc2[i0:i1]) < 1.e20]
        segments3 = [zip(xc3[i0:i1],yc3[i0:i1]) for a,i0,i1 in zip(coastsegarea[:-1],coastsegind[:-1],coastsegind[1:]) if a > area_thresh and max(xc3[i0:i1]) < 1.e20 and max(yc3[i0:i1]) < 1.e20]
        self.coastsegs = segments+segments2+segments3

        # same as above for country polygons.
        xc,yc = proj(N.array(cntrylons,'f'),N.array(cntrylats,'f'))
        xc2,yc2 = proj(N.array(cntrylons2,'f'),N.array(cntrylats,'f'))
        xc3,yc3 = proj(N.array(cntrylons2,'f'),N.array(cntrylats,'f'))
        segments = [zip(xc[i0:i1],yc[i0:i1]) for i0,i1 in zip(cntrysegind[:-1],cntrysegind[1:])]
        segments2 = [zip(xc2[i0:i1],yc2[i0:i1]) for i0,i1 in zip(cntrysegind[:-1],cntrysegind[1:]) if max(xc2[i0:i1]) < 1.e20 and max(yc2[i0:i1]) < 1.e20]
        segments3 = [zip(xc3[i0:i1],yc3[i0:i1]) for i0,i1 in zip(cntrysegind[:-1],cntrysegind[1:]) if max(xc3[i0:i1]) < 1.e20 and max(yc3[i0:i1]) < 1.e20]
        self.cntrysegs = segments+segments2+segments3

        # same as above for state polygons.
        xc,yc = proj(N.array(statelons,'f'),N.array(statelats,'f'))
        xc2,yc2 = proj(N.array(statelons2,'f'),N.array(statelats,'f'))
        xc3,yc3 = proj(N.array(statelons3,'f'),N.array(statelats,'f'))
        segments = [zip(xc[i0:i1],yc[i0:i1]) for i0,i1 in zip(statesegind[:-1],statesegind[1:])]
        segments2 = [zip(xc2[i0:i1],yc2[i0:i1]) for i0,i1 in zip(statesegind[:-1],statesegind[1:]) if max(xc2[i0:i1]) < 1.e20 and max(yc2[i0:i1]) < 1.e20]
        segments3 = [zip(xc3[i0:i1],yc3[i0:i1]) for i0,i1 in zip(statesegind[:-1],statesegind[1:]) if max(xc3[i0:i1]) < 1.e20 and max(yc3[i0:i1]) < 1.e20]
        self.statesegs = segments+segments2+segments3

        # store coast polygons for filling.
        # special treatment (kludge) for Antarctica.
        self.coastpolygons=[]
        if projection == 'merc': 
            xsp,ysp = proj(0.,-89.9) # s. pole coordinates.
            xa,ya = proj(0.,-68.0) # edge of antarctica.
        for seg in self.coastsegs:
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
            elif projection == 'merc':
                if x[-1] == 0.000 and y[-1] < ya: # close antarctica
                    x.append(0.)
                    y.append(ysp)
                    x.insert(0,360.)
                    y.insert(0,ysp)
                if x[-1] == 360.000 and y[-1] < ya: 
                    x.append(360.)
                    y.append(ysp)
                    x.insert(0,720.)
                    y.insert(0,ysp)
                if x[-1] == -360.000 and y[-1] < ya: 
                    x.append(-360.)
                    y.append(ysp)
                    x.insert(0,0.)
                    y.insert(0,ysp)
            self.coastpolygons.append((x,y))

    def __call__(self,x,y,inverse=False):
        """
 Calling a Basemap class instance with the arguments lon, lat will
 convert lon/lat (in degrees) to x/y native map projection 
 coordinates (in meters).  If optional keyword 'inverse' is
 True (default is False), the inverse transformation from x/y
 to lon/lat is performed.

 For cylindrical equidistant projection ('cyl'), this
 does nothing (i.e. x,y == lon,lat).

 For mercator projection ('merc'), x == lon, but y has units
 of meters.

 lon,lat can be either scalar floats or N arrays.
        """
        return self.projtran(x,y,inverse=inverse)
 
    def makegrid(self,nx,ny):
        """
 return arrays of shape (ny,nx) containing lon,lat coordinates of
 an equally spaced native projection grid.
        """
        return self.projtran.makegrid(nx,ny)

    def fillcontinents(self,ax,color=0.8):
        """
 Fill continents.

 ax - current axis instance.
 color - color to fill continents (default gray).
        """
        # define corners of map domain.
        p1 = (self.llcrnrx,self.llcrnry); p2 = (self.urcrnrx,self.urcrnry)
        p3 = (self.llcrnrx,self.urcrnry); p4 = (self.urcrnrx,self.llcrnry)
        for x,y in self.coastpolygons:
            xa = N.array(x,'f')
            ya = N.array(y,'f')
        # clip to map domain.
            xa = N.clip(xa, self.xmin, self.xmax)
            ya = N.clip(ya, self.ymin, self.ymax)
        # check to see if all four corners of domain in polygon (if so,
        # don't draw since it will just fill in the whole map).
            test1 = N.fabs(xa-self.xmax) < 10.
            test2 = N.fabs(xa-self.xmin) < 10.
            test3 = N.fabs(ya-self.ymax) < 10.
            test4 = N.fabs(ya-self.ymin) < 10.
            hasp1 = sum(test1*test3)
            hasp2 = sum(test2*test3)
            hasp4 = sum(test2*test4)
            hasp3 = sum(test1*test4)
            if not hasp1 or not hasp2 or not hasp3 or not hasp4:
                xy = zip(xa.tolist(),ya.tolist())
                poly = Polygon(xy,facecolor=color,edgecolor=color,linewidth=0)
                ax.add_patch(poly)

    def drawcoastlines(self,ax,linewidth=1.,color='k',antialiased=1):
        """
 Draw coastlines.

 ax - current axis instance.
 linewidth - coastline width (default 1.)
 color - coastline color (default black)
 antialiased - antialiasing switch for coastlines (default True).
        """
        coastlines = LineCollection(self.coastsegs,antialiaseds=(antialiased,))
        try:
            coastlines.set_color(color)
        except: # this was a typo that existed in matplotlib-0.71 and earlier
            coastlines.color(color)
        coastlines.set_linewidth(linewidth)
        ax.add_collection(coastlines)

    def drawcountries(self,ax,linewidth=0.5,color='k',antialiased=1):
        """
 Draw country boundaries.

 ax - current axis instance.
 linewidth - country boundary line width (default 0.5)
 color - country boundary line color (default black)
 antialiased - antialiasing switch for country boundaries (default True).
        """
        coastlines = LineCollection(self.cntrysegs,antialiaseds=(antialiased,))
        try:
            coastlines.set_color(color)
        except: # this was a typo that existed in matplotlib-0.71 and earlier
            coastlines.color(color)
        coastlines.set_linewidth(linewidth)
        ax.add_collection(coastlines)

    def drawstates(self,ax,linewidth=0.5,color='k',antialiased=1):
        """
 Draw state boundaries in Americas.

 ax - current axis instance.
 linewidth - state boundary line width (default 0.5)
 color - state boundary line color (default black)
 antialiased - antialiasing switch for state boundaries (default True).
        """
        coastlines = LineCollection(self.statesegs,antialiaseds=(antialiased,))
        try:
            coastlines.set_color(color)
        except: # this was a typo that existed in matplotlib-0.71 and earlier
            coastlines.color(color)
        coastlines.set_linewidth(linewidth)
        ax.add_collection(coastlines)

    def drawparallels(self,ax,circles,color='k',linewidth=1., \
                      linestyle='--',dashes=[1,1]):
        """
 draw parallels (latitude lines).

 ax - current axis instance.
 circles - list containing latitude values to draw (in degrees).
 color - color to draw parallels (default black).
 linewidth - line width for parallels (default 1.)
 linestyle - line style for parallels (default '--', i.e. dashed).
 dashes - dash pattern for parallels (default [1,1], i.e. 1 pixel on,
  1 pixel off).
        """
        if self.projection in ['merc','cyl']:
            lons = N.arange(self.llcrnrlon,self.urcrnrlon,1).astype('f')
        else:
            lons = N.arange(0,362,1).astype('f')
        # make sure 80 degree parallel is drawn if projection not merc or cyl
        try:
            circlesl = circles.tolist()
        except:
            circlesl = circles
        if self.projection not in ['merc','cyl'] and 80. not in circlesl: 
            circlesl.append(80.)
        xdelta = 0.1*(self.xmax-self.xmin)
        ydelta = 0.1*(self.ymax-self.ymin)
        for circ in circlesl:
            lats = circ*N.ones(len(lons),'f')
            x,y = self(lons,lats)
            # remove points outside domain.
            testx = N.logical_and(x>=self.xmin-xdelta,x<=self.xmax+xdelta)
            x = N.compress(testx, x)
            y = N.compress(testx, y)
            testy = N.logical_and(y>=self.ymin-ydelta,y<=self.ymax+ydelta)
            x = N.compress(testy, x)
            y = N.compress(testy, y)
            if len(x) > 1 and len(y) > 1:
                # split into separate line segments if necessary.
                # (not necessary for mercator or cylindrical).
                xd = (x[1:]-x[0:-1])**2
                yd = (y[1:]-y[0:-1])**2
                dist = N.sqrt(xd+yd)
                split = dist > 500000.
                if N.sum(split) and self.projection not in ['merc','cyl']:
                   ind = (N.compress(split,MLab.squeeze(split*N.indices(xd.shape)))+1).tolist()
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

    def drawmeridians(self,ax,meridians,color='k',linewidth=1., \
                      linestyle='--',dashes=[1,1]):
        """
 draw meridians (longitude lines).

 ax - current axis instance.
 meridians - list containing longitude values to draw (in degrees).
 color - color to draw meridians (default black).
 linewidth - line width for meridians (default 1.)
 linestyle - line style for meridians (default '--', i.e. dashed).
 dashes - dash pattern for meridians (default [1,1], i.e. 1 pixel on,
  1 pixel off).
        """
        if self.projection not in ['merc','cyl']:
            lats = N.arange(-80,81).astype('f')
        else:
            lats = N.arange(-90,91).astype('f')
        xdelta = 0.1*(self.xmax-self.xmin)
        ydelta = 0.1*(self.ymax-self.ymin)
        for merid in meridians:
            lons = merid*N.ones(len(lats),'f')
            x,y = self(lons,lats)
            # remove points outside domain.
            testx = N.logical_and(x>=self.xmin-xdelta,x<=self.xmax+xdelta)
            x = N.compress(testx, x)
            y = N.compress(testx, y)
            testy = N.logical_and(y>=self.ymin-ydelta,y<=self.ymax+ydelta)
            x = N.compress(testy, x)
            y = N.compress(testy, y)
            if len(x) > 1 and len(y) > 1:
                # split into separate line segments if necessary.
                # (not necessary for mercator or cylindrical).
                xd = (x[1:]-x[0:-1])**2
                yd = (y[1:]-y[0:-1])**2
                dist = N.sqrt(xd+yd)
                split = dist > 500000.
                if N.sum(split) and self.projection not in ['merc','cyl']:
                   ind = (N.compress(split,MLab.squeeze(split*N.indices(xd.shape)))+1).tolist()
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

def interp(datain,lonsin,latsin,lonsout,latsout,checkbounds=False,mode='nearest',cval=0.0,order=3):
    """
 dataout = interp(datain,lonsin,latsin,lonsout,latsout,mode='constant',cval=0.0,order=3)

 interpolate data (datain) on a rectilinear lat/lon grid (with lons=lonsin
 lats=latsin) to a grid with lons=lonsout, lats=latsout.

 datain is a rank-2 array with 1st dimension corresponding to longitude,
 2nd dimension latitude.

 lonsin, latsin are rank-1 Numeric arrays containing longitudes and latitudes
 of datain grid in increasing order (i.e. from Greenwich meridian eastward, and
 South Pole northward)

 lonsout, latsout are rank-2 Numeric arrays containing lons and lats out desired
 output grid (typically a native map projection grid).

 If checkbounds=True, values of lonsout and latsout are checked to see that
 they lie within the range specified by lonsin and latsing.  Default is
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
        if min(N.ravel(lonsout)) < min(lonsin) or \
           max(N.ravel(lonsout)) > max(lonsin) or \
           min(N.ravel(latsout)) < min(latsin) or \
           max(N.ravel(latsout)) > max(latsin):
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
        lonsoutflat = N.ravel(lonsout)
        latsoutflat = N.ravel(latsout)
        ix = N.searchsorted(lonsin,lonsoutflat)-1
        iy = N.searchsorted(latsin,latsoutflat)-1
        xcoords = N.zeros(ix.shape,'f')
        ycoords = N.zeros(iy.shape,'f')
        for n,i in enumerate(ix):
            if i < 0:
                xcoords[n] = -1 # outside of range on lonsin (lower end)
            elif i >= len(lonsin)-1:
                xcoords[n] = len(lonsin) # outside range on upper end.
            else:
                xcoords[n] = float(i)+(lonsoutflat[n]-lonsin[i])/(lonsin[i+1]-lonsin[i])
        xcoords = N.reshape(xcoords,lonsout.shape)
        for m,j in enumerate(iy):
            if j < 0:
                ycoords[m] = -1 # outside of range of latsin (on lower end)
            elif j >= len(latsin)-1:
                ycoords[m] = len(latsin) # outside range on upper end
            else:
                ycoords[m] = float(j)+(latsoutflat[m]-latsin[j])/(latsin[j+1]-latsin[j])
        ycoords = N.reshape(ycoords,latsout.shape)
    coords = [ycoords,xcoords]
    # interpolate to output grid using numarray.nd_image spline filter.
    if order:
        return nd_image.map_coordinates(datain,coords,mode=mode,cval=cval,order=order)
    else:
        # nearest neighbor interpolation if order=0.
        # uses index arrays, so first convert to numarray.
        datatmp = N.array(datain,datain.typecode())
        xi = N.around(xcoords).astype('i')
        yi = N.around(ycoords).astype('i')
        return datatmp[yi,xi]
