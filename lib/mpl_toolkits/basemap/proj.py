import numpy as np
import pyproj
import math

__version__ = '1.2.2'
_dg2rad = math.radians(1.)
_rad2dg = math.degrees(1.)

_cylproj = ['cyl','merc','mill','gall']
_pseudocyl = ['moll','robin','sinu','mbtfpq']

_upper_right_out_of_bounds = (
    'the upper right corner of the plot is not in the map projection region')

_lower_left_out_of_bounds = (
    'the lower left corner of the plot is not in the map projection region')


class Proj(object):
    """
    peforms cartographic transformations (converts from longitude,latitude
    to native map projection x,y coordinates and vice versa) using proj
    (http://proj.maptools.org/)
    Uses a pyrex generated C-interface to libproj.

    __init__ method sets up projection information.
    __call__ method compute transformations.
    See docstrings for __init__ and __call__ for details.

    Contact: Jeff Whitaker <jeffrey.s.whitaker@noaa.gov>
    """

    def __init__(self,projparams,llcrnrlon,llcrnrlat,
                      urcrnrlon,urcrnrlat,urcrnrislatlon=True):
        """
        initialize a Proj class instance.

        Input 'projparams' is a dictionary containing proj map
        projection control parameter key/value pairs.
        See the proj documentation (http://www.remotesensing.org/proj/)
        for details.

        llcrnrlon,llcrnrlat are lon and lat (in degrees) of lower
        left hand corner of projection region.

        urcrnrlon,urcrnrlat are lon and lat (in degrees) of upper
        right hand corner of projection region if urcrnrislatlon=True
        (default). Otherwise, urcrnrlon,urcrnrlat are x,y in projection
        coordinates (units meters), assuming the lower left corner is x=0,y=0.
        """
        self.projparams = projparams
        self.projection = projparams['proj']
        # rmajor is the semi-major axis.
        # rminor is the semi-minor axis.
        # esq is eccentricity squared.
        try:
            self.rmajor = projparams['a']
            self.rminor = projparams['b']
        except:
            try:
                self.rmajor = projparams['R']
            except:
                self.rmajor = projparams['bR_a']
            self.rminor = self.rmajor
        if self.rmajor == self.rminor:
            self.ellipsoid = False
        else:
            self.ellipsoid = True
        self.flattening = (self.rmajor-self.rminor)/self.rmajor
        self.esq = (self.rmajor**2 - self.rminor**2)/self.rmajor**2
        self.llcrnrlon = llcrnrlon
        self.llcrnrlat = llcrnrlat
        if self.projection not in ['ortho','geos','cyl'] + _pseudocyl:
            self._proj4 = pyproj.Proj(projparams)
            llcrnrx, llcrnry = self(llcrnrlon,llcrnrlat)
        elif self.projection == 'cyl':
            llcrnrx = llcrnrlon
            llcrnry = llcrnrlat
        elif self.projection == 'ortho':
            if (llcrnrlon == -180 and llcrnrlat == -90 and
                urcrnrlon == 180 and urcrnrlat == 90):
                self._fulldisk = True
                self._proj4 = pyproj.Proj(projparams)
                llcrnrx = -self.rmajor
                llcrnry = -self.rmajor
                urcrnrx = -llcrnrx
                urcrnry = -llcrnry
            else:
                self._fulldisk = False
                self._proj4 = pyproj.Proj(projparams)
                llcrnrx, llcrnry = self(llcrnrlon,llcrnrlat)
                if llcrnrx > 1.e20 or llcrnry > 1.e20:
                    raise ValueError(_lower_left_out_of_bounds)
        elif self.projection == 'geos':
            self._proj4 = pyproj.Proj(projparams)
            # find major and minor axes of ellipse defining map proj region.
            delta = 0.01
            lats = np.arange(0,90,delta)
            lon_0 = projparams['lon_0']
            lons = lon_0*np.ones(len(lats),'d')
            x, y = self._proj4(lons, lats)
            yi = (y > 1.e20).tolist()
            ny = yi.index(1)-1
            height = y[ny]
            lons = np.arange(lon_0,lon_0+90,delta)
            lats = np.zeros(len(lons),'d')
            x, y = self(lons, lats)
            xi = (x > 1.e20).tolist()
            nx = xi.index(1)-1
            width = x[nx]
            self._height = height
            self._width = width
            if (llcrnrlon == -180 and llcrnrlat == -90 and
                urcrnrlon == 180 and urcrnrlat == 90):
                self._fulldisk = True
                llcrnrx = -width
                llcrnry = -height
                urcrnrx = -llcrnrx
                urcrnry = -llcrnry
            else:
                self._fulldisk = False
                llcrnrx, llcrnry = self(llcrnrlon,llcrnrlat)
                if llcrnrx > 1.e20 or llcrnry > 1.e20:
                    raise ValueError(_lower_left_out_of_bounds)
        elif self.projection in _pseudocyl:
            self._proj4 = pyproj.Proj(projparams)
            xtmp,urcrnry = self(projparams['lon_0'],90.)
            urcrnrx,xtmp = self(projparams['lon_0']+180.,0)
            llcrnrx = -urcrnrx
            llcrnry = -urcrnry
        # compute x_0, y_0 so ll corner of domain is x=0,y=0.
        # note that for 'cyl' x,y == lon,lat
        self.projparams['x_0']=-llcrnrx
        self.projparams['y_0']=-llcrnry
        # reset with x_0, y_0.
        if self.projection != 'cyl':
            self._proj4 = pyproj.Proj(projparams)
            llcrnry = 0.
            llcrnrx = 0.
        else:
            llcrnrx = llcrnrlon
            llcrnry = llcrnrlat
        if urcrnrislatlon:
            self.urcrnrlon = urcrnrlon
            self.urcrnrlat = urcrnrlat
            if self.projection not in ['ortho','geos'] + _pseudocyl:
                urcrnrx,urcrnry = self(urcrnrlon,urcrnrlat)
            elif self.projection == 'ortho':
                if self._fulldisk:
                    urcrnrx = 2.*self.rmajor
                    urcrnry = 2.*self.rmajor
                else:
                    urcrnrx,urcrnry = self(urcrnrlon,urcrnrlat)
                    if urcrnrx > 1.e20 or urcrnry > 1.e20:
                        raise ValueError(_upper_right_out_of_bounds)
            elif self.projection == 'geos':
                if self._fulldisk:
                    urcrnrx = 2.*self._width
                    urcrnry = 2.*self._height
                else:
                    urcrnrx,urcrnry = self(urcrnrlon,urcrnrlat)
                    if urcrnrx > 1.e20 or urcrnry > 1.e20:
                        raise ValueError(_upper_right_out_of_bounds)
            elif self.projection in _pseudocyl:
                xtmp,urcrnry = self(projparams['lon_0'],90.)
                urcrnrx,xtmp = self(projparams['lon_0']+180.,0)
        else:
            urcrnrx = urcrnrlon
            urcrnry = urcrnrlat
            urcrnrlon, urcrnrlat = self(urcrnrx, urcrnry, inverse=True)
            self.urcrnrlon = urcrnrlon
            self.urcrnrlat = urcrnrlat
        # corners of domain.
        self.llcrnrx = llcrnrx
        self.llcrnry = llcrnry
        self.urcrnrx = urcrnrx
        self.urcrnry = urcrnry
        if urcrnrx > llcrnrx:
            self.xmin = llcrnrx
            self.xmax = urcrnrx
        else:
            self.xmax = llcrnrx
            self.xmin = urcrnrx
        if urcrnry > llcrnry:
            self.ymin = llcrnry
            self.ymax = urcrnry
        else:
            self.ymax = llcrnry
            self.ymin = urcrnry

    def __call__(self, *args, **kw):
        # x,y,inverse=False):
        """
        Calling a Proj class instance with the arguments lon, lat will
        convert lon/lat (in degrees) to x/y native map projection
        coordinates (in meters).  If optional keyword 'inverse' is
        True (default is False), the inverse transformation from x/y
        to lon/lat is performed.

        For cylindrical equidistant projection ('cyl'), this
        does nothing (i.e. x,y == lon,lat).

        lon,lat can be either scalar floats or N arrays.
        """
        if len(args) == 1:
            xy = args[0]
            onearray = True
        else:
            x,y = args
            onearray = False
        if self.projection == 'cyl': # for cyl x,y == lon,lat
            if onearray:
                return xy
            else:
                return x,y
        inverse = kw.get('inverse', False)
        if onearray:
            outxy = self._proj4(xy, inverse=inverse)
        else:
            outx,outy = self._proj4(x, y, inverse=inverse)
        if inverse:
            if self.projection in ['merc','mill','gall']:
                if self.projection == 'merc':
                    coslat = math.cos(math.radians(self.projparams['lat_ts']))
                    sinlat = math.sin(math.radians(self.projparams['lat_ts']))
                else:
                    coslat = 1.
                    sinlat = 0.
                # radius of curvature of the ellipse perpendicular to
                # the plane of the meridian.
                rcurv = self.rmajor*coslat/math.sqrt(1.-self.esq*sinlat**2)
                if onearray:
                    outxy[:,0] = _rad2dg*(xy[:,0]/rcurv) + self.llcrnrlon
                else:
                    try: # x a scalar or an array
                        outx = _rad2dg*(x/rcurv) + self.llcrnrlon
                    except: # x a sequence
                        outx = [_rad2dg*(xi/rcurv) + self.llcrnrlon for xi in x]
        else:
            if self.projection in ['merc','mill','gall']:
                if self.projection == 'merc':
                    coslat = math.cos(math.radians(self.projparams['lat_ts']))
                    sinlat = math.sin(math.radians(self.projparams['lat_ts']))
                else:
                    coslat = 1.
                    sinlat = 0.
                # radius of curvature of the ellipse perpendicular to
                # the plane of the meridian.
                rcurv = self.rmajor*coslat/math.sqrt(1.-self.esq*sinlat**2)
                if onearray:
                    outxy[:,0] = rcurv*_dg2rad*(xy[:,0]-self.llcrnrlon)
                else:
                    try: # x is a scalar or an array
                        outx = rcurv*_dg2rad*(x-self.llcrnrlon)
                    except: # x is a sequence.
                        outx = [rcurv*_dg2rad*(xi-self.llcrnrlon) for xi in x]
        if onearray:
            return outxy
        else:
            return outx, outy

    def makegrid(self,nx,ny,returnxy=False):
        """
        return arrays of shape (ny,nx) containing lon,lat coordinates of
        an equally spaced native projection grid.
        if returnxy=True, the x,y values of the grid are returned also.
        """
        dx = (self.urcrnrx-self.llcrnrx)/(nx-1)
        dy = (self.urcrnry-self.llcrnry)/(ny-1)
        x = self.llcrnrx+dx*np.indices((ny,nx),np.float32)[1,:,:]
        y = self.llcrnry+dy*np.indices((ny,nx),np.float32)[0,:,:]
        lons, lats = self(x, y, inverse=True)
        if returnxy:
            return lons, lats, x, y
        else:
            return lons, lats

    def makegrid3d(self,nx,ny,returnxy=False):
        """
        return array of shape (ny,nx, 2) containing lon,lat coordinates of
        an equally spaced native projection grid.
        if returnxy=True, the x,y values of the grid are returned also.
        """
        dx = (self.urcrnrx-self.llcrnrx)/(nx-1)
        dy = (self.urcrnry-self.llcrnry)/(ny-1)
        xy = np.empty((ny,nx,2), np.float64)
        xy[...,0] = self.llcrnrx+dx*np.indices((ny,nx),np.float32)[1,:,:]
        xy[...,1] = self.llcrnry+dy*np.indices((ny,nx),np.float32)[0,:,:]
        lonlat = self(xy, inverse=True)
        if returnxy:
            return lonlat, xy
        else:
            return lonlat


if __name__ == "__main__":

    params = {}
    params['proj'] = 'lcc'
    params['R'] = 6371200
    params['lat_1'] = 50
    params['lat_2'] = 50
    params['lon_0'] = -107
    nx = 349; ny = 277; dx = 32463.41; dy = dx
    awips221 = Proj(params,-145.5,1.0,(nx-1)*dx,(ny-1)*dy,urcrnrislatlon=False)
    # AWIPS grid 221 parameters
    # (from http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html)
    llcornerx, llcornery = awips221(-145.5,1.)
    # find 4 lon/lat corners of AWIPS grid 221.
    llcornerx = 0.; llcornery = 0.
    lrcornerx = dx*(nx-1); lrcornery = 0.
    ulcornerx = 0.; ulcornery = dy*(ny-1)
    urcornerx = dx*(nx-1); urcornery = dy*(ny-1)
    llcornerlon, llcornerlat = awips221(llcornerx, llcornery, inverse=True)
    lrcornerlon, lrcornerlat = awips221(lrcornerx, lrcornery, inverse=True)
    urcornerlon, urcornerlat = awips221(urcornerx, urcornery, inverse=True)
    ulcornerlon, ulcornerlat = awips221(ulcornerx, ulcornery, inverse=True)
    print '4 corners of AWIPS grid 221:'
    print llcornerlon, llcornerlat
    print lrcornerlon, lrcornerlat
    print urcornerlon, urcornerlat
    print ulcornerlon, ulcornerlat
    print 'from GRIB docs'
    print '(see http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html)'
    print '   -145.5  1.0'
    print '   -68.318 0.897'
    print '   -2.566 46.352'
    print '   148.639 46.635'
    # compute lons and lats for the whole AWIPS grid 221 (377x249).
    import time; t1 = time.clock()
    lons, lats = awips221.makegrid(nx,ny)
    t2 = time.clock()
    print 'compute lats/lons for all points on AWIPS 221 grid (%sx%s)' %(nx,ny)
    print 'max/min lons'
    print min(np.ravel(lons)),max(np.ravel(lons))
    print 'max/min lats'
    print min(np.ravel(lats)),max(np.ravel(lats))
    print 'took',t2-t1,'secs'
    print 'Same thing but with a single 3-D array'
    t1 = time.clock()
    lonlat, xy = awips221.makegrid3d(nx,ny, returnxy=True)
    t2 = time.clock()
    print 'took',t2-t1,'secs'

    assert (lons==lonlat[...,0]).all(), "The longitudes are different"
    assert (lats==lonlat[...,1]).all(), "The latitudes are different"
