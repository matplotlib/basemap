import pylab as N
import proj4, math

class Proj:
    """
 peforms cartographic transformations (converts from longitude,latitude
 to native map projection x,y coordinates and vice versa) using proj 
 (http://proj.maptools.org/)
 Uses a pyrex generated C-interface to libproj.
 
 __init__ method sets up projection information.
 __call__ method compute transformations.
 See docstrings for __init__ and __call__ for details.

 Version: 0.4.1 (20050509)
 Contact: Jeff Whitaker <jeffrey.s.whitaker@noaa.gov>
    """

    def __init__(self,projparams,llcrnrlon,llcrnrlat,urcrnrlon,urcrnrlat,urcrnrislatlon=True):
        """
 initialize a Proj class instance.

 Input 'projparams' is a dictionary containing proj map
 projection control parameter key/value pairs.
 See the proj documentation (http://www.remotesensing.org/proj/)
 for details.

 llcrnrlon,llcrnrlat are lon and lat (in degrees) of lower left hand corner
 of projection region.

 urcrnrlon,urcrnrlat are lon and lat (in degrees) of upper right hand corner
 of projection region if urcrnrislatlon=True (default). Otherwise, 
 urcrnrlon,urcrnrlat are x,y in projection coordinates (units meters), 
 assuming the lower left corner is x=0,y=0.
        """
        # set units to meters.
        if not projparams.has_key('units'):
            projparams['units']='m'
        elif projparams['units'] != 'm':
            print 'resetting units to meters ...'
            projparams['units']='m'
        self.projparams = projparams
        # make sure proj parameter specified.
        # (no other checking done in proj parameters)
        if 'proj' not in projparams.keys():
            raise KeyError, "need to specify proj parameter"
        if 'R' not in projparams.keys() and 'a' and 'b' not in projparams.keys():
            raise KeyError, "need to specify R (perfect sphere radius), or a and b (major and minor sphere radii)"
        # build proj string.
        self.proj4cmd = []
        for key,value in map(None,projparams.keys(),projparams.values()):
            if key == 'x_0' or key == 'y_0':
                # x_0 and y_0 are set based on llcrnrlon,llcrnrlat.
                pass
            else:
                self.proj4cmd.append('+'+key+"="+str(value)+' ')
        self.projection = projparams['proj']
        try:
            self.rsphere = projparams['R']
        except:
            self.rsphere = 0.5*(projparams['a']+projparams['b'])
        self.llcrnrlon = llcrnrlon
        self.llcrnrlat = llcrnrlat
        if self.projection != 'cyl':
            self._proj4 = proj4.Proj(''.join(self.proj4cmd))
        llcrnrx, llcrnry = self._fwd(llcrnrlon,llcrnrlat)
        # compute x_0, y_0 so ll corner of domain is x=0,y=0.
        # note that for 'cyl' x,y == lon,lat
        self.proj4cmd.append("+x_0="+str(-llcrnrx)+' ')
        self.projparams['x_0']=-llcrnrx
        self.proj4cmd.append("+y_0="+str(-llcrnry)+' ')
        self.projparams['y_0']=-llcrnry
        # reset with x_0, y_0. 
        if self.projection != 'cyl':
            self._proj4 = proj4.Proj(''.join(self.proj4cmd))
            llcrnry = 0.
            llcrnrx = 0.
        else:
            llcrnrx = llcrnrlon
            llcrnry = llcrnrlat
        if urcrnrislatlon:
            self.urcrnrlon = urcrnrlon
            self.urcrnrlat = urcrnrlat
            urcrnrx,urcrnry = self._fwd(urcrnrlon,urcrnrlat)
        else:
            urcrnrx = urcrnrlon
            urcrnry = urcrnrlat
            urcrnrlon, urcrnrlat = self._inv(urcrnrx, urcrnry)
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

    def _fwd(self,x,y):
        if self.projection == 'cyl':
            return x,y # for cyl, this does nothing.
        else:
            outx,outy = self._proj4.fwd(x,y)
            if self.projection in ['merc','mill']: # for merc, x = rsphere*cos(lat_ts)*deltalon
                if self.projection == 'merc':
                    coslat = math.cos(math.radians(self.projparams['lat_ts']))
                else:
                    coslat = 1.
                try: # x a sequence
                    outx = [self.rsphere*coslat*math.radians(xi-self.llcrnrlon) for xi in x]
                except: # x a scalar
                    outx = self.rsphere*coslat*math.radians(x-self.llcrnrlon)
            return outx,outy

    def _inv(self,x,y):
        if self.projection == 'cyl': # for cyl, does nothing
            return x,y
        else:
            outx,outy = self._proj4.inv(x,y)
            if self.projection in ['merc','mill']: # for merc, cos(lat_ts)*deltalon = x/rsphere

                if self.projection == 'merc':
                    coslat = math.cos(math.radians(self.projparams['lat_ts']))
                else:
                    coslat = 1.
                try: # x a sequence
                    outx = [(xi/(self.rsphere*coslat))*(180./math.pi) + self.llcrnrlon for xi in x]
                except: # x a scalar
                    outx = (x/(self.rsphere*coslat))*(180./math.pi) + self.llcrnrlon

            return outx,outy

    def __call__(self,lon,lat,inverse=False):
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
        if self.projection == 'cyl': # for cyl x,y == lon,lat
            return lon,lat
        # if inputs are numarray arrays, get shape and typecode.
        try:
            shapein = lon.shape
            lontypein = lon.typecode()
            lattypein = lat.typecode()
        except:
            shapein = False
        # make sure inputs have same shape.
        if shapein and lat.shape != shapein:
            raise ValueError, 'lon, lat must be scalars or numarrays with the same shape'
        if shapein:
            x = N.ravel(lon).tolist()
            y = N.ravel(lat).tolist()
            if inverse:
                outx, outy = self._inv(x,y)
            else:
                outx, outy = self._fwd(x,y)
            outx = N.reshape(N.array(outx,lontypein),shapein)
            outy = N.reshape(N.array(outy,lattypein),shapein)
        else:
            if inverse:
                outx,outy = self._inv(lon,lat)
            else:
                outx,outy = self._fwd(lon,lat)
        return outx,outy
    
    def makegrid(self,nx,ny,returnxy=False):
        """
 return arrays of shape (ny,nx) containing lon,lat coordinates of
 an equally spaced native projection grid.
 if returnxy=True, the x,y values of the grid are returned also.
        """
        dx = (self.urcrnrx-self.llcrnrx)/(nx-1)
        dy = (self.urcrnry-self.llcrnry)/(ny-1)  
        x = self.llcrnrx+dx*N.indices((ny,nx))[1,:,:]
        y = self.llcrnry+dy*N.indices((ny,nx))[0,:,:]
        lons, lats = self(x, y, inverse=True)
        if returnxy:
            return lons, lats, x, y
        else:
            return lons, lats

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
    print min(N.ravel(lons)),max(N.ravel(lons))
    print 'max/min lats'
    print min(N.ravel(lats)),max(N.ravel(lats))
    print 'took',t2-t1,'secs'
