"""
Pyrex wrapper to provide python interfaces to 
PROJ.4 (http://proj.maptools.org) functions.

Performs cartographic transformations (converts from longitude,latitude
to native map projection x,y coordinates and vice versa).

Example usage:

>>> from pyproj import Proj
>>> params = {}
>>> params['proj'] = 'utm'
>>> params['zone'] = 10
>>> p = Proj(params)
>>> x,y = p(-120.108, 34.36116666)
>>> print x,y
>>> print p(x,y,inverse=True)
765975.641091 3805993.13406
(-120.10799999995851, 34.361166659972767)

Input coordinates can be given as python arrays, sequences, scalars
or Numeric/numarray arrays. Optimized for Numeric/numarray arrays.

Download http://www.cdc.noaa.gov/people/jeffrey.s.whitaker/python/pyproj-1.2.tar.gz

See pyproj.Proj.__doc__ for more documentation.

Contact:  Jeffrey Whitaker <jeffrey.s.whitaker@noaa.gov
"""

# Make changes to this file, not the c-wrappers that Pyrex generates.

import math, copy, array

cdef double _rad2dg, _dg2rad
_dg2rad = math.radians(1.)
_rad2dg = math.degrees(1.)

cdef extern from "proj_api.h":
    ctypedef double *projPJ
    ctypedef struct projUV:
        double u
        double v
    projPJ pj_init_plus(char *)
    projUV pj_fwd(projUV, projPJ)
    projUV pj_inv(projUV, projPJ)

cdef extern from "Python.h":
  int PyObject_AsReadBuffer(object, void **rbuf, int *len)

cdef class Proj:
    """
 performs cartographic transformations (converts from longitude,latitude
 to native map projection x,y coordinates and vice versa) using proj 
 (http://proj.maptools.org/)

 A Proj class instance is initialized with a dictionary containing 
 proj map projection control parameter key/value pairs.
 See the proj documentation (http://www.remotesensing.org/proj/)
 for details.

 Calling a Proj class instance with the arguments lon, lat will
 convert lon/lat (in degrees) to x/y native map projection 
 coordinates (in meters).  If optional keyword 'inverse' is
 True (default is False), the inverse transformation from x/y
 to lon/lat is performed. Works with numarray or Numeric arrays,
 python arrays, sequences or scalars (fastest for arrays containing
 doubles).
    """

    cdef char *pjinitstring
    cdef double *projpj
    cdef object projparams

    def __new__(self, projparams):
        """
 initialize a Proj class instance.

 Input 'projparams' is a dictionary containing proj map
 projection control parameter key/value pairs.
 See the proj documentation (http://proj.maptools.org) for details.
        """
        cdef double *projpj
        # set units to meters.
        if not projparams.has_key('units'):
            projparams['units']='m'
        elif projparams['units'] != 'm':
            print 'resetting units to meters ...'
            projparams['units']='m'
        # make sure proj parameter specified.
        # (no other checking done in proj parameters)
        if 'proj' not in projparams.keys():
            raise KeyError, "need to specify proj parameter"
        pjargs = []
        for key,value in projparams.iteritems():
            pjargs.append('+'+key+"="+str(value)+' ')
        pjinitstring = ''.join(pjargs)
        self.projparams = projparams
        projpj = pj_init_plus(pjinitstring)
        self.projpj = projpj

    def __reduce__(self):
        """special method that allows projlib.Proj instance to be pickled"""
        return (self.__class__,(self.projparams,))

    def _fwd(self, lons, lats):
        """
 forward transformation - lons,lats to x,y.
        """
        cdef projUV projxyout, projlonlatin
        cdef int ndim, i, buflenx, bufleny
        cdef double u, v
        cdef double *lonsdata, *latsdata
        cdef void *londata, *latdata
        
        try:
            # if buffer api is supported, get pointer to data buffers.
            if PyObject_AsReadBuffer(lons, &londata, &buflenx) <> 0:
                raise RuntimeError
            if PyObject_AsReadBuffer(lats, &latdata, &bufleny) <> 0:
                raise RuntimeError
            hasbufapi= True
        except:
            hasbufapi = False
        
        if hasbufapi:
        # process data in buffer (for Numeric, numarray and python arrays).
            if buflenx != bufleny:
                raise RuntimeError("Buffer lengths not the same")
            ndim = buflenx/8

            lonsdata = <double *>londata
            latsdata = <double *>latdata

            for i from 0 <= i < ndim:
                projlonlatin.u = _dg2rad*lonsdata[i]
                projlonlatin.v = _dg2rad*latsdata[i]
                projxyout = pj_fwd(projlonlatin,self.projpj)
                lonsdata[i] = projxyout.u
                latsdata[i] = projxyout.v
            return lons, lats
        else:
                try: # inputs are sequences.
                    ndim = len(lons)
                    if len(lats) != ndim:
                        raise RuntimeError("sequences must have the same number of elements")
                    x = []; y = []
                    for i from 0 <= i < ndim:
                        projlonlatin.u = _dg2rad*lons[i]
                        projlonlatin.v = _dg2rad*lats[i]
                        projxyout = pj_fwd(projlonlatin,self.projpj)
                        x.append(projxyout.u)
                        y.append(projxyout.v)
                except: # inputs are scalars.
                    projlonlatin.u = lons*_dg2rad
                    projlonlatin.v = lats*_dg2rad
                    projxyout = pj_fwd(projlonlatin,self.projpj)
                    x = projxyout.u
                    y = projxyout.v
                return x,y

    def _inv(self, object x, object y):
        """
 inverse transformation - x,y to lons,lats
        """
        cdef projUV projxyin, projlonlatout
        cdef int ndim, i, buflenx, bufleny
        cdef double u, v
        cdef void *xdata, *ydata
        cdef double *xdatab, *ydatab

        try:
            # if buffer api is supported, get pointer to data buffers.
            if PyObject_AsReadBuffer(x, &xdata, &buflenx) <> 0:
                raise RuntimeError
            if PyObject_AsReadBuffer(y, &ydata, &bufleny) <> 0:
                raise RuntimeError
            hasbufapi= True
        except:
            hasbufapi = False
        
        if hasbufapi:
        # process data in buffer (for Numeric, numarray and python arrays).

            if buflenx != bufleny:
                raise RuntimeError("Buffer lengths not the same")
            ndim = buflenx/8

            xdatab = <double *>xdata
            ydatab = <double *>ydata

            for i from 0 <= i < ndim:
                projxyin.u = xdatab[i]
                projxyin.v = ydatab[i]
                projlonlatout = pj_inv(projxyin,self.projpj)
                xdatab[i] = _rad2dg*projlonlatout.u
                ydatab[i] = _rad2dg*projlonlatout.v
            return x,y
                
        else:
                try: # inputs are sequences.
                    ndim = len(x)
                    if len(y) != ndim:
                        raise RuntimeError("sequences must have the same number of elements")
                    lons = []; lats = []
                    for i from 0 <= i < ndim:
                        projxyin.u = x[i]
                        projxyin.v = y[i]
                        projlonlatout = pj_inv(projxyin, self.projpj)
                        lons.append(projlonlatout.u*_rad2dg)
                        lats.append(projlonlatout.v*_rad2dg)
                except: # inputs are scalars.
                    projxyin.u = x
                    projxyin.v = y
                    projlonlatout = pj_inv(projxyin, self.projpj)
                    lons = projlonlatout.u*_rad2dg
                    lats = projlonlatout.v*_rad2dg
                return lons, lats


    def __call__(self,lon,lat,inverse=False):
        """
 Calling a Proj class instance with the arguments lon, lat will
 convert lon/lat (in degrees) to x/y native map projection 
 coordinates (in meters).  If optional keyword 'inverse' is
 True (default is False), the inverse transformation from x/y
 to lon/lat is performed.

 Inputs should be doubles (they will be cast to doubles
 if they are not).

 Works with Numeric or numarray arrays, python sequences or scalars
 (fastest for arrays containing doubles).
        """
        try:
            # typecast Numeric/numarray arrays to double.
            if lon.typecode() != 'd':
                lon = lon.astype('d')
            if lat.typecode() != 'd':
                lat = lat.astype('d')
        except:
            # typecast regular python arrays to double
            try:
                if lon.typecode != 'd':
                    lon = array.array('d',lon)
                if lat.typecode != 'd':
                    lat = array.array('d',lat)
            except:
                pass
        # make copies of inputs.
        # (If buffer api is supported, the data buffer of
        # the inputs will be modified in place.)
        inx = copy.copy(lon); iny = copy.copy(lat)
        if inverse:
            outx, outy = self._inv(inx, iny)
        else:
            outx, outy = self._fwd(inx, iny)
        return outx,outy
