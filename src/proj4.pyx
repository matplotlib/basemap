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
765975.641091.4805993.13406
(-120.10799999995851, 34.361166659972767)

Input coordinates can be given as python arrays, sequences, scalars
or Numeric/numarray arrays. Optimized for objects that support
the Python buffer protocol (regular python, Numeric and numarray arrays).

Download http://www.cdc.noaa.gov/people/jeffrey.s.whitaker/python/pyproj-1.5.tar.gz

See pyproj.Proj.__doc__ for more documentation.

Contact:  Jeffrey Whitaker <jeffrey.s.whitaker@noaa.gov
""" 

# Make changes to this file, not the c-wrappers that Pyrex generates.

import math, array

cdef double _rad2dg, _dg2rad
cdef int _doublesize
_dg2rad = math.radians(1.)
_rad2dg = math.degrees(1.)
_doublesize = sizeof(double)
__version__ = 1.5

cdef extern from "proj_api.h":
    ctypedef double *projPJ
    ctypedef struct projUV:
        double u
        double v
    projPJ pj_init_plus(char *)
    projUV pj_fwd(projUV, projPJ)
    projUV pj_inv(projUV, projPJ)
    void pj_free(projPJ)

cdef extern from "Python.h":
  int PyObject_AsWriteBuffer(object, void **rbuf, int *len)
  int PyObject_CheckReadBuffer(object)

cdef class Proj:
    """
 performs cartographic transformations (converts from longitude,latitude
 to native map projection x,y coordinates and vice versa) using proj 
 (http://proj.maptools.org/)

 A Proj class instance is initialized with a dictionary containing 
 proj map projection control parameter key/value pairs.
 See http://www.remotesensing.org/geotiff/proj_list and the
 proj man page for details.

 Calling a Proj class instance with the arguments lon, lat will
 convert lon/lat (in degrees) to x/y native map projection 
 coordinates (in meters).  If optional keyword 'inverse' is
 True (default is False), the inverse transformation from x/y
 to lon/lat is performed. If optional keyword 'radians' is True
 (default is False) lon/lat are interpreted as radians instead
 of degrees. Works with numarray or Numeric arrays, python arrays,
 sequences or scalars (fastest for arrays containing doubles).
    """

    cdef double *projpj
    cdef object projparams
    cdef char *pjinitstring

    def __new__(self, projparams):
        """
 initialize a Proj class instance.

 Input 'projparams' is a dictionary containing proj map
 projection control parameter key/value pairs.
 See the proj documentation (http://proj.maptools.org) for details.
        """
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
        self.projparams = projparams
        pjinitstring = ''.join(pjargs)
        self.projpj = pj_init_plus(pjinitstring)

    def __dealloc__(self):
        """destroy projection definition"""
        pj_free(self.projpj)

    def __reduce__(self):
        """special method that allows projlib.Proj instance to be pickled"""
        return (self.__class__,(self.projparams,))

    def _fwd(self, object lons, object lats, radians=False):
        """
 forward transformation - lons,lats to x,y.
 if radians=True, lons/lats are radians instead of degrees.
        """
        cdef projUV projxyout, projlonlatin
        cdef int ndim, i, buflenx, bufleny
        cdef double u, v
        cdef double *lonsdata, *latsdata
        cdef void *londata, *latdata
        try:
            # if buffer api is supported, get pointer to data buffers.
            if PyObject_AsWriteBuffer(lons, &londata, &buflenx) <> 0:
                raise RuntimeError
            if PyObject_AsWriteBuffer(lats, &latdata, &bufleny) <> 0:
                raise RuntimeError
            hasbufapi= True
        except:
            hasbufapi = False
        if hasbufapi:
        # process data in buffer (for Numeric, numarray and python arrays).
            if buflenx != bufleny:
                raise RuntimeError("Buffer lengths not the same")
            ndim = buflenx/_doublesize
            lonsdata = <double *>londata
            latsdata = <double *>latdata
            if radians:
                for i from 0 <= i < ndim:
                    projlonlatin.u = lonsdata[i]
                    projlonlatin.v = latsdata[i]
                    projxyout = pj_fwd(projlonlatin,self.projpj)
                    lonsdata[i] = projxyout.u
                    latsdata[i] = projxyout.v
            else:
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
                    raise RuntimeError("Sequences must have the same number of elements")
                x = []; y = []
                if radians:
                    for i from 0 <= i < ndim:
                        projlonlatin.u = lons[i]
                        projlonlatin.v = lats[i]
                        projxyout = pj_fwd(projlonlatin,self.projpj)
                        x.append(projxyout.u)
                        y.append(projxyout.v)
                else:
                    for i from 0 <= i < ndim:
                        projlonlatin.u = _dg2rad*lons[i]
                        projlonlatin.v = _dg2rad*lats[i]
                        projxyout = pj_fwd(projlonlatin,self.projpj)
                        x.append(projxyout.u)
                        y.append(projxyout.v)
            except: # inputs are scalars.
                if radians:
                    projlonlatin.u = lons
                    projlonlatin.v = lats
                else:
                    projlonlatin.u = lons*_dg2rad
                    projlonlatin.v = lats*_dg2rad
                projxyout = pj_fwd(projlonlatin,self.projpj)
                x = projxyout.u
                y = projxyout.v
            return x,y

    def _inv(self, object x, object y, radians=False):
        """
 inverse transformation - x,y to lons,lats.
 if radians=True, lons/lats are radians instead of degrees.
        """
        cdef projUV projxyin, projlonlatout
        cdef int ndim, i, buflenx, bufleny
        cdef double u, v
        cdef void *xdata, *ydata
        cdef double *xdatab, *ydatab
        try:
            # if buffer api is supported, get pointer to data buffers.
            if PyObject_AsWriteBuffer(x, &xdata, &buflenx) <> 0:
                raise RuntimeError
            if PyObject_AsWriteBuffer(y, &ydata, &bufleny) <> 0:
                raise RuntimeError
            hasbufapi= True
        except:
            hasbufapi = False
        if hasbufapi:
        # process data in buffer (for Numeric, numarray and python arrays).
            if buflenx != bufleny:
                raise RuntimeError("Buffer lengths not the same")
            ndim = buflenx/_doublesize

            xdatab = <double *>xdata
            ydatab = <double *>ydata

            if radians:
                for i from 0 <= i < ndim:
                    projxyin.u = xdatab[i]
                    projxyin.v = ydatab[i]
                    projlonlatout = pj_inv(projxyin,self.projpj)
                    xdatab[i] = projlonlatout.u
                    ydatab[i] = projlonlatout.v
            else:
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
                    raise RuntimeError("Sequences must have the same number of elements")
                lons = []; lats = []
                if radians:
                    for i from 0 <= i < ndim:
                        projxyin.u = x[i]
                        projxyin.v = y[i]
                        projlonlatout = pj_inv(projxyin, self.projpj)
                        lons.append(projlonlatout.u)
                        lats.append(projlonlatout.v)
                else:
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
                if radians:
                    lons = projlonlatout.u
                    lats = projlonlatout.v
                else:
                    lons = projlonlatout.u*_rad2dg
                    lats = projlonlatout.v*_rad2dg
            return lons, lats


    def __call__(self,lon,lat,inverse=False,radians=False):
        """
 Calling a Proj class instance with the arguments lon, lat will
 convert lon/lat (in degrees) to x/y native map projection 
 coordinates (in meters).  If optional keyword 'inverse' is
 True (default is False), the inverse transformation from x/y
 to lon/lat is performed.  If optional keyword 'radians' is
 True (default is False) the units of lon/lat are radians instead
 of degrees.

 Inputs should be doubles (they will be cast to doubles
 if they are not, causing a slight performance hit).

 Works with Numeric or numarray arrays, python sequences or scalars
 (fastest for arrays containing doubles).
        """
        try:
            # typecast Numeric/numarray arrays to double, if necessary.
            if lon.typecode() != 'd':
                lon = lon.astype('d')
            if lat.typecode() != 'd':
                lat = lat.astype('d')
        except:
            # typecast regular python arrays to double, if necessary.
            try:
                if lon.typecode != 'd':
                    lon = array.array('d',lon)
                if lat.typecode != 'd':
                    lat = array.array('d',lat)
            except:
                pass
        # If the buffer API is supported, make copies of inputs.
        # This is necessary since the data buffer of the inputs
        # will be modified in place.  Raise an exception if
        # copies cannot be made.
        # Buffer API will be used if inputs are arrays (regular python, 
        # Numeric or numarray).
        if PyObject_CheckReadBuffer(lon) and PyObject_CheckReadBuffer(lon):
            try:
                # try to make copy using __copy__ method.
                inx = lon.__copy__(); iny = lat.__copy__()
            except:
                msg = """could not create copy of inputs.  
This could be because your are using regular python arrays with Python 2.3
(python arrays are not copy-able before python 2.4).  Try using lists
or Numeric/numarray arrays instead (the latter will be faster)."""
                raise RuntimeError, msg
            # call proj4 functions.
            if inverse:
                outx, outy = self._inv(inx, iny, radians=radians)
            else:
                outx, outy = self._fwd(inx, iny, radians=radians)
        # copy not needed if buffer API not supported.
        else:
            if inverse:
                outx, outy = self._inv(lon, lat, radians=radians)
            else:
                outx, outy = self._fwd(lon, lat, radians=radians)
        # all done.
        return outx,outy
