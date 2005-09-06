"""Pyrex code to provide python interfaces to PROJ.4 functions.
Make changes to this file, not the c-wrappers that Pyrex generates."""

import math

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

# get 32 bit integer type
cdef extern from "inttypes.h":
  ctypedef long int32_t

# Functions from numarray API
cdef extern from "numarray/libnumarray.h":
  ctypedef int32_t maybelong
  # numarray types.
  ctypedef enum NumarrayType:
    tAny
    tBool       
    tInt8
    tUInt8
    tInt16
    tUInt16
    tInt32
    tUInt32
    tInt64
    tUInt64
    tFloat32
    tFloat64
    tComplex32
    tComplex64
    tObject
    tDefault
    tLong
  cdef struct PyArray_Descr:
    int type_num # PyArray_TYPES
    int elsize   # bytes for 1 element
    char type    # One of "cb1silfdFD "  Object array not supported
    # function pointers omitted
  ctypedef class numarray._numarray._numarray [object PyArrayObject]:
     cdef char *data
     cdef int nd
     cdef maybelong *dimensions
     cdef maybelong *strides
     cdef object base
     cdef PyArray_Descr *descr
     cdef int flags
     cdef maybelong *_dimensions
     cdef maybelong *_strides
  # creates a new numeric object from a data pointer.
  _numarray NA_vNewArray(void *, NumarrayType, int, maybelong *)
  # The numarray initialization funtion
  double import_libnumarray()
    
# this function must be called to initialize the numarray C API.
import_libnumarray()

cdef class Proj:

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

    def fwd_array(self, _numarray lons, _numarray lats):
        """
 forward transformation - lons,lats to x,y.
 x, y, lons, lats are numarrays.
        """
        cdef projUV projxyout, projlonlatin
        cdef long ndim, i
        cdef double u, v
        cdef double *lonsdata, *latsdata
        ndim = 1
        for i from 0 <= i < lons.nd:
            ndim = ndim*lons.dimensions[i]
        lonsdata = <double *>lons.data
        latsdata = <double *>lats.data
        for i from 0 <= i < ndim:
            projlonlatin.u = _dg2rad*lonsdata[i]
            projlonlatin.v = _dg2rad*latsdata[i]
            projxyout = pj_fwd(projlonlatin,self.projpj)
            lonsdata[i] = projxyout.u
            latsdata[i] = projxyout.v
        x = NA_vNewArray(<void *>lonsdata, tFloat64, lons.nd, lons.dimensions)
        y = NA_vNewArray(<void *>latsdata, tFloat64, lats.nd, lats.dimensions)
        return x,y

    def fwd(self, lons, lats):
        """
 forward transformation - lons,lats to x,y.
 x, y, lons, lats can be scalars, sequences or numarrays.
        """
        cdef projUV projxyout, projlonlatin
        cdef double u, v
        try:
            shape = lons.shape
            isnumarray = True
        except:
            isnumarray = False
        if isnumarray:
            if lons.typecode() != 'd' or lats.typecode != 'd':
                lons = lons.astype('d'); lats = lats.astype('d')
            x,y = self.fwd_array(lons,lats)
        else:
            try: # inputs are lists
                ndim = len(lons)
                x = []; y = []
                for i from 0 <= i < ndim:
                    projlonlatin.u = _dg2rad*lons[i]
                    projlonlatin.v = _dg2rad*lats[i]
                    projxyout = pj_fwd(projlonlatin,self.projpj)
                    x.append(projxyout.u)
                    y.append(projxyout.v)
            except: # inputs are scalars
                projlonlatin.u = lons*_dg2rad
                projlonlatin.v = lats*_dg2rad
                projxyout = pj_fwd(projlonlatin,self.projpj)
                x = projxyout.u
                y = projxyout.v
        return x,y

    def inv_array(self, _numarray x, _numarray y):
        """
 inverse transformation - x,y to lons,lats
 x, y, lons, lats are numarrays.
        """
        cdef projUV projxyin, projlonlatout
        cdef long ndim, i
        cdef double u, v
        cdef double *xdata, *ydata
        ndim = 1
        for i from 0 <= i < x.nd:
            ndim = ndim*x.dimensions[i]
        xdata = <double *>x.data
        ydata = <double *>y.data
        for i from 0 <= i < ndim:
            projxyin.u = xdata[i]
            projxyin.v = ydata[i]
            projlonlatout = pj_inv(projxyin,self.projpj)
            xdata[i] = _rad2dg*projlonlatout.u
            ydata[i] = _rad2dg*projlonlatout.v
        lons = NA_vNewArray(<void *>xdata, tFloat64, x.nd, x.dimensions)
        lats = NA_vNewArray(<void *>ydata, tFloat64, y.nd, y.dimensions)
        return lons, lats

    def inv(self, x, y):
        """
 inverse transformation - x,y to lons,lats
 x, y, lons, lats can be scalars, sequences or numarrays.
        """
        cdef projUV projxyin, projlonlatout
        cdef double u, v
        try:
            shape = x.shape
            isnumarray = True
        except:
            isnumarray = False
        if isnumarray:
            if x.typecode() != 'd' or y.typecode != 'd':
                x = x.astype('d'); y = y.astype('d')
            lons, lats = self.inv_array(x,y)
        else:
            try: # inputs are lists
                ndim = len(x)
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
        return lons,lats
