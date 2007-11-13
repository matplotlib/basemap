import sys
import numpy

__version__ = "0.1"

# need some python C API functions for strings.
cdef extern from "Python.h":
    object PyString_FromString(char *)

# taken from numpy.pxi in numpy 1.0rc2.
cdef extern from "numpy/arrayobject.h":
    ctypedef int npy_intp 
    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef npy_intp *dimensions
        cdef npy_intp *strides
        cdef object base
#       cdef dtype descr
        cdef int flags
    npy_intp PyArray_SIZE(ndarray arr)
    npy_intp PyArray_ISCONTIGUOUS(ndarray arr)
    void import_array()

# Initialize numpy
import_array()

# GENERAL NOTES:
#
#	- Remember to call initGEOS() before any use of this library's
#	  functions, and call finishGEOS() when done.
#
#	- Currently you have to explicitly GEOSGeom_destroy() all
#	  GEOSGeom objects to avoid memory leaks, and to free()
#	  all returned char * (unless const). This might change
#	  before first release to ensure greater API stability.

cdef extern from "geos_c.h":
# Supported geometry type IDs
    cdef enum:
        GEOS_POINT
        GEOS_LINESTRING
        GEOS_LINEARRING
        GEOS_POLYGON
        GEOS_MULTIPOINT
        GEOS_MULTILINESTRING
        GEOS_MULTIPOLYGON
        GEOS_GEOMETRYCOLLECTION
    ctypedef struct GEOSGeom:
        pass
    ctypedef struct GEOSCoordSeq:
        pass
    ctypedef void (*GEOSMessageHandler)(char *fmt, char *list)
    void initGEOS(GEOSMessageHandler notice_function, GEOSMessageHandler error_function)
    void finishGEOS()
    GEOSCoordSeq *GEOSCoordSeq_create(unsigned int size, unsigned int dims)
    void GEOSCoordSeq_destroy(GEOSCoordSeq* s)
    int GEOSCoordSeq_setX(GEOSCoordSeq* s,unsigned int idx, double val)
    int GEOSCoordSeq_setY(GEOSCoordSeq* s,unsigned int idx, double val)
    int GEOSCoordSeq_getX(GEOSCoordSeq* s, unsigned int idx, double *val)
    int GEOSCoordSeq_getY(GEOSCoordSeq* s, unsigned int idx, double *val)
    GEOSGeom  *GEOSGeom_createPoint(GEOSCoordSeq* s)
    GEOSGeom  *GEOSGeom_createLineString(GEOSCoordSeq* s)
    GEOSGeom  *GEOSGeom_createPolygon(GEOSGeom* shell, GEOSGeom** holes, unsigned int nholes)
    GEOSGeom *GEOSGeom_createLinearRing(GEOSCoordSeq* s)
    void GEOSGeom_destroy(GEOSGeom* g)
# Topology operations - return NULL on exception.
    GEOSGeom  *GEOSIntersection(GEOSGeom* g1, GEOSGeom* g2)
# Binary/Unary predicate - return 2 on exception, 1 on true, 0 on false
    char  GEOSIntersects(GEOSGeom* g1, GEOSGeom* g2)
    char  GEOSWithin(GEOSGeom* g1, GEOSGeom* g2)
    char  GEOSContains(GEOSGeom* g1, GEOSGeom* g2)
    char  GEOSisEmpty(GEOSGeom* g1)
    char  GEOSisValid(GEOSGeom* g1)
    char  GEOSisSimple(GEOSGeom* g1)
    char  GEOSisRing(GEOSGeom* g1)
#  Geometry info
    char  *GEOSGeomType(GEOSGeom* g1)
    int GEOSGeomTypeId(GEOSGeom* g1)
# Functions: Return 0 on exception, 1 otherwise 
    int  GEOSArea(GEOSGeom* g1, double *area)
    int  GEOSLength(GEOSGeom* g1, double *length)
# returns -1 on error and 1 for non-multi geoms
    int  GEOSGetNumGeometries(GEOSGeom* g1)
# Return NULL on exception, Geometry must be a Collection.
# Returned object is a pointer to internal storage:
# it must NOT be destroyed directly.
    GEOSGeom  *GEOSGetGeometryN(GEOSGeom* g, int n)
    int  GEOSGetNumInteriorRings(GEOSGeom* g1)
# Return NULL on exception, Geometry must be a Polygon.
# Returned object is a pointer to internal storage:
# it must NOT be destroyed directly.
    GEOSGeom  *GEOSGetExteriorRing(GEOSGeom* g)
# Return NULL on exception.
# Geometry must be a LineString, LinearRing or Point.
    GEOSCoordSeq  *GEOSGeom_getCoordSeq(GEOSGeom* g)
    int GEOSCoordSeq_getSize(GEOSCoordSeq *s, unsigned int *size)

cdef void notice_h(char *fmt, char*msg):
    format = PyString_FromString(fmt)
    message = PyString_FromString(msg)
    try:
        warn_msg = format % message
    except:
        warn_msg = format
    sys.stdout.write('GEOS_NOTICE: %s\n' % warn_msg)

cdef void error_h(char *fmt, char*msg):
    format = PyString_FromString(fmt)
    message = PyString_FromString(msg)
    try:
        warn_msg = format % message
    except:
        warn_msg = format
    sys.stderr.write('GEOS_ERROR: %s\n' % warn_msg)

# intialize GEOS (parameters are notice and error function callbacks).
initGEOS(notice_h, error_h)

cdef class BaseGeometry:
    cdef GEOSGeom *_geom
    cdef unsigned int _npts

    def is_valid(self):
        cdef char valid
        valid = GEOSisValid(self._geom)
        if valid:
            return True
        else:
            return False

    def geom_type(self):
        return PyString_FromString(GEOSGeomType(self._geom))

    def within(self, BaseGeometry geom):
        cdef GEOSGeom *g1, *g2
        cdef char answer
        g1 = self._geom
        g2 = geom._geom
        answer =  GEOSWithin(g1, g2)
        if answer:
            return True
        else:
            return False

    def intersects(self, BaseGeometry geom):
        cdef GEOSGeom *g1, *g2
        cdef char answer
        g1 = self._geom
        g2 = geom._geom
        answer =  GEOSIntersects(g1, g2)
        if answer:
            return True
        else:
            return False

    def intersection(self, BaseGeometry geom):
        cdef GEOSGeom *g1, *g2, *g3, *gout
        cdef char answer
        cdef int numgeoms, i, typeid
        g1 = self._geom
        g2 = geom._geom
        g3 =  GEOSIntersection(g1, g2)
        typeid = GEOSGeomTypeId(g3)
        if typeid == GEOS_POLYGON:
            p = Polygon() # create an empty Polygon instance
            p = _add_geom(p,g3) # add geometry to it.
            # above should be faster than this ..
            #b = _get_coords(g3)
            #p = Polygon(b)
            pout = [p] # return a list with a single element
        elif typeid == GEOS_LINESTRING:
            p = LineString() # create an empty LineString instance
            p = _add_geom(p,g3) # add geometry to it.
            return [p]
        elif typeid == GEOS_MULTIPOLYGON:
            numgeoms = GEOSGetNumGeometries(g3)
            pout = []
            for i from 0 <= i < numgeoms:
                gout = GEOSGetGeometryN(g3, i)
                p = Polygon() # create an empty Polygon instance
                p = _add_geom(p,gout) # add geometry to it.
                pout.append(p)
        elif typeid == GEOS_MULTILINESTRING:
            numgeoms = GEOSGetNumGeometries(g3)
            pout = []
            for i from 0 <= i < numgeoms:
                gout = GEOSGetGeometryN(g3, i)
                p = LineString() # create an LineString instance
                p = _add_geom(p,gout) # add geometry to it.
                pout.append(p)
        else:
            raise NotImplementedError("intersections of type '%s' not yet implemented" % (type))
        return pout

    def get_coords(self):
        return _get_coords(self._geom)

    def __dealloc__(self):
        """destroy GEOS geometry"""
        GEOSGeom_destroy(self._geom)

    def __reduce__(self):
        """special method that allows geos instance to be pickled"""
        return (self.__class__,(self.get_coords(),))
    
cdef class Polygon(BaseGeometry):

    def __init__(self, ndarray b=None):
        cdef unsigned int M, m, n, i
        cdef double dx, dy
        cdef double *bbuffer
        cdef GEOSCoordSeq *cs
        cdef GEOSGeom *lr

        # just return an empty class
        if b is None:
            return

        # make sure data is contiguous.
        # if not, make a local copy.
        if not PyArray_ISCONTIGUOUS(b):
            b = b.copy()

        m = b.shape[0]

        # Add closing coordinates to sequence?
        if b[-1,0] != b[0,0] or b[-1,1] != b[0,1]:
            M = m + 1
        else:
            M = m
        self._npts = M

        # Create a coordinate sequence
        cs = GEOSCoordSeq_create(M, 2)

        # add to coordinate sequence
        bbuffer = <double *>b.data
        for i from 0 <= i < m:
            dx = bbuffer[2*i]
            dy = bbuffer[2*i+1]
            # Because of a bug in the GEOS C API, 
            # always set X before Y
            GEOSCoordSeq_setX(cs, i, dx)
            GEOSCoordSeq_setY(cs, i, dy)

        # Add closing coordinates to sequence?
        if M > m:
            dx = bbuffer[0]
            dy = bbuffer[1]
            GEOSCoordSeq_setX(cs, M-1, dx)
            GEOSCoordSeq_setY(cs, M-1, dy)

        # create LinearRing
        lr = GEOSGeom_createLinearRing(cs)

        # create Polygon from LinearRing (assuming no holes)
        self._geom = GEOSGeom_createPolygon(lr,NULL,0)


    def area(self):
        cdef double area
        GEOSArea(self._geom, &area)
        return area

cdef _add_geom(BaseGeometry p, GEOSGeom *geom):
    cdef GEOSCoordSeq *cs
    cdef GEOSGeom *lr
    cdef unsigned int M
    if GEOSGeomTypeId(geom) == GEOS_POLYGON:
        lr = GEOSGetExteriorRing(geom)
        cs = GEOSGeom_getCoordSeq(lr)
    else:
        cs = GEOSGeom_getCoordSeq(geom)
    GEOSCoordSeq_getSize(cs, &M)
    p._geom = geom
    p._npts = M
    return p

cdef _get_coords(GEOSGeom *geom):
    cdef GEOSCoordSeq *cs
    cdef GEOSGeom *lr
    cdef unsigned int i, M
    cdef double dx, dy
    cdef ndarray b
    cdef double *bbuffer
    if GEOSGeomTypeId(geom) == GEOS_POLYGON:
        lr = GEOSGetExteriorRing(geom)
        cs = GEOSGeom_getCoordSeq(lr)
    else:
        cs = GEOSGeom_getCoordSeq(geom)
    GEOSCoordSeq_getSize(cs, &M)
    b = numpy.empty((M,2), numpy.float64)
    bbuffer = <double *>b.data
    for i from 0 <= i < M:
        GEOSCoordSeq_getX(cs, i, &dx)
        GEOSCoordSeq_getY(cs, i, &dy)
        bbuffer[2*i] = dx
        bbuffer[2*i+1] = dy
    return b

cdef class LineString(BaseGeometry):
    def __init__(self, ndarray b=None):
        cdef double dx, dy
        cdef GEOSCoordSeq *cs
        cdef int i, M
        cdef double *bbuffer

        # just return an empty class
        if b is None:
            return

        # make sure data is contiguous.
        # if not, make a local copy.
        if not PyArray_ISCONTIGUOUS(b):
            b = b.copy()

        M = b.shape[0]
        self._npts = M

        # Create a coordinate sequence
        cs = GEOSCoordSeq_create(M, 2)

        # add to coordinate sequence
        bbuffer = <double *>b.data
        for i from 0 <= i < M:
            dx = bbuffer[2*i]
            dy = bbuffer[2*i+1]
            # Because of a bug in the GEOS C API, 
            # always set X before Y
            GEOSCoordSeq_setX(cs, i, dx)
            GEOSCoordSeq_setY(cs, i, dy)

        # create LineString
        self._geom = GEOSGeom_createLineString(cs)

cdef class Point(BaseGeometry):
    cdef public x,y
    def __init__(self, b=None):
        cdef double dx, dy
        cdef GEOSCoordSeq *cs
        # just return an empty class
        if b is None:
            return
        # Create a coordinate sequence
        cs = GEOSCoordSeq_create(1, 2)
        dx = b[0]; dy = b[1]
        GEOSCoordSeq_setX(cs, 0, dx)
        GEOSCoordSeq_setY(cs, 0, dy)
        self._geom = GEOSGeom_createPoint(cs)
        self._npts = 1
