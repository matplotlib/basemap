import sys
import numpy

__version__ = "0.3"

# need some python C API functions for strings.
cdef extern from "Python.h":
    object PyBytes_FromString(char *)

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
        GEOS_VERSION_MAJOR
    ctypedef struct GEOSGeom:
        pass
    ctypedef struct GEOSCoordSeq:
        pass
    ctypedef void (*GEOSMessageHandler)(char *fmt, char *list)
    char *GEOSversion()
    void initGEOS(GEOSMessageHandler notice_function, GEOSMessageHandler error_function)
    void finishGEOS()
    GEOSCoordSeq *GEOSCoordSeq_create(unsigned int size, unsigned int dims)
    void GEOSCoordSeq_destroy(GEOSCoordSeq* s)
    int GEOSCoordSeq_setX(GEOSCoordSeq* s,unsigned int idx, double val)
    int GEOSCoordSeq_setY(GEOSCoordSeq* s,unsigned int idx, double val)
    int GEOSCoordSeq_getX(GEOSCoordSeq* s, unsigned int idx, double *val)
    int GEOSCoordSeq_getY(GEOSCoordSeq* s, unsigned int idx, double *val)
    GEOSGeom *GEOSUnion(GEOSGeom* g1, GEOSGeom* g2)
    GEOSGeom *GEOSUnaryUnion(GEOSGeom* g1)
    GEOSGeom *GEOSEnvelope(GEOSGeom* g1)
    GEOSGeom *GEOSConvexHull(GEOSGeom* g1)
    GEOSGeom *GEOSGeom_createPoint(GEOSCoordSeq* s)
    GEOSGeom *GEOSGeom_createLineString(GEOSCoordSeq* s)
    GEOSGeom *GEOSGeom_createPolygon(GEOSGeom* shell, GEOSGeom** holes, unsigned int nholes)
    GEOSGeom  *GEOSGeom_createLinearRing(GEOSCoordSeq* s)
    void GEOSGeom_destroy(GEOSGeom* g)
# Topology operations - return NULL on exception.
    GEOSGeom *GEOSIntersection(GEOSGeom* g1, GEOSGeom* g2)
    GEOSGeom *GEOSSimplify(GEOSGeom* g1, double tolerance)
    GEOSGeom *GEOSBuffer(GEOSGeom* g1, double width, int quadsegs)
    GEOSGeom *GEOSTopologyPreserveSimplify(GEOSGeom* g1, double tolerance)
# Binary/Unary predicate - return 2 on exception, 1 on true, 0 on false
    char GEOSIntersects(GEOSGeom* g1, GEOSGeom* g2)
    char GEOSWithin(GEOSGeom* g1, GEOSGeom* g2)
    char GEOSContains(GEOSGeom* g1, GEOSGeom* g2)
    char GEOSisEmpty(GEOSGeom* g1)
    char GEOSisValid(GEOSGeom* g1)
    char GEOSisSimple(GEOSGeom* g1)
    char GEOSisRing(GEOSGeom* g1)
#  Geometry info
    char *GEOSGeomType(GEOSGeom* g1)
    int GEOSGeomTypeId(GEOSGeom* g1)
# Functions: Return 0 on exception, 1 otherwise 
    int GEOSArea(GEOSGeom* g1, double *area)
    int GEOSLength(GEOSGeom* g1, double *length)
# returns -1 on error and 1 for non-multi geoms
    int GEOSGetNumGeometries(GEOSGeom* g1)
# Return NULL on exception, Geometry must be a Collection.
# Returned object is a pointer to internal storage:
# it must NOT be destroyed directly.
    GEOSGeom  *GEOSGetGeometryN(GEOSGeom* g, int n)
    int GEOSGetNumInteriorRings(GEOSGeom* g1)
# Return NULL on exception, Geometry must be a Polygon.
# Returned object is a pointer to internal storage:
# it must NOT be destroyed directly.
    GEOSGeom  *GEOSGetExteriorRing(GEOSGeom* g)
# Return NULL on exception.
# Geometry must be a LineString, LinearRing or Point.
    GEOSCoordSeq  *GEOSGeom_getCoordSeq(GEOSGeom* g)
    int GEOSCoordSeq_getSize(GEOSCoordSeq *s, unsigned int *size)

cdef void notice_h(char *fmt, char*msg):
    pass
    #format = PyBytes_FromString(fmt)
    #message = PyBytes_FromString(msg)
    #try:
    #    warn_msg = format % message
    #except:
    #    warn_msg = format
    #sys.stdout.write('GEOS_NOTICE: %s\n' % warn_msg)

cdef void error_h(char *fmt, char*msg):
    format = PyBytes_FromString(fmt)
    message = PyBytes_FromString(msg)
    try:
        warn_msg = format % message
    except:
        warn_msg = format
    sys.stderr.write('GEOS_ERROR: %s\n' % warn_msg)

# check library version
cdef geos_version():
    return PyBytes_FromString(GEOSversion())
__geos_version__ = geos_version() # module variable.
__geos_major_version__ = GEOS_VERSION_MAJOR
#if __geos_version__ != "2.2.3-CAPI-1.1.1":
#     raise ValueError('version 2.2.3 of the geos library is required')

# intialize GEOS (parameters are notice and error function callbacks).
initGEOS(notice_h, error_h)

cdef class BaseGeometry:
    cdef GEOSGeom *_geom
    cdef unsigned int _npts
    cdef public object boundary

    def is_valid(self):
        cdef char valid
        valid = GEOSisValid(self._geom)
        if valid:
            return True
        else:
            return False

    def geom_type(self):
        return PyBytes_FromString(GEOSGeomType(self._geom))

    def within(self, BaseGeometry geom):
        cdef GEOSGeom *g1
        cdef GEOSGeom *g2
        cdef char answer
        g1 = self._geom
        g2 = geom._geom
        answer =  GEOSWithin(g1, g2)
        if answer:
            return True
        else:
            return False

    def union(self, BaseGeometry geom):
        cdef GEOSGeom *g1
        cdef GEOSGeom *g2
        cdef GEOSGeom *g3
        cdef GEOSGeom *gout
        cdef int numgeoms, i, typeid
        g1 = self._geom
        g2 = geom._geom
        g3 = GEOSUnion(g1, g2)
        typeid = GEOSGeomTypeId(g3)
        if typeid == GEOS_POLYGON:
            b = _get_coords(g3)
            p = Polygon(b)
        elif typeid == GEOS_LINESTRING:
            b = _get_coords(g3)
            p = LineString(b)
        # for multi-geom structures, just return first one.
        elif typeid == GEOS_MULTIPOLYGON:
            numgeoms = GEOSGetNumGeometries(g3)
            gout = GEOSGetGeometryN(g3, 0)
            b = _get_coords(gout)
            p = Polygon(b)
        elif typeid == GEOS_MULTILINESTRING:
            numgeoms = GEOSGetNumGeometries(g3)
            gout = GEOSGetGeometryN(g3, 0)
            b = _get_coords(gout)
            p = LineString(b)
        else:
            type = PyBytes_FromString(GEOSGeomType(g3))
            raise NotImplementedError("unions of type '%s' not yet implemented" % (type))
        GEOSGeom_destroy(g3)
        return p

    def simplify(self, tol):
        cdef GEOSGeom *g1
        cdef GEOSGeom *g3
        cdef GEOSGeom *gout
        cdef double tolerance
        cdef int numgeoms, i, typeid
        g1 = self._geom
        tolerance = tol
        g3 = GEOSSimplify(g1,tolerance)
        typeid = GEOSGeomTypeId(g3)
        if typeid == GEOS_POLYGON:
            b = _get_coords(g3)
            p = Polygon(b)
        elif typeid == GEOS_LINESTRING:
            b = _get_coords(g3)
            p = LineString(b)
        # for multi-geom structures, just return first one.
        elif typeid == GEOS_MULTIPOLYGON:
            numgeoms = GEOSGetNumGeometries(g3)
            gout = GEOSGetGeometryN(g3, 0)
            b = _get_coords(gout)
            p = Polygon(b)
        elif typeid == GEOS_MULTILINESTRING:
            numgeoms = GEOSGetNumGeometries(g3)
            gout = GEOSGetGeometryN(g3, 0)
            b = _get_coords(gout)
            p = LineString(b)
        else:
            type = PyBytes_FromString(GEOSGeomType(g3))
            raise NotImplementedError("intersections of type '%s' not yet implemented" % (type))
        GEOSGeom_destroy(g3)
        return p

    def fix(self):
        cdef GEOSGeom *g1
        cdef GEOSGeom *g3
        cdef GEOSGeom *gout
        cdef int numgeoms, i, typeid
        g1 = self._geom
        g3 = GEOSBuffer(g1, 0., 0)
        typeid = GEOSGeomTypeId(g3)
        if typeid == GEOS_POLYGON:
            b = _get_coords(g3)
            p = Polygon(b)
        elif typeid == GEOS_LINESTRING:
            b = _get_coords(g3)
            p = LineString(b)
        # for multi-geom structures, just return first one.
        elif typeid == GEOS_MULTIPOLYGON:
            numgeoms = GEOSGetNumGeometries(g3)
            gout = GEOSGetGeometryN(g3, 0)
            b = _get_coords(gout)
            p = Polygon(b)
        elif typeid == GEOS_MULTILINESTRING:
            numgeoms = GEOSGetNumGeometries(g3)
            gout = GEOSGetGeometryN(g3, 0)
            b = _get_coords(gout)
            p = LineString(b)
        else:
            type = PyBytes_FromString(GEOSGeomType(g3))
            raise NotImplementedError("intersections of type '%s' not yet implemented" % (type))
        GEOSGeom_destroy(g3)
        return p

    def intersects(self, BaseGeometry geom):
        cdef GEOSGeom *g1
        cdef GEOSGeom *g2
        cdef char answer
        g1 = self._geom
        g2 = geom._geom
        answer =  GEOSIntersects(g1, g2)
        if answer:
            return True
        else:
            return False

    def intersection(self, BaseGeometry geom):
        cdef GEOSGeom *g1
        cdef GEOSGeom *g2
        cdef GEOSGeom *g3
        cdef GEOSGeom *gout
        cdef char answer
        cdef int numgeoms, i, typeid
        g1 = self._geom
        g2 = geom._geom
        g3 =  GEOSIntersection(g1, g2)
        typeid = GEOSGeomTypeId(g3)
        if typeid == GEOS_POLYGON:
            b = _get_coords(g3)
            p = Polygon(b)
            pout = [p]
        elif typeid == GEOS_LINESTRING:
            b = _get_coords(g3)
            p = LineString(b)
            pout = [p]
        elif typeid == GEOS_MULTIPOLYGON:
            numgeoms = GEOSGetNumGeometries(g3)
            pout = []
            for i from 0 <= i < numgeoms:
                gout = GEOSGetGeometryN(g3, i)
                b = _get_coords(gout)
                p = Polygon(b)
                pout.append(p)
        elif typeid == GEOS_MULTILINESTRING:
            numgeoms = GEOSGetNumGeometries(g3)
            pout = []
            for i from 0 <= i < numgeoms:
                gout = GEOSGetGeometryN(g3, i)
                b = _get_coords(gout)
                p = LineString(b)
                pout.append(p)
        else:
            #type = PyBytes_FromString(GEOSGeomType(g3))
            #raise NotImplementedError("intersections of type '%s' not yet implemented" % (type))
            GEOSGeom_destroy(g3)
            return []
        GEOSGeom_destroy(g3)
        return pout

    def get_coords(self):
        return _get_coords(self._geom)

    def __dealloc__(self):
        """destroy GEOS geometry"""
        GEOSGeom_destroy(self._geom)

    def __reduce__(self):
        """special method that allows geos instance to be pickled"""
        return (self.__class__,(self.boundary,))
    
cdef class Polygon(BaseGeometry):

    def __init__(self, ndarray b):
        cdef unsigned int M, m, i
        cdef double dx, dy
        cdef double *bbuffer
        cdef GEOSCoordSeq *cs
        cdef GEOSGeom *lr



        # make sure data is contiguous.
        # if not, make a local copy.
        if not PyArray_ISCONTIGUOUS(b):
            b = b.copy()

        m = b.shape[0]
       
        # Add closing coordinates to sequence?
        if m > 0 and (b[-1,0] != b[0,0] or b[-1,1] != b[0,1]):
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
        self.boundary = b


    def area(self):
        cdef double area
        GEOSArea(self._geom, &area)
        return area

cdef class LineString(BaseGeometry):
    def __init__(self, ndarray b):
        cdef double dx, dy
        cdef GEOSCoordSeq *cs
        cdef int i, M
        cdef double *bbuffer

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
        self.boundary = b

cdef class Point(BaseGeometry):
    cdef public x,y
    def __init__(self, b):
        cdef double dx, dy
        cdef GEOSCoordSeq *cs
        # Create a coordinate sequence
        cs = GEOSCoordSeq_create(1, 2)
        dx = b[0]; dy = b[1]
        GEOSCoordSeq_setX(cs, 0, dx)
        GEOSCoordSeq_setY(cs, 0, dy)
        self._geom = GEOSGeom_createPoint(cs)
        self._npts = 1
        self.boundary = b

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
