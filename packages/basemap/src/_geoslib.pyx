import sys
import numpy
cimport numpy as cnp

__version__ = "1.4.1"


# Need some Python C-API functions for strings.
cdef extern from "Python.h":
    object PyBytes_FromString(char *)


# Taken from `numpy.pxi` in numpy 1.0rc2.
cdef extern from "numpy/arrayobject.h":
    ctypedef int npy_intp
    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef npy_intp *dimensions
        cdef npy_intp *strides
        cdef object base
        cdef int flags
    npy_intp PyArray_SIZE(ndarray arr)
    npy_intp PyArray_ISCONTIGUOUS(ndarray arr)


# Initialize numpy.
cnp.import_array()


# GENERAL NOTES:
# - Remember to call initGEOS() before any use of this library's
#   functions, and call finishGEOS() when done.
# - Currently you have to explicitly GEOSGeom_destroy() all
#   GEOSGeom objects to avoid memory leaks, and to free()
#   all returned char * (unless const). This might change
#   before first release to ensure greater API stability.
cdef extern from "geos_c.h":
    # Supported geometry type IDs.
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
    ctypedef struct GEOSGeometry:
        pass
    ctypedef struct GEOSCoordSequence:
        pass
    # Cython 3: Next ctypedef needs "noexcept" declaration unless
    # the compiler directive "legacy_implicit_noexcept" is used
    # ("noexcept" syntax supported since Cython 0.29.31).
    ctypedef void (*GEOSMessageHandler)(const char *fmt, ...)
    char *GEOSversion()
    void initGEOS(GEOSMessageHandler notice_function, GEOSMessageHandler error_function)
    void finishGEOS()
    GEOSCoordSequence *GEOSCoordSeq_create(unsigned int size, unsigned int dims)
    void GEOSCoordSeq_destroy(GEOSCoordSequence* s)
    int GEOSCoordSeq_setX(GEOSCoordSequence* s, unsigned int idx, double val)
    int GEOSCoordSeq_setY(GEOSCoordSequence* s, unsigned int idx, double val)
    int GEOSCoordSeq_getX(GEOSCoordSequence* s, unsigned int idx, double *val)
    int GEOSCoordSeq_getY(GEOSCoordSequence* s, unsigned int idx, double *val)
    GEOSGeometry *GEOSUnion(GEOSGeometry* g1, GEOSGeometry* g2)
    GEOSGeometry *GEOSUnaryUnion(GEOSGeometry* g1)
    GEOSGeometry *GEOSEnvelope(GEOSGeometry* g1)
    GEOSGeometry *GEOSConvexHull(GEOSGeometry* g1)
    GEOSGeometry *GEOSGeom_createPoint(GEOSCoordSequence* s)
    GEOSGeometry *GEOSGeom_createLineString(GEOSCoordSequence* s)
    GEOSGeometry *GEOSGeom_createPolygon(GEOSGeometry* shell, GEOSGeometry** holes, unsigned int nholes)
    GEOSGeometry *GEOSGeom_createLinearRing(GEOSCoordSequence* s)
    void GEOSGeom_destroy(GEOSGeometry* g)
    # Topology operations: Return NULL on exception.
    GEOSGeometry *GEOSIntersection(GEOSGeometry* g1, GEOSGeometry* g2)
    GEOSGeometry *GEOSSimplify(GEOSGeometry* g1, double tolerance)
    GEOSGeometry *GEOSBuffer(GEOSGeometry* g1, double width, int quadsegs)
    GEOSGeometry *GEOSTopologyPreserveSimplify(GEOSGeometry* g1, double tolerance)
    # Binary/Unary predicate: Return 2 on exception, 1 on true, 0 on false.
    char GEOSIntersects(GEOSGeometry* g1, GEOSGeometry* g2)
    char GEOSWithin(GEOSGeometry* g1, GEOSGeometry* g2)
    char GEOSContains(GEOSGeometry* g1, GEOSGeometry* g2)
    char GEOSisEmpty(GEOSGeometry* g1)
    char GEOSisValid(GEOSGeometry* g1)
    char GEOSisSimple(GEOSGeometry* g1)
    char GEOSisRing(GEOSGeometry* g1)
    # Geometry info.
    char *GEOSGeomType(GEOSGeometry* g1)
    int GEOSGeomTypeId(GEOSGeometry* g1)
    # Functions: Return 0 on exception, 1 otherwise.
    int GEOSArea(GEOSGeometry* g1, double *area)
    int GEOSLength(GEOSGeometry* g1, double *length)
    # Returns -1 on error and 1 for non-multi geoms.
    int GEOSGetNumGeometries(GEOSGeometry* g1)
    # Return NULL on exception, Geometry must be a Collection.
    # Returned object is a pointer to internal storage:
    # it must NOT be destroyed directly.
    GEOSGeometry *GEOSGetGeometryN(GEOSGeometry* g, int n)
    int GEOSGetNumInteriorRings(GEOSGeometry* g1)
    # Return NULL on exception, Geometry must be a Polygon.
    # Returned object is a pointer to internal storage:
    # it must NOT be destroyed directly.
    GEOSGeometry *GEOSGetExteriorRing(GEOSGeometry* g)
    # Return NULL on exception.
    # Geometry must be a LineString, LinearRing or Point.
    GEOSCoordSequence *GEOSGeom_getCoordSeq(const GEOSGeometry* g)
    int GEOSCoordSeq_getSize(const GEOSCoordSequence *s, unsigned int *size)


# Cython 3: Next cdef needs "noexcept" declaration unless
# the compiler directive "legacy_implicit_noexcept" is used
# ("noexcept" syntax supported since Cython 0.29.31).
cdef void notice_h(const char *fmt, ...):
    pass
    #format = PyBytes_FromString(fmt)
    #message = PyBytes_FromString(msg)
    #try:
    #    warn_msg = format % message
    #except:
    #    warn_msg = format
    #sys.stdout.write('GEOS_NOTICE: %s\n' % warn_msg)


# Cython 3: Next cdef needs "noexcept" declaration unless
# the compiler directive "legacy_implicit_noexcept" is used
# ("noexcept" syntax supported since Cython 0.29.31).
# FIXME: The type should be: error_h(const char *fmt, ...), but
# Cython does not currently support varargs functions.
cdef void error_h(const char *fmt, char*msg):
    format = PyBytes_FromString(fmt)
    message = PyBytes_FromString(msg)
    try:
        warn_msg = format % message
    except:
        warn_msg = format
    sys.stderr.write('GEOS_ERROR: %s\n' % warn_msg)


# Check library version.
cdef geos_version():
    return PyBytes_FromString(GEOSversion())

# Module variables.
__geos_version__ = geos_version()
__geos_major_version__ = GEOS_VERSION_MAJOR


# Initialize GEOS (parameters are notice and error function callbacks).
initGEOS(notice_h, <GEOSMessageHandler>error_h)


cdef class BaseGeometry:

    cdef GEOSGeometry *_geom
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
        cdef GEOSGeometry *g1
        cdef GEOSGeometry *g2
        cdef char answer
        g1 = self._geom
        g2 = geom._geom
        answer = GEOSWithin(g1, g2)
        if answer:
            return True
        else:
            return False

    def union(self, BaseGeometry geom):
        cdef GEOSGeometry *g1
        cdef GEOSGeometry *g2
        cdef GEOSGeometry *g3
        cdef const GEOSGeometry *gout
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
        # For multi-geom structures, just return first one.
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
        cdef GEOSGeometry *g1
        cdef GEOSGeometry *g3
        cdef const GEOSGeometry *gout
        cdef double tolerance
        cdef int numgeoms, i, typeid
        g1 = self._geom
        tolerance = tol
        g3 = GEOSSimplify(g1, tolerance)
        typeid = GEOSGeomTypeId(g3)
        if typeid == GEOS_POLYGON:
            b = _get_coords(g3)
            p = Polygon(b)
        elif typeid == GEOS_LINESTRING:
            b = _get_coords(g3)
            p = LineString(b)
        # For multi-geom structures, just return first one.
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
        cdef GEOSGeometry *g1
        cdef GEOSGeometry *g3
        cdef const GEOSGeometry *gout
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
        # For multi-geom structures, just return first one.
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
        cdef GEOSGeometry *g1
        cdef GEOSGeometry *g2
        cdef char answer
        g1 = self._geom
        g2 = geom._geom
        answer = GEOSIntersects(g1, g2)
        if answer:
            return True
        else:
            return False

    def intersection(self, BaseGeometry geom):
        cdef GEOSGeometry *g1
        cdef GEOSGeometry *g2
        cdef GEOSGeometry *g3
        cdef const GEOSGeometry *gout
        cdef char answer
        cdef int numgeoms, i, typeid
        g1 = self._geom
        g2 = geom._geom
        g3 = GEOSIntersection(g1, g2)
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
            for i in range(numgeoms):
                gout = GEOSGetGeometryN(g3, i)
                b = _get_coords(gout)
                p = Polygon(b)
                pout.append(p)
        elif typeid == GEOS_MULTILINESTRING:
            numgeoms = GEOSGetNumGeometries(g3)
            pout = []
            for i in range(numgeoms):
                gout = GEOSGetGeometryN(g3, i)
                b = _get_coords(gout)
                p = LineString(b)
                pout.append(p)
        elif typeid == GEOS_GEOMETRYCOLLECTION:
            numgeoms = GEOSGetNumGeometries(g3)
            pout = []
            for i in range(numgeoms):
                gout = GEOSGetGeometryN(g3, i)
                typeid = GEOSGeomTypeId(gout)
                if typeid == GEOS_POLYGON:
                    b = _get_coords(gout)
                    p = Polygon(b)
                    pout.append(p)
                elif typeid == GEOS_LINESTRING:
                    b = _get_coords(gout)
                    p = LineString(b)
                    pout.append(p)
                else:
                    # More cases might need to be handled here:
                    # - GEOS_MULTILINESTRING
                    # - GEOS_MULTIPOLYGON
                    # - GEOS_GEOMETRYCOLLECTION
                    # The underlying problem is the need of a generic
                    # converter from GEOSGeom pointers to `_geoslib`
                    # objects, since postprocessing `GeometryCollections`
                    # might need recursiveness.
                    pass
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
        return (self.__class__, (self.boundary,))


cdef class Polygon(BaseGeometry):

    def __init__(self, ndarray b):
        cdef unsigned int M, m, i
        cdef double dx, dy
        cdef double *bbuffer
        cdef GEOSCoordSequence *cs
        cdef GEOSGeometry *lr

        # Make sure data is contiguous. If not, make a local copy.
        if not PyArray_ISCONTIGUOUS(b):
            b = b.copy()

        m = b.shape[0]

        # Add closing coordinates to sequence?
        if m > 0 and (b[-1, 0] != b[0, 0] or b[-1, 1] != b[0, 1]):
            M = m + 1
        else:
            M = m
        self._npts = M

        # Create a coordinate sequence.
        cs = GEOSCoordSeq_create(M, 2)

        # Add to coordinate sequence.
        bbuffer = <double *>b.data
        for i in range(m):
            dx = bbuffer[2 * i]
            dy = bbuffer[2 * i + 1]
            # Because of a bug in the GEOS C API, always set X before Y.
            GEOSCoordSeq_setX(cs, i, dx)
            GEOSCoordSeq_setY(cs, i, dy)

        # Add closing coordinates to sequence?
        if M > m:
            dx = bbuffer[0]
            dy = bbuffer[1]
            GEOSCoordSeq_setX(cs, M - 1, dx)
            GEOSCoordSeq_setY(cs, M - 1, dy)

        # Create LinearRing.
        lr = GEOSGeom_createLinearRing(cs)

        # Create Polygon from LinearRing (assuming no holes).
        self._geom = GEOSGeom_createPolygon(lr, NULL, 0)
        self.boundary = b

    def area(self):
        cdef double area
        GEOSArea(self._geom, &area)
        return area


cdef class LineString(BaseGeometry):

    def __init__(self, ndarray b):

        cdef double dx, dy
        cdef GEOSCoordSequence *cs
        cdef int i, M
        cdef double *bbuffer

        # Make sure data is contiguous. If not, make a local copy.
        if not PyArray_ISCONTIGUOUS(b):
            b = b.copy()

        M = b.shape[0]
        self._npts = M

        # Create a coordinate sequence.
        cs = GEOSCoordSeq_create(M, 2)

        # Add to coordinate sequence.
        bbuffer = <double *>b.data
        for i in range(M):
            dx = bbuffer[2 * i]
            dy = bbuffer[2 * i + 1]
            # Because of a bug in the GEOS C API, always set X before Y.
            GEOSCoordSeq_setX(cs, i, dx)
            GEOSCoordSeq_setY(cs, i, dy)

        # Create LineString.
        self._geom = GEOSGeom_createLineString(cs)
        self.boundary = b


cdef class Point(BaseGeometry):

    cdef public x, y
    def __init__(self, b):

        cdef double dx, dy
        cdef GEOSCoordSequence *cs

        # Create a coordinate sequence.
        cs = GEOSCoordSeq_create(1, 2)
        dx = b[0]
        dy = b[1]
        GEOSCoordSeq_setX(cs, 0, dx)
        GEOSCoordSeq_setY(cs, 0, dy)
        self._geom = GEOSGeom_createPoint(cs)
        self._npts = 1
        self.boundary = b


cdef _get_coords(const GEOSGeometry *geom):

    cdef const GEOSCoordSequence *cs
    cdef const GEOSGeometry *lr
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
    b = numpy.empty((M, 2), numpy.float64)
    bbuffer = <double *>b.data
    for i in range(M):
        GEOSCoordSeq_getX(cs, i, &dx)
        GEOSCoordSeq_getY(cs, i, &dy)
        bbuffer[2 * i] = dx
        bbuffer[2 * i + 1] = dy
    return b
