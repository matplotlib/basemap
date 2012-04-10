// Test Suite for C-API LineString and Point functions

#include <tut.hpp>
// geos
#include <geos_c.h>
// std
#include <cstdarg>
#include <cstdio>
#include <cstdlib>

namespace tut
{
    //
    // Test Group
    //

    // Common data used in test cases.
    struct test_capilinestringpoint_data
    {
        GEOSGeometry* geom1_;

        static void notice(const char *fmt, ...)
        {
            std::fprintf( stdout, "NOTICE: ");

            va_list ap;
            va_start(ap, fmt);
            std::vfprintf(stdout, fmt, ap);
            va_end(ap);

            std::fprintf(stdout, "\n");
        }

        test_capilinestringpoint_data()
            : geom1_(0)
        {
            initGEOS(notice, notice);
        }

        ~test_capilinestringpoint_data()
        {
            GEOSGeom_destroy(geom1_);
            geom1_ = 0;
            finishGEOS();
        }

    };

    typedef test_group<test_capilinestringpoint_data> group;
    typedef group::object object;

    group test_capilinestringpoint_group("capi::GEOSLineStringPoint");

    //
    // Test Cases
    //

    template<>
    template<>
    void object::test<1>()
    {
        geom1_ = GEOSGeomFromWKT("LINESTRING(0 0, 5 5, 10 10)");
        GEOSGeometry *geom2;
        double x, y;
        ensure( 0 != geom1_ );

        char const r1 = GEOSisClosed(geom1_);

        ensure_equals(r1, 0);

        geom2 = GEOSGeomGetPointN(geom1_, 0);
        GEOSGeomGetX(geom2, &x);
        GEOSGeomGetY(geom2, &y);

        ensure_equals(x, 0);
        ensure_equals(y, 0);

        GEOSGeom_destroy(geom2);

        geom2 = GEOSGeomGetStartPoint(geom1_);
        GEOSGeomGetX(geom2, &x);
        GEOSGeomGetY(geom2, &y);

        ensure_equals(x, 0);
        ensure_equals(y, 0);

        GEOSGeom_destroy(geom2);

        geom2 = GEOSGeomGetEndPoint(geom1_);
        GEOSGeomGetX(geom2, &x);
        GEOSGeomGetY(geom2, &y);

        ensure_equals(x, 10);
        ensure_equals(y, 10);

        GEOSGeom_destroy(geom2);
    }

    template<>
    template<>
    void object::test<2>()
    {
        geom1_ = GEOSGeomFromWKT("LINESTRING(0 0, 5 5, 10 10)");
        double length;
        ensure( 0 != geom1_ );

        GEOSGeomGetLength(geom1_, &length);
        ensure(length != 0.0);
    }

    template<>
    template<>
    void object::test<3>()
    {
        geom1_ = GEOSGeomFromWKT("LINESTRING(0 0, 5 5, 10 10)");
        int points;
        ensure( 0 != geom1_ );

        points = GEOSGeomGetNumPoints(geom1_);
        ensure_equals(points, 3);
    }
} // namespace tut
