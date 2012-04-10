/**********************************************************************
 * $Id: UniqueCoordinateArrayFilter.h 2958 2010-03-29 11:29:40Z mloskot $
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 * Copyright (C) 2006 Refractions Research Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation. 
 * See the COPYING file for more information.
 *
 **********************************************************************/

#ifndef GEOS_UTIL_UNIQUECOORDINATEARRAYFILTER_H
#define GEOS_UTIL_UNIQUECOORDINATEARRAYFILTER_H

#include <geos/export.h>
#include <cassert>
#include <set>
#include <vector>

#include <geos/geom/CoordinateFilter.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/Coordinate.h>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4251) // warning C4251: needs to have dll-interface to be used by clients of class
#endif

namespace geos {
namespace util { // geos::util

/*
 *  A CoordinateFilter that fills a vector of Coordinate const pointers.
 *  The set of coordinates contains no duplicate points.
 *
 *  Last port: util/UniqueCoordinateArrayFilter.java rev. 1.17
 */
class GEOS_DLL UniqueCoordinateArrayFilter: public geom::CoordinateFilter
{
public:
	/**
	 * Constructs a CoordinateArrayFilter.
	 *
	 * @param target The destination set. 
	 */
	UniqueCoordinateArrayFilter(geom::Coordinate::ConstVect &target)
		: pts(target)
	{}

	/**
	 * Destructor.
	 * Virtual dctor promises appropriate behaviour when someone will
	 * delete a derived-class object via a base-class pointer.
	 * http://www.parashift.com/c++-faq-lite/virtual-functions.html#faq-20.7
	 */
	virtual ~UniqueCoordinateArrayFilter() {}

	/**
	 * Performs a filtering operation with or on coord in "read-only" mode.
	 * @param coord The "read-only" Coordinate to which
	 * 				the filter is applied.
	 */
	virtual void filter_ro(const geom::Coordinate *coord)
	{
		if ( uniqPts.insert(coord).second )
		{
			pts.push_back(coord);
		}
    }

private:
	geom::Coordinate::ConstVect &pts;	// target set reference
	geom::Coordinate::ConstSet uniqPts; 	// unique points set

    // Declare type as noncopyable
    UniqueCoordinateArrayFilter(const UniqueCoordinateArrayFilter& other);
    UniqueCoordinateArrayFilter& operator=(const UniqueCoordinateArrayFilter& rhs);
};

} // namespace geos::util
} // namespace geos

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // GEOS_UTIL_UNIQUECOORDINATEARRAYFILTER_H

/**********************************************************************
 * $Log$
 * Revision 1.4  2006/06/19 23:33:03  strk
 * Don't *require* CoordinateFilters to define both read-only and read-write methods.
 *
 * Revision 1.3  2006/06/12 10:10:39  strk
 * Fixed getGeometryN() to take size_t rather then int, changed unsigned int parameters to size_t.
 *
 * Revision 1.2  2006/04/10 09:21:23  mloskot
 * Added new test for UniqueCoordinateArrayFilter class. Small fixes related to signed/unsigned comparison.
 *
 * Revision 1.1  2006/03/09 16:46:49  strk
 * geos::geom namespace definition, first pass at headers split
 *
 **********************************************************************/
