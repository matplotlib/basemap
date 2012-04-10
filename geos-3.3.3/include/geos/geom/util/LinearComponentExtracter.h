/**********************************************************************
 * $Id: LinearComponentExtracter.h 2772 2009-12-03 19:30:54Z mloskot $
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

#ifndef GEOS_GEOM_UTIL_LINEARCOMPONENTEXTRACTER_H
#define GEOS_GEOM_UTIL_LINEARCOMPONENTEXTRACTER_H


#include <geos/export.h>
#include <vector>

#include <geos/geom/GeometryComponentFilter.h>
#include <geos/geom/Geometry.h> // to be removed when we have the .inl
#include <geos/geom/LineString.h> // to be removed when we have the .inl
//#include <geos/platform.h>

namespace geos {
namespace geom { // geos.geom
namespace util { // geos.geom.util

/**
 * Extracts all the 1-dimensional (LineString) components from a Geometry.
 */
class GEOS_DLL LinearComponentExtracter: public GeometryComponentFilter {

private:

	LineString::ConstVect &comps;

    // Declare type as noncopyable
    LinearComponentExtracter(const LinearComponentExtracter& other);
    LinearComponentExtracter& operator=(const LinearComponentExtracter& rhs);

public:
	/**
	 * Push the linear components from a single geometry into
	 * the provided vector.
	 * If more than one geometry is to be processed, it is more
	 * efficient to create a single LinearComponentExtracterFilter instance
	 * and pass it to multiple geometries.
	 */
	static void getLines(const Geometry &geom, std::vector<const LineString*> &ret)
	{
		LinearComponentExtracter lce(ret);
		geom.apply_ro(&lce);
	}

	/**
	 * Constructs a LinearComponentExtracterFilter with a list in which
	 * to store LineStrings found.
	 */
	LinearComponentExtracter(std::vector<const LineString*> &newComps)
		:
		comps(newComps)
		{}

	void filter_rw(Geometry *geom)
	{
if ( const LineString *ls=dynamic_cast<const LineString *>(geom) )
		comps.push_back(ls);
	}

	void filter_ro(const Geometry *geom)
	{
if ( const LineString *ls=dynamic_cast<const LineString *>(geom) )
		comps.push_back(ls);
	}

};

} // namespace geos.geom.util
} // namespace geos.geom
} // namespace geos

#endif

/**********************************************************************
 * $Log$
 * Revision 1.1  2006/03/09 16:46:49  strk
 * geos::geom namespace definition, first pass at headers split
 *
 **********************************************************************/
