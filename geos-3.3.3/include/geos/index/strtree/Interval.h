/**********************************************************************
 * $Id: Interval.h 2556 2009-06-06 22:22:28Z strk $
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2006 Refractions Research Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation. 
 * See the COPYING file for more information.
 *
 **********************************************************************/

#ifndef GEOS_INDEX_STRTREE_INTERVAL_H
#define GEOS_INDEX_STRTREE_INTERVAL_H

#include <geos/export.h>

namespace geos {
namespace index { // geos::index
namespace strtree { // geos::index::strtree

/// A contiguous portion of 1D-space. Used internally by SIRtree.
//
/// @see SIRtree
///
class GEOS_DLL Interval {
public:
	Interval(Interval *other);
	Interval(double newMin, double newMax);
	double getCentre();
	Interval* expandToInclude(Interval *other);
	bool intersects(Interval *other);
	bool equals(void *o);
private:
	double imin;
	double imax;
};


} // namespace geos::index::strtree
} // namespace geos::index
} // namespace geos

#endif // GEOS_INDEX_STRTREE_INTERVAL_H

/**********************************************************************
 * $Log$
 * Revision 1.1  2006/03/21 10:47:34  strk
 * indexStrtree.h split
 *
 **********************************************************************/

