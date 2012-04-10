/**********************************************************************
 * $Id: SingleSidedBufferResultMatcher.h 2809 2009-12-06 01:05:24Z mloskot $
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2009  Sandro Santilli <strk@keybit.net>
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation. 
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: original work
 *
 **********************************************************************/

#ifndef XMLTESTER_SINGLESIDEDBUFFERRESULTMATCHER_H
#define XMLTESTER_SINGLESIDEDBUFFERRESULTMATCHER_H

// Forward declarations
namespace geos {
	namespace geom {
		class Geometry;
	}
}

namespace geos {
namespace xmltester {

class SingleSidedBufferResultMatcher
{
public:
	bool isBufferResultMatch(const geom::Geometry& actualBuffer,
	                         const geom::Geometry& expectedBuffer,
	                         double distance);

private:

	static double MAX_HAUSDORFF_DISTANCE_FACTOR;

	/*
	 * The minimum distance tolerance which will be used.
	 * This is required because densified vertices do no lie
	 * precisely on their parent segment.
	 */
	static double MIN_DISTANCE_TOLERANCE;

	bool isBoundaryHausdorffDistanceInTolerance(
                        const geom::Geometry& actualBuffer,
	                const geom::Geometry& expectedBuffer,
	                double distance);
};

} // namespace geos::xmltester
} // namespace geos

#endif // XMLTESTER_SINGLESIDEDBUFFERRESULTMATCHER_H
