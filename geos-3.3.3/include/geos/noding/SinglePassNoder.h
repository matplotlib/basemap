/**********************************************************************
 * $Id: SinglePassNoder.h 2556 2009-06-06 22:22:28Z strk $
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2006      Refractions Research Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation. 
 * See the COPYING file for more information.
 *
 **********************************************************************/

#ifndef GEOS_NODING_SINGLEPASSNODER_H
#define GEOS_NODING_SINGLEPASSNODER_H

#include <geos/export.h>

#include <vector>

#include <geos/inline.h>

#include <geos/noding/Noder.h>

// Forward declarations
namespace geos {
	namespace noding {
		class SegmentString;
		class SegmentIntersector;
	}
}

namespace geos {
namespace noding { // geos.noding


/**
 * Base class for {@link Noder}s which make a single
 * pass to find intersections.
 * This allows using a custom {@link SegmentIntersector}
 * (which for instance may simply identify intersections, rather than
 * insert them).
 *
 * Last port: noding/SinglePassNoder.java rev. 1.3 (JTS-1.7)
 *
 * TODO: Noder inheritance (that's just an interface!)
 *
 */
class GEOS_DLL SinglePassNoder : public Noder { // implements Noder

protected:

	/// Externally owned
	SegmentIntersector* segInt;

public:

	SinglePassNoder(SegmentIntersector* nSegInt=NULL): segInt(nSegInt) {}

	virtual ~SinglePassNoder() {}

	/**
	 * Sets the SegmentIntersector to use with this noder.
	 * A SegmentIntersector will normally add intersection nodes
	 * to the input segment strings, but it may not - it may
	 * simply record the presence of intersections.
	 * However, some Noders may require that intersections be added.
	 *
	 * @param newSegInt
	 */
	virtual void setSegmentIntersector(SegmentIntersector* newSegInt) {
	  segInt = newSegInt;
	}

	/**
	 * Computes the noding for a collection of {@link SegmentString}s.
	 *
	 * @param segStrings a collection of {@link SegmentString}s to node
	 */
	virtual void computeNodes(std::vector<SegmentString*>* segStrings)=0;

	/**
	 * Returns a {@link Collection} of fully noded {@link SegmentStrings}.
	 * The SegmentStrings have the same context as their parent.
	 *
	 * @return a Collection of SegmentStrings
	 */
	virtual std::vector<SegmentString*>* getNodedSubstrings() const=0;
	
};

} // namespace geos.noding
} // namespace geos

//#ifdef GEOS_INLINE
//# include "geos/noding/SinglePassNoder.inl"
//#endif

#endif // GEOS_NODING_SINGLEPASSNODER_H

/**********************************************************************
 * $Log$
 * Revision 1.2  2006/03/24 09:52:41  strk
 * USE_INLINE => GEOS_INLINE
 *
 * Revision 1.1  2006/03/09 16:46:49  strk
 * geos::geom namespace definition, first pass at headers split
 *
 **********************************************************************/

