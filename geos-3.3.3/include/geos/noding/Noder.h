/**********************************************************************
 * $Id: Noder.h 2556 2009-06-06 22:22:28Z strk $
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

#ifndef GEOS_NODING_NODER_H
#define GEOS_NODING_NODER_H

#include <geos/export.h>

#include <vector>
#include <iostream>

#include <geos/inline.h>

// Forward declarations
namespace geos {
	namespace noding {
		class SegmentString;
	}
}

namespace geos {
namespace noding { // geos.noding


/** \brief
 * Computes all intersections between segments in a set of SegmentString.
 *
 * Intersections found are represented as {@link SegmentNode}s and added to the
 * {@link SegmentString}s in which they occur.
 * As a final step in the noding a new set of segment strings split
 * at the nodes may be returned.
 *
 * Last port: noding/Noder.java rev. 1.8 (JTS-1.7)
 *
 * TODO: this was really an interface, we should avoid making it a Base class
 *
 */
class GEOS_DLL Noder {
public:
	/** \brief
	 * Computes the noding for a collection of {@link SegmentString}s.
	 *
	 * Some Noders may add all these nodes to the input SegmentStrings;
	 * others may only add some or none at all.
	 *
	 * @param segStrings a collection of {@link SegmentString}s to node
	 */
	virtual void computeNodes(std::vector<SegmentString*>* segStrings)=0;

	/** \brief
	 * Returns a {@link Collection} of fully noded {@link SegmentStrings}.
	 * The SegmentStrings have the same context as their parent.
	 *
	 * @return a newly allocated std::vector of const SegmentStrings.
	 *         Caller is responsible to delete it
	 */
	virtual std::vector<SegmentString*>* getNodedSubstrings() const=0;

	virtual ~Noder() {}

protected:
	Noder(){};
};

} // namespace geos.noding
} // namespace geos

//#ifdef GEOS_INLINE
//# include "geos/noding/Noder.inl"
//#endif

#endif // GEOS_NODING_NODER_H

/**********************************************************************
 * $Log$
 * Revision 1.2  2006/03/24 09:52:41  strk
 * USE_INLINE => GEOS_INLINE
 *
 * Revision 1.1  2006/03/09 16:46:49  strk
 * geos::geom namespace definition, first pass at headers split
 *
 **********************************************************************/

