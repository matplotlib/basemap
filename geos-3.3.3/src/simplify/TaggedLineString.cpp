/**********************************************************************
 * $Id: TaggedLineString.cpp 3574 2012-03-22 08:34:59Z strk $
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2006 Refractions Research Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Licence as published
 * by the Free Software Foundation. 
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: simplify/TaggedLineString.java rev. 1.2 (JTS-1.7.1)
 *
 **********************************************************************/

#include <geos/simplify/TaggedLineString.h>
#include <geos/simplify/TaggedLineSegment.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/LineString.h>
#include <geos/geom/Geometry.h> // for auto_ptr destructor 
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/CoordinateSequenceFactory.h>

#include <cassert>
#include <memory>

#ifndef GEOS_DEBUG
#define GEOS_DEBUG 0
#endif

#ifdef GEOS_DEBUG
#include <iostream>
#endif

using namespace geos::geom;
using namespace std;

namespace geos {
namespace simplify { // geos::simplify

/*public*/
TaggedLineString::TaggedLineString(const geom::LineString* nParentLine,
			std::size_t nMinimumSize)
	:
	parentLine(nParentLine),
	minimumSize(nMinimumSize)
{
	init();
}

/*public*/
TaggedLineString::~TaggedLineString()
{
#if GEOS_DEBUG
	cerr << "TaggedLineString[" << this << "] destructor"
	     << endl;
#endif

	for (std::size_t i=0, n=segs.size(); i<n; i++)
		delete segs[i];

	for (std::size_t i=0, n=resultSegs.size(); i<n; i++)
		delete resultSegs[i];
}

/*private*/
void
TaggedLineString::init()
{
	assert(parentLine);
	const CoordinateSequence* pts = parentLine->getCoordinatesRO();

#if GEOS_DEBUG
	cerr << "TaggedLineString[" << this << "] pts.size() " << pts->size()
	     << endl;
#endif

	if ( pts->size() )
	{

		segs.reserve(pts->size()-1);

		for (std::size_t i=0, n=pts->size()-1; i<n; i++)
		{
			TaggedLineSegment* seg = new TaggedLineSegment(
					pts->getAt(i),
					pts->getAt(i+1),
					parentLine, i);

			segs.push_back(seg);
		}

	}

#if GEOS_DEBUG
	cerr << "TaggedLineString[" << this << "] segs.size " << segs.size()
	    << endl;
	cerr << "TaggedLineString[" << this << "] resultSegs.size " << resultSegs.size()
	    << endl;
#endif
}

/*public*/
std::size_t
TaggedLineString::getMinimumSize() const
{
	return minimumSize;
}

/*public*/
const geom::LineString* 
TaggedLineString::getParent() const
{
	return parentLine;
}

/*public*/
const CoordinateSequence*
TaggedLineString::getParentCoordinates() const
{
	assert(parentLine);
	return parentLine->getCoordinatesRO();
}

/*public*/
CoordinateSequence::AutoPtr
TaggedLineString::getResultCoordinates() const
{

#if GEOS_DEBUG
	cerr << __FUNCTION__ << " resultSegs.size: "
	     << resultSegs.size() << endl;
#endif

	CoordVectPtr pts = extractCoordinates(resultSegs);

#if GEOS_DEBUG
	cerr << __FUNCTION__ << " extracted Coords.size: "
	     << pts->size() << endl;
#endif


	CoordVect* v = pts.release();
	return CoordinateSequence::AutoPtr(parentLine->getFactory()->getCoordinateSequenceFactory()->create(v));

}

/*private static*/
TaggedLineString::CoordVectPtr
TaggedLineString::extractCoordinates(
		const std::vector<TaggedLineSegment*>& segs)
{
	CoordVectPtr pts(new CoordVect());

#if GEOS_DEBUG
	cerr << __FUNCTION__ << " segs.size: " << segs.size() << endl;
#endif

	std::size_t i=0, size=segs.size();

	if ( size ) {
		for (; i<size; i++)
		{
			TaggedLineSegment* seg = segs[i];
			assert(seg);
			pts->push_back(seg->p0);
		}

		// add last point
		pts->push_back(segs[size-1]->p1);
	}

	return pts;
}

/*public*/
std::size_t
TaggedLineString::getResultSize() const
{
	unsigned resultSegsSize = resultSegs.size();
	return resultSegsSize == 0 ? 0 : resultSegsSize + 1;
}

/*public*/
TaggedLineSegment*
TaggedLineString::getSegment(std::size_t i) 
{
	return segs[i];
}

/*public*/
const TaggedLineSegment*
TaggedLineString::getSegment(std::size_t i) const
{
	return segs[i];
}

/*public*/
vector<TaggedLineSegment*>&
TaggedLineString::getSegments()
{
	assert(0);
	return segs;
}

/*public*/
const vector<TaggedLineSegment*>&
TaggedLineString::getSegments() const
{
	return segs;
}

/*public*/
auto_ptr<Geometry>
TaggedLineString::asLineString() const
{
	return parentLine->getFactory()->createLineString(
			getResultCoordinates());
}

/*public*/
auto_ptr<Geometry>
TaggedLineString::asLinearRing() const
{
	return parentLine->getFactory()->createLinearRing(
			getResultCoordinates());
}

/*public*/
void
TaggedLineString::addToResult(auto_ptr<TaggedLineSegment> seg)
{
#if GEOS_DEBUG
	cerr << "TaggedLineString[" << this << "] adding "
	     << " seg " << seg.get() << " to result"
	     << endl;
#endif
	resultSegs.push_back(seg.release());
#if GEOS_DEBUG
	cerr << "TaggedLineString[" << this << "] adding "
	     << " seg " << seg.get() << " to result"
	     << endl;
#endif
}

} // namespace geos::simplify
} // namespace geos

/**********************************************************************
 * $Log$
 * Revision 1.5  2006/06/12 11:29:24  strk
 * unsigned int => size_t
 *
 * Revision 1.4  2006/04/13 21:52:35  strk
 * Many debugging lines and assertions added. Fixed bug in TaggedLineString class.
 *
 * Revision 1.3  2006/04/13 09:21:46  mloskot
 * Removed definition of copy ctor and assignment operator for TaggedLineString class.
 * According to following rule: Declaring, but not defining, private copy operations has
 * the effect of "turning off" copying for the class.
 *
 * Revision 1.2  2006/04/12 15:20:37  strk
 * LineSegmentIndex class
 *
 * Revision 1.1  2006/04/12 14:22:12  strk
 * Initial implementation of TaggedLineSegment and TaggedLineString classes
 *
 **********************************************************************/
