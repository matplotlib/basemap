/**********************************************************************
 * $Id: EdgeEnd.cpp 2755 2009-11-30 17:29:48Z mloskot $
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2005-2006 Refractions Research Inc.
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation. 
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: geomgraph/EdgeEnd.java rev. 1.6 (JTS-1.10)
 *
 **********************************************************************/

#include <geos/geomgraph/EdgeEnd.h>
#include <geos/geomgraph/Node.h> // for assertions 
#include <geos/algorithm/CGAlgorithms.h>
#include <geos/geomgraph/Label.h>
#include <geos/geomgraph/Quadrant.h>
#include <geos/geom/Coordinate.h>

#include <typeinfo>
#include <cmath>
#include <sstream>
#include <iostream>
#include <string>
#include <cassert>

using namespace geos::geom;

namespace geos {
namespace geomgraph { // geos.geomgraph

using namespace geos::algorithm;

/*public*/
EdgeEnd::~EdgeEnd()
{
	delete label;
}

/*public*/
EdgeEnd::EdgeEnd()
	:
	edge(NULL),
	label(NULL),
	node(NULL),
	dx(0.0),
	dy(0.0),
	quadrant(0)
{
}

/*protected*/
EdgeEnd::EdgeEnd(Edge* newEdge)
	:
	edge(newEdge),
	label(NULL),
	node(NULL),
	dx(0.0),
	dy(0.0),
	quadrant(0)
{
}

/*public*/
EdgeEnd::EdgeEnd(Edge* newEdge, const Coordinate& newP0,
		const Coordinate& newP1, Label* newLabel)
	:
	edge(newEdge),
	label(newLabel),
	node(NULL),
	dx(0.0),
	dy(0.0),
	quadrant(0)
{
	init(newP0, newP1);
}

/*public*/
void
EdgeEnd::init(const Coordinate& newP0, const Coordinate& newP1)
{
	p0=newP0;
	p1=newP1;
	dx=p1.x-p0.x;
	dy=p1.y-p0.y;
	quadrant=Quadrant::quadrant(dx,dy);

	// "EdgeEnd with identical endpoints found");
	assert(!(dx == 0 && dy == 0));
}

/*public*/
Label*
EdgeEnd::getLabel()
{
	return label;
}

/*public*/
Coordinate&
EdgeEnd::getCoordinate()
{
	return p0;
}

/*public*/
Coordinate&
EdgeEnd::getDirectedCoordinate()
{
	return p1;
}

/*public*/
int
EdgeEnd::getQuadrant()
{
	return quadrant;
}

/*public*/
double
EdgeEnd::getDx()
{
	return dx;
}

/*public*/
double
EdgeEnd::getDy()
{
	return dy;
}

/*public*/
void
EdgeEnd::setNode(Node* newNode)
{
	node=newNode;
	assert(node->getCoordinate().equals2D(p0));
}

/*public*/
Node*
EdgeEnd::getNode()
{
	return node;
}

/*public*/
int
EdgeEnd::compareTo(const EdgeEnd* e) const
{
	return compareDirection(e);
}

/*public*/
int
EdgeEnd::compareDirection(const EdgeEnd* e) const
{
	assert(e);
	if (dx == e->dx && dy == e->dy)
		return 0;

	// if the rays are in different quadrants,
	// determining the ordering is trivial
	if (quadrant>e->quadrant) return 1;
	if (quadrant<e->quadrant) return -1;

	// vectors are in the same quadrant - check relative
	// orientation of direction vectors
	// this is > e if it is CCW of e
	return CGAlgorithms::computeOrientation(e->p0, e->p1, p1);
}

/*public*/
void
EdgeEnd::computeLabel(const algorithm::BoundaryNodeRule& /*boundaryNodeRule*/)
{
	// subclasses should override this if they are using labels
}

/*public*/
std::string
EdgeEnd::print()
{
	std::ostringstream s;

	s<<*this;

	return s.str();
}

std::ostream&
operator<< (std::ostream& os, const EdgeEnd& ee)
{
	os << "EdgeEnd: ";
	os << ee.p0;
	os << " - ";
	os << ee.p1;
	os << " ";
	os << ee.quadrant << ":" << std::atan2(ee.dy, ee.dx);
	os << "  ";
	os << *(ee.label);

	return os;
}


} // namespace geos.geomgraph
} // namespace geos

/**********************************************************************
 * $Log$
 * Revision 1.19  2006/06/14 14:32:20  strk
 * EdgeEnd::getEdge() made non-virtual and inlined.
 *
 * Revision 1.18  2006/04/08 13:05:49  strk
 * Added assertion
 *
 * Revision 1.17  2006/04/06 12:47:31  strk
 * Fixed output function
 *
 * Revision 1.16  2006/04/06 09:39:55  strk
 * Added operator<<
 *
 * Revision 1.15  2006/04/03 17:05:22  strk
 * Assertion checking, port info, cleanups
 *
 * Revision 1.14  2006/03/15 17:16:29  strk
 * streamlined headers inclusion
 *
 * Revision 1.13  2006/03/09 16:46:47  strk
 * geos::geom namespace definition, first pass at headers split
 *
 * Revision 1.12  2006/03/06 19:40:46  strk
 * geos::util namespace. New GeometryCollection::iterator interface, many cleanups.
 *
 * Revision 1.11  2006/03/03 10:46:21  strk
 * Removed 'using namespace' from headers, added missing headers in .cpp files, removed useless includes in headers (bug#46)
 *
 * Revision 1.10  2006/02/28 14:34:04  strk
 * Added many assertions and debugging output hunting for a bug in BufferOp
 *
 * Revision 1.9  2006/02/19 19:46:49  strk
 * Packages <-> namespaces mapping for most GEOS internal code (uncomplete, but working). Dir-level libs for index/ subdirs.
 **********************************************************************/
