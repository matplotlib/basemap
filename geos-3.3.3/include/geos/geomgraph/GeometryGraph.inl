/**********************************************************************
 * $Id: GeometryGraph.inl 2557 2009-06-08 09:30:55Z strk $
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2005-2006 Refractions Research Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation. 
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: geomgraph/GeometryGraph.java rev. 1.9 (JTS-1.10)
 *
 **********************************************************************/

#ifndef GEOS_GEOMGRAPH_GEOMETRYGRAPH_INL
#define GEOS_GEOMGRAPH_GEOMETRYGRAPH_INL

#include <geos/geomgraph/GeometryGraph.h>

namespace geos {
namespace geomgraph { // geos::geomgraph

INLINE index::SegmentIntersector*
GeometryGraph::computeSelfNodes(
		algorithm::LineIntersector& li,
		bool computeRingSelfNodes)
{
	return computeSelfNodes(&li, computeRingSelfNodes);
}

INLINE void
GeometryGraph::getBoundaryNodes(std::vector<Node*>&bdyNodes)
{
	nodes->getBoundaryNodes(argIndex, bdyNodes);
}

INLINE const geom::Geometry*
GeometryGraph::getGeometry()
{
	return parentGeom;
}

INLINE 
GeometryGraph::~GeometryGraph()
{
}

} // namespace geos::geomgraph
} // namespace geos

#endif // GEOS_GEOMGRAPH_GEOMETRYGRAPH_INL
