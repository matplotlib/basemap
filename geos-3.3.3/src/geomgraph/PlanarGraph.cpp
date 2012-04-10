/**********************************************************************
 * $Id: PlanarGraph.cpp 3281 2011-04-18 15:08:24Z strk $
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
 * Last port: geomgraph/PlanarGraph.java rev. 1.6 (JTS-1.10)
 *
 **********************************************************************/

#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/Location.h>

#include <geos/geomgraph/PlanarGraph.h>
#include <geos/geomgraph/Node.h>
#include <geos/geomgraph/NodeFactory.h>
#include <geos/geomgraph/Edge.h>
#include <geos/geomgraph/EdgeEndStar.h>
#include <geos/geomgraph/DirectedEdge.h>
#include <geos/geomgraph/DirectedEdgeStar.h>
#include <geos/geomgraph/NodeMap.h>
#include <geos/geomgraph/Quadrant.h>

#include <geos/algorithm/CGAlgorithms.h>

#include <vector>
#include <sstream>
#include <string>
#include <cassert>


#ifndef GEOS_DEBUG
#define GEOS_DEBUG 0
#endif

using namespace std;
using namespace geos::algorithm;
using namespace geos::geom;

namespace geos {
namespace geomgraph { // geos.geomgraph

/*public*/
PlanarGraph::PlanarGraph(const NodeFactory &nodeFact)
	:
	edges(new vector<Edge*>()),
	nodes(new NodeMap(nodeFact)),
	edgeEndList(new vector<EdgeEnd*>())
{
}

/*public*/
PlanarGraph::PlanarGraph()
	:
	edges(new vector<Edge*>()),
	nodes(new NodeMap(NodeFactory::instance())),
	edgeEndList(new vector<EdgeEnd*>())
{
}

/*public*/
PlanarGraph::~PlanarGraph()
{
	delete nodes;
#if 1 // FIXME: PlanarGraph should *not* own edges!
	for(size_t i=0, n=edges->size(); i<n; i++) {
		delete (*edges)[i];
	}
#endif
	delete edges;

	for(size_t i=0, n=edgeEndList->size(); i<n; i++) {
		delete (*edgeEndList)[i];
	}
	delete edgeEndList;
}

/*public*/
vector<Edge*>::iterator
PlanarGraph::getEdgeIterator()
{
	assert(edges);
	return edges->begin();
}

/*public*/
vector<EdgeEnd*> *
PlanarGraph::getEdgeEnds()
{
	return edgeEndList;
}

/*public*/
bool
PlanarGraph::isBoundaryNode(int geomIndex, const Coordinate& coord)
{
	assert(nodes);

	Node *node=nodes->find(coord);
	if (node==NULL) return false;

	Label *label=node->getLabel();
	if (label!=NULL && label->getLocation(geomIndex)==Location::BOUNDARY)
			return true;

	return false;
}

/*protected*/
void
PlanarGraph::insertEdge(Edge *e)
{
	assert(e);
	assert(edges);
	edges->push_back(e);
}

/*public*/
void
PlanarGraph::add(EdgeEnd* e)
{

	assert(e);
	assert(nodes);
	nodes->add(e);

	assert(edgeEndList);
	edgeEndList->push_back(e);
}

/*public*/
NodeMap::iterator
PlanarGraph::getNodeIterator()
{
	assert(nodes);
	return nodes->begin();
}

/*public*/
void
PlanarGraph::getNodes(vector<Node*>& values)
{
	assert(nodes);
	NodeMap::iterator it=nodes->nodeMap.begin();
	while(it!=nodes->nodeMap.end()) {
		assert(it->second);
		values.push_back(it->second);
		it++;
	}
}

// arg cannot be const, NodeMap::addNode will
// occasionally label-merge first arg.
/*public*/
Node*
PlanarGraph::addNode(Node *node)
{
	assert(nodes);
#if GEOS_DEBUG
	cerr << "PlanarGraph::addNode(Node * " << *node 
		<< ")" << endl;
#endif
	return nodes->addNode(node);
}

/*public*/
Node*
PlanarGraph::addNode(const Coordinate& coord)
{
#if GEOS_DEBUG
	cerr << "PlanarGraph::addNode(Coordinate& "
		<< coord << ")" << endl;
#endif
	return nodes->addNode(coord);
}

/*public*/
Node*
PlanarGraph::find(Coordinate& coord)
{
	assert(nodes);
	return nodes->find(coord);
}

/*public*/
void
PlanarGraph::addEdges(const vector<Edge*>& edgesToAdd)
{
	// create all the nodes for the edges
	for (vector<Edge*>::const_iterator it=edgesToAdd.begin(),
		endIt=edgesToAdd.end(); it!=endIt; ++it)
	{
		Edge *e=*it;
		assert(e);
		edges->push_back(e);

		// PlanarGraph destructor will delete all DirectedEdges 
		// in edgeEndList, which is where these are added
		// by the ::add(EdgeEnd) call
		DirectedEdge *de1=new DirectedEdge(e, true);
		DirectedEdge *de2=new DirectedEdge(e, false);

		de1->setSym(de2);
		de2->setSym(de1);
		add(de1);
		add(de2);
	}
}

/*public static*/
void
PlanarGraph::linkResultDirectedEdges()
{
#if GEOS_DEBUG
	cerr<<"PlanarGraph::linkResultDirectedEdges called"<<endl;
#endif
	NodeMap::iterator nodeit=nodes->nodeMap.begin();
	for (;nodeit!=nodes->nodeMap.end();nodeit++) {
		Node *node=nodeit->second;
		assert(node);

		EdgeEndStar* ees=node->getEdges();
		assert(ees);
		assert(dynamic_cast<DirectedEdgeStar*>(ees));
		DirectedEdgeStar* des = static_cast<DirectedEdgeStar*>(ees);

		// this might throw an exception
		des->linkResultDirectedEdges();
	}
}

/*
 * Link the DirectedEdges at the nodes of the graph.
 * This allows clients to link only a subset of nodes in the graph, for
 * efficiency (because they know that only a subset is of interest).
 */
void
PlanarGraph::linkAllDirectedEdges()
{
#if GEOS_DEBUG
	cerr<<"PlanarGraph::linkAllDirectedEdges called"<<endl;
#endif
	NodeMap::iterator nodeit=nodes->nodeMap.begin();
	for (;nodeit!=nodes->nodeMap.end();nodeit++)
	{
		Node *node=nodeit->second;
		assert(node);

		EdgeEndStar *ees=node->getEdges();
		assert(ees);

		// Unespected non-DirectedEdgeStar in node
		assert(dynamic_cast<DirectedEdgeStar *>(ees));
		DirectedEdgeStar *des=static_cast<DirectedEdgeStar *>(ees);

		des->linkAllDirectedEdges();
	}
}

/*public*/
EdgeEnd*
PlanarGraph::findEdgeEnd(Edge *e)
{
	vector<EdgeEnd*>* eev=getEdgeEnds();
	assert(eev);

	for (vector<EdgeEnd*>::iterator i=eev->begin(), iEnd=eev->end();
			i != iEnd;
			++i)
	{
		EdgeEnd *ee=*i;
		assert(ee);

		// should test using values rather then pointers ?
		if (ee->getEdge()==e) return ee;
	}
	return NULL;
}

/*public*/
Edge*
PlanarGraph::findEdge(const Coordinate& p0, const Coordinate& p1)
{
	for (size_t i=0, n=edges->size(); i<n; ++i)
	{
		Edge *e=(*edges)[i];
		assert(e);

		const CoordinateSequence* eCoord=e->getCoordinates();
		assert(eCoord);

		if (p0==eCoord->getAt(0) && p1==eCoord->getAt(1))
			return e;
	}
	return NULL;
}

/*public*/
Edge*
PlanarGraph::findEdgeInSameDirection(const Coordinate& p0,
		const Coordinate& p1)
{
	for(size_t i=0, n=edges->size(); i<n; i++)
	{
		Edge *e=(*edges)[i];
		assert(e);

		const CoordinateSequence* eCoord=e->getCoordinates();
		assert(eCoord);

		size_t nCoords=eCoord->size();
		assert(nCoords>1);

		if (matchInSameDirection(p0, p1,
				eCoord->getAt(0),
				eCoord->getAt(1)))
		{
			return e;
		}

		if (matchInSameDirection(p0, p1,
				eCoord->getAt(nCoords-1),
				eCoord->getAt(nCoords-2)))
		{
			return e;
		}
	}

	return NULL;
}

/*private*/
bool
PlanarGraph::matchInSameDirection(const Coordinate& p0, const Coordinate& p1,
		const Coordinate& ep0, const Coordinate& ep1)
{
	if (!(p0==ep0))
		return false;

	if (CGAlgorithms::computeOrientation(p0,p1,ep1)==CGAlgorithms::COLLINEAR
		&& Quadrant::quadrant(p0,p1)==Quadrant::quadrant(ep0,ep1))
			return true;
	return false;
}

string
PlanarGraph::printEdges()
{
	
    std::ostringstream oss;
    oss << "Edges: ";
	for(size_t i=0, n=edges->size(); i<n; ++i)
	{
        Edge *e=(*edges)[i];
		oss << "edge " << i << ":\n" << e->print() << e->eiList.print();
	}
    return oss.str();
}

NodeMap*
PlanarGraph::getNodeMap()
{
	return nodes;
}

} // namespace geos.geomgraph
} // namespace geos

/**********************************************************************
 * $Log$
 * Revision 1.33  2006/06/12 11:29:23  strk
 * unsigned int => size_t
 *
 * Revision 1.32  2006/05/04 12:17:27  strk
 * Added comment about management of newly created DirectedEdges
 *
 * Revision 1.31  2006/04/27 15:07:15  strk
 * use output operators in debugging lines
 *
 * Revision 1.30  2006/04/07 16:52:20  strk
 * Port info, doxygen comments, assertion checking
 *
 * Revision 1.29  2006/03/15 17:16:29  strk
 * streamlined headers inclusion
 *
 * Revision 1.28  2006/03/14 15:46:53  strk
 * Added PlanarGraph::getNodes(vector&) func, to reduce useless heap allocations
 *
 * Revision 1.27  2006/03/09 16:46:47  strk
 * geos::geom namespace definition, first pass at headers split
 *
 * Revision 1.26  2006/03/06 19:40:46  strk
 * geos::util namespace. New GeometryCollection::iterator interface, many cleanups.
 *
 * Revision 1.25  2006/03/03 10:46:21  strk
 * Removed 'using namespace' from headers, added missing headers in .cpp files, removed useless includes in headers (bug#46)
 *
 * Revision 1.24  2006/03/02 16:21:26  strk
 * geos::precision namespace added
 *
 * Revision 1.23  2006/03/02 12:12:00  strk
 * Renamed DEBUG macros to GEOS_DEBUG, all wrapped in #ifndef block to allow global override (bug#43)
 *
 **********************************************************************/

