/**********************************************************************
 * $Id: Edge.cpp 1820 2006-09-06 16:54:23Z mloskot $
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2006 Refractions Research Inc.
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Licence as published
 * by the Free Software Foundation. 
 * See the COPYING file for more information.
 *
 **********************************************************************/

#include <geos/planargraph/Edge.h>
#include <geos/planargraph/DirectedEdge.h>
#include <geos/planargraph/Node.h>

#include <iostream>

namespace geos {
namespace planargraph {

/*
 * Initializes this Edge's two DirectedEdges, and for each DirectedEdge:
 *  sets the Edge,
 *  sets the symmetric DirectedEdge, and
 *  adds this Edge to its from-Node.
 */
void
Edge::setDirectedEdges(DirectedEdge *de0, DirectedEdge *de1)
{
	dirEdge.push_back(de0);
	dirEdge.push_back(de1);
	de0->setEdge(this);
	de1->setEdge(this);
	de0->setSym(de1);
	de1->setSym(de0);
	de0->getFromNode()->addOutEdge(de0);
	de1->getFromNode()->addOutEdge(de1);
}

/**
* Returns one of the DirectedEdges associated with this Edge.
* @param i 0 or 1
*/
DirectedEdge *
Edge::getDirEdge(int i)
{
	return dirEdge[i];
}

/*
 * Returns the DirectedEdge that starts from the given node, or null if the
 * node is not one of the two nodes associated with this Edge.
 */
DirectedEdge *
Edge::getDirEdge(Node *fromNode)
{
	if (dirEdge[0]->getFromNode()==fromNode) return dirEdge[0];
	if (dirEdge[1]->getFromNode()==fromNode) return dirEdge[1];
	// node not found
	// possibly should throw an exception here?
	return NULL;
}

/**
* If <code>node</code> is one of the two nodes associated with this Edge,
* returns the other node; otherwise returns null.
*/
Node*
Edge::getOppositeNode(Node *node)
{
	if (dirEdge[0]->getFromNode()==node) return dirEdge[0]->getToNode();
	if (dirEdge[1]->getFromNode()==node) return dirEdge[1]->getToNode();
	// node not found
	// possibly should throw an exception here?
	return NULL;
}

std::ostream& operator<<(std::ostream& os, const Edge& n) {
	os << "Edge ";
	if ( n.isMarked() ) os << " Marked ";
	if ( n.isVisited() ) os << " Visited ";
	return os;
}

} // namespace planargraph
} // namespace geos

/**********************************************************************
 * $Log$
 * Revision 1.1  2006/03/21 21:42:54  strk
 * planargraph.h header split, planargraph:: classes renamed to match JTS symbols
 *
 **********************************************************************/

