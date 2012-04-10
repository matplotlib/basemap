/**********************************************************************
 * $Id: Node.cpp 3404 2011-07-05 09:43:51Z strk $
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2006 Refractions Research Inc.
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation. 
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: index/quadtree/Node.java rev 1.8 (JTS-1.10)
 *
 **********************************************************************/

#include <geos/index/quadtree/Node.h> 
#include <geos/index/quadtree/Key.h> 
#include <geos/geom/Envelope.h>

#include <string>
#include <sstream>
#include <cassert>

#ifndef GEOS_DEBUG
#define GEOS_DEBUG 0
#endif

#if GEOS_DEBUG
#include <iostream>
#endif

using namespace std;
using namespace geos::geom;

namespace geos {
namespace index { // geos.index
namespace quadtree { // geos.index.quadtree

/* public static */
std::auto_ptr<Node>
Node::createNode(const Envelope& env)
{
	Key key(env);

	std::auto_ptr<Envelope> nenv ( new Envelope(key.getEnvelope()) );
	std::auto_ptr<Node> node (
		new Node(nenv, key.getLevel())
	);
	return node;
}

/* static public */
std::auto_ptr<Node>
Node::createExpanded(std::auto_ptr<Node> node, const Envelope& addEnv)
{
	Envelope expandEnv(addEnv);
	if ( node.get() ) // should this be asserted ?
	{
		expandEnv.expandToInclude(node->getEnvelope());
	}

#if GEOS_DEBUG
	cerr<<"Node::createExpanded computed "<<expandEnv.toString()<<endl;
#endif

	std::auto_ptr<Node> largerNode = createNode(expandEnv);
	if ( node.get() ) // should this be asserted ?
	{
		largerNode->insertNode(node);
	}

	return largerNode;
}

/*public*/
Node*
Node::getNode(const Envelope *searchEnv)
{
	int subnodeIndex = getSubnodeIndex(searchEnv, centre);
	// if subquadIndex is -1 searchEnv is not contained in a subquad
	if (subnodeIndex != -1)
	{
		// create the quad if it does not exist
		Node *node = getSubnode(subnodeIndex);
		// recursively search the found/created quad
		return node->getNode(searchEnv);
	}
	else
	{
		return this;
	}
}

/*public*/
NodeBase*
Node::find(const Envelope *searchEnv)
{
	int subnodeIndex=getSubnodeIndex(searchEnv, centre);
	if (subnodeIndex==-1)
		return this;
	if (subnode[subnodeIndex]!=NULL) {
		// query lies in subquad, so search it
		Node *node=subnode[subnodeIndex];
		return node->find(searchEnv);
	}
	// no existing subquad, so return this one anyway
	return this;
}

void
Node::insertNode(std::auto_ptr<Node> node)
{
	assert( env->contains(node->getEnvelope()) );

	int index = getSubnodeIndex(node->getEnvelope(), centre);
	assert(index >= 0);

	if (node->level == level-1)
	{
		// We take ownership of node 
		delete subnode[index];
		subnode[index] = node.release();

		//System.out.println("inserted");
	}
	else
	{
		// the quad is not a direct child, so make a new child
		// quad to contain it and recursively insert the quad
		std::auto_ptr<Node> childNode ( createSubnode(index) );

		// childNode takes ownership of node
		childNode->insertNode(node);

		// We take ownership of childNode 
		delete subnode[index];
		subnode[index] = childNode.release();
	}
}

Node*
Node::getSubnode(int index)
{
	assert(index >=0 && index < 4);
	if (subnode[index] == NULL)
	{
		subnode[index] = createSubnode(index).release();
	}
	return subnode[index];
}

std::auto_ptr<Node>
Node::createSubnode(int index)
{
	// create a new subquad in the appropriate quadrant
	double minx=0.0;
	double maxx=0.0;
	double miny=0.0;
	double maxy=0.0;

	switch (index) {
		case 0:
			minx=env->getMinX();
			maxx=centre.x;
			miny=env->getMinY();
			maxy=centre.y;
			break;
		case 1:
			minx=centre.x;
			maxx=env->getMaxX();
			miny=env->getMinY();
			maxy=centre.y;
			break;
		case 2:
			minx=env->getMinX();
			maxx=centre.x;
			miny=centre.y;
			maxy=env->getMaxY();
			break;
		case 3:
			minx=centre.x;
			maxx=env->getMaxX();
			miny=centre.y;
			maxy=env->getMaxY();
			break;
	}
	std::auto_ptr<Envelope> sqEnv ( new Envelope(minx,maxx,miny,maxy) );
	std::auto_ptr<Node> node ( new Node(sqEnv, level-1) );
	return node;
}

string
Node::toString() const
{
	ostringstream os;
	os <<"L"<<level<<" "<<env->toString()<<" Ctr["<<centre.toString()<<"]";
	os <<" "+NodeBase::toString();
	return os.str();
}


} // namespace geos.index.quadtree
} // namespace geos.index
} // namespace geos

/**********************************************************************
 * $Log$
 * Revision 1.2  2006/03/23 13:31:58  strk
 * Fixed to allow build with GEOS_DEBUG
 *
 * Revision 1.1  2006/03/22 14:28:53  strk
 * Filenames renamed to match class names (matching JTS)
 *
 * Revision 1.17  2006/03/22 12:22:50  strk
 * indexQuadtree.h split
 *
 **********************************************************************/

