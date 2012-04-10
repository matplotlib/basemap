/**********************************************************************
 * $Id: Key.cpp 2407 2009-04-25 00:23:06Z strk $
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2009  Sandro Santilli <strk@keybit.net>
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
 * Last port: index/quadtree/Key.java rev 1.8 (JTS-1.10)
 *
 **********************************************************************/

#include <geos/index/quadtree/Key.h>
#include <geos/index/quadtree/DoubleBits.h>
#include <geos/geom/Envelope.h>
#include <geos/geom/Coordinate.h>

#include <cmath>

#ifndef GEOS_DEBUG
#define GEOS_DEBUG 0
#endif

#ifdef GEOS_DEBUG
#include <iostream>
#endif

using namespace geos::geom;

namespace geos {
namespace index { // geos.index
namespace quadtree { // geos.index.quadtree

/* static public */
int
Key::computeQuadLevel(const Envelope& env)
{
	double dx = env.getWidth();
	double dy = env.getHeight();
	double dMax = dx > dy ? dx : dy;
	int level=DoubleBits::exponent(dMax)+1;
#if GEOS_DEBUG
	std::cerr<<"Maxdelta:"<<dMax<<" exponent:"<<(level-1)<<std::endl;
#endif
	return level;
}

Key::Key(const Envelope& itemEnv)
	:
	pt(),
	level(0),
	env()
{
	computeKey(itemEnv);
}

Key::~Key()
{
}

const Coordinate&
Key::getPoint() const
{
	return pt;
}

int
Key::getLevel() const
{
	return level;
}

const Envelope&
Key::getEnvelope() const
{
	return env;
}

Coordinate*
Key::getCentre() const
{
	return new Coordinate(
			( env.getMinX() + env.getMaxX() ) / 2,
			( env.getMinY() + env.getMaxY() ) / 2
		);
}

/*public*/
void
Key::computeKey(const Envelope& itemEnv)
{
	level=computeQuadLevel(itemEnv);
	env.init(); // reset to null 
	computeKey(level, itemEnv);
	// MD - would be nice to have a non-iterative form of this algorithm
	while (!env.contains(itemEnv)) {
		level+=1;
		computeKey(level, itemEnv);
	}
#if GEOS_DEBUG
	std::cerr<<"Key::computeKey:"<<std::endl;
	std::cerr<<" itemEnv: "<<itemEnv.toString()<<std::endl;
	std::cerr<<"  keyEnv: "<<env.toString()<<std::endl;
	std::cerr<<"  keyLvl: "<<level<<std::endl;

#endif
}

void
Key::computeKey(int level, const Envelope& itemEnv)
{
	double quadSize=DoubleBits::powerOf2(level);
	//double quadSize=pow2.power(level);
	pt.x = std::floor(itemEnv.getMinX()/quadSize)*quadSize;
	pt.y = std::floor(itemEnv.getMinY()/quadSize)*quadSize;
	env.init(pt.x, pt.x + quadSize, pt.y, pt.y + quadSize);
}

} // namespace geos.index.quadtree
} // namespace geos.index
} // namespace geos

/**********************************************************************
 * $Log$
 * Revision 1.1  2006/03/22 14:28:52  strk
 * Filenames renamed to match class names (matching JTS)
 *
 * Revision 1.13  2006/03/22 12:22:50  strk
 * indexQuadtree.h split
 *
 **********************************************************************/
