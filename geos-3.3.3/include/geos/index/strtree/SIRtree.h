/**********************************************************************
 * $Id: SIRtree.h 2961 2010-03-29 12:17:37Z mloskot $
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2006 Refractions Research Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation. 
 * See the COPYING file for more information.
 *
 **********************************************************************/

#ifndef GEOS_INDEX_STRTREE_SIRTREE_H
#define GEOS_INDEX_STRTREE_SIRTREE_H

#include <geos/export.h>

#include <geos/index/strtree/AbstractSTRtree.h> // for inheritance
#include <geos/index/strtree/Interval.h> // for inline

#include <vector>
#include <memory>

namespace geos {
namespace index { // geos::index
namespace strtree { // geos::index::strtree

/** \brief
 * One-dimensional version of an STR-packed R-tree.
 *
 * SIR stands for "Sort-Interval-Recursive".
 *
 * STR-packed R-trees are described in:
 * P. Rigaux, Michel Scholl and Agnes Voisard. Spatial Databases With
 * Application To GIS. Morgan Kaufmann, San Francisco, 2002.
 *
 * @see STRtree
 */
class GEOS_DLL SIRtree: public AbstractSTRtree {
using AbstractSTRtree::insert;
using AbstractSTRtree::query;

public:

	/** \brief
	 * Constructs an SIRtree with the default node capacity.
	 */
	SIRtree();

	/** \brief
	 * Constructs an SIRtree with the given maximum number of child nodes
	 * that a node may have
	 */
	SIRtree(std::size_t nodeCapacity);

	virtual ~SIRtree();

	void insert(double x1, double x2, void* item);

	/**
	 * Returns items whose bounds intersect the given bounds.
	 * @param x1 possibly equal to x2
	 */
	std::vector<void*>* query(double x1, double x2)
	{
		std::vector<void*>* results = new std::vector<void*>();
		Interval interval(std::min(x1, x2), std::max(x1, x2));
		AbstractSTRtree::query(&interval, *results);
		return results;
	}

	/**
	 * Returns items whose bounds intersect the given value.
	 */
	std::vector<void*>* query(double x) { return query(x,x); }


protected:

	class SIRIntersectsOp:public AbstractSTRtree::IntersectsOp {
	public:
		bool intersects(const void* aBounds, const void* bBounds);
	};

	/** \brief
	 * Sorts the childBoundables then divides them into groups of size M, where
	 * M is the node capacity.
	 */
	std::auto_ptr<BoundableList> createParentBoundables(
			BoundableList* childBoundables, int newLevel);

	AbstractNode* createNode(int level);

	IntersectsOp* getIntersectsOp() {return intersectsOp;};

	std::auto_ptr<BoundableList> sortBoundables(const BoundableList* input);

private:

	IntersectsOp* intersectsOp;
};
	

} // namespace geos::index::strtree
} // namespace geos::index
} // namespace geos

#endif // GEOS_INDEX_STRTREE_SIRTREE_H

/**********************************************************************
 * $Log$
 * Revision 1.2  2006/06/12 10:49:43  strk
 * unsigned int => size_t
 *
 * Revision 1.1  2006/03/21 10:47:34  strk
 * indexStrtree.h split
 *
 **********************************************************************/

