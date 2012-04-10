/**********************************************************************
 * $Id: BufferBuilder.h 3301 2011-04-27 09:42:31Z strk $
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2009-2011 Sandro Santilli <strk@keybit.net>
 * Copyright (C) 2008-2010 Safe Software Inc.
 * Copyright (C) 2006-2007 Refractions Research Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation. 
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: operation/buffer/BufferBuilder.java r378 (JTS-1.12)
 *
 **********************************************************************/

#ifndef GEOS_OP_BUFFER_BUFFERBUILDER_H
#define GEOS_OP_BUFFER_BUFFERBUILDER_H

#include <geos/export.h>

#include <vector>

#include <geos/operation/buffer/BufferOp.h> // for inlines (BufferOp enums)
#include <geos/operation/buffer/OffsetCurveBuilder.h> // for inline (OffsetCurveBuilder enums)
#include <geos/geomgraph/EdgeList.h> // for composition

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4251) // warning C4251: needs to have dll-interface to be used by clients of class
#endif

// Forward declarations
namespace geos {
	namespace geom {
		class PrecisionModel;
		class Geometry;
		class GeometryFactory;
	}
	namespace algorithm {
		class CGAlgorithms;
		class LineIntersector;
	}
	namespace noding {
		class Noder;
		class SegmentString;
		class IntersectionAdder;
	}
	namespace geomgraph {
		class Edge;
		class Label;
		class PlanarGraph;
	}
	namespace operation {
		namespace buffer {
			class BufferSubgraph;
		}
		namespace overlay {
			class PolygonBuilder;
		}
	}
}

namespace geos {
namespace operation { // geos.operation
namespace buffer { // geos.operation.buffer

/**
 *
 * \brief
 * Builds the buffer geometry for a given input geometry and precision model.
 *
 * Allows setting the level of approximation for circular arcs,
 * and the precision model in which to carry out the computation.
 * 
 * When computing buffers in floating point double-precision
 * it can happen that the process of iterated noding can fail to converge
 * (terminate).
 *
 * In this case a TopologyException will be thrown.
 * Retrying the computation in a fixed precision
 * can produce more robust results.
 *
 */
class GEOS_DLL BufferBuilder {

public:
	/**
	 * Creates a new BufferBuilder
	 *
	 * @param nBufParams buffer parameters, this object will
	 *                   keep a reference to the passed parameters
	 *                   so caller must make sure the object is
	 *                   kept alive for the whole lifetime of
	 *                   the buffer builder.
	 */
	BufferBuilder(const BufferParameters& nBufParams)
		:
		bufParams(nBufParams),
		workingPrecisionModel(NULL),
		li(NULL),
		intersectionAdder(NULL),
		workingNoder(NULL),
		geomFact(NULL),
		edgeList()
	{}

	~BufferBuilder();


	/**
	 * Sets the precision model to use during the curve computation
	 * and noding,
	 * if it is different to the precision model of the Geometry.
	 * If the precision model is less than the precision of the
	 * Geometry precision model,
	 * the Geometry must have previously been rounded to that precision.
	 *
	 * @param pm the precision model to use
	 */
	void setWorkingPrecisionModel(const geom::PrecisionModel *pm) {
		workingPrecisionModel=pm;
	}

	/**
	 * Sets the {@link noding::Noder} to use during noding.
	 * This allows choosing fast but non-robust noding, or slower
	 * but robust noding.
	 *
	 * @param noder the noder to use
	 */
	void setNoder(noding::Noder* newNoder) { workingNoder = newNoder; }

	geom::Geometry* buffer(const geom::Geometry *g, double distance);
		// throw (GEOSException);

	/// Not in JTS: this is a GEOS extension
	geom::Geometry* bufferLineSingleSided( const geom::Geometry* g,
	                                double distance, bool leftSide ) ;
		// throw (GEOSException);

private:
	/**
	 * Compute the change in depth as an edge is crossed from R to L
	 */
	static int depthDelta(geomgraph::Label *label);

	const BufferParameters& bufParams; 

	const geom::PrecisionModel* workingPrecisionModel;

	algorithm::LineIntersector* li;

	noding::IntersectionAdder* intersectionAdder;

	noding::Noder* workingNoder;

	const geom::GeometryFactory* geomFact;

	geomgraph::EdgeList edgeList;

	std::vector<geomgraph::Label *> newLabels;

	void computeNodedEdges(std::vector<noding::SegmentString*>& bufSegStr,
			const geom::PrecisionModel *precisionModel);
			// throw(GEOSException);

	/**
	 * Inserted edges are checked to see if an identical edge already
	 * exists.
	 * If so, the edge is not inserted, but its label is merged
	 * with the existing edge.
	 *
	 * The function takes responsability of releasing the Edge parameter
	 * memory when appropriate.
	 */
	void insertUniqueEdge(geomgraph::Edge *e);

	void createSubgraphs(geomgraph::PlanarGraph *graph,
			std::vector<BufferSubgraph*>& list);

	/**
	 * Completes the building of the input subgraphs by
	 * depth-labelling them,
	 * and adds them to the PolygonBuilder.
	 * The subgraph list must be sorted in rightmost-coordinate order.
	 *
	 * @param subgraphList the subgraphs to build
	 * @param polyBuilder the PolygonBuilder which will build
	 *        the final polygons
	 */
	void buildSubgraphs(const std::vector<BufferSubgraph*>& subgraphList,
			overlay::PolygonBuilder& polyBuilder);

	/// \brief
	/// Return the externally-set noding::Noder OR a newly created
	/// one using the given precisionModel.
	//
	/// NOTE: if an externally-set noding::Noder is available no
	/// check is performed to ensure it will use the
	/// given PrecisionModel
	///
	noding::Noder* getNoder(const geom::PrecisionModel* precisionModel);


	/**
	 * Gets the standard result for an empty buffer.
	 * Since buffer always returns a polygonal result,
	 * this is chosen to be an empty polygon.
	 *
	 * @return the empty result geometry, transferring ownership to caller.
	 */
	geom::Geometry* createEmptyResultGeometry() const;
    
    // Declare type as noncopyable
    BufferBuilder(const BufferBuilder& other);
    BufferBuilder& operator=(const BufferBuilder& rhs);
};

} // namespace geos::operation::buffer
} // namespace geos::operation
} // namespace geos

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // ndef GEOS_OP_BUFFER_BUFFERBUILDER_H

/**********************************************************************
 * $Log$
 * Revision 1.2  2006/03/14 16:08:21  strk
 * changed buildSubgraphs signature to use refs rather then pointers, made it const-correct. Reduced heap allocations in createSubgraphs()
 *
 * Revision 1.1  2006/03/14 00:19:40  strk
 * opBuffer.h split, streamlined headers in some (not all) files in operation/buffer/
 *
 **********************************************************************/

