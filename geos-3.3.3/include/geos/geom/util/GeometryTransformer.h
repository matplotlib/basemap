/**********************************************************************
 * $Id: GeometryTransformer.h 3179 2011-02-03 19:59:23Z strk $
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2011 Sandro Santilli <strk@keybit.net>
 * Copyright (C) 2006 Refractions Research Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation. 
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: geom/util/GeometryTransformer.java r320 (JTS-1.12)
 *
 **********************************************************************/

#ifndef GEOS_GEOM_UTIL_GEOMETRYTRANSFORMER_H
#define GEOS_GEOM_UTIL_GEOMETRYTRANSFORMER_H


#include <geos/export.h>
#include <geos/geom/Coordinate.h> // destructor visibility for vector
#include <geos/geom/Geometry.h> // destructor visibility for auto_ptr
#include <geos/geom/CoordinateSequence.h> // destructor visibility for auto_ptr

#include <memory> // for auto_ptr
#include <vector>

// Forward declarations
namespace geos {
	namespace geom {
		class Geometry;
		class GeometryFactory;
		class Point;
		class LinearRing;
		class LineString;
		class Polygon;
		class MultiPoint;
		class MultiPolygon;
		class MultiLineString;
		class GeometryCollection;
		namespace util {
			//class GeometryEditorOperation;
		}
	}
}


namespace geos {
namespace geom { // geos.geom
namespace util { // geos.geom.util

/** \brief
 * A framework for processes which transform an input {@link Geometry} into
 * an output {@link Geometry}, possibly changing its structure and type(s).
 *
 * This class is a framework for implementing subclasses
 * which perform transformations on
 * various different Geometry subclasses.
 * It provides an easy way of applying specific transformations
 * to given geometry types, while allowing unhandled types to be simply copied.
 * Also, the framework ensures that if subcomponents change type
 * the parent geometries types change appropriately to maintain valid structure.
 * Subclasses will override whichever <code>transformX</code> methods
 * they need to to handle particular Geometry types.
 * 
 * A typically usage would be a transformation that may transform Polygons into
 * Polygons, LineStrings
 * or Points.  This class would likely need to override the
 * {@link transformMultiPolygon} method to ensure that if input Polygons
 * change type the result is a GeometryCollection,
 * not a MultiPolygon
 * 
 * The default behaviour of this class is to simply recursively transform
 * each Geometry component into an identical object by copying.
 *
 * Note that all <code>transformX</code> methods may return <code>null</code>,
 * to avoid creating empty geometry objects. This will be handled correctly
 * by the transformer.
 * The {@link transform} method itself will always
 * return a geometry object.
 *
 * @see GeometryEditor
 *
 * Possible extensions:
 * getParent() method to return immediate parent e.g. of LinearRings in Polygons
 *
 */
class GEOS_DLL GeometryTransformer {

public:

	GeometryTransformer();

	virtual ~GeometryTransformer();

	std::auto_ptr<Geometry> transform(const Geometry* nInputGeom);

protected:

	const GeometryFactory* factory;

	/**
	 * Convenience method which provides standard way of
	 * creating a {@link CoordinateSequence}
	 *
	 * @param coords the coordinate array to copy
	 * @return a coordinate sequence for the array
	 *
	 * [final]
	 */
	CoordinateSequence::AutoPtr createCoordinateSequence(
			std::auto_ptr< std::vector<Coordinate> > coords);

	virtual CoordinateSequence::AutoPtr transformCoordinates(
			const CoordinateSequence* coords,
			const Geometry* parent);

	virtual Geometry::AutoPtr transformPoint(
			const Point* geom,
			const Geometry* parent);

	virtual Geometry::AutoPtr transformMultiPoint(
			const MultiPoint* geom,
			const Geometry* parent);

	virtual Geometry::AutoPtr transformLinearRing(
			const LinearRing* geom,
			const Geometry* parent);

	virtual Geometry::AutoPtr transformLineString(
			const LineString* geom,
			const Geometry* parent);

	virtual Geometry::AutoPtr transformMultiLineString(
			const MultiLineString* geom,
			const Geometry* parent);

	virtual Geometry::AutoPtr transformPolygon(
			const Polygon* geom,
			const Geometry* parent);

	virtual Geometry::AutoPtr transformMultiPolygon(
			const MultiPolygon* geom,
			const Geometry* parent);

	virtual Geometry::AutoPtr transformGeometryCollection(
			const GeometryCollection* geom,
			const Geometry* parent);

private:

	const Geometry* inputGeom;

	// these could eventually be exposed to clients
	/**
	 * <code>true</code> if empty geometries should not be included in the result
	 */
	bool pruneEmptyGeometry;

	/**
	 * <code>true</code> if a homogenous collection result
	 * from a {@link GeometryCollection} should still
	 * be a general GeometryCollection
	 */
	bool preserveGeometryCollectionType;

	/**
	 * <code>true</code> if the output from a collection argument should still be a collection
	 */
	bool preserveCollections;

	/**
	 * <code>true</code> if the type of the input should be preserved
	 */
	bool preserveType;

    // Declare type as noncopyable
    GeometryTransformer(const GeometryTransformer& other);
    GeometryTransformer& operator=(const GeometryTransformer& rhs);
};


} // namespace geos.geom.util
} // namespace geos.geom
} // namespace geos

//#ifdef GEOS_INLINE
//# include "geos/geom/util/GeometryTransformer.inl"
//#endif

#endif // GEOS_GEOM_UTIL_GEOMETRYTRANSFORMER_H

/**********************************************************************
 * $Log$
 * Revision 1.4  2006/06/19 21:20:22  strk
 * updated port info
 *
 * Revision 1.3  2006/04/13 14:25:17  strk
 * TopologyPreservingSimplifier initial port
 *
 * Revision 1.2  2006/04/11 12:56:06  strk
 * used typedef for auto_ptr<CoordinateSequence>
 *
 * Revision 1.1  2006/04/11 12:21:49  strk
 * GeometryTransformer class ported
 *
 **********************************************************************/
