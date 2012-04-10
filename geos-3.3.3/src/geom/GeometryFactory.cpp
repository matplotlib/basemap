/**********************************************************************
 * $Id: GeometryFactory.cpp 3181 2011-02-04 08:50:42Z strk $
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2011 Sandro Santilli <strk@keybit.net>
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 * Copyright (C) 2005 Refractions Research Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation. 
 * See the COPYING file for more information.
 *
 **********************************************************************
 *
 * Last port: geom/GeometryFactory.java r320 (JTS-1.12)
 *
 **********************************************************************/

#include <geos/geom/Coordinate.h>
#include <geos/geom/CoordinateArraySequenceFactory.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Point.h>
#include <geos/geom/LineString.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/MultiPoint.h>
#include <geos/geom/MultiLineString.h>
#include <geos/geom/MultiPolygon.h>
#include <geos/geom/GeometryCollection.h>
#include <geos/geom/PrecisionModel.h>
#include <geos/geom/Envelope.h>
#include <geos/util/IllegalArgumentException.h>

#include <cassert>
#include <vector>
#include <typeinfo>
#include <cmath>

#ifndef GEOS_DEBUG
#define GEOS_DEBUG 0
#endif

#ifdef GEOS_DEBUG
#include <iostream>
#endif

#ifndef GEOS_INLINE
# include <geos/geom/GeometryFactory.inl>
#endif

using namespace std;

namespace geos {
namespace geom { // geos::geom

//namespace { 
//	class gfCoordinateOperation: public CoordinateOperation {
//	using CoordinateOperation::edit;
//	public:
//		virtual CoordinateSequence* edit(const CoordinateSequence *coordinates,
//				const Geometry *geometry);
//	};
//}



/*public*/
GeometryFactory::GeometryFactory()
	:
	precisionModel(new PrecisionModel()),
	SRID(0),
	coordinateListFactory(CoordinateArraySequenceFactory::instance())
{
#if GEOS_DEBUG
	std::cerr << "GEOS_DEBUG: GeometryFactory["<<this<<"]::GeometryFactory()" << std::endl;
	std::cerr << "\tcreate PrecisionModel["<<precisionModel<<"]" << std::endl;
#endif
}

/*public*/
GeometryFactory::GeometryFactory(const PrecisionModel* pm, int newSRID,
		CoordinateSequenceFactory* nCoordinateSequenceFactory)
	:
	SRID(newSRID)
{
#if GEOS_DEBUG
	std::cerr << "GEOS_DEBUG: GeometryFactory["<<this<<"]::GeometryFactory(PrecisionModel["<<pm<<"], SRID)" << std::endl;
#endif
	if ( ! pm ) {
		precisionModel=new PrecisionModel();
	} else {
		precisionModel=new PrecisionModel(*pm);
	}

	if ( ! nCoordinateSequenceFactory ) {
		coordinateListFactory=CoordinateArraySequenceFactory::instance();
	} else {
		coordinateListFactory=nCoordinateSequenceFactory;
	}
}

/*public*/
GeometryFactory::GeometryFactory(
		CoordinateSequenceFactory* nCoordinateSequenceFactory)
	:
	precisionModel(new PrecisionModel()),
	SRID(0)
{
#if GEOS_DEBUG
	std::cerr << "GEOS_DEBUG: GeometryFactory["<<this<<"]::GeometryFactory(CoordinateSequenceFactory["<<nCoordinateSequenceFactory<<"])" << std::endl;
#endif
	if ( ! nCoordinateSequenceFactory ) {
		coordinateListFactory=CoordinateArraySequenceFactory::instance();
	} else {
		coordinateListFactory=nCoordinateSequenceFactory;
	}
}

/*public*/
GeometryFactory::GeometryFactory(const PrecisionModel *pm)
	:
	SRID(0),
	coordinateListFactory(CoordinateArraySequenceFactory::instance())
{
#if GEOS_DEBUG
	std::cerr << "GEOS_DEBUG: GeometryFactory["<<this<<"]::GeometryFactory(PrecisionModel["<<pm<<"])" << std::endl;
#endif
	if ( ! pm ) {
		precisionModel=new PrecisionModel();
	} else {
		precisionModel=new PrecisionModel(*pm);
	}
}

/*public*/
GeometryFactory::GeometryFactory(const PrecisionModel* pm, int newSRID)
	:
	SRID(newSRID),
	coordinateListFactory(CoordinateArraySequenceFactory::instance())
{
#if GEOS_DEBUG
	std::cerr << "GEOS_DEBUG: GeometryFactory["<<this<<"]::GeometryFactory(PrecisionModel["<<pm<<"], SRID)" << std::endl;
#endif
	if ( ! pm ) {
		precisionModel=new PrecisionModel();
	} else {
		precisionModel=new PrecisionModel(*pm);
	}
}

/*public*/
GeometryFactory::GeometryFactory(const GeometryFactory &gf)
{
	assert(gf.precisionModel);
	precisionModel=new PrecisionModel(*(gf.precisionModel));
	SRID=gf.SRID;
	coordinateListFactory=gf.coordinateListFactory;
}

/*public*/
GeometryFactory::~GeometryFactory(){
#if GEOS_DEBUG
	std::cerr << "GEOS_DEBUG: GeometryFactory["<<this<<"]::~GeometryFactory()" << std::endl;
#endif
	delete precisionModel;
}
  
/*public*/
Point*
GeometryFactory::createPointFromInternalCoord(const Coordinate* coord,
		const Geometry *exemplar) const
{
	assert(coord);
	Coordinate newcoord = *coord;
	exemplar->getPrecisionModel()->makePrecise(&newcoord);
	return exemplar->getFactory()->createPoint(newcoord);
}


/*public*/
Geometry*
GeometryFactory::toGeometry(const Envelope* envelope) const
{
	Coordinate coord;

	if (envelope->isNull()) {
		return createPoint();
	}
	if (envelope->getMinX()==envelope->getMaxX() && envelope->getMinY()==envelope->getMaxY()) {
		coord.x = envelope->getMinX();
		coord.y = envelope->getMinY();
		return createPoint(coord);
	}
	CoordinateSequence *cl=CoordinateArraySequenceFactory::instance()->
        create((size_t) 0, 2);
	coord.x = envelope->getMinX();
	coord.y = envelope->getMinY();
	cl->add(coord);
	coord.x = envelope->getMaxX();
	coord.y = envelope->getMinY();
	cl->add(coord);
	coord.x = envelope->getMaxX();
	coord.y = envelope->getMaxY();
	cl->add(coord);
	coord.x = envelope->getMinX();
	coord.y = envelope->getMaxY();
	cl->add(coord);
	coord.x = envelope->getMinX();
	coord.y = envelope->getMinY();
	cl->add(coord);

	Polygon *p = createPolygon(createLinearRing(cl), NULL);
	return p;
}

/*public*/
const PrecisionModel*
GeometryFactory::getPrecisionModel() const
{
	return precisionModel;
}

/*public*/
Point*
GeometryFactory::createPoint() const
{
	return new Point(NULL, this);
}

/*public*/
Point*
GeometryFactory::createPoint(const Coordinate& coordinate) const
{
	if (coordinate.isNull()) {
		return createPoint();
	} else {
		std::size_t dim = ISNAN(coordinate.z) ? 2 : 3;
		CoordinateSequence *cl = coordinateListFactory->create(new vector<Coordinate>(1, coordinate), dim);
		//cl->setAt(coordinate, 0);
		Point *ret = createPoint(cl);
		return ret;
	}
}

/*public*/
Point*
GeometryFactory::createPoint(CoordinateSequence *newCoords) const
{
	return new Point(newCoords,this);
}

/*public*/
Point*
GeometryFactory::createPoint(const CoordinateSequence &fromCoords) const
{
	CoordinateSequence *newCoords = fromCoords.clone();
	Point *g = NULL;
	try {
		g = new Point(newCoords,this); 
	} catch (...) {
		delete newCoords;
		throw;
	}
	return g;

}

/*public*/
MultiLineString*
GeometryFactory::createMultiLineString() const
{
	return new MultiLineString(NULL,this);
}

/*public*/
MultiLineString*
GeometryFactory::createMultiLineString(vector<Geometry *> *newLines)
	const
{
	return new MultiLineString(newLines,this);
}

/*public*/
MultiLineString*
GeometryFactory::createMultiLineString(const vector<Geometry *> &fromLines)
	const
{
	vector<Geometry *>*newGeoms = new vector<Geometry *>(fromLines.size());
	for (size_t i=0; i<fromLines.size(); i++)
	{
		const LineString *line = dynamic_cast<const LineString *>(fromLines[i]);
		if ( ! line ) throw geos::util::IllegalArgumentException("createMultiLineString called with a vector containing non-LineStrings");
		(*newGeoms)[i] = new LineString(*line);
	}
	MultiLineString *g = NULL;
	try {
		g = new MultiLineString(newGeoms,this);
	} catch (...) {
		for (size_t i=0; i<newGeoms->size(); i++) {
			delete (*newGeoms)[i];
		}
		delete newGeoms;
		throw;
	}
	return g;
}

/*public*/
GeometryCollection*
GeometryFactory::createGeometryCollection() const
{
	return new GeometryCollection(NULL,this);
}

/*public*/
Geometry*
GeometryFactory::createEmptyGeometry() const
{
	return new GeometryCollection(NULL,this);
}

/*public*/
GeometryCollection*
GeometryFactory::createGeometryCollection(vector<Geometry *> *newGeoms) const
{
	return new GeometryCollection(newGeoms,this);
}

/*public*/
GeometryCollection*
GeometryFactory::createGeometryCollection(const vector<Geometry *> &fromGeoms) const
{
	vector<Geometry *> *newGeoms = new vector<Geometry *>(fromGeoms.size());
	for (size_t i=0; i<fromGeoms.size(); i++) {
		(*newGeoms)[i] = fromGeoms[i]->clone();
	}
	GeometryCollection *g = NULL;
	try {
		g = new GeometryCollection(newGeoms,this);
	} catch (...) {
		for (size_t i=0; i<newGeoms->size(); i++) {
			delete (*newGeoms)[i];
		}
		delete newGeoms;
		throw;
	}
	return g;
}

/*public*/
MultiPolygon*
GeometryFactory::createMultiPolygon() const
{
	return new MultiPolygon(NULL,this);
}

/*public*/
MultiPolygon*
GeometryFactory::createMultiPolygon(vector<Geometry *> *newPolys) const
{
	return new MultiPolygon(newPolys,this);
}

/*public*/
MultiPolygon*
GeometryFactory::createMultiPolygon(const vector<Geometry *> &fromPolys) const
{
	vector<Geometry *>*newGeoms = new vector<Geometry *>(fromPolys.size());
	for (size_t i=0; i<fromPolys.size(); i++)
	{
		(*newGeoms)[i] = fromPolys[i]->clone();
	}
	MultiPolygon *g = NULL;
	try {
		g = new MultiPolygon(newGeoms,this);
	} catch (...) {
		for (size_t i=0; i<newGeoms->size(); i++) {
			delete (*newGeoms)[i];
		}
		delete newGeoms;
		throw;
	}
	return g;
}

/*public*/
LinearRing*
GeometryFactory::createLinearRing() const
{
	return new LinearRing(NULL,this);
}

/*public*/
LinearRing*
GeometryFactory::createLinearRing(CoordinateSequence* newCoords) const
{
	return new LinearRing(newCoords,this);
}

/*public*/
Geometry::AutoPtr
GeometryFactory::createLinearRing(CoordinateSequence::AutoPtr newCoords) const
{
	return Geometry::AutoPtr(new LinearRing(newCoords, this));
}

/*public*/
LinearRing*
GeometryFactory::createLinearRing(const CoordinateSequence& fromCoords) const
{
	CoordinateSequence *newCoords = fromCoords.clone();
	LinearRing *g = NULL;
	// construction failure will delete newCoords
	g = new LinearRing(newCoords, this);
	return g;
}

/*public*/
MultiPoint*
GeometryFactory::createMultiPoint(vector<Geometry *> *newPoints) const
{
	return new MultiPoint(newPoints,this);
}

/*public*/
MultiPoint*
GeometryFactory::createMultiPoint(const vector<Geometry *> &fromPoints) const
{
	vector<Geometry *>*newGeoms = new vector<Geometry *>(fromPoints.size());
	for (size_t i=0; i<fromPoints.size(); i++)
	{
		(*newGeoms)[i] = fromPoints[i]->clone();
	}

	MultiPoint *g = NULL;
	try {
		g = new MultiPoint(newGeoms,this);
	} catch (...) {
		for (size_t i=0; i<newGeoms->size(); i++) {
			delete (*newGeoms)[i];
		}
		delete newGeoms;
		throw;
	}
	return g;
}

/*public*/
MultiPoint*
GeometryFactory::createMultiPoint() const
{
	return new MultiPoint(NULL, this);
}

/*public*/
MultiPoint*
GeometryFactory::createMultiPoint(const CoordinateSequence &fromCoords) const
{
	size_t npts=fromCoords.getSize();
	vector<Geometry *> *pts=new vector<Geometry *>;
	pts->reserve(npts);
	for (size_t i=0; i<npts; ++i) {
		Point *pt=createPoint(fromCoords.getAt(i));
		pts->push_back(pt);
	}
	MultiPoint *mp = NULL;
	try {
		mp = createMultiPoint(pts);
	} catch (...) {
		for (size_t i=0; i<npts; ++i) delete (*pts)[i];
		delete pts;
		throw;
	}
	return mp;
}

/*public*/
MultiPoint*
GeometryFactory::createMultiPoint(const std::vector<Coordinate> &fromCoords) const
{
	size_t npts=fromCoords.size();
	vector<Geometry *> *pts=new vector<Geometry *>;
	pts->reserve(npts);
	for (size_t i=0; i<npts; ++i) {
		Point *pt=createPoint(fromCoords[i]);
		pts->push_back(pt);
	}
	MultiPoint *mp = NULL;
	try {
		mp = createMultiPoint(pts);
	} catch (...) {
		for (size_t i=0; i<npts; ++i) delete (*pts)[i];
		delete pts;
		throw;
	}
	return mp;
}

/*public*/
Polygon*
GeometryFactory::createPolygon() const
{
	return new Polygon(NULL, NULL, this);
}

/*public*/
Polygon*
GeometryFactory::createPolygon(LinearRing *shell, vector<Geometry *> *holes)
	const
{
	return new Polygon(shell, holes, this);
}

/*public*/
Polygon*
GeometryFactory::createPolygon(const LinearRing &shell, const vector<Geometry *> &holes)
	const
{
	LinearRing *newRing = dynamic_cast<LinearRing *>(shell.clone());
	vector<Geometry *>*newHoles = new vector<Geometry *>(holes.size());
	for (size_t i=0; i<holes.size(); i++)
	{
		(*newHoles)[i] = holes[i]->clone();
	}
	Polygon *g = NULL;
	try {
		g = new Polygon(newRing, newHoles, this);
	} catch (...) {
		delete newRing;
		for (size_t i=0; i<holes.size(); i++)
			delete (*newHoles)[i];
		delete newHoles;
		throw;
	}
	return g;
}

/*public*/
LineString *
GeometryFactory::createLineString() const
{
	return new LineString(NULL, this);
}

/*public*/
std::auto_ptr<LineString>
GeometryFactory::createLineString(const LineString& ls) const
{
	return std::auto_ptr<LineString>(new LineString(ls));
}

/*public*/
LineString*
GeometryFactory::createLineString(CoordinateSequence *newCoords)
	const
{
	return new LineString(newCoords, this);
}

/*public*/
Geometry::AutoPtr
GeometryFactory::createLineString(CoordinateSequence::AutoPtr newCoords)
	const
{
	return Geometry::AutoPtr(new LineString(newCoords, this));
}

/*public*/
LineString*
GeometryFactory::createLineString(const CoordinateSequence &fromCoords)
	const
{
	CoordinateSequence *newCoords = fromCoords.clone();
	LineString *g = NULL;
	// construction failure will delete newCoords
	g = new LineString(newCoords, this);
	return g;
}

/*public*/
Geometry*
GeometryFactory::buildGeometry(vector<Geometry *> *newGeoms) const
{
	string geomClass("NULL");
	bool isHeterogeneous=false;
	bool hasGeometryCollection=false;
    
	for (size_t i=0, n=newGeoms->size(); i<n; ++i)
	{
		Geometry* geom = (*newGeoms)[i];
		string partClass(typeid(*geom).name());
		if (geomClass=="NULL")
		{
			geomClass=partClass;
		}
		else if (geomClass!=partClass)
		{
			isHeterogeneous = true;
		}
		if ( dynamic_cast<GeometryCollection*>(geom) )
		{
			hasGeometryCollection=true;
		}
	}

	// for the empty geometry, return an empty GeometryCollection
	if (geomClass=="NULL")
	{
		// we do not need the vector anymore
		delete newGeoms;
		return createGeometryCollection();
	}
	if (isHeterogeneous || hasGeometryCollection)
	{
		return createGeometryCollection(newGeoms);
	}

	// At this point we know the collection is not hetereogenous.
	// Determine the type of the result from the first Geometry in the
	// list. This should always return a geometry, since otherwise
	// an empty collection would have already been returned
	Geometry *geom0=(*newGeoms)[0];
	bool isCollection=newGeoms->size()>1;
	if (isCollection)
	{
		if (typeid(*geom0)==typeid(Polygon)) {
			return createMultiPolygon(newGeoms);
		} else if (typeid(*geom0)==typeid(LineString)) {
			return createMultiLineString(newGeoms);
		} else if (typeid(*geom0)==typeid(LinearRing)) {
			return createMultiLineString(newGeoms);
		} else if (typeid(*geom0)==typeid(Point)) {
			return createMultiPoint(newGeoms);
		} else {
			return createGeometryCollection(newGeoms);
		}
	}

	// since this is not a collection we can delete vector
	delete newGeoms;
	return geom0;
}

/*public*/
Geometry*
GeometryFactory::buildGeometry(const vector<Geometry *> &fromGeoms) const
{
	string geomClass("NULL");
	bool isHeterogeneous=false;
	bool isCollection=fromGeoms.size()>1;
	size_t i;
    
	for (i=0; i<fromGeoms.size(); i++) {
		string partClass(typeid(*fromGeoms[i]).name());
		if (geomClass=="NULL") {
			geomClass=partClass;
		} else if (geomClass!=partClass) {
			isHeterogeneous = true;
		}
	}

	// for the empty geometry, return an empty GeometryCollection
	if (geomClass=="NULL") {
		return createGeometryCollection();
	}
	if (isHeterogeneous) {
		return createGeometryCollection(fromGeoms);
	}

	// At this point we know the collection is not hetereogenous.
	// Determine the type of the result from the first Geometry in the
	// list. This should always return a geometry, since otherwise
	// an empty collection would have already been returned
	Geometry *geom0=fromGeoms[0];
	if (isCollection) {
		if (typeid(*geom0)==typeid(Polygon)) {
			return createMultiPolygon(fromGeoms);
		} else if (typeid(*geom0)==typeid(LineString)) {
			return createMultiLineString(fromGeoms);
		} else if (typeid(*geom0)==typeid(LinearRing)) {
			return createMultiLineString(fromGeoms);
		} else if (typeid(*geom0)==typeid(Point)) {
			return createMultiPoint(fromGeoms);
		}
		assert(0); // buildGeomtry encountered an unkwnon geometry type
	}

	return geom0->clone();
}

/*public*/
Geometry*
GeometryFactory::createGeometry(const Geometry *g) const
{
	// could this be cached to make this more efficient? Or maybe it isn't enough overhead to bother
	return g->clone();
	//GeometryEditor *editor=new GeometryEditor(this);
	//gfCoordinateOperation *coordOp = new gfCoordinateOperation();
	//Geometry *ret=editor->edit(g, coordOp);
	//delete coordOp;
	//delete editor;
	//return ret;
}

/*public*/
void
GeometryFactory::destroyGeometry(Geometry *g) const
{
	delete g;
}

/*public static*/
const GeometryFactory*
GeometryFactory::getDefaultInstance() 
{
	static GeometryFactory defInstance;
	return &defInstance;
}

} // namespace geos::geom
} // namespace geos

/**********************************************************************
 * $Log$
 * Revision 1.71  2006/07/08 00:33:54  strk
 *         * configure.in: incremented CAPI minor version, to avoid                        falling behind any future version from the 2.2. branch.
 *         * source/geom/Geometry.cpp, source/geom/GeometryFactory.cpp,
 *         source/geomgraph/EdgeRing.cpp,
 *         source/headers/geos/geom/Geometry.h,
 *         source/headers/geos/geom/GeometryFactory.h,
 *         source/headers/geos/geom/GeometryFactory.inl,
 *         source/headers/geos/geomgraph/EdgeRing.h:
 *         updated doxygen comments (sync with JTS head).
 *         * source/headers/geos/platform.h.in: include <inttypes.h>
 *         rather then <stdint.h>
 *
 * Revision 1.70  2006/06/19 21:17:23  strk
 * port info and doxygen dox.
 *
 * Revision 1.69  2006/06/12 10:10:39  strk
 * Fixed getGeometryN() to take size_t rather then int, changed unsigned int parameters to size_t.
 *
 * Revision 1.68  2006/04/28 11:56:52  strk
 * * source/geom/GeometryFactory.cpp, source/headers/geos/geom/GeometryFactory.h: added LineString copy constructor.
 * * source/geom/Polygon.cpp: fixed getBoundary method to always return a geometry composed by LineStrings (not LinearRings)
 *
 * Revision 1.67  2006/04/11 11:16:25  strk
 * Added LineString and LinearRing constructors by auto_ptr
 *
 * Revision 1.66  2006/04/10 13:09:47  strk
 * Added GeometryFactory::defaultInstance()
 * Made Geometry::INTERNAL_GEOMETRY_FACTORY an alias for it
 * removed last deletion from Unload::Release class
 *
 * Revision 1.65  2006/04/06 12:33:04  strk
 * More debugging lines
 *
 * Revision 1.64  2006/03/31 17:51:24  strk
 * A few assertion checking, comments cleanup, use of initialization lists
 * in constructors, handled NULL parameters.
 *
 * Revision 1.63  2006/03/24 09:52:41  strk
 * USE_INLINE => GEOS_INLINE
 *
 * Revision 1.62  2006/03/20 10:11:50  strk
 * Bug #67 - Debugging helpers in GeometryFactory class
 *
 * Revision 1.61  2006/03/09 16:46:47  strk
 * geos::geom namespace definition, first pass at headers split
 *
 * Revision 1.60  2006/03/06 19:40:46  strk
 * geos::util namespace. New GeometryCollection::iterator interface, many cleanups.
 *
 * Revision 1.59  2006/03/03 10:46:21  strk
 * Removed 'using namespace' from headers, added missing headers in .cpp files, removed useless includes in headers (bug#46)
 *
 * Revision 1.58  2006/03/01 18:37:08  strk
 * Geometry::createPointFromInternalCoord dropped (it's a duplication of GeometryFactory::createPointFromInternalCoord).
 * Fixed bugs in InteriorPoint* and getCentroid() inserted by previous commits.
 *
 * Revision 1.57  2006/02/23 23:17:52  strk
 * - Coordinate::nullCoordinate made private
 * - Simplified Coordinate inline definitions
 * - LMGeometryComponentFilter definition moved to LineMerger.cpp file
 * - Misc cleanups
 *
 * Revision 1.56  2006/02/09 15:52:47  strk
 * GEOSException derived from std::exception; always thrown and cought by const ref.
 *
 **********************************************************************/

