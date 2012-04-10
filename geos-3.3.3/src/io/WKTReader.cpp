/**********************************************************************
 * $Id: WKTReader.cpp 3113 2010-10-15 14:55:22Z mloskot $
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
 * Last port: io/WKTReader.java rev. 1.1 (JTS-1.7)
 *
 **********************************************************************/

#include <geos/io/WKTReader.h>
#include <geos/io/StringTokenizer.h>
#include <geos/io/ParseException.h>
#include <geos/io/CLocalizer.h>
#include <geos/geom/GeometryFactory.h>
#include <geos/geom/Coordinate.h>
#include <geos/geom/Point.h>
#include <geos/geom/LinearRing.h>
#include <geos/geom/LineString.h>
#include <geos/geom/Polygon.h>
#include <geos/geom/MultiPoint.h>
#include <geos/geom/MultiLineString.h>
#include <geos/geom/MultiPolygon.h>
#include <geos/geom/CoordinateSequenceFactory.h>
#include <geos/geom/CoordinateSequence.h>
#include <geos/geom/PrecisionModel.h>
#include <geos/inline.h>

#include <sstream>
#include <string>
#include <cassert>

#ifndef GEOS_DEBUG
#define GEOS_DEBUG 0
#endif

#ifdef GEOS_DEBUG
#include <iostream>
#endif

#ifndef GEOS_INLINE
#include <geos/io/WKTReader.inl>
#endif

using namespace std;
using namespace geos::geom;

namespace geos {
namespace io { // geos.io

Geometry *
WKTReader::read(const string &wellKnownText)
{
	//auto_ptr<StringTokenizer> tokenizer(new StringTokenizer(wellKnownText));
        CLocalizer clocale;
	StringTokenizer tokenizer(wellKnownText);
	Geometry *g=NULL;
	g=readGeometryTaggedText(&tokenizer);
	return g;
}

CoordinateSequence*
WKTReader::getCoordinates(StringTokenizer *tokenizer)
{
	size_t dim;
	string nextToken=getNextEmptyOrOpener(tokenizer);
	if (nextToken=="EMPTY") {
		return geometryFactory->getCoordinateSequenceFactory()->create(NULL);
		//new CoordinateArraySequence(); 
	}

	Coordinate coord;
	getPreciseCoordinate(tokenizer, coord, dim);

	CoordinateSequence *coordinates = \
            geometryFactory->getCoordinateSequenceFactory()->create((size_t)0,dim);
	coordinates->add(coord);
	try {
		nextToken=getNextCloserOrComma(tokenizer);
		while (nextToken==",") {
			getPreciseCoordinate(tokenizer, coord, dim );
			coordinates->add(coord);
			nextToken=getNextCloserOrComma(tokenizer);
		}
	} catch (...) {
		delete coordinates;
		throw;
	}

	return coordinates;
}

void
WKTReader::getPreciseCoordinate(StringTokenizer *tokenizer, 
                                Coordinate& coord,
                                size_t &dim )
{
	coord.x=getNextNumber(tokenizer);
	coord.y=getNextNumber(tokenizer);
	if (isNumberNext(tokenizer)) {
		coord.z=getNextNumber(tokenizer);
		dim = 3;
        
        // If there is a fourth value (M) read and discard it.
        if (isNumberNext(tokenizer)) 
            getNextNumber(tokenizer);

	} else {
		coord.z=DoubleNotANumber;
		dim = 2;
	}
	precisionModel->makePrecise(coord);
}

bool
WKTReader::isNumberNext(StringTokenizer *tokenizer)
{
	return tokenizer->peekNextToken()==StringTokenizer::TT_NUMBER;
}

double
WKTReader::getNextNumber(StringTokenizer *tokenizer)
{
	int type=tokenizer->nextToken();
	switch(type){
		case StringTokenizer::TT_EOF:
			throw  ParseException("Expected number but encountered end of stream");
		case StringTokenizer::TT_EOL:
			throw  ParseException("Expected number but encountered end of line");
		case StringTokenizer::TT_NUMBER:
			return tokenizer->getNVal();
		case StringTokenizer::TT_WORD:
			throw  ParseException("Expected number but encountered word",tokenizer->getSVal());
		case '(':
			throw  ParseException("Expected number but encountered '('");
		case ')':
			throw  ParseException("Expected number but encountered ')'");
		case ',':
			throw  ParseException("Expected number but encountered ','");
	}
	assert(0); // Encountered unexpected StreamTokenizer type
	return 0;
}

string
WKTReader::getNextEmptyOrOpener(StringTokenizer *tokenizer)
{
	string nextWord=getNextWord(tokenizer);

    // Skip the Z, M or ZM of an SF1.2 3/4 dim coordinate. 
    if (nextWord == "Z" || nextWord == "M" || nextWord == "ZM" )
        nextWord = getNextWord(tokenizer);

	if (nextWord=="EMPTY" || nextWord=="(") {
		return nextWord;
	}
	throw  ParseException("Expected 'Z', 'M', 'ZM', 'EMPTY' or '(' but encountered ",nextWord);
}

string
WKTReader::getNextCloserOrComma(StringTokenizer *tokenizer)
{
	string nextWord=getNextWord(tokenizer);
	if (nextWord=="," || nextWord==")") {
		return nextWord;
	}
	throw  ParseException("Expected ')' or ',' but encountered",nextWord);
}

string
WKTReader::getNextCloser(StringTokenizer *tokenizer)
{
	string nextWord=getNextWord(tokenizer);
	if (nextWord==")") {
		return nextWord;
	}
	throw  ParseException("Expected ')' but encountered",nextWord);
}

string
WKTReader::getNextWord(StringTokenizer *tokenizer)
{
	int type=tokenizer->nextToken();
	switch(type){
		case StringTokenizer::TT_EOF:
			throw  ParseException("Expected word but encountered end of stream");
		case StringTokenizer::TT_EOL:
			throw  ParseException("Expected word but encountered end of line");
		case StringTokenizer::TT_NUMBER:
			throw  ParseException("Expected word but encountered number", tokenizer->getNVal());
		case StringTokenizer::TT_WORD:
        {
            string word = tokenizer->getSVal();
            int i = word.size();

            while( --i >= 0 )
            {
                word[i] = static_cast<char>(toupper(word[i]));
            }
			return word;
        }
		case '(':
			return "(";
		case ')':
			return ")";
		case ',':
			return ",";
	}
	assert(0);
	//throw  ParseException("Encountered unexpected StreamTokenizer type");
	return "";
}

Geometry*
WKTReader::readGeometryTaggedText(StringTokenizer *tokenizer)
{
	string type = getNextWord(tokenizer);
	if (type=="POINT") {
		return readPointText(tokenizer);
	} else if (type=="LINESTRING") {
		return readLineStringText(tokenizer);
	} else if (type=="LINEARRING") {
		return readLinearRingText(tokenizer);
	} else if (type=="POLYGON") {
		return readPolygonText(tokenizer);
	} else if (type=="MULTIPOINT") {
		return readMultiPointText(tokenizer);
	} else if (type=="MULTILINESTRING") {
		return readMultiLineStringText(tokenizer);
	} else if (type=="MULTIPOLYGON") {
		return readMultiPolygonText(tokenizer);
	} else if (type=="GEOMETRYCOLLECTION") {
		return readGeometryCollectionText(tokenizer);
	}
	throw  ParseException("Unknown type",type);
}

Point*
WKTReader::readPointText(StringTokenizer *tokenizer)
{
	size_t dim;
	string nextToken=getNextEmptyOrOpener(tokenizer);
	if (nextToken=="EMPTY") {
		return geometryFactory->createPoint(Coordinate::getNull());
	}

	Coordinate coord;
	getPreciseCoordinate(tokenizer, coord, dim);
	getNextCloser(tokenizer);

	return geometryFactory->createPoint(coord);
}

LineString* WKTReader::readLineStringText(StringTokenizer *tokenizer) {
	CoordinateSequence *coords = getCoordinates(tokenizer);
	LineString *ret = geometryFactory->createLineString(coords);
	return ret;
}

LinearRing* WKTReader::readLinearRingText(StringTokenizer *tokenizer) {
	CoordinateSequence *coords = getCoordinates(tokenizer);
	LinearRing *ret;
	ret = geometryFactory->createLinearRing(coords);
	return ret;
}

MultiPoint*
WKTReader::readMultiPointText(StringTokenizer *tokenizer)
{
	string nextToken=getNextEmptyOrOpener(tokenizer);
	if (nextToken=="EMPTY") {
		return geometryFactory->createMultiPoint();
	}

	int tok = tokenizer->peekNextToken();

	if ( tok == StringTokenizer::TT_NUMBER )
	{
		size_t dim;

		// Try to parse deprecated form "MULTIPOINT(0 0, 1 1)"
		const CoordinateSequenceFactory* csf = \
			geometryFactory->getCoordinateSequenceFactory();
		CoordinateSequence *coords = csf->create(NULL);
		try {
			do {
				Coordinate coord;
				getPreciseCoordinate(tokenizer, coord, dim);
				coords->add(coord);
				nextToken=getNextCloserOrComma(tokenizer);
			} while(nextToken == ",");

			MultiPoint *ret = geometryFactory->createMultiPoint(*coords);
			delete coords;
			return ret;
		} catch (...) {
			delete coords;
			throw;
		}
	}

	else if ( tok == '(' )
	{
		// Try to parse correct form "MULTIPOINT((0 0), (1 1))"
		vector<Geometry *> *points=new vector<Geometry *>();
		try {
			do {
				Point *point=readPointText(tokenizer);
				points->push_back(point);
				nextToken=getNextCloserOrComma(tokenizer);
			} while(nextToken == ",");
			return geometryFactory->createMultiPoint(points);
		} catch (...) {
			// clean up 
			for (size_t i=0; i<points->size(); i++)
			{
				delete (*points)[i];
			}
			delete points;
			throw;
		}
	}

	else 
	{
		stringstream err;
		err << "Unexpected token: ";
		switch (tok)
		{
			case StringTokenizer::TT_WORD:
				err << "WORD " << tokenizer->getSVal();
				break;
			case StringTokenizer::TT_NUMBER:
				err << "NUMBER " << tokenizer->getNVal();
				break;
			case StringTokenizer::TT_EOF:
			case StringTokenizer::TT_EOL:
				err << "EOF or EOL";
				break;
			case '(':
				err << "(";
				break;
			case ')':
				err << ")";
				break;
			case ',':
				err << ",";
				break;
			default:
				err << "??";
				break;
		}
		err << endl;
		throw ParseException(err.str());
	}
}

Polygon*
WKTReader::readPolygonText(StringTokenizer *tokenizer)
{
	Polygon *poly=NULL;
	LinearRing *shell=NULL;
	string nextToken=getNextEmptyOrOpener(tokenizer);
	if (nextToken=="EMPTY") {
		return geometryFactory->createPolygon(NULL,NULL);
	}

	vector<Geometry *> *holes=new vector<Geometry *>();
	try {
		shell=readLinearRingText(tokenizer);
		nextToken=getNextCloserOrComma(tokenizer);
		while(nextToken==",") {
			LinearRing *hole=readLinearRingText(tokenizer);
			holes->push_back(hole);
			nextToken=getNextCloserOrComma(tokenizer);
		}
		poly = geometryFactory->createPolygon(shell,holes);
	} catch (...) {
		for (unsigned int i=0; i<holes->size(); i++)
			delete (*holes)[i];
		delete holes;
		delete shell;
		throw;
	}
	return poly;
}

MultiLineString* WKTReader::readMultiLineStringText(StringTokenizer *tokenizer) {
	string nextToken=getNextEmptyOrOpener(tokenizer);
	if (nextToken=="EMPTY") {
		return geometryFactory->createMultiLineString(NULL);
	}
	vector<Geometry *> *lineStrings=new vector<Geometry *>();
	LineString *lineString=readLineStringText(tokenizer);
	lineStrings->push_back(lineString);
	nextToken=getNextCloserOrComma(tokenizer);
	while(nextToken==",") {
		LineString *lineString=readLineStringText(tokenizer);
		lineStrings->push_back(lineString);
		nextToken=getNextCloserOrComma(tokenizer);
	}
	MultiLineString *ret = geometryFactory->createMultiLineString(lineStrings);
	//for (int i=0; i<lineStrings->size(); i++) delete (*lineStrings)[i];
	//delete lineStrings;
	return ret;
}

MultiPolygon* WKTReader::readMultiPolygonText(StringTokenizer *tokenizer) {
	string nextToken=getNextEmptyOrOpener(tokenizer);
	if (nextToken=="EMPTY") {
		return geometryFactory->createMultiPolygon(NULL);
	}
	vector<Geometry *> *polygons=new vector<Geometry *>();
	Polygon *polygon=readPolygonText(tokenizer);
	polygons->push_back(polygon);
	nextToken=getNextCloserOrComma(tokenizer);
	while(nextToken==",") {
		Polygon *polygon=readPolygonText(tokenizer);
		polygons->push_back(polygon);
		nextToken=getNextCloserOrComma(tokenizer);
	}
	MultiPolygon *ret = geometryFactory->createMultiPolygon(polygons);
	//for (int i=0; i<polygons->size(); i++) delete (*polygons)[i];
	//delete polygons;
	return ret;
}

GeometryCollection* WKTReader::readGeometryCollectionText(StringTokenizer *tokenizer) {
	string nextToken=getNextEmptyOrOpener(tokenizer);
	if (nextToken=="EMPTY") {
		return geometryFactory->createGeometryCollection(NULL);
	}
	vector<Geometry *> *geoms=new vector<Geometry *>();
	Geometry *geom;
	geom=readGeometryTaggedText(tokenizer);
	geoms->push_back(geom);
	nextToken=getNextCloserOrComma(tokenizer);
	while(nextToken==",") {
		geom=readGeometryTaggedText(tokenizer);
		geoms->push_back(geom);
		nextToken=getNextCloserOrComma(tokenizer);
	}
	GeometryCollection *ret = geometryFactory->createGeometryCollection(geoms);
	//for (int i=0; i<geoms->size(); i++) delete (*geoms)[i];
	//delete geoms;
	return ret;
}

} // namespace geos.io
} // namespace geos

/**********************************************************************
 * $Log$
 * Revision 1.42  2006/04/26 16:35:40  strk
 * Had WKTReader accept correct form for MultiPoint
 *
 * Revision 1.41  2006/04/10 12:05:35  strk
 * Added inline-replicator implementation files to make sure
 * functions in .inl files are still available out-of-line.
 * A side effect is this should fix MingW build.
 *
 * Revision 1.40  2006/03/24 09:52:41  strk
 * USE_INLINE => GEOS_INLINE
 *
 * Revision 1.39  2006/03/22 16:58:35  strk
 * Removed (almost) all inclusions of geom.h.
 * Removed obsoleted .cpp files.
 * Fixed a bug in WKTReader not using the provided CoordinateSequence
 * implementation, optimized out some memory allocations.
 *
 **********************************************************************/
