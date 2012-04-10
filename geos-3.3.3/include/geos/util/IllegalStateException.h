/**********************************************************************
 * $Id$
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2011 Sandro Santilli <strk@keybit.net>
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation. 
 * See the COPYING file for more information.
 *
 **********************************************************************/

#ifndef GEOS_UTIL_ILLEGALSTATEEXCEPTION_H
#define GEOS_UTIL_ILLEGALSTATEEXCEPTION_H

#include <geos/export.h>
#include <string>

#include <geos/util/GEOSException.h>

namespace geos {
namespace util { // geos::util

/// Indicates an illegal state
class GEOS_DLL IllegalStateException: public GEOSException {
public:
	IllegalStateException()
		:
		GEOSException("IllegalStateException", "")
	{}

	IllegalStateException(const std::string& msg)
		:
		GEOSException("IllegalStateException", msg)
	{}

	~IllegalStateException() throw() {};
};

} // namespace geos::util
} // namespace geos


#endif // GEOS_UTIL_ILLEGALSTATEEXCEPTION_H
