/**********************************************************************
 * $Id: SweepLineEvent.cpp 1820 2006-09-06 16:54:23Z mloskot $
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
 **********************************************************************/

#include <geos/index/sweepline/SweepLineEvent.h>

namespace geos {
namespace index { // geos.index
namespace sweepline { // geos.index.sweepline

SweepLineEvent::SweepLineEvent(double x, SweepLineEvent *newInsertEvent,
		SweepLineInterval *newSweepInt)
	:
	xValue(x),
	eventType(SweepLineEvent::INSERT_EVENT),
	insertEvent(newInsertEvent),
	sweepInt(newSweepInt)
{
	if (insertEvent!=0)
		eventType=SweepLineEvent::DELETE_EVENT;
}

bool
SweepLineEvent::isInsert()
{
	return insertEvent==0;
}

bool
SweepLineEvent::isDelete()
{
	return insertEvent!=0;
}

SweepLineEvent*
SweepLineEvent::getInsertEvent()
{
	return insertEvent;
}

int
SweepLineEvent::getDeleteEventIndex()
{
	return deleteEventIndex;
}

void
SweepLineEvent::setDeleteEventIndex(int newDeleteEventIndex)
{
	deleteEventIndex=newDeleteEventIndex;
}

SweepLineInterval*
SweepLineEvent::getInterval()
{
	return sweepInt;
}

int
SweepLineEvent::compareTo(const SweepLineEvent *pe) const
{
	if (xValue<pe->xValue) return -1;
	if (xValue>pe->xValue) return 1;
	if (eventType<pe->eventType) return -1;
	if (eventType>pe->eventType) return 1;
	return 0;
}

#if 0
int
SweepLineEvent::compareTo(void *o) const
{
	SweepLineEvent *pe=(SweepLineEvent*) o;
	if (xValue<pe->xValue) return -1;
	if (xValue>pe->xValue) return 1;
	if (eventType<pe->eventType) return -1;
	if (eventType>pe->eventType) return 1;
	return 0;
}
#endif // 0

bool
SweepLineEventLessThen::operator() (const SweepLineEvent *first, const SweepLineEvent *second) const
{
	if (first->compareTo(second)<0)
		return true;
	else
		return false;
}


} // namespace geos.index.sweepline
} // namespace geos.index
} // namespace geos

/**********************************************************************
 * $Log$
 * Revision 1.1  2006/03/21 10:01:30  strk
 * indexSweepline.h header split
 *
 * Revision 1.8  2006/03/02 14:53:44  strk
 * SweepLineEvent::DELETE=>DELETE_EVENT, INSERT=>INSERT_EVENT (#45)
 *
 * Revision 1.7  2006/02/20 10:14:18  strk
 * - namespaces geos::index::*
 * - Doxygen documentation cleanup
 *
 * Revision 1.6  2004/07/02 13:28:27  strk
 * Fixed all #include lines to reflect headers layout change.
 * Added client application build tips in README.
 *
 * Revision 1.5  2003/11/07 01:23:42  pramsey
 * Add standard CVS headers licence notices and copyrights to all cpp and h
 * files.
 *
 *
 **********************************************************************/

