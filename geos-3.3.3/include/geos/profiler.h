/**********************************************************************
 * $Id: profiler.h 3255 2011-03-01 17:56:10Z mloskot $
 *
 * GEOS - Geometry Engine Open Source
 * http://geos.refractions.net
 *
 * Copyright (C) 2001-2002 Vivid Solutions Inc.
 *
 * This is free software; you can redistribute and/or modify it under
 * the terms of the GNU Lesser General Public Licence as published
 * by the Free Software Foundation. 
 * See the COPYING file for more information.
 *
 **********************************************************************/

#ifndef GEOS_PROFILER_H
#define GEOS_PROFILER_H

#include <geos/export.h>

/* For MingW builds with __STRICT_ANSI__ (-ansi) */
#if defined(__MINGW32__)
/* Allow us to check for presence of gettimeofday in MingW */ 
#include <config.h>

#include <sys/time.h>
extern "C" {
  extern _CRTIMP void __cdecl	_tzset (void);
  __MINGW_IMPORT int	_daylight;
  __MINGW_IMPORT long	_timezone;
  __MINGW_IMPORT char 	*_tzname[2];
}
#endif
 
#if defined(_MSC_VER) || defined(__MINGW32__) && !defined(HAVE_GETTIMEOFDAY)
#include <geos/timeval.h>
#else
#include <sys/time.h>
#endif

#include <map>
#include <memory>
#include <iostream>
#include <string>
#include <vector>

#ifndef PROFILE
#define PROFILE 0
#endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4251) // warning C4251: needs to have dll-interface to be used by clients of class
#endif

namespace geos {
namespace util {


/*
 * \class Profile utils.h geos.h
 *
 * \brief Profile statistics
 */
class GEOS_DLL Profile {
public:
	/** \brief Create a named profile */
	Profile(std::string name);

	/** \brief Destructor */
	~Profile();

	/** \brief start a new timer */
	void start() {
		gettimeofday(&starttime, NULL);
	}

	/** \brief stop current timer */
	void stop()
	{
		gettimeofday(&stoptime, NULL);
		double elapsed = 1000000*(stoptime.tv_sec-starttime.tv_sec)+
			(stoptime.tv_usec-starttime.tv_usec);

		timings.push_back(elapsed);
		totaltime += elapsed;
		if ( timings.size() == 1 ) max = min = elapsed;
		else
		{
			if ( elapsed > max ) max = elapsed;
			if ( elapsed < min ) min = elapsed;
		}
		avg = totaltime / timings.size();
	}

	/** \brief Return Max stored timing */
	double getMax() const;

	/** \brief Return Min stored timing */
	double getMin() const;

	/** \brief Return total timing */
	double getTot() const;

	/** \brief Return average timing */
	double getAvg() const;

	/** \brief Return number of timings */
	size_t getNumTimings() const;

	/** \brief Profile name */
	std::string name;


private:

	/* \brief current start and stop times */
	struct timeval starttime, stoptime;

	/* \brief actual times */
	std::vector<double> timings;

	/* \brief total time */
	double totaltime;

	/* \brief max time */
	double max;

	/* \brief max time */
	double min;

	/* \brief max time */
	double avg;

};

/*
 * \class Profiler utils.h geos.h
 *
 * \brief Profiling class
 *
 */
class GEOS_DLL Profiler {

public:

	Profiler();
	~Profiler();

	/**
	 * \brief
	 * Return the singleton instance of the
	 * profiler.
	 */
	static Profiler *instance(void);

	/**
	 * \brief
	 * Start timer for named task. The task is
	 * created if does not exist.
	 */
	void start(std::string name);

	/**
	 * \brief
	 * Stop timer for named task. 
	 * Elapsed time is registered in the given task.
	 */
	void stop(std::string name);

	/** \brief get Profile of named task */
	Profile *get(std::string name);

	std::map<std::string, Profile *> profs;
};


/** \brief Return a string representing the Profile */
std::ostream& operator<< (std::ostream& os, const Profile&);

/** \brief Return a string representing the Profiler */
std::ostream& operator<< (std::ostream& os, const Profiler&);

} // namespace geos::util
} // namespace geos

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif // ndef GEOS_PROFILER_H

/**********************************************************************
 * $Log$
 * Revision 1.8  2006/06/12 11:29:23  strk
 * unsigned int => size_t
 *
 * Revision 1.7  2006/03/06 19:40:47  strk
 * geos::util namespace. New GeometryCollection::iterator interface, many cleanups.
 *
 * Revision 1.6  2006/03/03 10:46:21  strk
 * Removed 'using namespace' from headers, added missing headers in .cpp files, removed useless includes in headers (bug#46)
 *
 * Revision 1.5  2005/02/01 14:18:04  strk
 * Made profiler start/stop inline
 *
 * Revision 1.4  2004/12/03 16:21:07  frank
 * dont try for sys/time.h with MSVC
 *
 * Revision 1.3  2004/11/30 16:44:16  strk
 * Added gettimeofday implementation for win32, curtesy of Wu Yongwei.
 *
 * Revision 1.2  2004/11/04 08:49:13  strk
 * Unlinked new documentation.
 *
 * Revision 1.1  2004/11/01 16:43:04  strk
 * Added Profiler code.
 * Temporarly patched a bug in DoubleBits (must check drawbacks).
 * Various cleanups and speedups.
 *
 **********************************************************************/
