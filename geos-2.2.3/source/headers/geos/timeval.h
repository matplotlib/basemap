/*
 * timeval.h    1.3 2003/01/14
 *
 * Defines gettimeofday, timeval, etc. for Win32
 *
 * By Wu Yongwei
 *
 */

#ifndef _TIMEVAL_H
#define _TIMEVAL_H

#ifdef _WIN32

#define WIN32_LEAN_AND_MEAN
#include <winsock2.h>
#include <time.h>

#if defined(_MSC_VER) || defined(__BORLANDC__)
#define EPOCHFILETIME (116444736000000000i64)
#else
#define EPOCHFILETIME (116444736000000000LL)
#endif

#else  /* _WIN32 */

#include <sys/time.h>

#endif /* _WIN32 */

#endif /* _TIMEVAL_H */
