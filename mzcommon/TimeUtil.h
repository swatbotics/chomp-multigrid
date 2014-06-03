/*
* Copyright (c) 2008-2014, Matt Zucker
*
* This file is provided under the following "BSD-style" License:
*
* Redistribution and use in source and binary forms, with or
* without modification, are permitted provided that the following
* conditions are met:
*
* * Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
*
* * Redistributions in binary form must reproduce the above
* copyright notice, this list of conditions and the following
* disclaimer in the documentation and/or other materials provided
* with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
* CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
* USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*/

#ifndef _TIMEUTIL_H_
#define _TIMEUTIL_H_

#include <iostream>
#include <math.h>

#ifdef _WIN32
// TODO: win32 time of day
#else
#include <sys/time.h>
#include <sys/types.h>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdint.h>
#endif

enum { NSEC_PER_SEC = 1000000000 };



#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

class Duration;

/** Encapsulate an absolute time reference. Represented as nanoseconds
 *  since the Unix epoch.
 */
class TimeStamp {
public:

  /** 64-bit integer holding nanoseconds since the Unix epoch.
   *  If this is equal to uint64_t(-1), then the TimeStamp is considered
   *  invalid.
   */
  uint64_t nsec;

  /** Construct the zero TimeStamp. */
  TimeStamp(): nsec(0) {}

  /** Construct with given number of nanoseconds. */
  explicit TimeStamp(uint64_t n): nsec(n) {}

  /** Construct with seconds and nanoseconds. */
  TimeStamp(unsigned long s, unsigned long n) {
    nsec = (uint64_t)s*(uint64_t)NSEC_PER_SEC + (uint64_t)n;
  }
  
  /** Add a Duration to a TimeStamp to get another TimeStamp. */
  TimeStamp& operator+=(const Duration& d);

  /** Subtrack a Duration from a TimeStamp to get another TimeStamp. */
  TimeStamp& operator-=(const Duration& d);

  /** Return true if nsec == uint64_t(-1). */
  bool valid() const {
    return (nsec != (uint64_t)-1);
  }

  /** Construct the invalid TimeStamp. */
  static TimeStamp invalid() { 
    return TimeStamp((uint64_t)-1);
  }

  /** Return the modification time of the given file or directory.
   *  @param path The path of the file or directory.
   *  @param bad  The TimeStamp to return if the file was not found.
   */
  static TimeStamp mtime(const char* path, const TimeStamp bad=TimeStamp::invalid()) {

    struct stat sb;

    if (stat(path, &sb) != 0) {
      return bad;
    }
    
    TimeStamp rval;
    rval.nsec = (uint64_t)sb.st_mtime * (uint64_t)NSEC_PER_SEC;

    return rval;

  }

  /** Return the current time as a TimeStamp. */
  static TimeStamp now() {
    struct timeval tv;
#ifdef _WIN32
    union {
      LONG_LONG ns100; /*time since 1 Jan 1601 in 100ns units */
      FILETIME ft;
    } now;
    GetSystemTimeAsFileTime (&now.ft);
    tv.tv_usec = (long) ((now.ns100 / 10LL) % 1000000LL);
    tv.tv_sec = (long) ((now.ns100 - 116444736000000000LL) / 10000000LL);
#else
    gettimeofday(&tv, 0);
#endif
    return fromTimeval(tv);
  }
      
  /** Construct from double storing seconds since the Unix epoch. */
  static TimeStamp fromDouble(double d) {
    if (d < 0) { return TimeStamp(); }
    return TimeStamp((uint64_t) floor(d * double(NSEC_PER_SEC)));
  }

  /** Convert to double storing seconds since the Unix
   *  epoch. Converting back may result in loss of precision. */
  double toDouble() const {
    return nsec / double(NSEC_PER_SEC);
  }

  /** Construct a TimeStamp from a Unix timespec (used by
   *  pthread_cond_timedwait and other functions). */
  static TimeStamp fromTimespec(const struct timespec& ts) {
    TimeStamp rval;
    rval.nsec = (uint64_t)ts.tv_sec * (uint64_t)NSEC_PER_SEC +
      ts.tv_nsec;
    return rval;
  }

  /** Convert a TimeStamp to a Unix timespec (used by
   *  pthread_cond_timedwait and other functions). */
  struct timespec toTimespec() const {
    struct timespec ts;
    ts.tv_sec = nsec / (uint64_t)NSEC_PER_SEC;
    ts.tv_nsec = nsec % (uint64_t)NSEC_PER_SEC;
    return ts;
  }

  /** Construct a TimeStamp from a Unix timeval (used by
   *  select and other functions). */
  static TimeStamp fromTimeval(const struct timeval& tv) {
    TimeStamp rval;
    rval.nsec = (uint64_t)tv.tv_sec * (uint64_t)NSEC_PER_SEC +
      (uint64_t)tv.tv_usec * (uint64_t)1000;
    return rval;
  }

  /** Convert a TimeStamp to a Unix timeval (used by
   *  select and other functions). */
  struct timeval toTimeval() const {
    struct timeval rval;
    rval.tv_sec = nsec / (uint64_t)NSEC_PER_SEC;
    rval.tv_usec = (nsec % (uint64_t)NSEC_PER_SEC) / 1000;
    return rval;
  }

};

/** Encapsulate a relative time interval, represented as signed number
 *  of nanoseconds.
 */
class Duration {
public:

  /** Signed 64-bit integer representing number of nanoseconds. */
  int64_t nsec;

  /** Construct a zero duration. */
  Duration(): nsec(0) {}

  /** Construct with given number of nanoseconds. */
  explicit Duration(int64_t n): nsec(n) {}

  /** Cast to void in order to allow boolean checking of non-zero
   *  durations. */
  operator void*() const {
    return (void*)(nsec != 0);
  }

  /** Construct from seconds and nanoseconds. */
  Duration(long s, long n) {
    nsec = (int64_t)s*(int64_t)NSEC_PER_SEC + (int64_t)n;
  }

  /** Addition-assignment. */
  Duration& operator+=(const Duration& d);

  /** Subtractinon-assignment. */
  Duration& operator-=(const Duration& d);

  /** Integer multiplication-assignment. */
  Duration& operator*=(int s);

  /** Double-precision multiplication-assignment. May result in loss of precision. */
  Duration& operator*=(double s);

  /** Double-precision division-assignment. May result in loss of precision. */
  Duration& operator/=(double s);

  /** Construct a Duration from a Unix timeval (used by
   *  select and other functions). */
  static Duration fromTimeval(const struct timeval& tv) {
    Duration rval;
    rval.nsec = (int64_t)tv.tv_sec*(int64_t)NSEC_PER_SEC + (int64_t)tv.tv_usec * (int64_t)1000;
    return rval;
  }

  /** Convert a Duration to a Unix timeval (used by select and other
   *  functions). */
  struct timeval toTimeval() const {
    struct timeval tv;
    tv.tv_sec = nsec / (int64_t)NSEC_PER_SEC;
    tv.tv_usec = (nsec % (int64_t)NSEC_PER_SEC) / (int64_t)1000;
    return tv;
  }


  /** Convert a Duration to a Unix timespec (used by
   *  pthread_cond_timedwait and other functions). */
  struct timespec toTimespec() const {
    struct timespec ts;
    ts.tv_sec = nsec / (int64_t)NSEC_PER_SEC;
    ts.tv_nsec = nsec % (int64_t)NSEC_PER_SEC;
    return ts;
  }

  /** Construct from double storing seconds. */
  static Duration fromDouble(double d) {
    return Duration((int64_t) trunc(d * double(NSEC_PER_SEC)));
  }

  /** Convert to double storing seconds. */
  double toDouble() const {
    return nsec / double(NSEC_PER_SEC);
  }

  /** Return a duration of one microsecond. */
  static Duration microsecond() {
    return Duration(1000);
  }

  /** Return a duration of one millisecond. */
  static Duration millisecond() {
    return Duration(1000000);
  }

  /** Return a duration of one second. */
  static Duration second() {
    return Duration(1000000000);
  }

};

//////////////////////////////////////////////////////////////////////

/** @relates TimeStamp
 *  Subtract two TimeStamps to get a Duration.
 */
inline Duration operator-(const TimeStamp& t1, const TimeStamp& t2) {
  return Duration(t1.nsec - t2.nsec);
}

/** @relates Duration
 *  Subtract two Durations.
 */
inline Duration operator-(const Duration& d1, const Duration& d2) {
  return Duration(d1.nsec - d2.nsec);
}

/** @relates Duration
 *  Add two Durations.
 */
inline Duration operator+(const Duration& d1, const Duration& d2) {
  return Duration(d1.nsec + d2.nsec);
}

/** @relates TimeStamp
 *  Add a Duration to a TimeStamp to get a new TimeStamp.
 */
inline TimeStamp operator+(const TimeStamp& t, const Duration& d) {
  return TimeStamp(t.nsec + d.nsec);
}

/** @relates TimeStamp
 *  Subtract a Duration from a TimeStamp to get a new TimeStamp.
 */
inline TimeStamp operator-(const TimeStamp& t, const Duration& d) {
  return TimeStamp(t.nsec - d.nsec);
}

/** @relates Duration
 *  Multiply a Duration by double-precision float. May result in loss of precision.
 */
inline Duration operator*(const Duration& d, double s) {
  return Duration::fromDouble(d.toDouble()*s);
}

/** @relates Duration
 *  Multiply a Duration by single-precision float. May result in loss of precision.
 */
inline Duration operator*(const Duration& d, float s) {
  return Duration::fromDouble(d.toDouble()*s);
}

/** @relates Duration
 *  Multiply a Duration by double-precision float. May result in loss of precision.
 */
inline Duration operator*(double s, const Duration& d) {
  return Duration::fromDouble(d.toDouble()*s);
}

/** @relates Duration
 *  Multiply a Duration by single-precision float. May result in loss of precision.
 */
inline Duration operator*(float s, const Duration& d) {
  return Duration::fromDouble(d.toDouble()*s);
}

/** @relates Duration
 *  Multiply a Duration by signed integer. 
 */
inline Duration operator*(int s, const Duration& d) {
  return Duration(d.nsec * int64_t(s));
}

/** @relates Duration
 *  Multiply a Duration by signed integer. 
 */
inline Duration operator*(const Duration& d, int s) {
  return Duration(d.nsec * int64_t(s));
}

/** @relates Duration
 *  Multiply a Duration by unsigned integer. 
 */
inline Duration operator*(unsigned int s, const Duration& d) {
  return Duration(d.nsec * int64_t(s));
}

/** @relates Duration
 *  Multiply a Duration by unsigned integer. 
 */
inline Duration operator*(const Duration& d, unsigned int s) {
  return Duration(d.nsec * int64_t(s));
}

/** @relates Duration
 *  Divide a Duration by double-precision float. May result in loss of precision.
 */
inline Duration operator/(const Duration& d, double s) {
  return Duration::fromDouble(d.toDouble()/s);
}


/** @relates Duration
 *  Subtract two Durations. 
 */
inline Duration operator-(const Duration& d) {
  return Duration(-(d.nsec));
}

//////////////////////////////////////////////////////////////////////

/** @relates Duration
 *  Stream output for duration. 
 */
inline std::ostream& operator<<(std::ostream& ostr, const Duration& d) {
  ostr << d.toDouble() << "s";
  return ostr;
}

/** @relates TimeStamp
 *  Stream output for TimeStamp.
 */
inline std::ostream& operator<<(std::ostream& ostr, const TimeStamp& t) {
  time_t clock = time_t(t.toDouble());
  struct tm tm;
  localtime_r(&clock, &tm);
  char buf[1024];
  strftime(buf, 1024, "%a %b %e %H:%M:%S %Z %Y", &tm);
  ostr << buf;
  return ostr;
}

//////////////////////////////////////////////////////////////////////

/** @relates TimeStamp
 *  Compare TimeStamps.
 */
inline bool operator<(const TimeStamp& t1, const TimeStamp& t2) {
  return (t1.nsec < t2.nsec);
}

/** @relates TimeStamp
 *  Compare TimeStamps.
 */
inline bool operator>(const TimeStamp& t1, const TimeStamp& t2) {
  return (t1.nsec > t2.nsec);
}

/** @relates TimeStamp
 *  Compare TimeStamps.
 */
inline bool operator<=(const TimeStamp& t1, const TimeStamp& t2) {
  return (t1.nsec <= t2.nsec);
}

/** @relates TimeStamp
 *  Compare TimeStamps.
 */
inline bool operator>=(const TimeStamp& t1, const TimeStamp& t2) {
  return (t1.nsec >= t2.nsec);
}

/** @relates TimeStamp
 *  Compare TimeStamps.
 */
inline bool operator==(const TimeStamp& t1, const TimeStamp& t2) {
  return (t1.nsec == t2.nsec);
}

/** @relates TimeStamp
 *  Compare TimeStamps.
 */
inline bool operator!=(const TimeStamp& t1, const TimeStamp& t2) {
  return (t1.nsec != t2.nsec);
}

//////////////////////////////////////////////////////////////////////

/** @relates Duration
 *  Compare Durations.
 */
inline bool operator<(const Duration& d1, const Duration& d2) {
  return (d1.nsec < d2.nsec);
}

/** @relates Duration
 *  Compare Durations.
 */
inline bool operator>(const Duration& d1, const Duration& d2) {
  return (d1.nsec > d2.nsec);
}

/** @relates Duration
 *  Compare Durations.
 */
inline bool operator<=(const Duration& d1, const Duration& d2) {
  return (d1.nsec <= d2.nsec);
}

/** @relates Duration
 *  Compare Durations.
 */
inline bool operator>=(const Duration& d1, const Duration& d2) {
  return (d1.nsec >= d2.nsec);
}

/** @relates Duration
 *  Compare Durations.
 */
inline bool operator==(const Duration& d1, const Duration& d2) {
  return (d1.nsec == d2.nsec);
}

/** @relates Duration
 *  Compare Durations.
 */
inline bool operator!=(const Duration& d1, const Duration& d2) {
  return (d1.nsec != d2.nsec);
}

/** @relates Duration
 *  Integer division of two durations.
 */
inline int64_t operator/(const Duration& d1, const Duration& d2) {
  return d1.nsec / d2.nsec;
}

/** @relates Duration
 *  Compute modular division on Durations.
 */
inline Duration operator%(const Duration& d1, const Duration& d2) {
  return Duration(d1.nsec % d2.nsec);
}

//////////////////////////////////////////////////////////////////////

inline Duration& Duration::operator+=(const Duration& d) {
  return (*this = *this + d);
}

inline Duration& Duration::operator-=(const Duration& d) {
  return (*this = *this - d);
}

inline Duration& Duration::operator*=(int s) {
  return (*this = *this * s);
}

inline Duration& Duration::operator*=(double s) {
  return (*this = *this * s);
}

inline Duration& Duration::operator/=(double s) {
  return (*this = *this / s);
}

inline TimeStamp& TimeStamp::operator+=(const Duration& d) {
  return (*this = *this + d);
}

inline TimeStamp& TimeStamp::operator-=(const Duration& d) {
  return (*this = *this - d);
}

#endif
