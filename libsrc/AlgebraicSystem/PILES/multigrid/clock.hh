#ifndef FILE_CLOCK
#define FILE_CLOCK

/**************************************************************************/
/* File:   clock.hh                                                       */
/* Author: Wolfram Muehlhuber                                             */
/* Date:   July 1999                                                      */
/**************************************************************************/


#include <sys/time.h>
#include <unistd.h>

/*
   Local clock class: needs POSIX compliance, returns the amount of CPU time
   used since the creation of the object. This class is needed for long term
   simulations, as the system function clock() wraps around after a rather
   short amount of time (usually 36 min)
*/


class Clock 
{
  /// real time when the clock object was created
  timeval starttime;

public:
  /// Constructor
  Clock ();
  /// Destructor
  virtual ~Clock ();

  /// get time in seconds since the object was created, this is wall clock time
  double GetTime () const;
};


class CPUClock 
{
  /// CPU time when the clock object was created (user and system time)
  timeval userstarttime;
  timeval sysstarttime;

public:
  /// Constructor
  CPUClock ();
  /// Destructor
  virtual ~CPUClock ();

  /// get CPU time in seconds since the object was created (user and system time)
  double GetTime () const;
  /// get user time in seconds since the object was created
  double GetUserTime () const;
  /// get system time in seconds since the object was created
  double GetSystemTime () const;
};


#endif // FILE_CLOCK
