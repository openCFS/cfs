#include <iostream.h>

#include "clock.hh"
#include <sys/resource.h>

// the following define is needed for SUN
#ifndef MSEC_PER_SEC
#define MSEC_PER_SEC 1e6
#endif

/* ********************************** Clock ******************************** */

Clock :: Clock ()
{
  starttime.tv_sec = 0;    // time in seconds
  starttime.tv_usec = 0;   // time in nano seconds

  if (gettimeofday (& starttime, NULL))
    cerr << "Clock :: Clock : can not initialize timer" << endl;
}


Clock :: ~Clock ()
{
  ;
}


double Clock :: GetTime () const
{
  timeval timenow;

  if (gettimeofday (& timenow, NULL))
    {
      cerr << "Clock :: GetTime : can not get current time" << endl;
      return 0;
    }

  timenow.tv_sec -= starttime.tv_sec;
  timenow.tv_usec -= starttime.tv_usec;

  return timenow.tv_sec + ( (double) timenow.tv_usec) / MSEC_PER_SEC;
}


/* ******************************** CPUClock ******************************* */

CPUClock :: CPUClock ()
{
  rusage rusage; // resource sructure

  if (getrusage (RUSAGE_SELF, & rusage))
    {
      cerr << "CPUClock :: CPUClock : can not initialize timer" << endl;
      userstarttime.tv_sec = 0;   // time in seconds
      userstarttime.tv_usec = 0;  // time in micro seconds
      sysstarttime.tv_sec = 0;   
      sysstarttime.tv_usec = 0;
    }
  else
    {
      userstarttime = rusage.ru_utime;
      sysstarttime = rusage.ru_stime;
    }
}


CPUClock :: ~CPUClock ()
{
  ;
}


double CPUClock :: GetTime () const
{
  rusage rusage;

  if (getrusage (RUSAGE_SELF, & rusage))
    {
      cerr << "CPUClock :: GetTime : can not get current time" << endl;
      return 0;
    }

  rusage.ru_utime.tv_sec -= userstarttime.tv_sec;
  rusage.ru_utime.tv_usec -= userstarttime.tv_usec;
  rusage.ru_stime.tv_sec -= sysstarttime.tv_sec;
  rusage.ru_stime.tv_usec -= sysstarttime.tv_usec;

  return rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec +
    ( rusage.ru_utime.tv_usec + rusage.ru_stime.tv_usec ) / 1e6;
}


double CPUClock :: GetUserTime () const
{
  rusage rusage;

  if (getrusage (RUSAGE_SELF, & rusage))
    {
      cerr << "CPUClock :: GetUserTime : can not get current time" << endl;
      return 0;
    }

  rusage.ru_utime.tv_sec -= userstarttime.tv_sec;
  rusage.ru_utime.tv_usec -= userstarttime.tv_usec;

  return rusage.ru_utime.tv_sec + rusage.ru_utime.tv_usec / 1e6;
}


double CPUClock :: GetSystemTime () const
{
  rusage rusage;

  if (getrusage (RUSAGE_SELF, & rusage))
    {
      cerr << "CPUClock :: GetSystemTime : can not get current time" << endl;
      return 0;
    }

  rusage.ru_stime.tv_sec -= sysstarttime.tv_sec;
  rusage.ru_stime.tv_usec -= sysstarttime.tv_usec;

  return rusage.ru_stime.tv_sec + rusage.ru_stime.tv_usec / 1e6;
}
