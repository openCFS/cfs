// - C++ -
/***************************************************************************
    File        : Chronometer.h
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Tue Sep 10 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef CHRONOMETER_H
#define CHRONOMETER_H


/*
 * The Class "Time"
 * -----------------
 *
 * Description  :   This class is based on the system call times() and is
 *                  designed to measure cpu and wallclock time of function
 *                  calls.
 *
 * Date:            4.11.1993
 *
 */

#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <string.h>


namespace grd {

class Chronometer {
public:
  // Constructors
  Chronometer() : startime(0),endtime(0) { }

  // Destructor
  virtual ~Chronometer() {}

  // Member Functions
  void start(void) {
      startime = times(&before);
  }
  void stop(void) {
      endtime = times(&after);
  }
  void print(char* title);

private:
  struct tms before, after;
  clock_t utime, stime;
  clock_t startime, endtime;
};

} // namespace grd
#endif /* CHRONOMETER_H */