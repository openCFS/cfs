/***************************************************************************
    File        : Chronometer.cpp
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Tue Sep 10 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#include "Chronometer.h"

namespace grd {

// Description
//
void
Chronometer::print(char* title)
{
  double cpu_time, wallclok_time, clock_tick;

  cpu_time = ((double)(after.tms_utime - before.tms_utime))/((double)CLK_TCK);
  wallclok_time = ((double) (endtime - startime)) / ((double) CLK_TCK);
  clock_tick    = (endtime - startime);

  cerr << title << '\n';
  cerr << "---------------------------------------------------------------------" << '\n'
       << "|    CPU time (sec)   |  Wall clock time (sec) |    Clock Ticks     |" << '\n'
       << "---------------------------------------------------------------------" << '\n'
       << setw(21) << '|'
       << setprecision(6)
       << cpu_time
       << setw(24) << '|'
       << setprecision(6)
       << wallclok_time
       << setw(20) << '|'
       << setprecision(6)
       << clock_tick
       << '|'
       << '\n';
}

}