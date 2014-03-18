#include <stdio.h>
#include <stdlib.h>

#include <ilupack_config.h>

#ifdef WIN32
#include <windows.h>
#else
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#include <sys/resource.h>
#else
#ifdef HAVE_SYS_TIMES_H
#include <sys/times.h>
#include <unistd.h>
#else
#error "No timing headers found!"
#endif
#endif // HAVE_SYS_TIME_H
#endif // WIN32

#include <ilupack.h>
#include <ilupackmacros.h>

float evaluate_time( float *user, float *sys  )
{
#ifdef WIN32
    FILETIME creation_time, exit_time, kernel_time, user_time;
#else
#ifdef HAVE_SYS_TIME_H
#else
#ifdef HAVE_SYS_TIMES_H
  /*
     code written by Georg Fuss, September 2003
  */
    
    clock_t ttt, time_u, time_s;
    long    div;
    struct tms buf;
#endif
#endif
#endif

    float t0, t1;

#ifdef WIN32
    // http://www.philosophicalgeek.com/2009/01/03/determine-cpu-usage-of-current-process-c-and-c/
    // The implementation has been taken from the etime function of the G95 compiler
    // www.g95.org
    GetProcessTimes(GetCurrentProcess(), &creation_time, &exit_time,
                    &kernel_time, &user_time);

    t0 = 429.4967296 * ((unsigned) user_time.dwHighDateTime)
        + 1.0e-7 * ((unsigned) user_time.dwLowDateTime);

    t1 = 429.4967296 * ((unsigned) kernel_time.dwHighDateTime)
        + 1.0e-7 * ((unsigned) kernel_time.dwLowDateTime);
#else

#ifdef HAVE_SYS_TIME_H
    // The implementation has been taken from the etime function of the G95 compiler
    // www.g95.org
    struct rusage ru;

    if (getrusage(RUSAGE_SELF, &ru) < 0) {
        t0 = -1.0;
        t1 = -1.0;

    } else {
        t0 = ru.ru_utime.tv_sec + 1.0e-6 * ru.ru_utime.tv_usec;
        t1 = ru.ru_stime.tv_sec + 1.0e-6 * ru.ru_stime.tv_usec;
    }
#else
#ifdef HAVE_SYS_TIMES_H
    // Old implementation from ILUPACK
    div = sysconf( _SC_CLK_TCK );
    
    ttt = times( &buf );
    time_u = buf.tms_utime;
    time_s = buf.tms_stime;

    t0 = (float)time_u / (float)div;
    t1  = (float)time_s / (float)div;
#endif    
#endif // HAVE_SYS_TIME_H

#endif // WIN32

    if (user != NULL) {
        *user = t0;
    }

    if (sys != NULL) {
        *sys = t1;
    }

    return t0+t1;

}
