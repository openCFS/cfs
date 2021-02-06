#ifndef TIMER_HH_
#define TIMER_HH_

#include <ctime>
#include <string>
#include <boost/chrono.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace CoupledField
{

/** The Timer class is based on the code timer.hh from Ken Wilder (http://sites.google.com/site/jivsoft/Home/timer)
 * It implements a timer object which sums up time intervals.
 *
 * As it can be that a timed function calls another function which wants to start the same timer again, we keep the
 * current level, otherwise a stop would of the nested function would stop the timer and the calling functions does not expect
 * this. This can be controlled via nesting. Basically it shall be better to work w/o nesting
 * and  you can use Timer to detect and fix it.
 *
 * A convenient way to create a timer is via ParamNode::AsTimer().
 *
 * performance.py parses the .info.xml and presents the timer statistics
 */
class Timer
{
 public:
  /** 'running' is initially 0 als long as we don't set start_immediately
       A Timer needs to be explicitly started vi Start() or ResetStart() if not start_immediately is set
       @param sub means that performance.py shall not use the time for the 'not_measured' calculation because it is just a finer
       time for an action already measured. */
  Timer(const std::string& name = "", bool sub = false, bool start_immediately = false);

  /** Start a timer. Handles nesting by throwing an exception it already running and not nesting.
   * Increments calls only when really starting the timer.
   * @return false if was already running (nothing is changed), true for a real start */
  bool Start();

  /** Stop the timer, can be restarted via Start. Decrements running and stops only when 0.
   * Throws and exception if not running and not nesting set.
   * @return false if it was not running before */
  bool Stop();

  /** Turn the timer off and start it again from 0.
   * @return false if it was not running */
  bool ResetStart();

  bool IsRunning() const { return running > 0; }

  /** current nesting level. 0 is not running */
  int GetLevel() const { assert(running >= 0); return running; }

  /** is nesting allowed? */
  bool IsNesting() const { return nesting; }

  void SetNesting(bool val) { nesting = val; }

  void SetLabel(const std::string& name) {
    label_ = name;
  }

  /** to set sub when created as ParamNode::AsTimer() */
  void SetSub() {
    sub_ = true;
  }

  /** The number of Start() calls since construction or the last ResetStart()+1.
   * Only real starting of timers are counted
   * ElapseTime() / GetCalls() is the average runtime. */
  int GetCalls() const { return calls_; }

  /** elapsed seconds. */
  double GetWallTime() const;

  /** elapsed CPU seconds, might be higher than wall time if some parallelization
   * (e.g. OpenMP MKL) is used */
  double GetCPUTime() const;

  /** This writes writes a xml element with time attributes for ParamNode::ToFile().
   * @param name is used for "<name wall='xx' .../>
   *        This is the name of the param element in in->Get("timer")->SetValue(timer) */
  std::string ToXMLFormat(const std::string& name) const;
  
  /** static function; converts a boost time interval
   *  to a formatted string suitable for command line output 
   *  
   *  Usage:
   *  1. ptime start_time = microsec_clock::local_time();
   *    or
   *     ptime start_time = second_clock::local_time(); 
   *  2. cout << GetTimeString(second (or microsec)_clock::local_time() - start_time);*/
  static const std::string GetTimeString(const boost::posix_time::time_duration period);

  static const std::string GetTimeString(double seconds);

  /** Prints the time in human readable format to specified stream
   *  @param (in) stream output stream
   */
  void PrintTime(std::ostream & stream);

 private:
  /** The number of Start() calls since construction or the last ResetStart()+1 */
  int calls_  = 0;

  /** our nested running counter. */
  int running = 0;

  /** is nesting allowed? Set in constructor to experiment w/o recompiling all cfs */
  bool nesting;

  /** this is for debug purpose */
  int max_nesting_ = 0;

  /** clock is cpu time */
  clock_t start_clock = 0;
  clock_t sum_clock = 0;

  /** wall time: boost time, as time_t is only in seconds */
  boost::posix_time::ptime start_time_;
  double sum_time_ = 0;

  std::string label_;
  /** sub and generic counters shall not be summed up by performance.py */
  bool sub_ = false;

}; // class timer


} // end of namespace



#endif /* TIMER_HH_ */
