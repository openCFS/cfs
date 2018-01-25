#ifndef TIMER_HH_
#define TIMER_HH_

#include <ctime>
#include <string>
#include "boost/date_time/posix_time/posix_time.hpp"
#include "General/defs.hh"

namespace CoupledField
{

/** The Timer class is based on the code timer.hh from Ken Wilder (http://sites.google.com/site/jivsoft/Home/timer)
 * It implements a timer object which sums up time intervals.
 */
class Timer
{
 public:
  /** 'running' is initially false.
       A Timer needs to be explicitly started vi Start() or ResetStart()
       @param sub means that performance.py shall not use the time for the 'not_measured' calculation because it is just a finer
       time for an action already measued. */
  Timer(const std::string& name = "", bool sub = false);

  bool IsRunning() const { return running; }

  /** Start or resume a timer.
    If it is already running, let it continue running. */
  void Start();

  /** Turn the timer off and start it again from 0 */
  void ResetStart();

  /** Stop the timer, ca be restarted via Start */
  void Stop();

  void SetLabel(const std::string& name) {
    label_ = name;
  }

  /** The number of Start() calls since construction or the last ResetStart()+1.
   * ElapseTime() / GetCalls() is the average runtime. */
  int GetCalls() const { return calls_; }

  /** elapsed seconds. */
  Double GetWallTime() const;

  /** elapsed CPU seconds, might be higher than wall time if some parallelization
   * (e.g. OpenMP MKL) is used */
  Double GetCPUTime() const;

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

  /** Prints the time in human readable format to specified stream
   *  @param (in) stream output stream
   */
  void PrintTime(std::ostream & stream);

 private:
  /** The number of Start() calls since construction or the last ResetStart()+1 */
  int    calls_;
  bool   running;
  clock_t start_clock;
  time_t  start_time;
  Double sum_time;
  clock_t sum_clock;
  std::string label_;
  bool sub_;

}; // class timer


} // end of namespace



#endif /* TIMER_HH_ */
