#ifndef TIMER_HH_
#define TIMER_HH_

#include <ctime>
#include <string>
#include <chrono>
#include <cassert>
#include <boost/shared_ptr.hpp>

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

  /** parent is not stored, just it's id_ is read. */
  Timer(const boost::shared_ptr<Timer>& parent);
  Timer(const Timer* parent);
  Timer(int parent);

  /** To set 'sub' when created as ParamNode::AsTimer().
   * A sub-timer is not counted by performance.py. E. the parent "assemble" has sub timers 
   * which are not counted, otherwise the sum of timers would be larger than the total runtime. 
   * @param parent of parent available it's id. But still sets to sub-timer */
  void SetSub(const boost::shared_ptr<Timer>& parent);
  void SetSub(int parent = -1); 

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
  
  /** Prints the time in human readable format to specified stream */
  void PrintTime(std::ostream& stream);

  /** convert a std::chrono::steady_clock time interval
   *  to a formatted string suitable for command line output. 
   *  With hour, minutes and second handling. 
   *  On macOS with ARM64 the shortest time period is 4.2e-8 seconds. 
   *  
   *  Usage if no Timer object is used:
   *  1. auto start_time = std::chrono::steady_clock::now();
   *  2. cout << Duration(std::chrono::steady_clock::now() - start_time);*/
  static const std::string Duration(double seconds);
  static const std::string Duration(const std::chrono::milliseconds period);
  static const std::string Duration(const std::chrono::nanoseconds period);

  /** today in format 11-Feb-26 */
  static const std::string Today() { 
    return Timer::TimeStamp("%d-%b-%y"); 
  }
  
  /** easy human readable time stamp in format 11-Feb-26 01:30:39 */
  static const std::string TimeStamp() {
    return Timer::TimeStamp("%d-%b-%y %H:%M:%S"); 
  }

  /** easy human readable time stamp in format 2026-Feb-11 01:30:39 */
  static const std::string TimeStampYYYYmmDD() {
    return Timer::TimeStamp("%Y-%b-%d %H:%M:%S"); 
  }

  /** generic format for now(). 
   * E.g. "%Y-%b-%d %H:%M:%S" for "YYYY-mmm-DD HH:MM:SS"
   *      "%d-%b-%y %H:%M:%S" for "11-Feb-26 02:03:39"*/
  static const std::string TimeStamp(const std::string& format);

  /** generic format for time_point */
  static const std::string TimeStamp(std::chrono::system_clock::time_point, const std::string& format = "%Y-%b-%d %H:%M:%S");

  private:

  /** thread safe id for the timer starting with 0 */
  static int GenerateID();

  /** our unique timer id for xml output and reconstruct parent / sub-timer graphs */
  int id_ = -1; 

  /** parent id if sub with known parent. Partly redundant with sub_.
   * @see sub_ */
  int parent_ = -1;

  /** To help performance.py that we are a sub timer and our duration shall not be summed up.
   * @see parent_ */
  bool sub_ = false;

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
  std::chrono::time_point<std::chrono::steady_clock> start_time_;
  double sum_time_ = 0;

  std::string label_;

}; // class timer


} // end of namespace



#endif /* TIMER_HH_ */
