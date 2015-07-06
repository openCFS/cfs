#include "Utils/Timer.hh"
#include <sstream>

using namespace CoupledField;
using namespace boost;
using posix_time::ptime;
using posix_time::second_clock;
using posix_time::microsec_clock;
using posix_time::time_duration;
using std::string;

Timer::Timer() :
  calls_(0),
  running(false),
  start_clock(0),
  start_time(0),
  sum_time(0),
  sum_clock(0)
{
}


void Timer::Start()
{
  // Return immediately if the timer is already running
  if (running) return;

  calls_++;

  // Set timer status to running and set the start time
  running = true;
  start_clock = clock();
  start_time = time(0);

}

void Timer::ResetStart()
{
  // Set timer status to running, reset accumulated time, and set start time
  running = true;
  sum_time    = 0;
  sum_clock   = 0;
  start_clock = clock();
  start_time  = time(0);

  calls_ = 1;
}

void Timer::Stop()
{
  // Compute accumulated running time and set timer status to not running
  if (running)
  {
    sum_time  += time(NULL) - start_time;
    sum_clock += clock() - start_clock;
  }
  running = false;
}


double Timer::GetWallTime() const
{
  return sum_time + (running ? time(NULL) - start_time : 0);
}

double Timer::GetCPUTime() const
{
  clock_t total = sum_clock + (running ? clock() - start_clock : 0);
  return total / (1.0 * CLOCKS_PER_SEC);
}


string Timer::ToXMLFormat(const string& name) const
{
  std::ostringstream os;

  os << "<" << name;
  os << " wall=\"" << GetWallTime() << "\"";
  os << " cpu=\"" << GetCPUTime() << "\"";
  os << " calls=\"" << calls_ << "\"";
  os << "/>";

  return os.str();
}

const string Timer::GetTimeString(const time_duration period)
{ 
  if(date_time::micro != period.resolution() &&
     date_time::sec   != period.resolution())
  {
    return "NA";
  }
  
  string time_output(to_simple_string(period));
  string suffix(" h");
  int max_length(11); // for micro, e. g. 01:02:03.04
  if(date_time::sec == period.resolution()) max_length -= 3;
  if(time_output.substr(0, 2) == "00") // remove 0 hours
  {
    time_output = time_output.substr(3);
    max_length -= 3;
    suffix = " m";
  }
  else // we have positive hours, so we need the whole string
    return time_output.substr(0, max_length) + suffix; // cut off microseconds
  
  if(time_output.substr(0, 2) == "00") // remove 0 minutes
  {
    time_output = time_output.substr(3);
    max_length -= 3;
    suffix = " s";
  }
  return time_output.substr(0, max_length) + suffix; // cut off microseconds
}

void Timer::PrintTime(std::ostream & stream){
  if(this->running){
    stream << "Timer is still running! ";
    stream << ">> Current time: wall clock: '";
  }else{
    stream << ">> Elapsed time: wall clock: '";
  }


  const int walltime((int) this->GetWallTime());
  const int cputime((int) this->GetCPUTime());

  if(walltime > 120)
  {
    const int wallmin((int) (walltime / 60.0));
    const int cpumin((int) (cputime / 60.0));
    if(wallmin > 60)
    {
      stream << wallmin / 60 << "h " << (wallmin % 60)
           << "m' CPU time: '" << cpumin / 60 << "h " << (cpumin % 60) << "m'";
    }
    else
    {
      stream << wallmin << "m " << (walltime % 60)
           << "s' CPU time: '" << cpumin << "m " << (cputime % 60) << "s'";
    }
  }
  else
  {
    stream << walltime << "s' CPU time: '"
         << cputime << "s'";
  }

  stream << std::endl << std::endl;
}


