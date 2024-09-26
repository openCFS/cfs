#include "Utils/Timer.hh"
#include "General/Exception.hh"
#include "DataInOut/ProgramOptions.hh"
#include <sstream>

using namespace CoupledField;
using std::string;

Timer::Timer(const std::string& name, bool sub, bool start_immediately)
{
  label_ = name;
  sub_ = sub;
  nesting = true; // false means very clean code. Set global or local to true to be more robust
  if(start_immediately)
    Start();
}


bool Timer::Start()
{
  if(IsRunning() && !nesting)
    throw Exception("attempt to start already running timer " + label_ + " with nesting off");

  calls_++;

  // Set timer status to running and set the start time
  running++;
  max_nesting_ = std::max(running, max_nesting_);

  if(running > 1)
    return false;

  start_clock = clock();

  // https://stackoverflow.com/questions/6734375/get-current-time-in-milliseconds-using-c-and-boost
  start_time_ = boost::posix_time::microsec_clock::local_time();

  return true;
}

bool Timer::ResetStart()
{
  // Set timer status to running, reset accumulated time, and set start time
  bool was_running = IsRunning();
  running = 1;

  start_clock = clock();
  sum_clock   = 0;

  start_time_ = boost::posix_time::microsec_clock::local_time();
  sum_time_ =   0;

  calls_ = 1;
  return !was_running;
}

bool Timer::Stop()
{
  if(!IsRunning() && !nesting)
    throw Exception("attempt to stop not running timer with nesting off");

  if(running > 0) // don't make negative
    running--;

  // if nesting is still to high, we don't stop
  if(IsRunning())
    return false;

  // Compute accumulated running time and set timer status to not running
  sum_clock += clock() - start_clock;

  boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
  boost::posix_time::time_duration delta_t = now - start_time_;
  sum_time_ += delta_t.total_milliseconds() / 1000.0;

  return true;
}


double Timer::GetWallTime() const
{
  double wall = sum_time_;
  if(running) {
    boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
    boost::posix_time::time_duration delta_t = now - start_time_;
    wall += delta_t.total_milliseconds() / 1000.0;
  }
  return wall;
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
  if (label_ != "")
    os << " label=\""<< label_ << "\"";
  os << " wall-clock=\"" << GetTimeString(GetWallTime()) << "\"";
  os << " wall=\"" << GetWallTime() << "\"";
  os << " cpu=\"" << GetCPUTime() << "\"";
  os << " calls=\"" << calls_ << "\"";
  if(progOpts->DoDetailedInfo())
    os << " max_nesting=\"" << max_nesting_ << "\"";
  if (sub_)
    os << " sub=\"true\"";
  os << "/>";

  return os.str();
}

const string Timer::GetTimeString(const boost::posix_time::time_duration period)
{ 
  if(boost::date_time::micro != period.resolution() &&
      boost::date_time::sec   != period.resolution())
  {
    return "NA";
  }
  
  string time_output(to_simple_string(period));
  string suffix(" h");
  int max_length(11); // for micro, e. g. 01:02:03.04
  if(boost::date_time::sec == period.resolution()) max_length -= 3;
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

const std::string Timer::GetTimeString(double seconds)
{
  int h = (int) (seconds/3600);
  int m = (int) ((seconds - (h*3600))/60);
  double s = seconds - (h*3600) - (m*60);

  std::stringstream ss;
  if(h > 0)
    ss << h << "h ";
  if(m > 0)
    ss << m << "m ";
  ss.precision(4);
  ss << s << "s";
  return ss.str();
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


