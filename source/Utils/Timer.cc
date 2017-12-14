#include "Utils/Timer.hh"
#include <sstream>

using namespace CoupledField;
using std::string;

Timer::Timer(const std::string& name, bool sub) :
  calls_(0),
  running(false),
  start_clock(0),
  sum_clock(0),
  label_(name),
  sub_(sub)
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

  // https://stackoverflow.com/questions/6734375/get-current-time-in-milliseconds-using-c-and-boost
  start_time_ = boost::posix_time::microsec_clock::local_time();

}

void Timer::ResetStart()
{
  // Set timer status to running, reset accumulated time, and set start time
  running = true;

  start_clock = clock();
  sum_clock   = 0;

  start_time_ = boost::posix_time::microsec_clock::local_time();
  sum_time_ =   0;

  calls_ = 1;
}

void Timer::Stop()
{
  // Compute accumulated running time and set timer status to not running
  if(running)
  {
    sum_clock += clock() - start_clock;

    boost::posix_time::ptime now = boost::posix_time::microsec_clock::local_time();
    boost::posix_time::time_duration delta_t = now - start_time_;
    sum_time_ += delta_t.total_milliseconds() / 1000.0;
  }
  running = false;
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
  os << " wall=\"" << GetWallTime() << "\"";
  os << " cpu=\"" << GetCPUTime() << "\"";
  os << " calls=\"" << calls_ << "\"";
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


