#include "Utils/Timer.hh"
#include "General/Exception.hh"
#include "DataInOut/ProgramOptions.hh"
#include <sstream>
#include <iomanip>

using namespace CoupledField;
using std::string;

Timer::Timer(const std::string& name, bool sub, bool start_immediately)
{
  id_ = GenerateID();
  label_ = name;
  sub_ = sub;
  nesting = true; // false means very clean code. Set global or local to true to be more robust
  if(start_immediately)
    Start();
}

Timer::Timer(int parent) : Timer("")
{
  if(parent >= 0)
    SetSub(parent);
}

Timer::Timer(const boost::shared_ptr<Timer>& parent) : Timer("")
{
  if(parent)
    SetSub(parent->id_);
}

Timer::Timer(const Timer* parent) : Timer("")
{
  if(parent != nullptr)
    SetSub(parent->id_);
} 

void Timer::SetSub(const boost::shared_ptr<Timer>& parent)
{
  if(parent)
    SetSub(parent->id_);
  else
    sub_ = true; 
}

void Timer::SetSub(int parent)
{
  parent_ = parent;
  sub_ = true;
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
  start_time_ = std::chrono::steady_clock::now();

  return true;
}

bool Timer::ResetStart()
{
  // Set timer status to running, reset accumulated time, and set start time
  bool was_running = IsRunning();
  running = 1;

  start_clock = clock();
  sum_clock   = 0;

  start_time_ = std::chrono::steady_clock::now();
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

  auto now = std::chrono::steady_clock::now();
  auto delta_t = now - start_time_;
  sum_time_ += std::chrono::duration<double>(delta_t).count();

  return true;
}


double Timer::GetWallTime() const
{
  double wall = sum_time_;
  if(running) {
    auto now = std::chrono::steady_clock::now();
    auto delta_t = now - start_time_;
    wall += std::chrono::duration<double>(delta_t).count();
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
  os << " id=\"" << id_ << "\"";
  if (label_ != "")
    os << " label=\""<< label_ << "\"";
  os << " wall-clock=\"" << Duration(GetWallTime()) << "\"";
  os << " wall=\"" << GetWallTime() << "\"";
  os << " cpu=\"" << GetCPUTime() << "\"";
  os << " calls=\"" << calls_ << "\"";
  if(progOpts->DoDetailedInfo())
    os << " max_nesting=\"" << max_nesting_ << "\"";
  if(sub_)
    os << " sub=\"true\"";
  if(parent_ >= 0)
    os << " parent=\"" << parent_ << "\"";
  os << "/>";

  return os.str();
}

const std::string Duration(const std::chrono::milliseconds period)
{
  double d = std::chrono::duration<double>(period).count();
  return Timer::Duration(d);
}

const string Timer::Duration(const std::chrono::nanoseconds period)
{
  double d = std::chrono::duration<double>(period).count();
  int h = (int) (d / 3600);
  int m = (int) ((d - (h*3600))/60);
  double s = d - h*3600 - m*60;
  
  std::stringstream ss;
  if(h > 0)
    ss << h << ":";
  if(m > 0 || h > 0)
    ss << std::setfill('0') << std::setw(2) << m << ":";
  
  ss << std::fixed << std::setprecision(2) << std::setfill('0') << std::setw(5) << s;
  if(h > 0)
    ss << " h";
  else if(m > 0)
    ss << " m";
  else
    ss << " s";

  return ss.str();
}

const std::string Timer::Duration(double seconds)
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

const std::string Timer::TimeStamp(const string& format)
{
  return Timer::TimeStamp(std::chrono::system_clock::now(), format);
}

const std::string Timer::TimeStamp(std::chrono::system_clock::time_point time_point, const std::string& format)
{
  std::time_t in_time_t = std::chrono::system_clock::to_time_t(time_point);
  std::tm tm_buf;
  
#ifdef _WIN32
  // Windows: thread-safe localtime_s
  if (localtime_s(&tm_buf, &in_time_t) != 0)
    throw Exception("localtime_s failed in Timer::TimeStamp");
#else
  // POSIX: thread-safe localtime_r
  if (localtime_r(&in_time_t, &tm_buf) == nullptr)
    throw Exception("localtime_r failed in Timer::TimeStamp");
#endif
  
  std::stringstream ss;
  ss << std::put_time(&tm_buf, format.c_str());
  return ss.str();
}




#include <atomic>

int Timer::GenerateID()
{
  static std::atomic<int> next_id = 0;
  return next_id++;
}
