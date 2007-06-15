#include "Utils/VecStat.hh"
#include "Utils/tools.hh"
#include "General/exception.hh"

using namespace CoupledField;

void VecStat::Calc(const double* data, unsigned int size)
{
  initialized_ = true;
  
  min_ = +1e19;
  max_ = -1e19;

  if(size == 0) EXCEPTION("empty vector");
  if(data == NULL) EXCEPTION("NULL data");
  
  for(unsigned int i = 0; i < size; i++)
  {
    min_ = std::min(min_, *(data + i));
    max_ = std::max(max_, *(data + i));
  }

  // this could be clearly faster if we reimplement.
  // But as long as valgrind does not meassure it 
  // significantly we avoid copy & past mess 
  normL2_ = NormL2(data, size);
  average_ = Average(data, size);
  stdDev_ = StandardDeviation(data, size);
}

std::string VecStat::ToString()
{
  std::ostringstream os;
  
  os << "avg=" << average_ << " stdDev=" << stdDev_ << " L2=" << normL2_ << " min=" 
     << min_ << " max=" << max_;
  return os.str();
}
      
double VecStat::GetAverage()
{
  if(!initialized_) EXCEPTION("not initialized");
  return average_;
}

double VecStat::GetMin()
{
  if(!initialized_) EXCEPTION("not initialized");
  return min_;
}

double VecStat::GetMax()
{
  if(!initialized_) EXCEPTION("not initialized");
  return max_;
}

double VecStat::GetNormL2()
{
  if(!initialized_) EXCEPTION("not initialized");
  return normL2_;
}

double VecStat::GetStdDev()
{
  if(!initialized_) EXCEPTION("not initialized");
  return stdDev_;
}
