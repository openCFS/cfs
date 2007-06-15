#ifndef VEC_STAT_HH_
#define VEC_STAT_HH_

#include "Utils/StdVector.hh"

namespace CoupledField
{
  /** This class does some statistics about vectors. Norms, ... */
  class VecStat
  {
    public:
      /** empty constructor for class attribute */
      VecStat()
      {
        initialized_ = false;
      }
    
      VecStat(StdVector<double> vec)
      {
        Calc(vec.GetPointer(), vec.GetSize());
      }
    
      VecStat(const double* data, unsigned int length)
      {
        Calc(data, length);
      }

      /** this actually does the job and sets all attributes (incl. initialized_) */
      void Calc(const double* data, unsigned int length);

      /** Prints a shoer summery! */
      std::string ToString();
      
      /** This is a version for logging */
      static std::string ToString(const double* data, unsigned int length)
      {
        VecStat vs(data, length);
        return vs.ToString();
      }
      
      double GetAverage();
      double GetMin();
      double GetMax();
      double GetStdDev();
      double GetNormL2();
      
      bool IsInitialized() { return initialized_; }
      
    private:
      

      double average_;
      double min_;
      double max_;
      double stdDev_;
      double normL2_;
      
      bool initialized_;  
  };
}

#endif /* VEC_STAT_HH_*/

