#ifndef FILE_TIMEERROR_2001
#define FILE_TIMEERROR_2001

#include "tools.hh"

namespace CoupledField
{

//! Class where we implement time error estimator
class TimeErrorEstimator
{
public:
  //!
  TimeErrorEstimator();

  //! Deconstructor
  virtual ~TimeErrorEstimator();

  //! 
  virtual Boolean TestError();

  //!
  void CalcError();

  //! 
  void ChangeStep();
 
};

} // end of namespace

#endif
