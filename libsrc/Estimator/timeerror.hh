#ifndef FILE_TIMEERROR_2001
#define FILE_TIMEERROR_2001

#include "tools.hh"
#include "basepde.hh"

namespace CoupledField
{

//! Class where we implement time error estimator
class TimeErrorEstimator
{
public:
  //! constructor with pointer to PDE
  TimeErrorEstimator(BasePDE* aptpde);

  //! Deconstructor
  virtual ~TimeErrorEstimator(){;}

  //!  return true, if error is more than tolerance
  virtual Boolean TestError() { return FALSE; }

  //! 
  virtual void ChangeStep(Double &){;}

protected:
  
   //!
   BasePDE* ptPDE_;

   //!
   Double relativeerror; 
};

} // end of namespace

#endif
