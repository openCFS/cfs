#ifndef FILE_ACOUSTICTIMEERROR_2001
#define FILE_ACOUSTICTIMEERROR_2001

#include "Utils/tools.hh"
#include "PDE/acoustic2dPDE.hh"
#include "timeerror.hh"

namespace CoupledField
{

//! Class where we implement time error estimator
class AcousticTimeErrorEstimator: virtual public TimeErrorEstimator
{
public:
  //! constructor with pointer to PDE
  AcousticTimeErrorEstimator(BasePDE * aptpde);

  //! Deconstructor
  virtual ~AcousticTimeErrorEstimator(){;}

  //!  return true, if error is more than tolerance
  virtual Boolean TestError() { return FALSE; }

  //! 
  virtual void ChangeStep(Double &);

private:
  
  //!
  void CalcError(){;}

};

} // end of namespace

#endif
