#ifndef FILE_ACOUSTICTIMEERROR_2001
#define FILE_ACOUSTICTIMEERROR_2001

#include "tools.hh"
#include "acoustic2dPDE.hh"
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
  virtual ~AcousticTimeErrorEstimator();

  //!  return true, if error is more than tolerance
  Boolean TestError(const Double dt);

  //! 
  void ChangeStep(Double &);

private:
  
  //!
  void CalcError(const Double dt);

  //!
  Double relativeerror_;

  //!
  Vector<Double> thirddersol_;

  //!
  Double tol_, beta_, theta_;

  //!
  Integer numrepeat_, counter_;
};

inline AcousticTimeErrorEstimator::~AcousticTimeErrorEstimator()
{
#ifdef TRACE
  (*trace) << "entering AcousticTimeErrorEstimation::~AcousticTimeErrorEstimator" << std::endl;
#endif 
;
}

} // end of namespace

#endif
