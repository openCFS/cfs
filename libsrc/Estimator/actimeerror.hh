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

  //!
  void CalcThirdDer();  

  //!
  void CalcError(const Double dt);

private:
  //!
  Double maxdt_,mindt_;

  //!
  Double maxrelativeerror_, maxnormsol_, maxnorml2sol_;

  //!
  Vector<Double> thirddersoldt_;

  //!
  Double tol_, beta_, theta_;

  //!
  Integer numrepeat_, counter_;

  //!
  Boolean Calc3DerFromEquation_;

  //!
  Double err1_,err2_;  

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
