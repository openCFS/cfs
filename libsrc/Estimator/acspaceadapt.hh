#ifndef FILE_ACOUSTICSPACEERROR_2002
#define FILE_ACOUSTICSPACEERROR_2002

#include <vector>

#include "tools.hh"
#include "acoustic2dPDE.hh"
#include "spaceerror.hh"

namespace CoupledField
{

//! Class where we implement space error estimator
class AcousticSpaceErrorEstimator:virtual public SpaceErrorEstimator
{
public:
  //! constructor with pointer to PDE
  AcousticSpaceErrorEstimator(BasePDE * , Grid * );

  //! Deconstructor
  virtual ~AcousticSpaceErrorEstimator();

  //!  return true, if error is more than tolerance
  Boolean TestError();

  //!
  void CalcError();

  //!
  void RefineMesh();

  //!
  Boolean TestLocError(grd::Element * t);

private:
  //!
  Double tol_, theta_;

  //!
  Double error_, maxenergynormsol_;

  //!
  void DefineRefinedElems(const Integer level, std::vector<Integer> & elems2ref);

  //!
  Double CalcLocError(const Integer iElem);  

};

inline AcousticSpaceErrorEstimator::~AcousticSpaceErrorEstimator()
{
#ifdef TRACE
  (*trace) << "entering AcousticSpaceErrorEstimation::~AcousticSpaceErrorEstimator" << std::endl;
#endif 
;
}

// class SetRefFlagAcoustic
// {
// public:

//   SetRefFlag(AcousticSpaceErrorEstimator * aptASEE){ ptSpaceError_=apASEE;}
//   ~SetRefFlag(){;}

//  void operator() (grd::Element * t)
//  {
//   t->markForRefinement();
//  }

// };

} // end of namespace

#endif
