#ifndef FILE_ELECSTERRORESTIMATOR_2002
#define FILE_ELECSTERRORESTIMATOR_2002

#include <vector>

#include "Utils/tools.hh"
#include "PDE/elecst2dPDE.hh"
#include "spaceerror.hh"

namespace CoupledField
{

//! Class where we implement space error estimator
class Elecst2dErrEstimator:virtual public SpaceErrorEstimator
{
public:
  //! constructor with pointer to PDE
  Elecst2dErrEstimator(BasePDE * , Grid * );

  //! Deconstructor
  virtual ~Elecst2dErrEstimator();

  //!  return true, if error is more than tolerance
  Boolean TestError();

  //!
  void CalcError();

  //!
  void RefineMesh();

  //!
#ifdef ADAPTGRID
  Boolean TestLocError(grd::Element * t);
#endif

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

inline Elecst2dErrEstimator::~Elecst2dErrEstimator()
{
#ifdef TRACE
  (*trace) << "entering Elecst2dErrEstimation::~Elecst2dErrEstimator" << std::endl;
#endif 
;
}

} // end of namespace

#endif
