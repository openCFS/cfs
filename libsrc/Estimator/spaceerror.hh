#ifndef FILE_SPACEERROR_2002
#define FILE_SPACEERROR_2002

#include "tools.hh"
#include "spaceerror.hh"

namespace CoupledField
{

class SpaceErrorEstimator
{
public:
  //! constructor with pointer to PDE
  SpaceErrorEstimator(BasePDE * aptPDE, Grid * aptGrid)
  { ptPDE_=aptPDE; ptGrid_=aptGrid;}

  //! Deconstructor
  virtual ~SpaceErrorEstimator(){;}

  //!  return true, if error is more than tolerance
  virtual Boolean TestError()=0;

  //!
  virtual void CalcError()=0;

  //!
  virtual void RefineMesh()=0;

protected:
  //!
  BasePDE * ptPDE_;  

  //!
  Grid * ptGrid_;  

};

} // end of namespace

#endif
