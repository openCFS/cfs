#ifndef FILE_SPACEERROR_2002
#define FILE_SPACEERROR_2002

#include "tools.hh"
#include "basepde.hh"

#ifdef ADAPTGRID
#include "Vertex.h"
#include "Edge.h"
#include "Element.h"
#endif

namespace CoupledField
{

/// class for space error estimation
template<Integer dim>
class SpaceErrorEstimator
{
public:
  //! constructor with pointer to PDE
  SpaceErrorEstimator(Grid * aptGrid);

  //! Deconstructor
  virtual ~SpaceErrorEstimator();

  //! Init pointer to PDE
  virtual void Init(BasePDE * aptPDE){ ptPDE_=aptPDE;}

  //!
  void KellyError(std::vector<std::string> &sbdoms, const Integer level);

  //! calculation of error map
 void CalcErrorMap(const Vector<Double> * sol, std::vector<std::string> & subdoms,
 Grid * ptgrid, Vector<Double> & errorMap, Vector<Double> & gradSPRElemL2norm);

  //! recovery procedure
  void RecoveryProcedure4ElemsPatch(const std::vector<Elem*> &Elems, Grid * ptgrid, const Vector<Double> & sol, const Integer aComponent, Vector<Double> &result,std::vector<Integer>&locations);

  //! recovery procedure for 3d
  void RecoveryProcedure4ElemsPatch3d(const std::vector<Elem*> &Elems, Grid * ptgrid, const Vector<Double> & sol, const Integer aComponent, Vector<Double> &result,std::vector<Integer>&locations);

  //! only for 1 element
  void RecoveryProcedure4Elem(Elem * aptElem, Grid * ptgrid, const Vector<Double> & sol, const Integer aComponent, Vector<Double> &result);

protected:
  //!
  BasePDE * ptPDE_;  

  //!
  Grid * ptGrid_;  

  //! calculation error for one element
  void CalcErrorForElem(const Elem* elem, const Vector<Double>* SPRgrad,  Double & error, Double & normGradSPR, const Vector<Double>* sol, Grid * ptgrid);

private:  

  //! auxiliary function for calculation of Kelly error
  Double SpaceErrorEstimator::IntegralOverRegularFace_KellyError();
  
  //! calculation of RHS for the recovery procedure
  void ComputeRHS4RecoverySol(Grid * ptgrid,  const std::vector<std::string> & subdoms, const Integer as_sysid, AbstractAlgebraicSys * ptalgsys, const Vector<Double>&sol);

  //! auxialary procedure for calculating number of nodes in patch
  /* direct calculating is implemented */
  Integer CalcNumberOfNodesInPatch(const std::vector<Elem*> & patch);

};

template <Integer dim>
inline SpaceErrorEstimator<dim>::~SpaceErrorEstimator()
{
#ifdef TRACE
  (*trace) << "entering SpaceErrorEstimator::~SpaceErrorEstimator()" << std::endl;
#endif 
;
}

#ifdef __GNUC__
template class SpaceErrorEstimator<2>;
template class SpaceErrorEstimator<3>;
#endif

} // end of namespace

#endif
