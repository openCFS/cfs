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

//! class for space error estimation
template<Integer dim>
class SpaceErrorEstimator
{
public:
  //! constructor with pointer to PDE
  /*!
    \param aptGrid pointer to the Grid
  */
  SpaceErrorEstimator(Grid * aptGrid);

  //! Deconstructor
  virtual ~SpaceErrorEstimator();

  //! Init pointer to PDE
  /*!
    \param aptPDE pointer to the class of the PDE
  */
  virtual void Init(BasePDE * aptPDE){ ptPDE_=aptPDE;}

  //! kelly error estimation
  /*!
    \param sbdoms vector with names of subdomains, on which results are calculated
    \param level level of the Grid
  */
  void KellyError(std::vector<std::string> &sbdoms, const Integer level);

  //! calculation of error map
  /*!
  \param sol solution 
  \param subdoms vector with names of subdomains, on which we do calculation
  \param ptgrid pointer to the Grid
  \param errorMap out: vector with calculated error for each element. \f $||grad_{\SPR} - grad_{FEM}||_L2$
   \param  gradSPRElemL2norm out: vector with calculated L2-norm of the SPR grad for each element.\f $||grad_{SPR}||_L2$
  */
 void CalcErrorMap(const Vector<Double> * sol, std::vector<std::string> & subdoms,
 Grid * ptgrid, Vector<Double> & errorMap, Vector<Double> & gradSPRElemL2norm);

  //! recovery procedure for the elements-patch
  /*!
    \param Elems vector with information about elements
    \param ptgrid pointer to the Grid
     \param sol solution 
      \param aComponent number of component of gradient 
     \param result (output) result
     \param locations (output) global numbers of nodes in patch. needful for the vector with result
  */
  void RecoveryProcedure4ElemsPatch(const std::vector<Elem*> &Elems, Grid * ptgrid, const Vector<Double> & sol, const Integer aComponent, Vector<Double> &result,std::vector<Integer>&locations);

  //! recovery procedure for the elements-patch: 3d
   /*!
    \param Elems vector with information about elements
    \param ptgrid pointer to the Grid
     \param sol solution 
     \param aComponent number of component of gradient 
     \param result (output) result
     \param locations (output) global numbers of nodes in patch. needful for the vector with result
  */
  void RecoveryProcedure4ElemsPatch3d(const std::vector<Elem*> &Elems, Grid * ptgrid, const Vector<Double> & sol, const Integer aComponent, Vector<Double> &result,std::vector<Integer>&locations);

  //! only for 1 element. this function only for testing
  void RecoveryProcedure4Elem(Elem * aptElem, Grid * ptgrid, const Vector<Double> & sol, const Integer aComponent, Vector<Double> &result);

protected:
  //!
  BasePDE * ptPDE_;  

  //!
  Grid * ptGrid_;  

  //! calculation error for one element
  /*!
    \param elem pointer to the element
    \param SPRgrad (output) returned calculated SPR gradient 
    \param error (output) retuned value of the error
    \param normGradSPR (output) L2-norm of SPR gradient
    \param sol (input) solution
    \param ptgrid pointer to Grid
  */
  void CalcErrorForElem(const Elem* elem, const Vector<Double>* SPRgrad,  Double & error,
  Double & normGradSPR, const Vector<Double>* sol, Grid * ptgrid);

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
