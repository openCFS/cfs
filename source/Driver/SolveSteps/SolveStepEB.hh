// ================================================================================================
/*!
 * description: Currently only for non-dynamical systems (only non-time dependent stiffness parts in the bilinear form)      
 * date     Jul 15, 2022
 * author   kroppert
 */
//================================================================================================

#ifndef FILE_SOLVESTEPEB
#define FILE_SOLVESTEPEB

#include "StdSolveStep.hh"
#include <Domain/CoefFunction/CoefFunctionMaterialModel.hh>

namespace CoupledField
{

  //! Base class for solution of a single step: Electrostatic-PDE

  class SolveStepEB : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepEB(StdPDE& apde);
    SolveStepEB(StdPDE& apde, UInt is_pseudo_time_stepping);

    //! Destructor
    virtual ~SolveStepEB();

    //! base method for solving one transient step 
    void SolveStepTrans();

    //! solves for one linear transient step 
    //void StepTransLin(){}

    void StepTransNonLin();

    Double ExactLineSearch(SBM_Vector& solIncrement, SBM_Vector& actSol);
    Double ExactLineSearch(SBM_Vector& solIncrement, SBM_Vector& actSol, SBM_Vector& Linform_nm1);

    Double InexactLineSearch(SBM_Vector& solIncrement, SBM_Vector& actSol);

    Double LineSearchArmijo(SBM_Vector& solIncrement, SBM_Vector& actSol);

    //! Original discrete line search (searches eta in {0.1, 0.2, ..., 1.0})
    void LineSearchHeavy(SBM_Vector& solIncrement, SBM_Vector& actSol, Double& etaLineSearch);
    void LineSearchHeavy(SBM_Vector& solIncrement, SBM_Vector& actSol, Double& etaLineSearch, SBM_Vector& Linform_nm1);

    double GetLineSearchDerivativeFunctionValue(SBM_Vector& solIncrement, SBM_Vector& actSol, Double eta);
    double GetLineSearchDerivativeFunctionValue(SBM_Vector& solIncrement, SBM_Vector& actSol, Double eta, SBM_Vector& Linform_nm1);

    double BrentMethod(SBM_Vector& solIncrement, SBM_Vector& actSol, Double a, Double b);
    double BrentMethod(SBM_Vector& solIncrement, SBM_Vector& actSol, Double a, Double b, SBM_Vector& Linform_nm1);

  protected:

    UInt pseudo_time_stepping_ = 0;

  private:

    // Coefficient function for material model
    std::map<RegionIdType, shared_ptr<CoefFunctionMaterialModel<Complex>> >  matModelCoefm_;
  };
} // end of namespace

#endif
