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

    Double InexactLineSearch(SBM_Vector& solIncrement, SBM_Vector& actSol);

    Double LineSearchArmijo(SBM_Vector& solIncrement, SBM_Vector& actSol);

    double GetLineSearchDerivativeFunctionValue(SBM_Vector& solIncrement, SBM_Vector& actSol, Double eta);

    double BrentMethod(SBM_Vector& solIncrement, SBM_Vector& actSol, Double a, Double b);

  protected:

    UInt pseudo_time_stepping_ = 0;

  private:

    // Coefficient function for material model
    std::map<RegionIdType, shared_ptr<CoefFunctionMaterialModel<Complex>> >  matModelCoefm_;
  };
} // end of namespace

#endif
