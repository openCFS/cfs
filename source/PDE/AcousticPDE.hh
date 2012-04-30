#ifndef ACOUSTICPDE_HH
#define ACOUSTICPDE_HH

#include "SinglePDE.hh"

namespace CoupledField{

  // forward class declaration
  class BaseResult;
  class ResultHandler;
  class LinearFormContext;
  class ElecForceOp;
  class BaseBDBInt;
  class ResultFunctor;
  class CoefFunctionMulti;

  class AcousticPDE : public SinglePDE{

  public:
    //!  Constructor.
    /*!
      \param aGrid pointer to grid
    */
    AcousticPDE( Grid* aGrid, PtrParamNode paramNode );

    virtual ~AcousticPDE(){};

    //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> >
    CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode );

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators( );

    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //! Initialize all the nodes by this value
    void SetInitialCondition();

    // ======================================================
    // COUPLING SECTION
    // ======================================================

    //! Initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling) {
      EXCEPTION("Coupling not supported for acousticPDE.");
    };

    //! Calculate coupling terms
    void CalcOutputCoupling() {
      EXCEPTION("Coupling not supported for acousticPDE.");
    };

    //! Returns if PDE can compute the quantity
    bool HasOutput(SolutionType output) { return false;};

    //! Return acoustic formulation. Can either be pressure or potential.
    SolutionType GetFormulation() const { return formulation_; }

    //! Indicate that acoustic PDE takes part in coupling to mechanics.
    void SetMechanicCoupling() {
      isMechCoupled_ = true;
    }

  protected:

    //!  Define available primary results
    void DefinePrimaryResults();
    
    //! Define available postprocessing results
    void DefinePostProcResults();

    //! Init the time stepping
    void InitTimeStepping();

    //! create feFunction for meanFluidMech velocity
    void CreateMeanFlowFunction(StdVector<std::string> dofNames);

  private:

    //! stores if the Acoustic PDE is in potential or pressure form
    SolutionType formulation_;

    //! Coefficient function for the flow field

    //! This coefficient function describes the flow field. As this 
    //! is in general different for each region and will most likely
    //! not be given in a close form, it is described by a CoefFunctionMulti.
    shared_ptr<CoefFunctionMulti> meanFlowCoef_;
    
    //! stores if the Acoustic PDE is coupled to mechanics
    bool isMechCoupled_;
  };

}

#endif
