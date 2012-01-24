#ifndef ACOUSTICPDE_HH
#define ACOUSTICPDE_HH

#include "SinglePDE.hh"

namespace CoupledField{

  // forward class declaration
  class BaseResult;
  class ResultHandler;
  class linElecInt;
  class LinearFormContext;
  class ElecForceOp;
  class BaseBDBInt;
  class ResultFunctor;

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

    // ======================================================
    // POSTPROCESSING SECTION
    // ======================================================

    //! perform postprocessing on data
    void CalcResults( shared_ptr<BaseResult> result );

  protected:

    //!  Define available postprocessing results
    void DefinePrimaryResults();

    //! Init the time stepping
    void InitTimeStepping();

  private:

    //! stores if the Acoustic PDE is in potential or pressure form
    SolutionType formulation_;


  };

}

#endif
