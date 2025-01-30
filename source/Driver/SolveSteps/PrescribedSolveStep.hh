#ifndef FILE_PRESCRIBEDSOLVESTEP
#define FILE_PRESCRIBEDSOLVESTEP

#include "StdSolveStep.hh"

namespace CoupledField
{

  //! Solve step that does NOT assemble or solve anything.
  //!
  //! Instead of building and solving a linear system, it writes a prescribed external field
  //! directly into the fe-function solution vector(s) via BaseFeFunction::ApplyExternalData()
  //! (the same node-wise injection path used by <initialValues><initialField>). The external
  //! data source(s) must have been registered on the fe-function(s) beforehand
  //! (PDE side: BaseFeFunction::AddExternalDataSource).
  //!
  //! Used by SmoothPDE in "prescribed displacement" mode, where the whole-domain deformation
  //! is provided externally (e.g. computed by OpenFOAM and read via preCICE) and openCFS only
  //! needs to make it available internally as SMOOTH_DISPLACEMENT - without solving.
  //!
  //! It derives from StdSolveStep purely to reuse its construction (wiring of feFunctions_ and
  //! the solution vectors); all assembling/solving methods are overridden to the pure injection.
  class PrescribedSolveStep : public StdSolveStep
  {
  public:
    //! Constructor
    PrescribedSolveStep(StdPDE& apde);

    //! Destructor
    virtual ~PrescribedSolveStep();

    //----------------------- STATIC -----------------------
    virtual void PreStepStatic() override {}
    virtual void SolveStepStatic() override { ApplyPrescribed(); }
    virtual void PostStepStatic() override {}

    //----------------------- TRANSIENT -----------------------
    //! No time-integration setup is needed: the field is set directly each step.
    virtual void InitTimeStepping() override {}
    virtual void PreStepTrans() override {}
    virtual void SolveStepTrans() override { ApplyPrescribed(); }
    virtual void PostStepTrans() override {}

    //! re-inject the (refreshed) external field, see BaseSolveStep::RefreshPrescribed()
    virtual void RefreshPrescribed() override { ApplyPrescribed(); }

    //----------------------- HARMONIC (unsupported) -----------------------
    virtual void PreStepHarmonic() override {}
    virtual void SolveStepHarmonic() override
    {
      EXCEPTION("PrescribedSolveStep: harmonic analysis is not supported for a prescribed field.");
    }
    virtual void PostStepHarmonic() override {}

  private:
    //! Write the registered external data source(s) directly into the solution vector(s).
    //! No assembly, no linear solve. Re-evaluated every call, so time-dependent sources
    //! (analytic expressions, files, or a preCICE-fed grid) are picked up automatically.
    void ApplyPrescribed();
  };

} // end of namespace

#endif
