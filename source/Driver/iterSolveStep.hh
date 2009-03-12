// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ITERSOLVESTEP
#define FILE_ITERSOLVESTEP

#include "PDE/basePDE.hh"
#include "Driver/baseSolveStep.hh"
#include "MatVec/SingleVector.hh"

namespace CoupledField
{
  // forward class declarations
  class IterCoupledPDE;
  class PDECoupling;

  //! Derived class for step-wise solving of iterative coupled StdPDEs
  class IterSolveStep : public BaseSolveStep
  {

  public:

    //! Constructor
    IterSolveStep(IterCoupledPDE& apde);

    //! Destructor
    virtual ~IterSolveStep();


    //----------------------- STATIC---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepStatic()  {;};
 
    /** base method for solving one static step
     * @partam step is an optinal information for exoport linsys only */
    virtual void SolveStepStatic(const std::string& comment = "");

    //! routine for acttions after the SolveStep-method
    virtual void PostStepStatic() {;}


    //----------------------- TRANSIENT---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepTrans() {;};


    //! base method for solving one transient step 
    virtual void SolveStepTrans();
    
    //! routine for actions after the SolveStep-method
    virtual void PostStepTrans() {;};

    //----------------------- HARMONIC---------------------------------------
    
    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepHarmonic(  ) {;};


    //!  base method for solving one harmonic step 
    virtual void SolveStepHarmonic(const std::string& comment = "");


    //!  routine for actions after the SolveStep-method
    virtual void PostStepHarmonic(  ) {;};


    //----------------------- SET/GET METHODS--------------------------------
    
    //! Set actual time
    virtual void SetActTime( const Double actTime );

    //! Set actual frequency
    virtual void SetActFreq( const Double actFreq );

    //! Set actual time / frequency step
    virtual void SetActStep( const UInt actStep );

    //! Set the current time step value
    void SetTimeStep( Double dt );

    //! Set number of time steps
    virtual void SetNumTimeSteps (UInt numTimeStep);
    
    //! Set restart time / frequency step
    virtual void SetStartStep (const UInt startStep);

  protected:

    //! calculates the norm of a vector
    Double CalcNorm(NormType normtype, SingleVector & val, SingleVector & oldval); 
    

    //! reference to PDE
    IterCoupledPDE &rPDE_;

    //! reference to coupling
    StdVector<PDECoupling*> & rCouplings_;

    //! analysis type of all iteratively coupled PDEs is retrieved
    BasePDE::AnalysisType actAnalysisType_;

  };

} // end of namespace

#endif

