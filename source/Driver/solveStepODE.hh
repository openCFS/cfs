// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SOLVESTEP_ODE_HH
#define FILE_SOLVESTEP_ODE_HH

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step of an ODE, wrapped in a PDE

  class SolveStepODE : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepODE(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepODE();


    //----------------------- STATIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    void PreStepStatic() {
      Error( "An ODE can not be solved statically!", __FILE__, __LINE__ );

    }

    //! base method for solving one static step 
    void SolveStepStatic() {
      Error( "An ODE can not be solved statically!", __FILE__, __LINE__ );
    }

    //! solves for one nonlinear static step 
    void StepStaticNonLin() {
      Error( "An ODE can not be solved statically!", __FILE__, __LINE__ );
    }
  
    //! routine for acttions after the SolveStep-method 
    void PostStepStatic() {
      Error( "An ODE can not be solved statically!", __FILE__, __LINE__ );
 }

    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    void PreStepTrans() {; };

    //! base method for solving one transient step 
    void SolveStepTrans();

    //! solves for one linear transient step 
    void StepTransLin() {;};

    //! routine for actions after the SolveStep-method
     void PostStepTrans() {;};

    //----------------------- HARMONIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepHarmonic() {
      Error( "An ODE can not be solved harmonically!", __FILE__, __LINE__ );
    }

    //!  base method for solving one harmonic step 
    virtual void SolveStepHarmonic() {
      Error( "An ODE can not be solved harmonically!", __FILE__, __LINE__ );
    }
    
    //! solves for one linear frequency step 
    virtual void StepHarmonicLin() {
      Error( "An ODE can not be solved harmonically!", __FILE__, __LINE__ );
    }

    //! solves for one nonlinear frequency step 
    virtual void StepHarmonicNonLin() {
      Error( "An ODE can not be solved harmonically!", __FILE__, __LINE__ );
    }

    //!  routine for actions after the SolveStep-method
    virtual void PostStepHarmonic() {
      Error( "An ODE can not be solved harmonically!", __FILE__, __LINE__ );
    }

  private:

  };

} // end of namespace

#endif

