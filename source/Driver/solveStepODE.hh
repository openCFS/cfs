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
      EXCEPTION( "An ODE can not be solved statically!" );

    }

    //! base method for solving one static step 
    void SolveStepStatic() {
      EXCEPTION( "An ODE can not be solved statically!" );
    }

    //! solves for one nonlinear static step 
    void StepStaticNonLin() {
      EXCEPTION( "An ODE can not be solved statically!" );
    }
  
    //! routine for acttions after the SolveStep-method 
    void PostStepStatic() {
      EXCEPTION( "An ODE can not be solved statically!" );
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
      EXCEPTION( "An ODE can not be solved harmonically!" );
    }

    //!  base method for solving one harmonic step 
    virtual void SolveStepHarmonic() {
      EXCEPTION( "An ODE can not be solved harmonically!" );
    }
    
    //! solves for one linear frequency step 
    virtual void StepHarmonicLin() {
      EXCEPTION( "An ODE can not be solved harmonically!" );
    }

    //! solves for one nonlinear frequency step 
    virtual void StepHarmonicNonLin() {
      EXCEPTION( "An ODE can not be solved harmonically!" );
    }

    //!  routine for actions after the SolveStep-method
    virtual void PostStepHarmonic() {
      EXCEPTION( "An ODE can not be solved harmonically!" );
    }

  private:

  };

} // end of namespace

#endif

