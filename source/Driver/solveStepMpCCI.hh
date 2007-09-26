// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SOLVESTEPMPCCI
#define FILE_SOLVESTEPMPCCI

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: MpCCI-PDE

  class SolveStepMpCCI : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepMpCCI(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepMpCCI();


    //----------------------- STATIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    void PreStepStatic( );


    //! solves for one linear static step 
    void StepStaticLin( ) {};


    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    void PreStepTrans( ) {};
    
    //! base method for solving one transient step 
    void SolveStepTrans( ) {};

    //! solves for one linear transient step 
    //! \param reset true: perfrom new assembly, etc
    void StepTransLin( ) {};

    //! routine for actions after the SolveStep-method
    void PostStepTrans();
    

  private:


  };

} // end of namespace

#endif

