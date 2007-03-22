// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SOLVESTEPSMOOTH
#define FILE_SOLVESTEPSMOOTH

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step

  class SolveStepSmooth : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepSmooth(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepSmooth();

    //----------------------- STATIC---------------------------------------

    //! solves for one nonlinear static step 
    //! \param reset true: perfrom new assembly, etc
    void StepStaticNonLin( const bool reset );

    //----------------------- TRANSIENT---------------------------------------  
    //! base method for solving one transient step 
    //! \param reset true: perfrom new assembly, etc
    void SolveStepTrans( const bool updatesysmat )
    {SolveStepStatic(updatesysmat);};

    //! routine for initilizations befor execution the SolveStep-method
    //! \param reset true: perfrom new assembly, etc
    void PreStepTrans( const bool reset )
    {PreStepStatic(reset);};
  
    //! routine for actions after the SolveStep-method
    void PostStepTrans()
    {PostStepStatic();};
  
    //! solves for one nonlinear transient step 
    //! \param reset true: perfrom new assembly, etc
    void StepTransNonLin( const bool reset )
    {StepStaticNonLin(reset);};


  private:


  };

} // end of namespace

#endif

