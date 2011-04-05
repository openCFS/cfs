// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SOLVE_STEP_MAG_EDGE
#define FILE_SOLVE_STEP_MAG_EDGE

#include "stdSolveStep.hh"
#include "General/Enum.hh"

namespace CoupledField
{

  class MagEdgePDE;

  //! Base class for solution of a single step: Magnetic-Edge PDE  
  class SolveStepMagEdge : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepMagEdge( StdPDE& apde );

    //! Destructor
    virtual ~SolveStepMagEdge();

    //! Set solver strategy
    void SetSolStrategy( SolStrategyType strategy );
    
    //----------------------- STATIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    void PreStepStatic( );

    //! base method for solving one static step 
    void SolveStepStatic(InfoNode* analysis_id);

    /** @see SolveStepStatic() */ 
    void StepStaticLin(InfoNode* analysis_id);

    //! solves for one nonlinear static step 
    void StepStaticNonLin(InfoNode* analysis_id);
    
    
    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    void PreStepTrans( )
    {PreStepTrans();};

    //! base method for solving one transient step 
    void SolveStepTrans(InfoNode* analysis_id)
    {SolveStepTrans(analysis_id);};

    //! solves for one linear transient step 
    void StepTransLin(InfoNode* analysis_id)
    {StepTransLin(analysis_id);};

    //! routine for actions after the SolveStep-method
     void PostStepTrans()
    {PostStepTrans();};

  private:
     
     //! Explicit pointer to MagEdgePDE
     MagEdgePDE * magEdgePDE_;

     //! Solution strategy
     SolStrategyType strategy_;
  };

} // end of namespace

#endif

