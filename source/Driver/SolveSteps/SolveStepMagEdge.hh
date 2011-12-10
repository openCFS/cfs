#ifndef FILE_SOLVE_STEP_MAG_EDGE
#define FILE_SOLVE_STEP_MAG_EDGE

#include "StdSolveStep.hh"
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

    //----------------------- STATIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    void PreStepStatic( );

    //! base method for solving one static step 
    void SolveStepStatic(PtrParamNode analysis_id);

    /** @see SolveStepStatic() */ 
    void StepStaticLin(PtrParamNode analysis_id);

    //! solves for one nonlinear static step 
    void StepStaticNonLin(PtrParamNode analysis_id);
    
    
    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    void PreStepTrans( )
    {PreStepTrans();};

    //! base method for solving one transient step 
    void SolveStepTrans(PtrParamNode analysis_id)
    {SolveStepTrans(analysis_id);};

    //! solves for one linear transient step 
    void StepTransLin(PtrParamNode analysis_id)
    {StepTransLin(analysis_id);};

    //! routine for actions after the SolveStep-method
     void PostStepTrans()
    {PostStepTrans();};

  private:
     
     //! Explicit pointer to MagEdgePDE
     MagEdgePDE * magEdgePDE_;

     //! log file for time
     std::ofstream timeFile_;
  };

} // end of namespace

#endif

