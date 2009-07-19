#ifndef FILE_SOLVESTEPSMOOTH
#define FILE_SOLVESTEPSMOOTH

#include "OLAS/algsys/basesystem.hh"

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Electrostatic-PDE

  class SolveStepSmooth : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepSmooth(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepSmooth();


    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    //void PreStepTrans();

    //! base method for solving one transient step 
    //void SolveStepTrans();

//     //! solves for one linear transient step 
//     void StepTransLin( )
//     {StepStaticLin();};

    //! routine for actions after the SolveStep-method
    //void PostStepTrans();

    //! solves for one nonlinear transient step 
    void StepTransNonLin(InfoNode* analysis_id);

  private:
  };

} // end of namespace

#endif

