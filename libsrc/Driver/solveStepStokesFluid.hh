#ifndef FILE_SOLVESTEPSTOKESFLUID
#define FILE_SOLVESTEPSTOKESFLUID

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Magnetics

  class SolveStepStokesFluid : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepStokesFluid(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepStokesFluid();

    //----------------------- STATIC---------------------------------------
    //! solves for one nonlinear static step 
    void StepStaticNonLin();

    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    void PreStepTrans( )
    {PreStepStatic();};

    //! base method for solving one transient step 
    //! \param reset true: perfrom new assembly, etc
    void SolveStepTrans()
    {SolveStepStatic();};

    //! solves for one linear transient step 
    void StepTransLin( )
    {StepStaticLin();};

    //! routine for actions after the SolveStep-method
     void PostStepTrans()
    {PostStepStatic();};


    //! solves for one nonlinear transient step 
    void StepTransNonLin()
    {StepStaticNonLin();};

    //! compute nonlinear part of RHS
    void AddNonLinRHS();

  private:


  };

} // end of namespace

#endif

