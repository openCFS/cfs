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
    //! \param reset TRUE: perfrom new assembly, etc
    void StepStaticNonLin( const Boolean reset );

    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    void PreStepTrans( const Boolean reset )
    {PreStepStatic(reset);};

    //! base method for solving one transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    void SolveStepTrans( const Boolean reset )
    {SolveStepStatic(reset);};

    //! solves for one linear transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    void StepTransLin( const Boolean reset )
    {StepStaticLin(reset);};

    //! routine for actions after the SolveStep-method
     void PostStepTrans()
    {PostStepStatic();};


    //! solves for one nonlinear transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    void StepTransNonLin( const Boolean reset )
    {StepStaticNonLin(reset);};

    //! compute nonlinear part of RHS
    void AddNonLinRHS();

  private:


  };

} // end of namespace

#endif

