#ifndef FILE_SOLVESTEPMAG
#define FILE_SOLVESTEPMAG

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Magnetics

  class SolveStepMag :  public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepMag(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepMag();

    //----------------------- STATIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    void PreStepStatic( const Boolean reset );
    
    //! solves for one linear static step 
    //! \param reset TRUE: perfrom new assembly, etc
    void StepStaticNonLin( const Boolean reset );
    
    //! routine for acttions after the SolveStep-method
    void PostStepStatic();


    //----------------------- TRANSIENT---------------------------------------
    //! solves for one nonlinear transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    void StepTransNonLin( const Boolean reset );


  private:


  };

} // end of namespace

#endif

