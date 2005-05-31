#ifndef FILE_SOLVESTEPMECH
#define FILE_SOLVESTEPMECH

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step

  class SolveStepMech : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepMech(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepMech();

    //----------------------- STATIC---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    void PreStepStatic( const Boolean reset );

    //! solves for one nonlinear static step 
    //! \param reset TRUE: perfrom new assembly, etc
    void StepStaticNonLin( const Boolean reset );

    //! routine for acttions after the SolveStep-method
    void PostStepStatic();


    //----------------------- TRANSIENT---------------------------------------
    
    //! solves for one linear transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    void StepTransNonLin( const Boolean reset );


  private:


  };

} // end of namespace

#endif

