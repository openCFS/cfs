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

    //! routine for initilizations befor execution the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    void PreStepStatic( const Boolean reset );

    //! solves for one nonlinear static step 
    //! \param reset TRUE: perfrom new assembly, etc
    void StepStaticNonLin( const Boolean reset );

    //! routine for acttions after the SolveStep-method
    void PostStepStatic();


    //----------------------- TRANSIENT---------------------------------------  
    //! base method for solving one transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    void SolveStepTrans( const Boolean updatesysmat )
    {SolveStepStatic(updatesysmat);};

    //! routine for initilizations befor execution the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    void PreStepTrans( const Boolean reset )
    {PreStepStatic(reset);};
  
    //! routine for actions after the SolveStep-method
    void PostStepTrans()
    {PostStepStatic();};
  
    //! solves for one nonlinear transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    void StepTransNonLin( const Boolean reset )
    {StepStaticNonLin(reset);};


  private:


  };

} // end of namespace

#endif

