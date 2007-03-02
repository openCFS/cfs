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
    void StepStaticLin( )
    { ENTER_FCN( "SolveStepMpCCI::StepStaticLin");};


    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
     void PreStepTrans( )
    { ENTER_FCN( "SolveStepMpCCI::PreStepTrans");};
    
    //! base method for solving one transient step 
    void SolveStepTrans( )
    { ENTER_FCN( "SolveStepMpCCI::SolveStepTrans"); };

    //! solves for one linear transient step 
    //! \param reset true: perfrom new assembly, etc
    void StepTransLin( )
    { ENTER_FCN( "SolveStepMpCCI::StepTransLin");};

    //! routine for actions after the SolveStep-method
    void PostStepTrans();


  private:


  };

} // end of namespace

#endif

