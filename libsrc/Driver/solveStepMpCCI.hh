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
    //! \param reset TRUE: perfrom new assembly, etc
    void PreStepStatic( const Boolean reset );

    //! routine for acttions after the SolveStep-method
    void PostStepStatic();
                                

    //! solves for one linear static step 
    //! \param reset TRUE: perfrom new assembly, etc
    void StepStaticLin( const Boolean reset )
    { ENTER_FCN( "SolveStepMpCCI::StepStaticLin");};


    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
     void PreStepTrans( const Boolean reset )
    { ENTER_FCN( "SolveStepMpCCI::PreStepTrans");};
    
    //! base method for solving one transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    void SolveStepTrans( const Boolean reset )
    { ENTER_FCN( "SolveStepMpCCI::SolveStepTrans"); };

    //! solves for one linear transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    void StepTransLin( const Boolean reset )
    { ENTER_FCN( "SolveStepMpCCI::StepTransLin");};

    //! routine for actions after the SolveStep-method
    void PostStepTrans();


  private:


  };

} // end of namespace

#endif

