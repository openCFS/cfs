#ifndef FILE_SOLVESTEP_ODE_HH
#define FILE_SOLVESTEP_ODE_HH

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step of an ODE, wrapped in a PDE

  class SolveStepODE : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepODE(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepODE();


    //----------------------- STATIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    void PreStepStatic( const Boolean reset ) {
      Error( "An ODE can not be solved statically!", __FILE__, __LINE__ );

    }

    //! base method for solving one static step 
    //!\param reset TRUE: perfrom new assembly, etc
    void SolveStepStatic( const Boolean reset ) {
      Error( "An ODE can not be solved statically!", __FILE__, __LINE__ );
    }

    //! solves for one nonlinear static step 
    //! \param reset TRUE: perfrom new assembly, etc
    void StepStaticNonLin( const Boolean reset ) {
      Error( "An ODE can not be solved statically!", __FILE__, __LINE__ );
    }
  
    //! routine for acttions after the SolveStep-method 
    void PostStepStatic() {
      Error( "An ODE can not be solved statically!", __FILE__, __LINE__ );
 }

    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    void PreStepTrans( const Boolean reset ) {; };

    //! base method for solving one transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    void SolveStepTrans( const Boolean reset );

    //! solves for one linear transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    void StepTransLin( const Boolean reset ) {;};

    //! routine for actions after the SolveStep-method
     void PostStepTrans() {;};

    //----------------------- HARMONIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void PreStepHarmonic( const Boolean reset ) {
      Error( "An ODE can not be solved harmonically!", __FILE__, __LINE__ );
    }

    //!  base method for solving one harmonic step 
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void SolveStepHarmonic( const Boolean reset ) {
      Error( "An ODE can not be solved harmonically!", __FILE__, __LINE__ );
    }
    
    //! solves for one linear frequency step 
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void StepHarmonicLin( const Boolean reset ) {
      Error( "An ODE can not be solved harmonically!", __FILE__, __LINE__ );
    }

    //! solves for one nonlinear frequency step 
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void StepHarmonicNonLin( const Boolean reset ) {
      Error( "An ODE can not be solved harmonically!", __FILE__, __LINE__ );
    }

    //!  routine for actions after the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void PostStepHarmonic( const Boolean reset ) {
      Error( "An ODE can not be solved harmonically!", __FILE__, __LINE__ );
    }

  private:

  };

} // end of namespace

#endif

