#ifndef FILE_ITERSOLVESTEP
#define FILE_ITERSOLVESTEP

#include "baseSolveStep.hh"
//#include "CoupledPDE/itercoupledpde.hh"
//#include "CoupledPDE/pdecoupling.hh"

namespace CoupledField
{
  // forward class declarations
  class IterCoupledPDE;
  class PDECoupling;

  //! Derived class for step-wise solving of iterative coupled StdPDEs
  class IterSolveStep : public BaseSolveStep
  {

  public:

    //! Constructor
    IterSolveStep(IterCoupledPDE& apde);

    //! Destructor
    virtual ~IterSolveStep();


    //----------------------- STATIC---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void PreStepStatic( const Boolean reset )  {;};
 
    //! base method for solving one static step 
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void SolveStepStatic( const Boolean reset );

    //! routine for acttions after the SolveStep-method
    virtual void PostStepStatic() {;}


    //----------------------- TRANSIENT---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void PreStepTrans( const Boolean reset ) {;};


    //! base method for solving one transient step 
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void SolveStepTrans( const Boolean updatesysmat );
    
    //! routine for actions after the SolveStep-method
    virtual void PostStepTrans() {;};

    //----------------------- HARMONIC---------------------------------------
    
    //! routine for initilizations befor execution the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void PreStepHarmonic( const Boolean reset ) {;};


    //!  base method for solving one harmonic step 
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void SolveStepHarmonic( const Boolean reset );


    //!  routine for actions after the SolveStep-method
    //! \param reset TRUE: perfrom new assembly, etc
    virtual void PostStepHarmonic( const Boolean reset ) {;};

    //! Set actual time
    virtual void SetActTime( const Double actTime );
    

    //! Set actual frequency
    virtual void SetActFreq( const Double actFreq );

    //! Set actual time / frequency step
    virtual void SetActStep( const UInt actStep );

    //! Set the current time step value
    void SetTimeStep( Double dt );


  protected:

    //! calculates the norm of a vector
    Double CalcNorm(NormType normtype, CFSVector & val, CFSVector & oldval); 
    

    //! reference to PDE
    IterCoupledPDE &rPDE_;

    //! reference to coupling
    StdVector<PDECoupling*> & rCouplings_;

  };

} // end of namespace

#endif

