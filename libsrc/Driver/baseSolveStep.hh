#ifndef FILE_BASESOLVESTEP
#define FILE_BASESOLVESTEP

#include "Utils/StdVector.hh"
#include "General/environment.hh"

namespace CoupledField
{

  //! Base class for solution of a single step

  class BaseSolveStep
  {

  public:

    

    //! Destructor
    virtual ~BaseSolveStep();


    //----------------------- STATIC---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepStatic( const Boolean reset) = 0;
    
    //! base method for solving one static step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepStatic( const Boolean reset) = 0;

    //! routine for acttions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
    */  
    virtual void PostStepStatic()  = 0;



    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepTrans( const Boolean reset) = 0;

    //! base method for solving one transient step 
    /*!
      \param kstep time step counter
      \param asteptime current time
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepTrans(const Boolean reset) = 0;

    //! routine for actions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
    */  
    virtual void PostStepTrans() = 0;
    
    //----------------------- HARMONIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param reset TRUE: perfrom new assembly, etc
    */   
    virtual void PreStepHarmonic(const Boolean reset) = 0;

    //!  base method for solving one harmonic step 
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepHarmonic(const Boolean reset) = 0;

    //!  routine for actions after the SolveStep-method
    /*!
      \param freqStep frequency step counter
      \param frequency current frequency
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void PostStepHarmonic(const Boolean reset) = 0;

    //----------------------- SET/ GET METHODS--------------------------------

    //! Set actual time
    virtual void SetActTime( const Double actTime ) {
      actTime_ = actTime;
    }

    //! Set actual frequency
    virtual void SetActFreq( const Double actFreq ) {
      actFreq_ = actFreq;
    }

    //! Set actual time / frequency step
    virtual void SetActStep( const UInt actStep ) {
      actStep_ = actStep;
    }

    //! Return actual time / frequency step
    virtual UInt GetActStep() {
      return actStep_;
    }

    //! Return actual frequency 
    virtual Double GetActFreq() {
      return actFreq_;
    }

    //! Return actual time
    virtual Double GetActTime() {
      return actTime_;
    }
    
    
    //! Set the current time step value
    virtual void SetTimeStep( Double dt ) = 0;

    

  protected:

    //! Constructor
    BaseSolveStep();


    //! Actual time / frequency step
    UInt actStep_;

    //! Actual time
    Double actTime_;

    //! Actual frequency 
    Double actFreq_;
    //!
  };


} // end of namespace

#endif

