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
       \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepStatic( const Boolean reset) = 0;
    
    //! base method for solving one static step 
    /*!
       \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepStatic( const Boolean reset) = 0;

    //! routine for acttions after the SolveStep-method
    /*!
       \param asteptime current time
    */  
    virtual void PostStepStatic()  = 0;



    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    /*!
       \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepTrans( const Boolean reset) = 0;

    //! base method for solving one transient step 
    /*!
       \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepTrans(const Boolean reset) = 0;

    //! base method for solving one transient step with slicing method 
    /*!
       \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepTrans4Slice(const Boolean reset)
    {Error("SolveStepTrans4Slice not implemented!",__FILE__,__LINE__);};

    //! routine for actions after the SolveStep-method
    /*!
      \param kstep time step counter
      \param asteptime current time
    */  
    virtual void PostStepTrans() = 0;
    
    //----------------------- HARMONIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    /*!
      \param reset TRUE: perfrom new assembly, etc
    */   
    virtual void PreStepHarmonic(const Boolean reset) = 0;

    //!  base method for solving one harmonic step 
    /*!
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


    //----------------------- SPECIAL FUNCTIONS ------------------------------
    //! transform solution and derivatives due to slicing technique
    virtual void TransformSol4Slice(UInt & nodeShift, UInt & shiftFactor, 
                                    const UInt flag)
    {Error("TransformSol4Slice not implemented", __FILE__, __LINE__); };

    //! save important nodes for the a post analysis
    //! \param shiftFactor number of nodes in x and y direction
    //! \param numShift how many shifts have been made
    //! \param nodeShift how many nodes have been shifted
    //! save solution of special nodes
    void SaveNodes(const UInt shiftFactor, const Double timeStep,
		   const UInt numShift, const Integer nodeShift, 
		   const UInt maxnumelemz_)
    {Error("TransformSol4Slice not implemented"); };



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

    //! Set number of time steps
    virtual void SetNumTimeSteps( UInt numTimeStep ) {
      ENTER_FCN( "BaseSolveStep::SetNumSteps") ;
      numTimeStep_ = numTimeStep;
    };
    
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

    //! number of time steps
    UInt numTimeStep_;
    //!
  };


} // end of namespace

#endif

