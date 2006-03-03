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
    virtual void PostStepStatic()  = 0;



    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations before execution of the SolveStep-method
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
    virtual void PostStepTrans() = 0;
    
    //----------------------- HARMONIC---------------------------------------
    //! routine for initilizations before execution of the SolveStep-method
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
      \param reset TRUE: perfrom new assembly, etc
    */
    virtual void PostStepHarmonic(const Boolean reset) = 0;

    //----------------------- HARMONIC ---------------------------------------
    //! Calculate the Eigenfrequencies of a generalized eigenvalue problem
    virtual UInt CalcEigenFrequencies( Vector<Double> & frequencies, 
                                       UInt numFreq, Double shift, Boolean shiftMode ) {
      Error( "Not implemented her!", __FILE__, __LINE__ );
      return 0;
    }
    
    //! Calculate the numMode-th eigenmode of a generalized eigenvalue problem.
    //! Therefore, previously CalcEigenFrequencies() has to be called.
    virtual void CalcEigenMode( UInt numMode ) {
      Error( "Not implemented her!", __FILE__, __LINE__ );
    }

    //----------------------- TRANSIENTHARMONIC------------------------------
    //! routine for initilizations before execution of the SolveStep-method
    /*!
       \param reset TRUE: perfrom new assembly, etc
    */  
    virtual void PreStepTransHarmonic( const Boolean reset)
    {Error("PreStepTransHarmonic not implemented!",__FILE__,__LINE__);};

    //! base method for solving one transient-harmonic coupled step 
    /*!
       \param reset TRUE: perfrom new assembly, etc
    */
    virtual void SolveStepTransHarmonic(const Boolean reset)
    {Error("SolveStepTransHarmonic not implemented!",__FILE__,__LINE__);};

    //! routine for actions after the SolveStep-method
    virtual void PostStepTransHarmonic()
    {Error("PostStepTransHarmonic not implemented!",__FILE__,__LINE__);};


    //----------------------- SPECIAL FUNCTIONS ------------------------------
    //! transform solution and derivatives due to slicing technique
    virtual void TransformSol4Slice(UInt & shiftFactor, UInt & nodeShift,
		UInt & elemgrid, Double &  meshsize, const UInt flag)
    {Error("TransformSol4Slice not implemented", __FILE__, __LINE__); };

    //! save important nodes for the a post analysis
    //! \param shiftFactor number of nodes in x and y direction
    //! \param numShift how many shifts have been made
    //! \param nodeShift how many nodes have been shifted
    //! save solution of special nodes
    virtual void SaveNodes(const UInt shiftFactor, const Double timeStep,
		   const UInt numShift, const Integer nodeShift, 
		   const UInt maxnumelemz_)
    {Error("TransformSol4Slice not implemented",__FILE__, __LINE__); };



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

    //! Set restart time / frequency step
    virtual void SetStartStep( const UInt startStep ) {
      startStep_ = startStep;
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

    //! Start time / frequency step
    UInt startStep_;

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

