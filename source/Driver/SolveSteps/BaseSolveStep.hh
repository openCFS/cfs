#ifndef FILE_BASESOLVESTEP
#define FILE_BASESOLVESTEP

#include "Utils/StdVector.hh"
#include "General/Environment.hh"
#include "Utils/tools.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{

  class BaseDriver;
  class InfoNode;
  class AdjointParameters;

  //! Base class for solution of a single step

  class BaseSolveStep
  {

  public:

    //! Destructor
    virtual ~BaseSolveStep();


    //----------------------- STATIC---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepStatic( ) = 0;
    
    /** base method for solving one static step 
     * @param analysis_id references the "base" analysis step. 
     *        is the the info/OLAS/process/step element and required the attribute
     *        "analysis_id" to be set!. In the non-lin case subelements are created. */
    virtual void SolveStepStatic(PtrParamNode analysis_id, AdjointParameters* adjointParams = NULL) = 0;

    //! routine for acttions after the SolveStep-method
    virtual void PostStepStatic()  = 0;

    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations before execution of the SolveStep-method
    virtual void PreStepTrans( ) = 0;

    //! routine for computing a predictor step
    // neede in case of FSI-Iterative-Coupling
    //virtual void PredictorStep() = 0;

    /** base method for solving one transient step
     * @param analysis_id @see SolveStepStatic() */
    virtual void SolveStepTrans(PtrParamNode analysis_id, AdjointParameters* adjointParams = NULL) = 0;

    //! base method for solving one transient step with slicing method
    virtual void SolveStepTrans4Slice()
    {EXCEPTION("SolveStepTrans4Slice not implemented!");};

    //! routine for actions after the SolveStep-method
    virtual void PostStepTrans() = 0;
    
    //! initialize timestepping special variables
    virtual void InitTimeStepping(){
      EXCEPTION("InitTimeStepping not implemented!");
    }

    //----------------------- HARMONIC---------------------------------------
    //! routine for initilizations before execution of the SolveStep-method
    virtual void PreStepHarmonic() = 0;

    /** base method for solving one harmonic step
     * @param analysis_id @see SolveStepStatic() */
    virtual void SolveStepHarmonic(PtrParamNode analysis_id) = 0;

    //!  routine for actions after the SolveStep-method
    virtual void PostStepHarmonic() = 0;

    //----------------------- HARMONIC ---------------------------------------
    //! Calculate the Eigenfrequencies of a generalized eigenvalue problem
    virtual UInt CalcEigenFrequencies( Vector<Double> & frequencies,
                                       Vector<Double> & errBounds,
                                       UInt numFreq, Double shift ) {
      EXCEPTION( "Not implemented her!");
      return 0;
    }

    //! Calculate the Eigenfrequencies of a quadratic eigenvalue problem
    virtual UInt CalcEigenFrequencies( Vector<Complex> & frequencies,
                                       Vector<Double> & errBounds,
                                       UInt numFreq, Double shift ) {
      EXCEPTION( "Not implemented here!" );
      return 0;
    }
    
    //! Calculate the numMode-th eigenmode of a generalized eigenvalue problem.
    //! Therefore, previously CalcEigenFrequencies() has to be called.
    virtual void CalcEigenMode( UInt numMode ) {
      EXCEPTION( "Not implemented here!" );
    }

    //----------------------- TRANSIENTHARMONIC------------------------------
    //! routine for initilizations before execution of the SolveStep-method
    virtual void PreStepTransHarmonic( )
    {EXCEPTION("PreStepTransHarmonic not implemented!");};

    //! base method for solving one transient-harmonic coupled step
    virtual void SolveStepTransHarmonic()
    {EXCEPTION("SolveStepTransHarmonic not implemented!");};

    //! routine for actions after the SolveStep-method
    virtual void PostStepTransHarmonic()
    {EXCEPTION("PostStepTransHarmonic not implemented!");};


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
      numTimeStep_ = numTimeStep;
    };
    
  protected:

    //! Constructor
    BaseSolveStep();

    BaseSolveStep(BaseDriver* driver);

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

