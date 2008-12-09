// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_BASESOLVESTEP
#define FILE_BASESOLVESTEP

#include "Utils/StdVector.hh"
#include "General/environment.hh"
#include "Utils/tools.hh"

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
    virtual void PreStepStatic( ) = 0;
    
    /** base method for solving one static step
     * @param comment will be used for filenames if export linear systems */
    virtual void SolveStepStatic(const std::string& comment = "") = 0;

    //! routine for acttions after the SolveStep-method
    virtual void PostStepStatic()  = 0;



    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations before execution of the SolveStep-method
    virtual void PreStepTrans( ) = 0;

    //! routine for computing a predictor step
    // neede in case of FSI-Iterative-Coupling
    //virtual void PredictorStep() = 0;

    //! base method for solving one transient step
    virtual void SolveStepTrans() = 0;

    //! base method for solving one transient step with slicing method
    virtual void SolveStepTrans4Slice()
    {Error("SolveStepTrans4Slice not implemented!",__FILE__,__LINE__);};

    //! routine for actions after the SolveStep-method
    virtual void PostStepTrans() = 0;
    
    //----------------------- HARMONIC---------------------------------------
    //! routine for initilizations before execution of the SolveStep-method
    virtual void PreStepHarmonic() = 0;

    /** base method for solving one harmonic step
     * @param comment - see SolveStepStatic() */
    virtual void SolveStepHarmonic(const std::string& comment = "") = 0;

    //!  routine for actions after the SolveStep-method
    virtual void PostStepHarmonic() = 0;

    //----------------------- HARMONIC ---------------------------------------
    //! Calculate the Eigenfrequencies of a generalized eigenvalue problem
    virtual UInt CalcEigenFrequencies( Vector<Double> & frequencies,
                                       Vector<Double> & errBounds,
                                       UInt numFreq, Double shift ) {
      Error( "Not implemented her!", __FILE__, __LINE__ );
      return 0;
    }

    //! Calculate the Eigenfrequencies of a quadratic eigenvalue problem
    virtual UInt CalcEigenFrequencies( Vector<Complex> & frequencies,
                                       Vector<Double> & errBounds,
                                       UInt numFreq, Double shift ) {
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
    virtual void PreStepTransHarmonic( )
    {Error("PreStepTransHarmonic not implemented!",__FILE__,__LINE__);};

    //! base method for solving one transient-harmonic coupled step
    virtual void SolveStepTransHarmonic()
    {Error("SolveStepTransHarmonic not implemented!",__FILE__,__LINE__);};

    //! routine for actions after the SolveStep-method
    virtual void PostStepTransHarmonic()
    {Error("PostStepTransHarmonic not implemented!",__FILE__,__LINE__);};


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

