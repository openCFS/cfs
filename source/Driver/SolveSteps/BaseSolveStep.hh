#ifndef FILE_BASESOLVESTEP
#define FILE_BASESOLVESTEP

#include "Utils/StdVector.hh"
#include "General/Environment.hh"
#include "Utils/ToolsFull.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{

  class BaseDriver;
  class BaseVector;

  //! Base class for solution of a single step

  class BaseSolveStep
  {

  public:

    //! Destructor
    virtual ~BaseSolveStep();


    //----------------------- STATIC---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    virtual void PreStepStatic( ) = 0;
    
    /** base method for solving one static step */
    virtual void SolveStepStatic() = 0;

    //! routine for acttions after the SolveStep-method
    virtual void PostStepStatic()  = 0;

    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations before execution of the SolveStep-method
    virtual void PreStepTrans( ) = 0;

    //! routine for computing a predictor step
    // neede in case of FSI-Iterative-Coupling
    //virtual void PredictorStep() = 0;

    /** base method for solving one transient step */
    virtual void SolveStepTrans() = 0;

    //! base method for solving one transient step with slicing method
    virtual void SolveStepTrans4Slice()
    {EXCEPTION("SolveStepTrans4Slice not implemented!");};

    //! routine for actions after the SolveStep-method
    virtual void PostStepTrans() = 0;

    /** Re-apply prescribed external fields (no-op for regular solve steps).
     * Called by TransientDriverPrecice after a coupling window CONVERGED and the
     * preCICE read buffers were refreshed to the converged sample, so that
     * prescribed solutions (e.g. SmoothPDE prescribedDisplacement) reflect the
     * converged data instead of the last-used iteration's read. */
    virtual void RefreshPrescribed() {}
    
    //! initialize timestepping special variables
    virtual void InitTimeStepping(){
      EXCEPTION("InitTimeStepping not implemented!");
    }

    virtual bool InitEigenvalueProblem(const bool isQuadratic=false){
        EXCEPTION("not implemented here");
    }

    //----------------------- HARMONIC---------------------------------------
    //! routine for initilizations before execution of the SolveStep-method
    virtual void PreStepHarmonic() = 0;

    /** base method for solving one harmonic step */
    virtual void SolveStepHarmonic() = 0;

    //!  routine for actions after the SolveStep-method
    virtual void PostStepHarmonic() = 0;

    //! same as GetSoltionVal and GetRHSVal but only in the
    //! multiharmonic case and it's triggered by the MultiHarmonicDriver
    //! in the SolveProblem() method
    virtual void GetSolutionValMultHarm(const UInt& h){
      EXCEPTION("Not implemented here");
    }
    virtual void GetRHSValMultHarm(const UInt& h){
      EXCEPTION("Not implemented here");
    }

    //----------------------- HARMONIC ---------------------------------------
    //! Calculate the Eigenfrequencies of a generalized eigenvalue problem
    virtual UInt CalcEigenFrequencies( Vector<Double>& frequencies, Vector<Double>& errBounds,
                                       UInt numFreq, double shift, bool sort) {
      EXCEPTION( "Not implemented her!");
      return 0;
    }

    /** Calculate the Eigenfrequencies of a quadratic eigenvalue problem
    @param bloch quadratic problem or bloch modes */
    virtual UInt CalcEigenFrequencies( Vector<Complex>& frequencies, Vector<Double>& errBounds,
                                       UInt numFreq, double shift, bool sort, bool bloch) {
      EXCEPTION( "Not implemented here!" );
      return 0;
    }

    //! Calculate the Eigenfrequencies in the interval [minVal,maxVal]
    virtual UInt CalcEigenFrequencies( Vector<Double>& frequencies, Vector<Double>& errBounds,
                                       Double minVal, Double maxVal) {
      EXCEPTION( "Not implemented her!");
      return 0;
    }

    virtual void CalcEigenValues(BaseVector &sol, BaseVector &err, Double minVal, Double maxVal){
        EXCEPTION( "Not implemented here!");
    }

    //! Call for the "solverDefined" eigenvalue selection.
    virtual void CalcEigenValues(BaseVector &sol, BaseVector &err){
        EXCEPTION( "Not implemented here!");
    }

    //! extract the numMode-th eigenmode of a generalized eigenvalue problem.
    //! Therefore, previously CalcEigenFrequencies() has to be called.
    virtual void GetEigenMode( UInt numMode, bool right=true ) {
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


    //----------------------- GetRidOfZeros-----------------------------------
    virtual void SetupGetRidOfZerosActive()
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

    //! Set multiharmonic frequency vector
    virtual void SetMultHarmonicFreq( const StdVector<Double> freqVec ) {
      multHarmFreqVec_ = freqVec;
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
    
    //! Set value for coupling counter
    void SetCouplingIter( UInt count ) {
      couplingIter_ = count;
    }
    
    //! if AMG is used, the auxiliary matrix only needs to be assembled once
    void SetAuxMat(bool set) {
    	auxSet_ = set;
    }

    //! set adjoint
    void SetAdjointSource() {
    	adjointSource_ = true;
    }

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

    //! Vector of considered frequencies (harmonics) for multiharmonic analysis
    StdVector<Double> multHarmFreqVec_;

    //! number of time steps
    UInt numTimeStep_;
    
    //! Counter for iterative coupling
    UInt couplingIter_;

    //! if AMG is used, is auxiliary built?
    bool auxSet_;

    //!
    bool adjointSource_;

  };


} // end of namespace

#endif

