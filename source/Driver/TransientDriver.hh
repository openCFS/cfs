#ifndef FILE_TRANSIENTDRIVER_2001
#define FILE_TRANSIENTDRIVER_2001

#include "SingleDriver.hh"

#include <memory>

namespace CoupledField {

  //! forward class declarations
  class Timer;
  class AdaptiveTimesteppingData;

  //! Class for transient simulations
  
  //! This class implements a time dependent problem with a fixed time step
  //! size "dt" and a fixed number of time steps.
  //! It defines the following muParser variables:
  //!   - t    : Current time (in s), always starting at t0
  //!   - dt   : Time step increment (in s)
  //!   - step : Current time step number (always starting at 1)
  //!   - t0   : Initial time (either 0 or accumulated time)
  class TransientDriver : virtual public SingleDriver {

  public:
    //! Constructor
    //! \param sequenceStep current step in multisequence simulation
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    TransientDriver( UInt sequenceStep,
                     bool isPartOfSequence,
                     shared_ptr<SimState> state, Domain* domain,
                     PtrParamNode paramNode, PtrParamNode infoNode );
  
    //! Default destructor
    virtual ~TransientDriver();

    //! Initialization method
    void Init( bool restart);

    //! main method, where time-stepping is implemented. it is for transient and static problem
    void SolveProblem();

    //! Return time increment
    Double GetDeltaT() { return firstdt_;}

    /** Return total number of time steps
     * @see BaseDriver::GetNumSteps() */
    UInt GetNumSteps() { return numstep_; }

    //! Return current time / frequency step of simulation
    UInt GetActStep( const std::string& pdename ) {
      return actTimeStep_;
    }
    
    //! Get total duration of current simulation (in s)
    Double GetDuration() { return firstdt_ * endStep_; }
    
    //! Set accumulated time (may be used as initial time value)
    void SetAccumulatedTime( Double time );

    //! \copydoc SingleDriver::SetToStepValue
    virtual void SetToStepValue(UInt stepNum, Double stepVal );
    
    //! Helper method which determines if an AnalyisType is complex.
    virtual bool IsComplex() { return false; };
    
    //! Static method being called in the case of a Ctr-C signal
    static void SignalHandler( int sig);

    //! Reads stepRejected, toleranceNotReachable, localError from atData_ after each solve; returns true if the step is accepted. Implements retry cap and PI anti-windup.
    bool adaptTimestep(int retryCount);
    

  protected:

    //! Read restart information
    void ReadRestart();
    
    //! Flag, if initial time starts at 0 or is accumulated 
    bool useAccumulatedTime_;

    //! Current simulation time (in s)
    Double actTime_;
    
    //! Accumulated time so far (including time of all previous MS steps)
    Double accTime_;
    
    //! Number of timesteps
    UInt numstep_;

    //! Initial time value used for "t0"
    Double initialTime_;
    
    //! current time step (always starting at initialTime)
    UInt actTimeStep_;
    
    //! Last time step (not necessarily the same as numstep_)
    UInt endStep_; 

    //! Delta t: increment of the time between two steps
    Double firstdt_;

    //! Current time step size (updated each step when adaptive is on)
    Double dt_;

    //! Adaptive only: physical end time of the last accepted step; each attempt solves at t = stepStartTime_ + dt_.
    Double stepStartTime_;

    // =======================================================================
    //  Adaptive timestepping related data
    // =======================================================================

    //! True when the <adaptiveTimeStepping> XML block is present.
    bool adaptiveEnabeled_;

    //! deltaTMin_& deltaTMax_ :  Determin Bounds for adaptivity
    Double deltaTMin_;
    Double deltaTMax_;

    //! tol_: Error Tolerance -> determins when a step has to be rerun
    Double tol_;

     //! Sigma: A factor[0-1], user defines step increase 
    Double sigma_;

    //! restartCount_ : tracks howmany restarts were done
    UInt restartCount_;

    //! Step-size controller type: 0=I, 1=PI (PI.3.4), 2=PID (H312).
    int controllerType_;

    // =======================================================================
    //  Restart related data
    // =======================================================================

    //! Flag, if analysis is restarted
    bool isRestarted_;
    
    //! Flag, if restart file is to be written
    bool writeRestart_;
    
    //! Time step to proceed from when performing restarted simulation
    UInt restartStep_;    
    
    // =======================================================================
    //  Timing estimation
    // =======================================================================

    //! Estimated time per step
    Double timePerStep_;
    
    //! Timer for estimating remaining runtime 
    shared_ptr<Timer> timer_;

    //! Target simulation end time, computed as firstdt_ * numstep_; loop exits once actTime_ >= this value.
    double simulationENDTime_;

    //! Flag, Is set when simulationENDTime_ is reached
    bool simulationEndTimeReached_;

    //! LTE error from the last accepted step; fed to MathParser as prevError for the PI controller's integral term.
    double prevLTEerror_;

    //! LTE error saved at the first rejection of a retry sequence; fed back to the PI controller when toleranceNotReachable is set to prevent integrator windup.
    double antiWindupError_;

    //! Cumulative count of all step rejections across the whole simulation run.
    int totalRetryCount_ = 0;

    //! Shared adaptive timestepping state; null when adaptive is disabled.
    shared_ptr<AdaptiveTimesteppingData> atData_;


  };

}

#endif // FILE_TRANSIENTDRIVER
