#ifndef FILE_TRANSIENTDRIVER_2001
#define FILE_TRANSIENTDRIVER_2001

#include "SingleDriver.hh"

#include <boost/shared_ptr.hpp>

namespace CoupledField {

  //! forward class declarations
  class Timer;

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

    void adaptTimestep();

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

    // =======================================================================
    //  Adaptive timestepping related data
    // =======================================================================

    //! adaptiveEnabeled_ : default false, if time stepping enabeled true
    bool adaptiveEnabeled_;

    //! AdaptiveTimestepping: determins Timestepping Scheme
    std::string adaptiveTimestepping_;

    //! deltaTMin_& deltaTMax_ :  Determin Bounds for adaptivity
    Double deltaTMin_;
    Double deltaTMax_;

    //! tol_: Error Tolerance -> determins when a step has to be rerun
    Double tol_;

     //! Sigma: A factor[0-1], user defines step increase 
    Double sigma_;

    //! restartCount_ : tracks howmany restarts were done
    UInt restartCount_;

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
    boost::shared_ptr<Timer> timer_;

  };

}

#endif // FILE_TRANSIENTDRIVER
