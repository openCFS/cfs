#ifndef FILE_TRANSIENTDRIVER_2001
#define FILE_TRANSIENTDRIVER_2001

#include "singleDriver.hh"

namespace CoupledField {

  //! driver for transient problems.it is derived from BaseDriver;
  class TransientDriver : virtual public SingleDriver {

  public:
    //! Constructor
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    //! \param driverTag tag for current driver section
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    TransientDriver( UInt stepOffset = 0,
                     Double timeOffset = 0.0,
                     std::string driverTag = "anyTag",
                     bool isPartOfSequence = false );
  
    //! Default destructor
    virtual ~TransientDriver();

    //! Initialization method
    void Init();

    //! main method, where time-stepping is implemented. it is for transient and static problem
    void SolveProblem();

    //! Return time step
    Double GetTimeStep() { return firstdt_;}

    //! Return current time / frequency step of simulation
    UInt GetActStep( const std::string& pdename ) {
      return actTimeStep_;
    }

  protected:

    //! Read restart information
    void ReadRestart();

    //! offset for first timestep (due to multiSequence )
    UInt stepOffset_;

    //! offset for first time (due to multiSequence)
    Double timeOffset_;

    //! Number of timesteps
    UInt numstep_;

    //! current time step
    UInt actTimeStep_;

    //! Delta t: increment of the time between two steps
    Double firstdt_;

    // =======================================================================
    //  Results related information
    // =======================================================================
    
    //! First time step to store results
    UInt isavebegin_;

    //! Last time step to store results
    UInt isaveincr_;

    //! Increment for storing results
    UInt isaveend_;
    
    // =======================================================================
    //  Restart related data
    // =======================================================================

    //! Number of steps before a restart file is stored
    UInt restartIncr_;

    //! Time step to proceed from when performing restarted simulation
    UInt restartStep_;

  };

}

#endif // FILE_TRANSIENTDRIVER
