// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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
    //! \param sequenceStep current step in multisequence simulation
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    TransientDriver( UInt stepOffset,
                     Double timeOffset,
                     UInt sequenceStep,
                     bool isPartOfSequence = false );
  
    //! Default destructor
    virtual ~TransientDriver();

    //! Initialization method
    void Init();

    //! main method, where time-stepping is implemented. it is for transient and static problem
    void SolveProblem();

    //! Return time increment
    Double GetDeltaT() { return firstdt_;}

    //! Return total number of time steps
    UInt GetNumSteps() { return numstep_; }

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
    //  Restart related data
    // =======================================================================

    //! Number of steps before a restart file is stored
    UInt restartIncr_;

    //! Time step to proceed from when performing restarted simulation
    UInt restartStep_;

  };

}

#endif // FILE_TRANSIENTDRIVER
