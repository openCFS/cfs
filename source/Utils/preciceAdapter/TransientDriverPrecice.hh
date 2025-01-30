#ifndef FILE_TRANSIENTDRIVERPRECICE
#define FILE_TRANSIENTDRIVERPRECICE

#include "Driver/SingleDriver.hh"
#include "Driver/TransientDriver.hh"

namespace CoupledField {

  //! forward class declarations
  class IPreciceAdapter;

  //! Class for transient simulations

  //! This class implements a time dependent problem with a fixed time step
  //! size "dt" and a fixed number of time steps.
  //! It defines the following muParser variables:
  //!   - t    : Current time (in s), always starting at t0
  //!   - dt   : Time step increment (in s)
  //!   - step : Current time step number (always starting at 1)
  //!   - t0   : Initial time (either 0 or accumulated time)
  class TransientDriverPrecice : virtual public TransientDriver {

  public:
    //! Constructor
    //! \param sequenceStep current step in multisequence simulation
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    TransientDriverPrecice( UInt sequenceStep,
                     bool isPartOfSequence,
                     shared_ptr<SimState> state, Domain* domain,
                     PtrParamNode paramNode, PtrParamNode infoNode );

    //! Default destructor
    virtual ~TransientDriverPrecice();

    void CheckpointingToTimestep(UInt stepNum);

    //! main method, where time-stepping is implemented. it is for transient and static problem
    void SolveProblem();

  protected:
    IPreciceAdapter* preciceAdapter_;
  };

}

#endif // FILE_TRANSIENTDRIVERPRECICE
