#ifndef FILE_TRANSIENTDRIVER_2001
#define FILE_TRANSIENTDRIVER_2001

#include "singleDriver.hh"

namespace CoupledField {

  //! driver for transient problems.it is derived from BaseDriver;
  class TransientDriver : virtual public SingleDriver {

  public:
    //! Constructor
    //! \param adomain pointer to class Domain
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    //! \param driverTag tag for current driver section
    //! \param isPartOfSequence true, if driver is part of  multiSequence
    TransientDriver( Domain * adomain,
                     Integer stepOffset = 0,
                     Double timeOffset = 0.0,
                     std::string driverTag = "anyTag",
                     Boolean isPartOfSequence = FALSE );
  
    //! Default destructor
    virtual ~TransientDriver();

    //! main method, where time-stepping is implemented. it is for transient and static problem
    virtual void SolveProblem();

  protected:

  private:
    //!
    Integer numstep_,isavebegin_,isaveincr_,isaveend_;

    //!
    Double firstdt_;
  };

}

#endif // FILE_TRANSIENTDRIVER
